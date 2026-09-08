#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Streaming implementation of ImputeAccure for large input files.

This version keeps only one input line and one sortable output chunk in memory.
It preserves the metrics and output columns of ImputeAccure.py, including the
final position sort and EWMA-based accuracy column.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import heapq
import os
import shutil
import sys
import tempfile
import time
from array import array
from pathlib import Path
from typing import Iterator, TextIO


EPS = sys.float_info.epsilon
ALPHA = 0.1
EWMA_BETA = 1.0 - ALPHA
ALLELE_MAP = {"0": 0, "1": 1, "a": 0, "A": 1}


def parse_bool(value: str) -> bool:
    if value in {"0", "False", "false"}:
        return False
    if value in {"1", "True", "true"}:
        return True
    raise argparse.ArgumentTypeError(
        "use one of 0, False, false, 1, True, or true"
    )


def read_list_file(path: str, cast=str) -> list:
    with open(path, "rt", encoding="utf-8") as handle:
        text = handle.read().strip()
    if not text:
        return []
    return [cast(item) for item in text.replace("\n", "").split(",") if item != ""]


def read_spec_file(path: str) -> list[str]:
    args: list[str] = []
    with open(path, "rt", encoding="utf-8") as spec_file:
        for raw_line in spec_file:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            opt, value = line.split(None, 1)
            args.extend([opt, value])
    return args


def read_panel_file(
    path: str,
    id_column: str,
    panel_column: str,
    unmatched_panel_name: str,
) -> dict[str, list[str]]:
    panels_by_id: dict[str, list[str]] = {}
    with open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"panel file has no header row: {path}")

        fieldnames = set(reader.fieldnames)
        if id_column not in fieldnames:
            raise ValueError(
                "panel file is missing ID column %r; available columns: %s"
                % (id_column, ", ".join(reader.fieldnames))
            )
        if panel_column not in fieldnames:
            if panel_column == "panel" and "Genotype Panel" in fieldnames:
                panel_column = "Genotype Panel"
            else:
                raise ValueError(
                    "panel file is missing panel column %r; available columns: %s"
                    % (panel_column, ", ".join(reader.fieldnames))
                )

        for row_number, row in enumerate(reader, 2):
            indiv_id = (row.get(id_column) or "").strip()
            if not indiv_id:
                continue
            panel = (row.get(panel_column) or "").strip() or unmatched_panel_name
            panels = panels_by_id.setdefault(indiv_id, [])
            if panel not in panels:
                panels.append(panel)
    return panels_by_id


def safe_filename_part(value: str) -> str:
    cleaned = []
    for char in value:
        if char.isalnum() or char in {".", "-", "_"}:
            cleaned.append(char)
        else:
            cleaned.append("_")
    safe = "".join(cleaned).strip("._")
    return safe or "panel"


def output_path_for_panel(
    output_path: Path,
    panel: str,
    panel_output_dir: Path | None,
) -> Path:
    output_dir = panel_output_dir if panel_output_dir is not None else output_path.parent
    suffix = output_path.suffix or ".accuracy"
    stem = output_path.stem if output_path.suffix else output_path.name
    return output_dir / f"{stem}.{safe_filename_part(panel)}{suffix}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Memory-efficient ImputeAccure analysis for very large files."
    )
    parser.add_argument("-f", "--specificationFile")
    parser.add_argument(
        "-i",
        "--input",
        help=(
            "input file; for genotype input this may be a comma-separated "
            "list of files"
        ),
    )
    parser.add_argument(
        "--input-files",
        nargs="+",
        help=(
            "multiple genotype input files with the same locus columns; "
            "genotype cells are merged by sample ID and locus"
        ),
    )
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "--input-format",
        choices=("probability", "haplotype", "genotype"),
        default="probability",
        help=(
            "probability: ImputeAccure marker-row input; haplotype: two "
            "haplotype lines per individual with ID plus allele columns; "
            "genotype: one row per individual with ID plus 0/1/2 genotype calls"
        ),
    )
    parser.add_argument("-l", "--leadingColumns", type=int, default=5)
    parser.add_argument("-c", "--calcThirdColumn", type=parse_bool, default=False)
    parser.add_argument("-p", "--pplToExclude")
    parser.add_argument("-m", "--markersToExclude")
    parser.add_argument("-n", "--columnNames")
    parser.add_argument(
        "--panel-file",
        help=(
            "CSV file mapping individual IDs to panels; supported for "
            "haplotype and genotype input formats"
        ),
    )
    parser.add_argument(
        "--panel-column",
        default="panel",
        help="panel column name in --panel-file; auto-detects 'Genotype Panel'",
    )
    parser.add_argument(
        "--id-column",
        default="ID",
        help="individual ID column name in --panel-file",
    )
    parser.add_argument(
        "--unmatched-panel-name",
        default="not_in_panel",
        help="group name for input IDs absent from --panel-file or with blank panel",
    )
    parser.add_argument(
        "--panel-output-dir",
        help="directory for per-panel accuracy outputs; defaults to output directory",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=250_000,
        help="number of per-marker result rows sorted in memory at once",
    )
    parser.add_argument(
        "--temp-dir",
        help="directory for temporary sorted chunks; defaults to the output directory",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100_000,
        help="print progress every N processed markers; 0 disables progress",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="keep temporary chunk directory for debugging",
    )
    return parser


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = build_parser()
    first_pass, _ = parser.parse_known_args(argv)
    if first_pass.specificationFile:
        spec_args = read_spec_file(first_pass.specificationFile)
        override_args: list[str] = []
        skip_next = False
        for index, item in enumerate(argv):
            if skip_next:
                skip_next = False
                continue
            if item in {"-f", "--specificationFile"}:
                skip_next = index + 1 < len(argv)
                continue
            override_args.append(item)
        argv = spec_args + override_args
    return parser.parse_args(argv)


def input_paths_from_args(args: argparse.Namespace) -> list[Path]:
    if args.input_files:
        raw_paths = args.input_files
    elif args.input:
        raw_paths = [
            item.strip()
            for item in args.input.split(",")
            if item.strip()
        ]
    else:
        raw_paths = []
    return [Path(path) for path in raw_paths]


def open_input(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")  # type: ignore[return-value]
    return open(path, "rt", encoding="utf-8")


def best_guess_index(p0: float, p1: float, p2: float) -> float:
    max_value = max(p0, p1, p2)
    winners_sum = 0.0
    winners_count = 0
    if p0 == max_value:
        winners_count += 1
    if p1 == max_value:
        winners_sum += 1.0
        winners_count += 1
    if p2 == max_value:
        winners_sum += 2.0
        winners_count += 1
    return winners_sum / winners_count


def add_best_guess_counts(
    p0: float, p1: float, p2: float, best_guess_mean: list[float]
) -> None:
    max_value = max(p0, p1, p2)
    winners = 0
    if p0 == max_value:
        winners += 1
    if p1 == max_value:
        winners += 1
    if p2 == max_value:
        winners += 1
    best_guess_mean[0] += 1.0 / winners if p0 == max_value else 0.0
    best_guess_mean[1] += 1.0 / winners if p1 == max_value else 0.0
    best_guess_mean[2] += 1.0 / winners if p2 == max_value else 0.0


def finalize_metrics(
    total_sum: float,
    f_a_total: float,
    q_mean: float,
    dosage_mean: list[float],
    best_guess_mean: list[float],
    e: float,
    e2: float,
    f: float,
    fe2: float,
    z: float,
    z2: float,
    ze: float,
) -> tuple[float, float, int, float, float, float, float, float]:
    n_indiv = int(round(total_sum, 0))
    f_a = f_a_total / (2.0 * n_indiv) if n_indiv > 0 else 0.0
    if f_a == 0.0:
        f_a = EPS
    elif f_a == 1.0:
        f_a = 1.0 - EPS

    q_hwe = -2.0 * f_a * (f_a - 1.0) * (3.0 * f_a**2 - 3.0 * f_a + 2.0)

    if n_indiv > 0:
        q_mean /= n_indiv
        dosage_mean = [value / n_indiv for value in dosage_mean]
        best_guess_mean = [value / n_indiv for value in best_guess_mean]

        iam_chance = 1.0 - q_mean * 3.0 / 2.0
        iam_hwe = 1.0 - q_mean / q_hwe
        similarity = sum(
            (best_guess_mean[i] ** 0.5) * (dosage_mean[i] ** 0.5)
            for i in range(3)
        )
        hiq = 1.0 - (1.0 - min(1.0, similarity)) ** 0.5
        maf = f_a * 100.0

        denom = 2.0 * f_a * (1.0 - f_a)
        r2_mach = ((e2 / n_indiv) - (e / n_indiv) ** 2) / denom if denom > 0 else -9

        beagle_denom = (f - (e**2 / n_indiv)) * (z2 - z**2 / n_indiv)
        if beagle_denom > 0:
            r2_beagle = ((ze - z * e / n_indiv) ** 2) / (
                f - (e**2 / n_indiv)
            ) / (z2 - z**2 / n_indiv)
        else:
            r2_beagle = -9

        impute_denom = 2.0 * n_indiv * f_a * (1.0 - f_a)
        info_impute = 1.0 - fe2 / impute_denom if impute_denom > 0 else -9
    else:
        iam_chance = iam_hwe = hiq = maf = info_impute = r2_mach = r2_beagle = -9

    return (
        round(iam_chance, 3),
        round(iam_hwe, 3),
        round(n_indiv, 0),
        round(maf, 1),
        round(hiq, 3),
        round(info_impute, 3),
        round(r2_mach, 3),
        round(r2_beagle, 3),
    )


def process_probabilities(
    values: list[str],
    calculate_third_column: bool,
    exclude_ppl: set[int],
) -> tuple[float, float, int, float, float, float, float, float]:
    stride = 2 if calculate_third_column else 3
    total_sum = 0.0
    f_a_total = 0.0
    q_mean = 0.0
    dosage_mean = [0.0, 0.0, 0.0]
    best_guess_mean = [0.0, 0.0, 0.0]
    e = e2 = f = fe2 = z = z2 = ze = 0.0

    sample_index = 0
    for offset in range(0, len(values), stride):
        if sample_index in exclude_ppl:
            p0 = p1 = p2 = 0.0
        elif calculate_third_column:
            p0 = float(values[offset])
            p1 = float(values[offset + 1])
            p2 = 1.0 - (p0 + p1)
        else:
            p0 = float(values[offset])
            p1 = float(values[offset + 1])
            p2 = float(values[offset + 2])

        if p0 < 0.0:
            p0 = p1 = p2 = 0.0

        row_sum = p0 + p1 + p2
        if row_sum != 1.0 and row_sum > 0.0:
            p0 /= row_sum
            p1 /= row_sum
            p2 /= row_sum
            row_sum = p0 + p1 + p2

        q_mean += p0 * (1.0 - p0) + p1 * (1.0 - p1) + p2 * (1.0 - p2)
        f_a_total += p1 + 2.0 * p2
        dosage_mean[0] += p0
        dosage_mean[1] += p1
        dosage_mean[2] += p2

        if row_sum > 0.0:
            add_best_guess_counts(p0, p1, p2, best_guess_mean)

        z_current = best_guess_index(p0, p1, p2)
        expected = p1 + 2.0 * p2
        expected_square = expected**2
        f_current = p1 + 4.0 * p2
        z += z_current
        z2 += z_current**2
        e += expected
        e2 += expected_square
        f += f_current
        fe2 += f_current - expected_square
        ze += z_current * expected
        total_sum += row_sum
        sample_index += 1

    return finalize_metrics(
        total_sum,
        f_a_total,
        q_mean,
        dosage_mean,
        best_guess_mean,
        e,
        e2,
        f,
        fe2,
        z,
        z2,
        ze,
    )


def allele_values(parts: list[str], n_loci: int, line_no: int) -> list[int]:
    if len(parts) != n_loci + 1:
        raise ValueError(
            "line %i has %i loci; expected %i"
            % (line_no, len(parts) - 1, n_loci)
        )

    values: list[int] = []
    for token in parts[1:]:
        allele = ALLELE_MAP.get(token)
        if allele is None:
            raise ValueError(
                "line %i has invalid allele %r; allowed values are 0, 1, a, A"
                % (line_no, token)
            )
        values.append(allele)
    return values


def genotype_values(parts: list[str], n_loci: int, line_no: int) -> list[int]:
    if len(parts) != n_loci + 1:
        raise ValueError(
            "line %i has %i loci; expected %i"
            % (line_no, len(parts) - 1, n_loci)
        )

    values: list[int] = []
    # Column 1 is the sample ID; every later column is a real locus.
    for token in parts[1:]:
        if token not in {"0", "1", "2", "9"}:
            raise ValueError(
                "line %i has invalid genotype %r; allowed values are 0, 1, 2, 9"
                % (line_no, token)
            )
        values.append(int(token))
    return values


def nonempty_lines(handle: TextIO) -> Iterator[tuple[int, str]]:
    for line_no, raw_line in enumerate(handle, 1):
        if raw_line.strip():
            yield line_no, raw_line


def combine_genotype_values(
    values_by_file: list[list[int]],
    line_no: int,
    sample_id: str,
) -> list[int]:
    combined: list[int] = []
    for locus_index, values in enumerate(zip(*values_by_file), 1):
        genotype = values[0]
        for next_genotype in values[1:]:
            if genotype == next_genotype:
                continue
            if genotype == 9 and next_genotype != 9:
                genotype = next_genotype
                continue
            if next_genotype == 9:
                continue
            # Conflicting known values keep the earlier input file's value.
        combined.append(genotype)
    return combined


def update_haplotype_counts(
    count0: array,
    count1: array,
    count2: array,
    paternal: list[int],
    maternal: list[int],
    active_loci: range | list[int],
) -> None:
    for locus_index in active_loci:
        genotype = paternal[locus_index] + maternal[locus_index]
        if genotype == 0:
            count0[locus_index] += 1
        elif genotype == 1:
            count1[locus_index] += 1
        else:
            count2[locus_index] += 1


def update_genotype_counts(
    count0: array,
    count1: array,
    count2: array,
    genotypes: list[int],
    active_loci: range | list[int],
    locus_offset: int = 0,
) -> None:
    for locus_index in active_loci:
        genotype = genotypes[locus_index]
        if genotype == 9:
            continue
        target_index = locus_offset + locus_index
        if genotype == 0:
            count0[target_index] += 1
        elif genotype == 1:
            count1[target_index] += 1
        else:
            count2[target_index] += 1


def metrics_from_haplotype_counts(
    count0: int,
    count1: int,
    count2: int,
) -> tuple[float, float, int, float, float, float, float, float]:
    total_sum = float(count0 + count1 + count2)
    f_a_total = float(count1 + 2 * count2)
    dosage_mean = [float(count0), float(count1), float(count2)]
    best_guess_mean = [float(count0), float(count1), float(count2)]
    e = f_a_total
    e2 = float(count1 + 4 * count2)
    f = e2
    z = e
    z2 = e2
    ze = e2
    return finalize_metrics(
        total_sum,
        f_a_total,
        0.0,
        dosage_mean,
        best_guess_mean,
        e,
        e2,
        f,
        0.0,
        z,
        z2,
        ze,
    )


def prepare_columns(column_names: str | None, leading_columns: int) -> list[str]:
    if column_names:
        columns = column_names.split(",")
    else:
        columns = []

    if len(columns) != leading_columns:
        if columns:
            print(
                "Number of column headers does not correspond to specified "
                "number of leading columns! Switching to default names!"
            )
        columns = [f"col{i}" for i in range(leading_columns)]

    if columns[1] != "SNP":
        columns[1] = "SNP"
        print("Renamed second column header to SNP!")
    if columns[2] != "position":
        columns[2] = "position"
        print("Renamed third column header to position!")
    return columns


def write_chunk(
    chunk: list[tuple[int, int, str]],
    temp_dir: Path,
    chunk_index: int,
) -> Path:
    chunk.sort(key=lambda item: (item[0], item[1]))
    chunk_path = temp_dir / f"chunk_{chunk_index:06d}.tmp"
    with open(chunk_path, "wt", encoding="utf-8") as handle:
        for position, row_index, row_text in chunk:
            handle.write(f"{position}\t{row_index}\t{row_text}\n")
    return chunk_path


def sorted_metric_rows(chunk_paths: list[Path]) -> Iterator[str]:
    files = [open(path, "rt", encoding="utf-8") for path in chunk_paths]
    try:
        heap: list[tuple[int, int, int, str]] = []
        for file_index, handle in enumerate(files):
            line = handle.readline()
            if not line:
                continue
            position, row_index, row_text = line.rstrip("\n").split("\t", 2)
            heap.append((int(position), int(row_index), file_index, row_text))
        heapq.heapify(heap)

        while heap:
            _position, _row_index, file_index, row_text = heapq.heappop(heap)
            yield row_text
            next_line = files[file_index].readline()
            if next_line:
                position, row_index, next_row_text = next_line.rstrip("\n").split(
                    "\t", 2
                )
                heapq.heappush(
                    heap, (int(position), int(row_index), file_index, next_row_text)
                )
    finally:
        for handle in files:
            handle.close()


def classify_accuracy(hwe_value: float, hiq_value: float, state: dict[str, float]) -> str:
    hwe_for_ewma = 0.47 if hwe_value == -9.0 else hwe_value
    hiq_for_ewma = 0.97 if hiq_value == -0.9 else hiq_value

    state["hwe_num"] = hwe_for_ewma + EWMA_BETA * state["hwe_num"]
    state["hiq_num"] = hiq_for_ewma + EWMA_BETA * state["hiq_num"]
    state["den"] = 1.0 + EWMA_BETA * state["den"]

    hwe_ewma = round(state["hwe_num"] / state["den"], 3)
    hiq_ewma = round(state["hiq_num"] / state["den"], 3)

    if hwe_ewma > 0.47 and hiq_ewma > 0.97:
        return "cold"
    if hwe_ewma < 0.47 / 2.0 and hiq_ewma < 0.97 / 2.0:
        return "very hot"
    if hwe_ewma < 0.47 and hiq_ewma < 0.97:
        return "hot"
    return "tepid"


def add_accuracy_column(row_text: str, leading_columns: int, state: dict[str, float]) -> str:
    fields = row_text.split(" ")
    hwe_value = float(fields[leading_columns + 3])
    hiq_value = float(fields[leading_columns + 4])
    accuracy = classify_accuracy(hwe_value, hiq_value, state)
    fields.insert(leading_columns + 5, accuracy)
    return " ".join(fields)


def format_metric_row(
    names: list[str],
    metrics: tuple[float, float, int, float, float, float, float, float],
) -> str:
    (
        iam_chance,
        iam_hwe,
        n,
        maf,
        hiq,
        info_impute,
        r2_mach,
        r2_beagle,
    ) = metrics
    return (
        " ".join(names)
        + f" {int(n)} {maf}% {iam_chance} {iam_hwe} {hiq} "
        + f"{info_impute} {r2_mach} {r2_beagle}"
    )


def make_count_arrays(n_loci: int) -> tuple[array, array, array]:
    return (
        array("Q", [0]) * n_loci,
        array("Q", [0]) * n_loci,
        array("Q", [0]) * n_loci,
    )


def write_count_accuracy(
    output_path: Path,
    header: str,
    count0: array,
    count1: array,
    count2: array,
    exclude_markers: set[str],
) -> tuple[int, int]:
    row_counter = 0
    skipped_markers = 0
    state = {"hwe_num": 0.0, "hiq_num": 0.0, "den": 0.0}
    with open(output_path, "wt", encoding="utf-8") as results:
        results.write("created at ")
        results.write(_dt.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        results.write("\n")
        results.write(header)

        for locus_index in range(len(count0)):
            marker = f"loci_{locus_index + 1}"
            if marker in exclude_markers:
                skipped_markers += 1
                continue
            metrics = metrics_from_haplotype_counts(
                count0[locus_index],
                count1[locus_index],
                count2[locus_index],
            )
            row_text = format_metric_row(
                [marker, marker, str(locus_index + 1)],
                metrics,
            )
            results.write(add_accuracy_column(row_text, 3, state))
            results.write("\n")
            row_counter += 1
    return row_counter, skipped_markers


def write_grouped_count_accuracy(
    output_path: Path,
    header: str,
    grouped_counts: dict[str, tuple[array, array, array]],
    exclude_markers: set[str],
    panel_output_dir: Path | None,
) -> tuple[int, int, dict[str, Path]]:
    total_rows = 0
    skipped_markers = 0
    outputs: dict[str, Path] = {}
    if panel_output_dir is not None:
        panel_output_dir.mkdir(parents=True, exist_ok=True)

    for panel in sorted(grouped_counts):
        count0, count1, count2 = grouped_counts[panel]
        panel_path = output_path_for_panel(output_path, panel, panel_output_dir)
        row_counter, skipped = write_count_accuracy(
            panel_path,
            header,
            count0,
            count1,
            count2,
            exclude_markers,
        )
        outputs[panel] = panel_path
        total_rows += row_counter
        skipped_markers = skipped

    return total_rows, skipped_markers, outputs


def get_group_counts(
    grouped_counts: dict[str, tuple[array, array, array]],
    panel: str,
    n_loci: int,
) -> tuple[array, array, array]:
    counts = grouped_counts.get(panel)
    if counts is None:
        counts = make_count_arrays(n_loci)
        grouped_counts[panel] = counts
    return counts


def process_haplotype_input(
    input_path: Path,
    output_path: Path,
    header: str,
    exclude_ppl: set[int],
    exclude_markers: set[str],
    progress_every: int,
    panel_by_id: dict[str, list[str]] | None,
    unmatched_panel_name: str,
    panel_output_dir: Path | None,
) -> tuple[int, int, int, dict[str, Path]]:
    n_loci: int | None = None
    active_loci: range | list[int] = range(0)
    count0: array | None = None
    count1: array | None = None
    count2: array | None = None
    grouped_counts: dict[str, tuple[array, array, array]] = {}
    paternal: list[int] | None = None
    paternal_id: str | None = None
    sample_index = 0
    tic = time.perf_counter()

    with open_input(input_path) as file:
        for line_no, raw_line in enumerate(file, 1):
            if not raw_line.strip():
                continue

            parts = raw_line.split()
            if len(parts) < 2:
                raise ValueError(f"line {line_no}: expected ID plus at least one allele")

            if n_loci is None:
                n_loci = len(parts) - 1
                if panel_by_id is None:
                    count0, count1, count2 = make_count_arrays(n_loci)
                if exclude_markers:
                    active_loci = [
                        index
                        for index in range(n_loci)
                        if f"loci_{index + 1}" not in exclude_markers
                    ]
                else:
                    active_loci = range(n_loci)
                print(f"Number of loci: {n_loci}")

            alleles = allele_values(parts, n_loci, line_no)
            if paternal is None:
                paternal = alleles
                paternal_id = parts[0]
                continue

            if paternal_id != parts[0]:
                raise ValueError(
                    "haplotype pair ending at line %i has mismatched IDs: %r and %r"
                    % (line_no, paternal_id, parts[0])
                )

            if sample_index not in exclude_ppl:
                if panel_by_id is None:
                    assert count0 is not None and count1 is not None and count2 is not None
                    target_counts_list = [(count0, count1, count2)]
                else:
                    panels = panel_by_id.get(paternal_id or "", [unmatched_panel_name])
                    target_counts_list = [
                        get_group_counts(grouped_counts, panel, n_loci)
                        for panel in panels
                    ]
                for target_counts in target_counts_list:
                    update_haplotype_counts(
                        target_counts[0],
                        target_counts[1],
                        target_counts[2],
                        paternal,
                        alleles,
                        active_loci,
                    )
            sample_index += 1
            paternal = None
            paternal_id = None

            if progress_every and sample_index % progress_every == 0:
                elapsed = time.perf_counter() - tic
                print(f"Processed {sample_index:,} individuals in {elapsed:.1f}s")

    if n_loci is None:
        raise ValueError("input file has no data lines")
    if paternal is not None:
        raise ValueError("file must have an even number of non-empty haplotype lines")

    invalid_exclusions = {
        value for value in exclude_ppl if value < 0 or value >= sample_index
    }
    if invalid_exclusions:
        print("Invalid ID(s) among excluded samples detected! They were ignored!")

    print(f"Number of samples: {sample_index}")
    if panel_by_id is None:
        assert count0 is not None and count1 is not None and count2 is not None
        row_counter, skipped_markers = write_count_accuracy(
            output_path,
            header,
            count0,
            count1,
            count2,
            exclude_markers,
        )
        return row_counter, skipped_markers, sample_index, {"all": output_path}

    row_counter, skipped_markers, outputs = write_grouped_count_accuracy(
        output_path,
        header,
        grouped_counts,
        exclude_markers,
        panel_output_dir,
    )
    return row_counter, skipped_markers, sample_index, outputs


def process_genotype_input(
    input_path: Path,
    output_path: Path,
    header: str,
    exclude_ppl: set[int],
    exclude_markers: set[str],
    progress_every: int,
    panel_by_id: dict[str, list[str]] | None,
    unmatched_panel_name: str,
    panel_output_dir: Path | None,
) -> tuple[int, int, int, dict[str, Path]]:
    return process_genotype_inputs(
        [input_path],
        output_path,
        header,
        exclude_ppl,
        exclude_markers,
        progress_every,
        panel_by_id,
        unmatched_panel_name,
        panel_output_dir,
    )


def process_genotype_inputs(
    input_paths: list[Path],
    output_path: Path,
    header: str,
    exclude_ppl: set[int],
    exclude_markers: set[str],
    progress_every: int,
    panel_by_id: dict[str, list[str]] | None,
    unmatched_panel_name: str,
    panel_output_dir: Path | None,
) -> tuple[int, int, int, dict[str, Path]]:
    n_loci: int | None = None
    active_loci: range | list[int] = range(0)
    count0, count1, count2 = make_count_arrays(0)
    grouped_counts: dict[str, tuple[array, array, array]] = {}
    combined_by_id: dict[str, list[int]] = {}
    sample_order: list[str] = []
    rows_read_by_file: list[int] = []
    duplicate_rows_by_file: list[int] = []
    tic = time.perf_counter()

    for file_index, input_path in enumerate(input_paths, 1):
        rows_read = 0
        duplicate_rows = 0
        seen_ids: set[str] = set()
        with open_input(input_path) as file:
            for line_no, raw_line in nonempty_lines(file):
                parts = raw_line.split()
                if len(parts) < 2:
                    raise ValueError(
                        "%s line %i: expected ID plus at least one genotype"
                        % (input_path, line_no)
                    )

                if n_loci is None:
                    n_loci = len(parts) - 1
                    count0.extend(array("Q", [0]) * n_loci)
                    count1.extend(array("Q", [0]) * n_loci)
                    count2.extend(array("Q", [0]) * n_loci)
                    if exclude_markers:
                        active_loci = [
                            index
                            for index in range(n_loci)
                            if f"loci_{index + 1}" not in exclude_markers
                        ]
                    else:
                        active_loci = range(n_loci)
                    print(f"Number of genotype input files: {len(input_paths)}")
                    print("ID column: 1; first locus column: 2")
                    print(f"Loci per file: {n_loci}")

                sample_id = parts[0]
                if sample_id in seen_ids:
                    duplicate_rows += 1
                    continue
                seen_ids.add(sample_id)

                genotypes = genotype_values(parts, n_loci, line_no)
                if sample_id in combined_by_id:
                    combined_by_id[sample_id] = combine_genotype_values(
                        [combined_by_id[sample_id], genotypes],
                        line_no,
                        sample_id,
                    )
                else:
                    combined_by_id[sample_id] = genotypes
                    sample_order.append(sample_id)

                rows_read += 1
                if progress_every and rows_read % progress_every == 0:
                    elapsed = time.perf_counter() - tic
                    print(
                        "Read %,d individuals from %s in %.1fs"
                        % (rows_read, input_path.name, elapsed)
                    )

        if rows_read == 0:
            raise ValueError(f"input file has no data lines: {input_path}")
        rows_read_by_file.append(rows_read)
        duplicate_rows_by_file.append(duplicate_rows)
        print(
            "Input %i/%i: %s; rows read: %i; duplicate rows ignored: %i"
            % (
                file_index,
                len(input_paths),
                input_path,
                rows_read,
                duplicate_rows,
            )
        )

    if n_loci is None:
        raise ValueError("genotype input files have no data lines")

    sample_count = len(sample_order)
    invalid_exclusions = {
        value for value in exclude_ppl if value < 0 or value >= sample_count
    }
    if invalid_exclusions:
        print("Invalid ID(s) among excluded samples detected! They were ignored!")

    for sample_index, sample_id in enumerate(sample_order):
        if sample_index in exclude_ppl:
            continue
        genotypes = combined_by_id[sample_id]
        if panel_by_id is None:
            target_counts_list = [(count0, count1, count2)]
        else:
            panels = panel_by_id.get(sample_id, [unmatched_panel_name])
            target_counts_list = [
                get_group_counts(grouped_counts, panel, n_loci)
                for panel in panels
            ]
        for target_counts in target_counts_list:
            update_genotype_counts(
                target_counts[0],
                target_counts[1],
                target_counts[2],
                genotypes,
                active_loci,
            )

    print(f"Number of loci: {n_loci}")
    print(
        "Rows read per file: "
        + ", ".join(str(value) for value in rows_read_by_file)
    )
    print(
        "Duplicate rows ignored per file: "
        + ", ".join(str(value) for value in duplicate_rows_by_file)
    )
    print(f"Unique samples after combining: {sample_count}")
    if panel_by_id is None:
        row_counter, skipped_markers = write_count_accuracy(
            output_path,
            header,
            count0,
            count1,
            count2,
            exclude_markers,
        )
        return row_counter, skipped_markers, sample_count, {"all": output_path}

    row_counter, skipped_markers, outputs = write_grouped_count_accuracy(
        output_path,
        header,
        grouped_counts,
        exclude_markers,
        panel_output_dir,
    )
    return row_counter, skipped_markers, sample_count, outputs


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    if args.input_format in {"haplotype", "genotype"}:
        args.leadingColumns = 3
        if not args.columnNames:
            args.columnNames = "SNP_no,SNP,position"

    if args.leadingColumns < 3:
        print("There needs to be at least 3 leading columns!", file=sys.stderr)
        return 1
    input_paths = input_paths_from_args(args)
    if not input_paths:
        print("Please specify an input file!", file=sys.stderr)
        return 1
    if len(input_paths) > 1 and args.input_format != "genotype":
        print(
            "Multiple input files are currently supported only with "
            "--input-format genotype",
            file=sys.stderr,
        )
        return 1
    if args.chunk_rows < 1:
        print("--chunk-rows must be at least 1", file=sys.stderr)
        return 1

    for input_path in input_paths:
        if not input_path.is_file():
            print(f"Input file DOES NOT exist: {input_path}", file=sys.stderr)
            return 1
    input_path = input_paths[0]

    output_path = Path(args.output) if args.output else input_path.with_suffix("")
    if not args.output:
        if input_path.suffix == ".gz":
            output_path = output_path.with_suffix(".accuracy")
        else:
            output_path = input_path.with_suffix(".accuracy")

    exclude_ppl = set()
    if args.pplToExclude:
        if not os.path.isfile(args.pplToExclude):
            print("File with samples to exclude DOES NOT exist!", file=sys.stderr)
            return 1
        exclude_ppl = {value - 1 for value in read_list_file(args.pplToExclude, int)}

    exclude_markers = set()
    if args.markersToExclude:
        if not os.path.isfile(args.markersToExclude):
            print("File with markers to exclude DOES NOT exist!", file=sys.stderr)
            return 1
        exclude_markers = set(read_list_file(args.markersToExclude, str))

    panel_by_id = None
    panel_output_dir = Path(args.panel_output_dir) if args.panel_output_dir else None
    if args.panel_file:
        if args.input_format not in {"haplotype", "genotype"}:
            print(
                "--panel-file is only supported with haplotype or genotype input",
                file=sys.stderr,
            )
            return 1
        if not os.path.isfile(args.panel_file):
            print(f"Panel CSV file DOES NOT exist: {args.panel_file}", file=sys.stderr)
            return 1
        try:
            panel_by_id = read_panel_file(
                args.panel_file,
                args.id_column,
                args.panel_column,
                args.unmatched_panel_name,
            )
        except ValueError as exc:
            print(f"Panel CSV error: {exc}", file=sys.stderr)
            return 1

    columns = prepare_columns(args.columnNames, args.leadingColumns)
    header = (
        " ".join(columns)
        + " N MAF Iam_chance Iam_hwe hiQ accuracy info_IMPUTE r2_MACH r2_BEAGLE\n"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_parent = Path(args.temp_dir) if args.temp_dir else output_path.parent
    temp_parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="imputeaccure_large_", dir=temp_parent))

    print("--- Specifications        ---------------------------")
    print("  input format:                     ", args.input_format)
    print("  leading columns:                  ", args.leadingColumns)
    print("  calculate third column/probability:", args.calcThirdColumn)
    print("  exclude marker:                   ", len(exclude_markers), "entries")
    print("  exclude probes:                   ", len(exclude_ppl), "entries")
    if panel_by_id is not None:
        print("  panel file IDs:                   ", len(panel_by_id))
        print(
            "  IDs in multiple panels:           ",
            sum(1 for panels in panel_by_id.values() if len(panels) > 1),
        )
    print("---- Data and Result file: -------------------------")
    if len(input_paths) == 1:
        print("  input-file:    ", input_path)
    else:
        print("  input-files:")
        for path in input_paths:
            print("    ", path)
    print("  output-file:   ", output_path)
    if panel_output_dir is not None:
        print("  panel-output-dir:", panel_output_dir)
    if args.input_format == "probability":
        print("  temp-dir:      ", temp_dir)
    print("--- ---------------------- -------------------------")

    if args.input_format in {"haplotype", "genotype"}:
        tic = time.perf_counter()
        try:
            if args.input_format == "haplotype":
                (
                    row_counter,
                    skipped_markers,
                    sample_count,
                    grouped_outputs,
                ) = process_haplotype_input(
                    input_path,
                    output_path,
                    header,
                    exclude_ppl,
                    exclude_markers,
                    args.progress_every,
                    panel_by_id,
                    args.unmatched_panel_name,
                    panel_output_dir,
                )
            else:
                (
                    row_counter,
                    skipped_markers,
                    sample_count,
                    grouped_outputs,
                ) = process_genotype_inputs(
                    input_paths,
                    output_path,
                    header,
                    exclude_ppl,
                    exclude_markers,
                    args.progress_every,
                    panel_by_id,
                    args.unmatched_panel_name,
                    panel_output_dir,
                )
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

        runtime = time.perf_counter() - tic
        print("\nRuntime: %f" % runtime)
        print("Markers written: %i" % row_counter)
        print("Markers skipped: %i" % skipped_markers)
        print("Samples read: %i" % sample_count)
        if panel_by_id is not None:
            print("Panel outputs:")
            for panel in sorted(grouped_outputs):
                print(f"  {panel}: {grouped_outputs[panel]}")
        if row_counter:
            print("Average runtime per marker: %f" % (runtime / row_counter))
        return 0

    row_counter = 0
    skipped_markers = 0
    chunk_index = 0
    chunk: list[tuple[int, int, str]] = []
    chunk_paths: list[Path] = []
    num_samples: int | None = None
    stride = 2 if args.calcThirdColumn else 3

    tic = time.perf_counter()
    try:
        with open_input(input_path) as file:
            for raw_row in file:
                if not raw_row.strip():
                    continue

                row_data = raw_row.split()
                if len(row_data) <= args.leadingColumns:
                    raise ValueError(f"row {row_counter + skipped_markers + 1} has no probability columns")

                if num_samples is None:
                    remainder = len(row_data) - args.leadingColumns
                    if remainder % stride != 0:
                        raise ValueError("Invalid column specifications!")
                    num_samples = remainder // stride
                    invalid_exclusions = {
                        value for value in exclude_ppl if value < 0 or value >= num_samples
                    }
                    if invalid_exclusions:
                        print(
                            "Invalid ID(s) among excluded samples detected! "
                            "They will be ignored!"
                        )
                        exclude_ppl -= invalid_exclusions
                    print(f"Number of samples: {num_samples}")
                elif len(row_data) != args.leadingColumns + num_samples * stride:
                    raise ValueError(
                        "row %i has %i columns; expected %i"
                        % (
                            row_counter + skipped_markers + 1,
                            len(row_data),
                            args.leadingColumns + num_samples * stride,
                        )
                    )

                if row_data[1] in exclude_markers:
                    skipped_markers += 1
                    continue

                names = row_data[: args.leadingColumns]
                values = row_data[args.leadingColumns :]
                metrics = process_probabilities(values, args.calcThirdColumn, exclude_ppl)
                row_text = format_metric_row(names, metrics)
                position = int(names[2])
                chunk.append((position, row_counter, row_text))
                row_counter += 1

                if len(chunk) >= args.chunk_rows:
                    chunk_paths.append(write_chunk(chunk, temp_dir, chunk_index))
                    chunk = []
                    chunk_index += 1

                if args.progress_every and row_counter % args.progress_every == 0:
                    elapsed = time.perf_counter() - tic
                    print(f"Processed {row_counter:,} markers in {elapsed:.1f}s")

        if chunk:
            chunk_paths.append(write_chunk(chunk, temp_dir, chunk_index))

        state = {"hwe_num": 0.0, "hiq_num": 0.0, "den": 0.0}
        with open(output_path, "wt", encoding="utf-8") as results:
            results.write("created at ")
            results.write(_dt.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            results.write("\n")
            results.write(header)
            for row_text in sorted_metric_rows(chunk_paths):
                results.write(add_accuracy_column(row_text, args.leadingColumns, state))
                results.write("\n")

    finally:
        if args.keep_temp:
            print(f"Temporary files kept in {temp_dir}")
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)

    runtime = time.perf_counter() - tic
    print("\nRuntime: %f" % runtime)
    print("Markers written: %i" % row_counter)
    print("Markers skipped: %i" % skipped_markers)
    if row_counter:
        print("Average runtime per row: %f" % (runtime / row_counter))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
