# Introduction

This repository contains the code that convert the [``AlphaPeel``](https://github.com/AlphaGenes/AlphaPeel) genotype probability output file to the input format of the [``IMPUTATIONQUALITY``](https://gitlab.gwdg.de/kolja.thormann1/imputationquality) [1].

The ``ImputeAccure.py`` is copied from the [``IMPUTATIONQUALITY`` repository](https://gitlab.gwdg.de/kolja.thormann1/imputationquality).

# Example

In the ``example\``:
- ``.geno_prob.txt`` is the example genotype probability output from ``AlphaPeel``
- ``converted_geno_prob.accuracy`` is the final output of ``run.sh``
- ``converted_geno_prob.txt`` is the ``IMPUTATIONQUALITY`` input file format version of ``.geno_prob.txt``
- ``params.txt`` is the parameter file for ``IMPUTATIONQUALITY``
- ``run.sh`` is the bash code to run the example

To run the example:

```bash
cd example
bash run.sh
```

# Citation

[1] Thormann, K.A., Tozzi, V., Starke, P. et al. ImputAccur: fast and user-friendly calculation of genotype-imputation accuracy-measures. BMC Bioinformatics 23, 316 (2022). https://doi.org/10.1186/s12859-022-04863-z


