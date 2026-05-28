# Introduction

This repository contains the code that convert the [``AlphaPeel``](https://github.com/AlphaGenes/AlphaPeel) genotype probability output file to the input format of the [``IMPUTATIONQUALITY``](https://gitlab.gwdg.de/kolja.thormann1/imputationquality) [1].

The ``ImputeAccure.py`` is modified from the [``IMPUTATIONQUALITY`` repository](https://gitlab.gwdg.de/kolja.thormann1/imputationquality).

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

# Assess the accuracy test results with visualization

After AlphaPeel accuracy test run, we would obtain a sets of outputs. 

You can use the existing code in ``run_all\`` to evaluate the outputs.

```bash
cd run_all
```

Install the requirements

```bash
pip install -r requirements.txt # written for Python 3.13, but may work for other versions of Python
```

Create a file called ``.env`` and put the absolute paths of your accuracy outputs and simulation file in it, such as 

```bash
DATA_DIR=~/AlphaPeel/tests/accuracy_tests/outputs
TRUE_DATA_DIR=~/AlphaPeel/tests/accuracy_tests/sim_for_alphapeel_accu_test
```

Run the bash file

```bash
bash run_list.txt
```

Then it would convert the genotype probabilities in AlphaPeel format from the AlphaPeel accuracy test directory to the input file format of ``ImputeAccure.py`` and run the program to obtain the metric evalutions of them. 

To plot the accuracies, you can use ``plot_separate.py`` to generate the plots for each file.

```bash
python plot_separate.py
```

Or, you can also compare multiple files together with the interactive plot by running 

```bash
python plot_interactive.py
```

Follow the output instruction to access the interactive plot.

If you want to host the interactive plot in a container, you can also use the existing docker file to generate a docker image.

# Citation

[1] Thormann, K.A., Tozzi, V., Starke, P. et al. ImputAccur: fast and user-friendly calculation of genotype-imputation accuracy-measures. BMC Bioinformatics 23, 316 (2022). https://doi.org/10.1186/s12859-022-04863-z


