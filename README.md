# Figures for approximation of musculoskeletal kinematics

The repository has the data and scripts for generating figures for "Functional description of musculoskeletal transformation using multidimensional polynomials" (2019) by Anton Sobinov, Matthew Boots, Valeriya Gritsenko, Rob Gaunt, Sergiy Yakovenko.

## Getting Started

This project's code is available on [GitHub](https://github.com/neurowired/musculoskeletal_approximation).

### Prerequisites

* [MATLAB 2017a+](https://www.mathworks.com/products/matlab.html) with built-in and configured for Microsoft Visual Studio compilation (`mex`) tool.

* [Python 2.7](https://www.python.org/download/releases/2.7/) In the description I will use `py -2` command to call Python.

* [NumPy and SciPy](https://www.scipy.org/scipylib/download.html) modules installed.

* [Scikit-learn](https://scikit-learn.org/stable/install.html)

* [GEPA](https://github.com/nishbo/gepa) - module with functions for approximation of differential systems.

## Structure of the repository

* "approximations_progress" has pickled Python structures for polynomial approximations.

* [arm_metaDOF.csv](arm_metaDOF.csv) and [arm_metaMuscle.csv](arm_metaMuscle.csv) are comma-delimited files with information about approxiamted muscles and DOFs that they cross. [extractMetaData.m](extractMetaData.m), [extractMetaMuscle.m](extractMetaMuscle.m) and [meta.py](meta.py) import these tables into MATLAB and Python.

* [external_measure_report.txt](external_measure_report.txt) has measurements of memory and generation time required for different methods.

* [accuracy_test_report.txt](accuracy_test_report.txt) has a report produced by accuracyTest.m on all muscles and a report by [params_and_aiccs.py](params_and_aiccs.py) with substituted values of the quality.

## Generating figures from the manuscript

### Figure 2

Run `msdFig2.m` in MATLAB environment.

### Figure 3

```
musculoskeletal_approximation>py -2 msd_fig3.py
```

### Figure 4

Run `msdFig4.m` in MATLAB environment.

### Figure 5

```
musculoskeletal_approximation>py -2 msd_fig5.py
```

### Figure 6

This also generates figures for Supplementary Figure 1.

```
musculoskeletal_approximation>py -2 msd_fig6.py
```

### Figure 7

```
musculoskeletal_approximation>py -2 msd_fig7.py
```

### Supplementary Figure 2

```
musculoskeletal_approximation>py -2 msd_supp_fig2.py
```

## Authors

* **Anton Sobinov** - *Code* - [nishbo](https://github.org/nishbo)
