# DASC

Detecting hidden batch factors through data adaptive adjustment for biological effects

## Overview

DASC is an R package used for identifying batches and classifying samples into different batches in a high dimensional gene expression dataset. The batch information can be further used as a covariate in conjunction with other variables of interest among standard bioinformatics analysis like differential expression analysis. 

## Installation

**_Dependencies_**

Before installing the package, please make sure the prerequisite packages described in `DESCRIPTION` file have been installed successfully. For example, [cvxcluster](https://github.com/echi/cvxclustr) needs to be installed from Eric Chi's github repository as it is not available through CRAN database.

To install the DASC package, start the terminal and enter:

On linux/MacOs
```bash
git clone https://github.com/HaidYi/DASC.git
R CMD build DASC
R CMD INSTALL DASC_*.tar.gz
```

## Usage

For concrete usages, please refer the *.html* vignettes in `inst/doc` file.

## Citation info

If you use DASC for your analysis, please cite our paper as here below.

```
@article{Yi2018Detecting,
    title={Detecting hidden batch factors through data-adaptive adjustment for biological effects},
    author={Haidong Yi* and
            Ayush T. Raman* and
            Han Zhang and
            Genevera I. Allen and
            Zhandong Liu},
    journal={Bioinformatics},
    volume={34},
    number={7},
    pages={1141--1147},
    year={2018},
}
```
