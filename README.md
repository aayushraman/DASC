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

## References:

1. Yi H*, Raman AT*, Zhang H, Allen GI, Liu Z. Detecting hidden batch factors through data-adaptive adjustment for biological effects. [Bioinformatics (2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6454417/) (PMID: 29617963; PMCID: PMC6454417)
2. Wang X, Yi H, Wang J, Liu Z, Yin Y, Zhang H. GDASC: a GPU parallel-based web server for detecting hidden batch factors. [Bioinformatics (2020)](https://academic.oup.com/bioinformatics/article/36/14/4211/5835275) (PMID: 32386292)
