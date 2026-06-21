
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TOSCA <a href="https://caravagnalab.github.io/TOSCA/"><img src="man/figures/logo.png" align="right" height="180" /></a>

<!-- badges: start -->
<!-- badges: end -->

TOSCA (Timing Of Somatic Events in Cancer), is a framework that
leverages longitudinal whole genome sequencing data to time clonal
expantions associated to drug resistance. It uses a Bayesian Poisson
process to model mutation accumulation along the phylogenetic branch
connecting the Early Common Ancestor (ECA) and the Most Recent Common
Ancestor (MRCA), where genetic resistance emerges. The framework infers
the timing of resistance-associated genetic events that alter the
intensity of the mutation accumulation process, such as copy-number
alterations and point mutations associated with a hypermutator
phenotype. Importantly, these inferred timings are mapped onto the
patient’s clinical history, enabling events to be interpreted in real
time rather than solely in pseudotime.

#### Citation

If you use `TOSCA`, please cite: \* *Timing the origin of genetic
therapy resistance in human cancers* Riccardo Bergamin, Alice Antonello,
Giulio Caravagna.

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/TOSCA/-yellow.svg)](https://caravagnalab.github.io/TOSCA)

## Installation

You can install the released version of from
[GitHub](https://github.com) with:

``` r
devtools::install_github("caravagnalab/TOSCA")
```

**Note:** TOSCA requires [cmdstanr](https://mc-stan.org/cmdstanr/) to
fit models. Please follow the [installation
instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
before use.
