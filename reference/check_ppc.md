# Posterior Predictive checks

Posterior Predictive checks

## Usage

``` r
check_ppc(x)
```

## Arguments

- x:

  TOSCA object

## Value

dataframe with columns: name, the name of the mutation cluster, pass:
whether the ppc test was passed, i.e. the real number of mutations falls
within 1 sd from the mean of the posterior predictive distribution
