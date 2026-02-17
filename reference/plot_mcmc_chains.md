# Plot the MCMC chains

Plot the MCMC chains

## Usage

``` r
plot_mcmc_chains(
  x,
  days_from = NULL,
  parameters = NULL,
  nuts_params = "divergent__"
)
```

## Arguments

- x:

  TOSCA obj

- days_from:

  date from which to compute the distance from for times

- parameters:

  parameters one desires to retrieve the posterior for

- nuts_params:

  annotations from bayes plot ("lp\_\_", "energy\_\_",
  "stepsize\_\_","n_leapfrog\_\_", "accept_stat\_\_", "treedepth\_\_")

## Value

chains plot
