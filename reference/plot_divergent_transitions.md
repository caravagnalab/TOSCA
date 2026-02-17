# Plot the sampling of parameters highlighting the iterations where divergent transitions occurred

Plot the sampling of parameters highlighting the iterations where
divergent transitions occurred

## Usage

``` r
plot_divergent_transitions(
  x,
  days_from = NULL,
  parameters = NULL,
  nuts_params = "divergent__"
)
```

## Arguments

- x:

  TOSCA object

- days_from:

  date from which to compute the distance from for times

- parameters:

  parameters one desires to retrieve the posterior for

- nuts_params:

  annotations from bayes plot ("lp\_\_", "energy\_\_",
  "stepsize\_\_","n_leapfrog\_\_", "accept_stat\_\_", "treedepth\_\_")

## Value

plot
