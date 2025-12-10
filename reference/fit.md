# Infer the timing of the events of interest fitting the appropriate TOSCA model

Infer the timing of the events of interest fitting the appropriate TOSCA
model

## Usage

``` r
fit(
  x,
  n_iterations = 10000,
  n_chains = 4,
  warm_up = 5000,
  cores = 4,
  model_name = "Driver",
  dormancy = F,
  verbose = F,
  initialisation = F,
  init_fun = NULL
)
```

## Arguments

- x:

  TOSCA object

- n_iterations:

  Number of iteration of each MCMC chain. The default is 10000.

- n_chains:

  A positive integer specifying the number of Markov chains. The default
  is 4.

- warm_up:

  A positive integer specifying the number of warmup (aka burnin)
  iterations per chain. The default is 5000.

- cores:

  Number of cores to use when executing the chains in parallel

- model_name:

  Name of the model: "CNA" or "Driver"

- dormancy:

  Boolean specifying if the model should include dormancy (available
  only for "CNA" model). The default is F

- verbose:

  Boolen specifying if the fit function should output the iterations
  progression. The default is F.

## Value

TOSCA object
