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
  adapt_delta = 0.95,
  stepsize = 0.01,
  seed = 5,
  model_name = "Driver",
  dormancy = F,
  verbose = F,
  initialisation = F,
  init_fun = NULL,
  max_mrca = NA,
  reg_dormancy = 0
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

- adapt_delta:

  target acceptance probability for the Metropolis step in HMC. The
  default is 0.95.

- seed:

  seed of the computation. Default is 5.

- model_name:

  Name of the model: "CNA" or "Driver"

- dormancy:

  Boolean specifying if the model should include dormancy (available
  only for "CNA" model). The default is F

- verbose:

  Boolen specifying if the fit function should output the iterations
  progression. The default is F.

- step_size:

  how far the particle moves in each leapfrog step along the trajectory.
  The default is 0.01.

## Value

TOSCA object
