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

- stepsize:

  how far the particle moves in each leapfrog step along the trajectory.
  The default is 0.01.

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

- initialisation:

  Boolean, default= F. If true, must specify an init_fun

- init_fun:

  Vector with initial values used if initialisation = T.

- max_mrca:

  maximun date for t_mrca

- reg_dormancy:

  Boolean (0 or 1). If reg_dormancy = 0, a likelihood regularisation is
  imposed on the timing of dormancy end

## Value

TOSCA object

## Examples

``` r
library(TOSCA)
library(dplyr)
library(ggplot2)
data("exampleData_CNA")
set.seed(123)
mutations = exampleData_CNA$Mutations
parameters = exampleData_CNA$Parameters
samples = exampleData_CNA$Samples
therapies = exampleData_CNA$Therapies
x = init(mutations=mutations, samples=samples, therapies=therapies, parameters=parameters)
#> Error in init(mutations = mutations, samples = samples, therapies = therapies,     parameters = parameters): Input must be a data.frame.
fit = TOSCA::fit(x, model_name='CNA',
                 n_iterations = 10000,
                 n_chains = 4,
                 warm_up = 5000)
#> Error: object 'x' not found
```
