fit = function(
    x,
    n_iterations = 50000,
    n_chains = 4,
    warm_up = 10000,
    max_treedepth = 10,
    cores = 4,
    adapt_delta = 0.99,
    stepsize = 0.01,
    model = 'cna'
    ){

  data = get_inference_data(x)
  if (model == 'cna'){
    file = '/stan/CNA.stan'
  }else{
    file = '~/Docs/GitHub/TOSCA/stan/Driver.stan'
  }

  fit = rstan::stan(
    file = file,
    data = data,
    control = list(
      adapt_delta = adapt_delta,
      stepsize = stepsize,
      max_treedepth = max_treedepth
    ),
    iter = n_iterations,
    chains = n_chains,
    warmup = warm_up,
    cores = cores
  )
  fit
}
