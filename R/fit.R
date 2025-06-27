fit = function(
    x,
    n_iterations = 10000,
    n_chains = 4,
    warm_up = 5000,
    max_treedepth = 10,
    cores = 4,
    adapt_delta = 0.99,
    stepsize = 0.01,
    model = 'standard'
    ){

  data = get_inference_data(x, model=model)
  code = generate_stan_code(x, model=model)

  fit = rstan::stan(
    model_code = code,
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
