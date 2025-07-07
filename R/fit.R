fit = function(
    x,
    n_iterations = 10000,
    n_chains = 4,
    warm_up = 5000,
    #max_treedepth = 10,
    cores = 4,
    adapt_delta = 0.99,
    stepsize = 0.01,
    model_name = 'Driver',
    fixed_pars = c()
    ){

  data = get_inference_data(x, model=model_name, fixed_pars=fixed_pars)
  model = get_model(model_name=model_name, fixed_pars=fixed_pars)

  fit <- model$sample(data = data,
                      iter_warmup = warm_up,
                      iter_sampling = n_iterations,
                      chains = n_chains
                      #parallel_chains = 8,
                      #output_dir = tmp_file_path
                      )

  # fit = rstan::stan(
  #   file = file,
  #   data = data,
  #   control = list(
  #     adapt_delta = adapt_delta,
  #     stepsize = stepsize,
  #     max_treedepth = max_treedepth
  #   ),
  #   iter = n_iterations,
  #   chains = n_chains,
  #   warmup = warm_up,
  #   cores = cores
  # )

  fit
}
