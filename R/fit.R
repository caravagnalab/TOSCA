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
    fixed_omega,
    fixed_mu,
    seed
    ){

  data = get_inference_data(x, model=model_name, fixed_omega = fixed_omega, fixed_mu = fixed_mu)
  model = get_model(model_name=model_name, fixed_omega = fixed_omega, fixed_mu = fixed_mu)

  fit <- model$sample(data = data,
                      iter_warmup = warm_up,
                      iter_sampling = n_iterations,
                      chains = n_chains,
                      parallel_chains = n_chains,
                      seed = seed
                      )

  x$tosca_fit = list(
    'posterior' = posterior::as_draws_df(fit$draws()),
    'summary' = fit$summary(),
    'diagnostic_summary' = fit$diagnostic_summary()
  )
  x
}
