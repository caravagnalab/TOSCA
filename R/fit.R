#' Infer the timing of the events of interest fitting the appropriate TOSCA model
#'
#' @param x TOSCA object
#' @param n_iterations Number of iteration of each MCMC chain. The default is 10000.
#' @param n_chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param warm_up A positive integer specifying the number of warmup (aka burnin) iterations per chain. The default is 5000.
#' @param cores Number of cores to use when executing the chains in parallel
#' @param model_name Name of the model: "CNA" or "Driver"
#' @param dormancy Boolean specifying if the model should include dormancy (available only for "CNA" model). The default is F
#' @param verbose Boolen specifying if the fit function should output the iterations progression. The default is F.
#'
#' @return TOSCA object
#' @export
#'
#' @examples

fit = function(
    x,
    n_iterations = 10000,
    n_chains = 4,
    warm_up = 5000,
    #max_treedepth = 10,
    cores = 4,
    #adapt_delta = 0.99,
    #stepsize = 0.01,
    model_name = 'Driver',
    dormancy=F,
    verbose = F,
    initialisation = F,
    init_fun = NULL
    ){

  data = TOSCA:::get_inference_data(x, model=model_name, dormancy = dormancy)
  model = TOSCA:::get_model(model_name=model_name, dormancy = dormancy)

  cat("\n--- Start Sampling ---\n")

  requireNamespace("cmdstanr")
  if (initialisation){
    fit <- model$sample(data = data,
                        iter_warmup = warm_up,
                        iter_sampling = n_iterations,
                        chains = n_chains,
                        show_messages = verbose,
                        init = init_fun
    )
  }else{
  fit <- model$sample(data = data,
                        iter_warmup = warm_up,
                        iter_sampling = n_iterations,
                        chains = n_chains,
                        show_messages = verbose
                        )
  }

  x$Fit = list(
    'posterior' = posterior::as_draws_df(fit$draws()),
    'summary' = fit$summary(),
    'diagnostic_summary' = fit$diagnostic_summary(),
    'model_info' = list('model_name'=model_name, 'dormancy'=dormancy)
  )

  cat("\n--- End Sampling ---\n")



  # ---- Print Diagnostics ----
  cat("\n--- Sampling Diagnostics ---\n")

  # Convergence: Check all Rhat < 1.01
  max_rhat <- max(x$Fit$summary$rhat, na.rm = TRUE)
  converged <- all(x$Fit$summary$rhat <= 1.01, na.rm = TRUE)
  cat(sprintf("Convergence (Rhat < 1.01): %s (max Rhat = %.3f)\n",
              if (converged) "\u2705 Yes" else "\u274C No", max_rhat))

  # Divergent transitions
  divergences <- x$Fit$diagnostic_summary$num_divergent
  total_div <- sum(divergences)
  total_samples <- n_iterations*n_chains
  div_pct <- 100 * total_div / total_samples
  cat(sprintf("Divergent transitions: %d / %d (%.2f%%)\n", total_div, total_samples, div_pct))

  # EBFMI check
  ebfmi <- x$Fit$diagnostic_summary$ebfmi
  low_ebfmi <- sum(ebfmi < 0.3)
  cat(sprintf("EBFMI < 0.3 in %d / %d chains\n", low_ebfmi, length(ebfmi)))
  #cat("--- End Diagnostics ---\n")

  cat("\n--- Posterior Predictive Checks ---\n")
  print.data.frame(TOSCA:::check_ppc(x))

  return(x)

}

