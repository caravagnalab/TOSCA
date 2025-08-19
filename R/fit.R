fit = function(
    x,
    n_iterations = 10000,
    n_chains = 4,
    warm_up = 5000,
    max_treedepth = 10,
    cores = 4,
    adapt_delta = 0.99,
    stepsize = 0.01,
    model_name = 'Driver',
    dormancy=F,
    verbose = F
    ){

  data = get_inference_data(x, model=model_name, dormancy = dormancy)
  model = get_model(model_name=model_name, dormancy = dormancy)

  cat("\n--- Start Sampling ---\n")

  fit <- model$sample(data = data,
                        iter_warmup = warm_up,
                        iter_sampling = n_iterations,
                        chains = n_chains,
                        show_messages = verbose
                        )

  x$tosca_fit = list(
    'posterior' = posterior::as_draws_df(fit$draws()),
    'summary' = fit$summary(),
    'diagnostic_summary' = fit$diagnostic_summary(),
    'model_info' = list('model_name'=model_name, 'dormancy'=dormancy)
  )

  cat("\n--- End Sampling ---\n")



  # ---- Print Diagnostics ----
  cat("\n--- Sampling Diagnostics ---\n")

  # Convergence: Check all Rhat < 1.01
  max_rhat <- max(x$tosca_fit$summary$rhat, na.rm = TRUE)
  converged <- all(x$tosca_fit$summary$rhat <= 1.01, na.rm = TRUE)
  cat(sprintf("Convergence (Rhat < 1.01): %s (max Rhat = %.3f)\n",
              if (converged) "✅ Yes" else "❌ No", max_rhat))

  # Divergent transitions
  divergences <- x$tosca_fit$diagnostic_summary$num_divergent
  total_div <- sum(divergences)
  total_samples <- n_iterations*n_chains
  div_pct <- 100 * total_div / total_samples
  cat(sprintf("Divergent transitions: %d / %d (%.2f%%)\n", total_div, total_samples, div_pct))

  # EBFMI check
  ebfmi <- x$tosca_fit$diagnostic_summary$ebfmi
  low_ebfmi <- sum(ebfmi < 0.3)
  cat(sprintf("EBFMI < 0.3 in %d / %d chains\n", low_ebfmi, length(ebfmi)))
  #cat("--- End Diagnostics ---\n")

  cat("\n--- Posterior Predictive Checks ---\n")
  print.data.frame(check_ppc(x))

  return(x)

  # Additional warnings :
  #   - posterior predictive checks
}

# fit = function(
#     x,
#     n_iterations = 10000,
#     n_chains = 4,
#     warm_up = 5000,
#     cores = 4,
#     adapt_delta = 0.99,
#     stepsize = 0.01,
#     model_name = 'Driver',
#     fixed_omega = FALSE,
#     fixed_mu = FALSE,
#     verbose = FALSE
# ) {
#   # Prepare data and model
#   data = get_inference_data(x, model = model_name, fixed_omega = fixed_omega, fixed_mu = fixed_mu)
#   model = get_model(model_name = model_name, fixed_omega = fixed_omega, fixed_mu = fixed_mu)
#
#   if (!verbose) {
#     # Setup message animation
#     keep_spinning <- TRUE
#     spinner <- parallel::mcparallel({
#       i <- 0
#       while (keep_spinning) {
#         dots <- paste0(rep(".", (i %% 4) + 1), collapse = "")
#         cat(sprintf("\rPerforming inference%s   ", dots))
#         flush.console()
#         Sys.sleep(0.5)
#         i <- i + 1
#       }
#     })
#
#     # Silence both stdout and stderr from cmdstanr
#     nullfile <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
#     stdout_con <- file(nullfile, open = "wt")
#     stderr_con <- file(nullfile, open = "wt")
#     sink(stdout_con, type = "output")
#     sink(stderr_con, type = "message")
#
#     # Ensure sinks are closed and spinner stopped
#     on.exit({
#       sink(type = "message"); sink(type = "output")
#       close(stdout_con); close(stderr_con)
#       keep_spinning <<- FALSE
#       parallel::mccollect(spinner)
#       cat("\rSampling complete.             \n")
#     }, add = TRUE)
#   } else {
#     message("Sampling in progress...")
#   }
#
#   # Run the sampling
#   fit <- model$sample(
#     data = data,
#     iter_warmup = warm_up,
#     iter_sampling = n_iterations,
#     chains = n_chains,
#     parallel_chains = cores,
#     adapt_delta = adapt_delta,
#     step_size = stepsize,
#     show_messages = verbose,
#     refresh = if (verbose) 100 else 0
#   )
#
#   # Extract posterior and diagnostics
#   posterior_draws <- posterior::as_draws_df(fit$draws())
#   summary_stats <- fit$summary()
#   diag_summary <- fit$diagnostic_summary()
#
#   x$tosca_fit <- list(
#     posterior = posterior_draws,
#     summary = summary_stats,
#     diagnostic_summary = diag_summary
#   )
#
#   # Convergence diagnostics
#   divergences <- sum(diag_summary$num_divergent__)
#   treedepth_hits <- sum(diag_summary$num_max_treedepth__)
#   low_ebfmi_chains <- sum(diag_summary$ebfmi < 0.2)
#
#   converged <- divergences == 0 && treedepth_hits == 0 && low_ebfmi_chains == 0
#
#   msg <- paste0(
#     "\nModel convergence status: ", if (converged) "✅ Converged" else "⚠️ Potential Issues",
#     "\n  - Divergent transitions: ", divergences,
#     "\n  - Max treedepth hit: ", treedepth_hits,
#     "\n  - Chains with low E-BFMI (< 0.2): ", low_ebfmi_chains
#   )
#   message(msg)
#
#   return(x)
# }







