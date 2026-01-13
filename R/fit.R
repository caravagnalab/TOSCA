#' Infer the timing of the events of interest fitting the appropriate TOSCA model
#'
#' @param x TOSCA object
#' @param n_iterations Number of iteration of each MCMC chain. The default is 10000.
#' @param n_chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param warm_up A positive integer specifying the number of warmup (aka burnin) iterations per chain. The default is 5000.
#' @param cores Number of cores to use when executing the chains in parallel
#' @param adapt_delta target acceptance probability for the Metropolis step in HMC. The default is 0.95.
#' @param stepsize how far the particle moves in each leapfrog step along the trajectory. The default is 0.01.
#' @param seed seed of the computation. Default is 5.
#' @param model_name Name of the model: "CNA" or "Driver"
#' @param dormancy Boolean specifying if the model should include dormancy (available only for "CNA" model). The default is F
#' @param verbose Boolen specifying if the fit function should output the iterations progression. The default is F.
#' @param initialisation Boolean, default= F. If true, must specify an init_fun
#' @param init_fun Vector with initial values used if initialisation = T.
#' @param max_mrca maximun date for t_mrca
#' @param reg_dormancy Boolean (0 or 1). If reg_dormancy = 0, a likelihood regularisation is imposed on the timing of dormancy end
#'
#' @return TOSCA object
#' @export
#'
#' @examples
#' library(TOSCA)
#' library(dplyr)
#' library(ggplot2)
#' data("exampleData_CNA")
#' set.seed(123)

#' mutations = exampleData_CNA$Mutations
#' parameters = exampleData_CNA$Parameters
#' samples = exampleData_CNA$Samples
#' therapies = exampleData_CNA$Therapies

#' x = init(mutations=mutations, samples=samples, therapies=therapies, parameters=parameters)
#' fit = TOSCA::fit(x, model_name='CNA',
#'                  n_iterations = 10000,
#'                  n_chains = 4,
#'                  warm_up = 5000)
fit = function(
    x,
    n_iterations = 10000,
    n_chains = 4,
    warm_up = 5000,
    #max_treedepth = 10,
    cores = 4,
    adapt_delta = 0.95,
    stepsize = 0.01,
    seed = 5,
    model_name = 'Driver',
    dormancy=F,
    verbose = F,
    initialisation = F,
    init_fun = NULL,
    max_mrca = NA,
    reg_dormancy = 0
    ){

  if (!is.na(max_mrca)) max_mrca = TOSCA:::convert_real_date(x, date = max_mrca)
  data = TOSCA:::get_inference_data(x, model=model_name, dormancy = dormancy, max_mrca=max_mrca, reg_dormancy=reg_dormancy)
  model = TOSCA:::get_model(model_name=model_name, dormancy = dormancy)

  cat("\n--- Start Sampling ---\n")

  requireNamespace("cmdstanr")
  if (initialisation){
    fit <- model$sample(data = data,
                        iter_warmup = warm_up,
                        iter_sampling = n_iterations,
                        chains = n_chains,
                        adapt_delta = adapt_delta,
                        step_size =  stepsize ,
                        show_messages = verbose,
                        seed = seed,
                        init = init_fun
    )
  }else{
  fit <- model$sample(data = data,
                        iter_warmup = warm_up,
                        iter_sampling = n_iterations,
                        adapt_delta = adapt_delta,
                        step_size =  stepsize,
                        chains = n_chains,
                        seed = seed,
                        show_messages = verbose
                        )
  }


  # posterior -> posterior as it came out of the stan model
  posterior = posterior::as_draws_df(fit$draws())

  # variables of interest
  timing_vars = posterior %>% dplyr::select(colnames(posterior)[!(grepl("_tr", colnames(posterior)))]) %>% dplyr::select(starts_with('t_')) %>% colnames()
  ppc_vars = posterior %>% dplyr::select(starts_with('m_')) %>% colnames()
  ind = x$Input$Parameters %>% dplyr::filter(Name == "alpha_th_step") %>% nrow()
  rates_vars = c("omega", paste0("mu_th_step[",1:ind,"]"))

  # timing_posteriors -> timings without the translated times, converted to dates
  timing_estimates = posterior %>% dplyr::select(colnames(posterior)[!(grepl("_tr", colnames(posterior)))]) %>%
    dplyr::select(timing_vars)
  timing_estimates = timing_estimates %>%
    apply(2, TOSCA:::convert_date_real, x=x) %>%
    dplyr::as_tibble()
  # non_timing = posterior %>%
  #   dplyr::select(!timing_vars)
  # posterior = cbind(timing_estimates, non_timing)
  # In summary dataframe

  # posterior_predictive_checks -> ppc with converted variable names
  ppc_posterior = posterior %>% select(ppc_vars)
  new_names = c("m_clock_primary_rep"="m_clock_primary",
                "m_clock_rep"="m_clock")
  convert_names = function(x, n){
    if (grepl("m_th_step",n) | grepl("m_alpha",n) | grepl("m_beta",n)){
      n1 = strsplit(strsplit(n, "\\[")[[1]][1], "_rep")[[1]][1]
      index = strsplit(strsplit(n, "\\[")[[1]][2], "\\]")[[1]][1]
      TOSCA:::get_original_mutation_name(x, n1, as.integer(index))
    }else{
      strsplit(n, "_rep")[[1]][1]
    }
  }
  new_names = sapply(colnames(ppc_posterior), convert_names, x=x)
  names(ppc_posterior) <- new_names[names(ppc_posterior)]

  # rates_posteriors -> rates posteriors with converted names
  rates_posterior = posterior %>% select(rates_vars)
  convert_names_par = function(x, n){
    if (grepl("mu_th_step",n)){
      index = strsplit(strsplit(n, "\\[")[[1]][2], "\\]")[[1]][1]
      paste0("mu_",TOSCA:::get_name_mu_step(x, as.integer(index)))
    }else{
      n
    }
  }
  new_names_par = sapply(colnames(rates_posterior), convert_names_par, x=x)
  names(rates_posterior) <- new_names_par[names(rates_posterior)]

  summary = fit$summary()
  # converted_summary = summary %>% dplyr::filter(variable %in% timing_vars) %>% dplyr::select(median, mean, q5, q95) %>%
  #   apply(2, TOSCA:::convert_date_real, x=x)
  # rest_of_cols = summary %>% dplyr::filter(variable %in% timing_vars) %>% dplyr::select(!c(median, mean, q5, q95))
  # converted_summary = cbind(converted_summary, rest_of_cols)
  # summary_2 = summary %>% dplyr::filter(!(variable %in% timing_vars) )
  # summary = rbind(converted_summary, summary_2)
  # summary <- summary[, c("variable", "mean", "median", "sd", "mad", "q5", "q95", "rhat",  "ess_bulk", "ess_tail")]
  #
  # posterior = dplyr::as_tibble(posterior)
  # summary = dplyr::as_tibble(summary) %>% dplyr::filter(variable %in% colnames(posterior)[!(grepl("_tr", colnames(posterior)))])

  x$Fit = list(
    'posteriors' = list(
      "stan_posterior" = posterior,
      "timing_posteriors" = timing_estimates,
      "posterior_predictive_checks" = ppc_posterior,
      "rates_posteriors" = rates_posterior
      ), # posterior::as_draws_df(fit$draws()),

    "summary" = summary,
    # 'summary' = list(
    #   "summary" = summary, # dplyr::as_tibble(summary),
    #   "timing_summary" = summary %>% dplyr::filter(variable %in% timing_vars),
    #   "rates_summary" = summary %>% dplyr::filter(variable %in% rates_vars)
    #   ), # fit$summary(),

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

  cat("\n--- Inference summary ---\n")
  print.data.frame(TOSCA::get_fit_summary(x))

  return(x)

}

