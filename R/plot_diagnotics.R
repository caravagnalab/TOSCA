#' Get mapping of input variables names and names used in the model
#'
#' @param x TOSCA obj
#'
#' @return dataframe with papping
get_variables_names = function(x){
  summary = x$Fit$summary

  # variables of interest
  timing_vars <- summary$variable[startsWith(summary$variable, "t_") &!grepl("_tr", summary$variable)]
  ppc_vars = summary$variable[startsWith(summary$variable, "m_")]
  ind = x$Input$Parameters %>% dplyr::filter(Name == "alpha_th_step") %>% nrow()
  rates_vars = c("omega")
  if (ind > 0){
    rates_vars = c(rates_vars, paste0("mu_th_step[",1:ind,"]"))
  }
  if (grepl("driver", x$Input$Parameters$Name) %>% sum() > 1){
    rates_vars = c(rates_vars, "mu_driver")
  }
  vars_of_interest = c(timing_vars, ppc_vars, rates_vars)

  convert_names_par = function(x, n){
    if (grepl("mu_th_step",n)){
      index = strsplit(strsplit(n, "\\[")[[1]][2], "\\]")[[1]][1]
      paste0("mu_",TOSCA:::get_name_mu_step(x, as.integer(index)))
    }else{
      n
    }
  }
  new_names_par = sapply(rates_vars, convert_names_par, x=x)

  convert_names = function(x, n){
    if (grepl("m_th_step",n) | grepl("m_alpha",n) | grepl("m_beta",n)){
      n1 = strsplit(strsplit(n, "\\[")[[1]][1], "_rep")[[1]][1]
      index = strsplit(strsplit(n, "\\[")[[1]][2], "\\]")[[1]][1]
      TOSCA:::get_original_mutation_name(x, n1, as.integer(index))
    }else{
      strsplit(n, "_rep")[[1]][1]
    }
  }
  new_names_ppc = sapply(ppc_vars, convert_names, x=x)

  timing_vars_df = data.frame("original_name"=timing_vars, "name_in_stan_model"=timing_vars, type = "time")
  parameters_df = data.frame("original_name"=new_names_par, "name_in_stan_model"=names(new_names_par), type="parameter")
  mutations_df = data.frame("original_name"=new_names_ppc, "name_in_stan_model"=names(new_names_ppc), type="mutation")
  rbind(timing_vars_df, parameters_df, mutations_df)
}

#' Get posterior draws from cmdstanr obj
#'
#' @param x TOSCA obj
#' @param days_from date from wich to compute the distance from for times
#' @param parameters parameters one desires to retrive the posterior for
#'
#' @return posterior draws
#' @export
get_cmdstanr_posterior = function(x, days_from=NULL, parameters=NULL){
  stan = x$Fit$posteriors$stan_fit
  posterior_cp <- as.array(stan$draws())
  all_pars_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters)) parameters = all_pars_df$original_name
  if (is.null(days_from)) days_from = (x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% pull(Date))[1]

  # Shift timing samples in dates or days from
  timing_pars = all_pars_df %>% dplyr::filter(type == "time")
  for (t in timing_pars$name_in_stan_model){
    #posterior_cp[,,t] = as.Date(TOSCA:::convert_date_real(posterior_cp[,,t], x))
    #if (!is.null(days_from)){
      posterior_cp[,,t] = as.Date(TOSCA:::convert_date_real(posterior_cp[,,t], x)) - as.Date(days_from)
    #}
  }


  # Convert names of parameters
  old_to_new <- all_pars_df$original_name
  names(old_to_new) = all_pars_df$name_in_stan_model
  param_names <- dimnames(posterior_cp)[[3]]
  # Replace only those that appear in the mapping
  param_names[param_names %in% names(old_to_new)] <-
    old_to_new[param_names[param_names %in% names(old_to_new)]]
  dimnames(posterior_cp)[[3]] <- param_names

  posterior_cp

  }

#' Plot the MCMC chains
#'
#' @param x TOSCA obj
#' @param days_from date from which to compute the distance from for times
#' @param parameters parameters one desires to retrieve the posterior for
#' @param nuts_params annotations from bayes plot ("lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__")
#'
#' @return chains plot
#' @export
plot_mcmc_chains = function(x, days_from=NULL, parameters=NULL, nuts_params = "divergent__"){

  posterior_cp = TOSCA:::get_cmdstanr_posterior(x, days_from=days_from, parameters=parameters)
  all_pars_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters)) parameters = all_pars_df$original_name

  stan = x$Fit$posteriors$stan_fit
  lp_cp <- bayesplot::log_posterior(stan)
  np_cp_div <- bayesplot::nuts_params(stan, pars = c(nuts_params)) # "lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__"

  if (!is.null(days_from)){
    date_pars=c()
  }else{
    date_pars= parameters[grepl("t_", parameters)]
  }
  chains = bayesplot::mcmc_trace(posterior_cp,pars = parameters,np = np_cp_div) +
                                         ggplot2::xlab("Post-warmup iteration")+
                                         TOSCA:::my_ggplot_theme()+
                                         # ggplot2::theme_bw()+
                                         ggplot2::theme(
                                           panel.grid.minor = ggplot2::element_blank(),
                                           strip.text = ggplot2::element_text(color = "black"),
                                           panel.border = ggplot2::element_blank(),
                                           strip.background = ggplot2::element_blank(),  # removes facet background
                                           panel.background = ggplot2::element_rect(fill = "white", color = NA)  # white panel, no border
                                         )

  # chains = list()
  #
  # for (p in 1:length(parameters)){
  #   chains[[p]] = bayesplot::mcmc_trace(posterior_cp,
  #                                       pars = parameters[p],
  #                                       np = np_cp_div) +
  #     ggplot2::xlab("Post-warmup iteration")+
  #     ggplot2::theme_bw()+
  #     ggplot2::theme(
  #       panel.grid.minor = ggplot2::element_blank(),
  #       panel.border = ggplot2::element_blank(),
  #       strip.background = ggplot2::element_blank(),  # removes facet background
  #       panel.background = ggplot2::element_rect(fill = "white", color = NA)  # white panel, no border
  #     )
  #   if (parameters[p] %in% date_pars){
  #     chains[[p]] = chains[[p]] #+ ggplot2::scale_y_continuous()
  #   }
  # }
  chains

}

#' Plot the sampling of parameters highlighting the iterations where divergent transitions occurred
#'
#' @param x TOSCA object
#' @param days_from date from which to compute the distance from for times
#' @param parameters parameters one desires to retrieve the posterior for
#' @param nuts_params annotations from bayes plot ("lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__")
#'
#' @return plot
#' @export
plot_divergent_transitions = function(x, days_from=NULL, parameters=NULL, nuts_params = "divergent__"){
  posterior_cp = TOSCA:::get_cmdstanr_posterior(x, days_from=days_from, parameters=parameters)
  all_pars_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters)) parameters = all_pars_df$original_name

  stan = x$Fit$posteriors$stan_fit
  lp_cp <- bayesplot::log_posterior(stan)
  np_cp_div <- bayesplot::nuts_params(stan, pars = c(nuts_params)) # "lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__"

  bayesplot::mcmc_parcoord(posterior_cp, np = np_cp_div, pars = parameters) +
    #ggplot2::xlab("Post-warmup iteration")+
    #ggplot2::theme_bw()+
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),  # removes facet background
      panel.background = ggplot2::element_rect(fill = "white", color = NA))
}

#' Plot univariate and multivariate distributions of selected parameters
#'
#' @param x TOSCA object
#' @param days_from date from which to compute the distance from for times
#' @param parameters parameters one desires to retrieve the posterior for
#' @param nuts_params annotations from bayes plot ("lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__")
#'
#' @return plot
#' @export
plot_pairs = function(x, days_from=NULL, parameters=NULL, nuts_params = "divergent__"){
  posterior_cp = TOSCA:::get_cmdstanr_posterior(x, days_from=days_from, parameters=parameters)
  all_pars_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters)) parameters = all_pars_df$original_name

  stan = x$Fit$posteriors$stan_fit
  lp_cp <- bayesplot::log_posterior(stan)
  np_cp_div <- bayesplot::nuts_params(stan, pars = c(nuts_params)) # "lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__"

  ggplot2::set_theme(TOSCA:::my_ggplot_theme()+ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                                                        panel.border = ggplot2::element_blank()))
  pairs = bayesplot::mcmc_pairs(posterior_cp, np = np_cp_div,
                     pars = parameters,
                     off_diag_args = list(size = 0.3, alpha = .5),
  )
  pairs
}

#' Plot Rhats of selected parameters
#'
#' @param x TOSCA object
#' @param parameters
#'
#' @return plot
#' @export
rhats_plot = function(x, parameters=NULL){
  stan = x$Fit$posteriors$stan_fit
  rhats <- bayesplot::rhat(stan)
  names_df = TOSCA:::get_variables_names(x)
  map = names_df$original_name
  names(map) = names_df$name_in_stan_model
  names(rhats)[names(rhats) %in% names(map)] <- map[names(rhats)[names(rhats) %in% names(map)]]
  if (is.null(parameters))  parameters = names_df$original_name
  rhats = rhats[names(rhats) %in% parameters]

  bayesplot::mcmc_rhat(rhats) + #ggplot2::xlim(1, 1.0005) +
    bayesplot::yaxis_text(hjust = 1)+
    #ggplot2::theme_bw()+
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
      )+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2))
}

#' Plot the Effective Sample Size of selected parameters
#'
#' @param x TOSCA object
#' @param parameters
#'
#' @return plot
#' @export
ess_plot = function(x, parameters=NULL){
  stan = x$Fit$posteriors$stan_fit
  ratios <- bayesplot::neff_ratio(stan)
  names_df = TOSCA:::get_variables_names(x)
  map = names_df$original_name
  names(map) = names_df$name_in_stan_model
  names(ratios)[names(ratios) %in% names(map)] <- map[names(ratios)[names(ratios) %in% names(map)]]
  if (is.null(parameters))  parameters = names_df$original_name
  ratios = ratios[names(ratios) %in% parameters]

  bayesplot::mcmc_neff(ratios, size = 2)+ bayesplot::yaxis_text(hjust = 1)+
    #ggplot2::theme_bw()+
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )+ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2))
}

#' Plot the autocorrelation of MCMC samples with progressive iterations
#'
#' @param x TOSCA object
#' @param parameters
#' @param days_from date from which to compute the distance from for times
#'
#' @return plot
#' @export
autocorrelation_plot = function(x, parameters=NULL, days_from=NULL){

  names_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters))  parameters = names_df$original_name

  posterior_cp = TOSCA:::get_cmdstanr_posterior(x, days_from=days_from, parameters=parameters)
  all_pars_df = TOSCA:::get_variables_names(x)
  if (is.null(parameters)) parameters = all_pars_df$original_name

  # stan = x$Fit$posteriors$stan_fit
  # lp_cp <- bayesplot::log_posterior(stan)
  # np_cp_div <- bayesplot::nuts_params(stan, pars = c(nuts_params)) # "lp__", "energy__", "stepsize__","n_leapfrog__", "accept_stat__", "treedepth__"


  bayesplot::mcmc_acf(posterior_cp,pars = parameters, lags = 10)+
  #ggplot2::theme_bw()+
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(color = "black"),
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()  # removes facet background
    #panel.background = element_rect(fill = "white", color = NA)  # white panel, no border
  )
}

#' Plot to quantiy the heaviness of the tails of the posterior
#'
#' @param x TOSCA object
#'
#' @return plot
#' @export
energy_plot = function(x){
  stan = x$Fit$posteriors$stan_fit
  np_cp_2 = bayesplot::nuts_params(stan, pars = c("energy__"))
  bayesplot::mcmc_nuts_energy(np_cp_2)+
    #ggplot2::theme_bw()+
    TOSCA:::my_ggplot_theme()+
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()  # removes facet background
      #panel.background = element_rect(fill = "white", color = NA)  # white panel, no border
    )+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 2))
}


