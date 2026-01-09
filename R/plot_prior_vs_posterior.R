#' Title
#'
#' @param x TOSCA object
#' @param index index of mu_step parameter
#'
#' @return name of the drug corresponding to the index
#'
get_name_mu_step = function(x, index){
  mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="alpha_th_step")
  mu_alpha = left_join(
    x$Input$Therapies %>% filter(Class=="Mutagenic"),
    mu_alpha %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
    dplyr::arrange(as.Date(Start))
  names = mu_alpha %>%
    group_by(Name) %>%
    slice_min(as.Date(Start, format="%Y-%m-%d"), with_ties = FALSE) %>%
    dplyr::arrange(as.Date(Start)) %>%
    pull(Name)
  names[index]
}

#' Prior vs Posterior distribution of inferred parameter
#'
#' @param x TOSCA object
#' @param parameter parameter name
#'
#' @return ggplot
#' @export
plot_prior_vs_posterior_single_parameter = function(x, parameter){

  posterior = TOSCA:::get_inferred_parameters(x) %>% dplyr::select(parameter)
  colnames(posterior) = "par"

  # if (!(parameter %in% colnames(posterior)) & parameter != "mrca") return(eplot())

  if (parameter=="omega"){
    alpha = TOSCA:::get_parameter(x, "omega_alpha")
    beta = TOSCA:::get_parameter(x, "omega_beta")
  }
  name_mu_step = ""
  if (grepl("mu_th_step", parameter)){
    index = as.integer(strsplit(strsplit(parameter, "\\[")[[1]][2], "\\]")[[1]][1])
    alpha = TOSCA:::get_mutation_rate(x, type = "th_step")[["alpha"]][index]
    beta = TOSCA:::get_mutation_rate(x, type = "th_step")[["beta"]][index]
    name_mu_step = TOSCA:::get_name_mu_step(x, index)
  }
  if (grepl("mu_driver", parameter)){
    alpha = TOSCA:::get_mutation_rate(x, type = "driver")[["alpha"]]
    beta = TOSCA:::get_mutation_rate(x, type = "driver")[["beta"]]
  }

  draws = stats::rgamma(1000000, alpha, beta)

  #if (prior == 'gamma') draws = stats::rgamma(1000000, alpha, beta)
  #if (prior == 'beta') draws = stats::rbeta(1000000, alpha, beta)

  #if (parameter == "mrca") { posterior = as.data.frame(posterior[["rho_mrca"]]) } else { posterior = as.data.frame(posterior[[parameter]]) }

  lb = min(min(draws), min(posterior$par))
  ub = max(max(draws), max(posterior$par))

  if (lb < 1e-10 & ub>2){
    density_post_vs_prior = ggplot2::ggplot() +
      TOSCA:::my_ggplot_theme() +
      ggplot2::geom_histogram(data= posterior, ggplot2::aes(x= par, y = ..density..),
                              bins = 100,
                              alpha= .6)+
      ggplot2::geom_vline(
        xintercept = mean(posterior$par),
        color = 'indianred3',
        linetype = 'dashed',
        size = 1
      )
  }else{
    density_post_vs_prior = ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = draws), alpha=.6, color= 'white', fill= 'grey') +
      # ggplot2::theme_bw() +
      TOSCA:::my_ggplot_theme() +
      ggplot2::geom_histogram(data= posterior, ggplot2::aes(x= par, y = ..density..),
                              bins = 100,
                              alpha= .6)+
      ggplot2::geom_vline(
        xintercept = mean(posterior$par),
        color = 'indianred3',
        linetype = 'dashed',
        size = 1
      )
  }

  if (parameter == "omega"){
    density_post_vs_prior = density_post_vs_prior +
      ggplot2::scale_x_continuous(name=bquote(.(rlang::sym(parameter))), breaks = scales::pretty_breaks(n=3) )
  }
  if (parameter == "mu_driver") {
    density_post_vs_prior = density_post_vs_prior +
      ggplot2::scale_x_continuous(name=bquote(.(rlang::sym(mu[.("driver")]))), breaks = scales::pretty_breaks(n=3) )
  }
  if (name_mu_step!=""){
    density_post_vs_prior = density_post_vs_prior +
      ggplot2::scale_x_continuous(name=bquote(mu[.(name_mu_step)]), breaks = scales::pretty_breaks(n=3) )
  }

  L = ggplot2::ggplot_build(density_post_vs_prior)$layout$panel_params[[1]]

  # density_post_vs_prior = density_post_vs_prior +
  #   ggrepel::geom_label_repel(
  #     ggplot2::aes(
  #       x = c(L$x.range[2] * .9, NA),
  #       y = c(L$y.range[2] * .9, NA),
  #       label = paste0("MAP ", round(mean(posterior$par), 3)),
  #     ),
  #     #ylim = c(L$y.range[2] * .9, NA),
  #     size = 2,
  #     fill = 'indianred3',
  #     color = 'white',
  #     nudge_y = 0,
  #     nudge_x = 0,
  #     show.legend = FALSE
  #   )

  return(density_post_vs_prior)
}

# get_prior_hyperparameters = function(x, name){
#   # name = mu_driver*, mu_th_step*, scale_th_cauchy*, omega*, mrca*, s*
#   if (name %in% x$Input$Parameters$Name){
#     x$Input$Parameters %>% dplyr::filter(Name==name) %>% dplyr::pull(Value) %>% as.double()
#   }else{
#     alpha= x$Input$Parameters %>% dplyr::filter(Name==paste0(name,'_alpha')) %>% dplyr::pull(Value) %>% as.double()
#     beta= x$Input$Parameters %>% dplyr::filter(Name==paste0(name,'_beta')) %>% dplyr::pull(Value) %>% as.double()
#     list('alpha'=alpha, 'beta'=beta)
#   }
# }

#' Produce collective plot of prior vs posterior for all inferred parameters
#'
#' @param x TOSCA obj
#'
#' @return plot
#' @export
#'
#' @examples
plot_prior_vs_posterior = function(x){
  parameters = c()
  posterior = TOSCA:::get_inferred_parameters(x)



  is_omega = grepl("omega", colnames(posterior)) %>% sum()
  # is_s = grepl("s", colnames(posterior)) %>% sum()
  is_mu_driver = grepl("mu_driver", colnames(posterior)) %>% sum()
  if (is_omega >= 1) parameters = c("omega", parameters)
  #if (is_s == 1) parameters = c("s", parameters)
  if (is_mu_driver == 1) parameters = c("mu_driver", parameters)

  n_mu_th = grepl("mu_th_step", colnames(posterior)) %>% sum()
  n_scales = grepl("scales_th_cauchy", colnames(posterior)) %>% sum()

  if (n_mu_th > 0){
    for (n in 1:n_mu_th){
      parameters = c(parameters, paste0("mu_th_step[",n,"]"))
    }
  }

  if (n_scales > 0){
    for (n in 1:n_scales){
      parameters = c(parameters, paste0("scales_th_cauchy[",n,"]"))
    }
  }

  parameter_plots = lapply(parameters, function(p){
    # print(p)
    TOSCA:::plot_prior_vs_posterior_single_parameter(x, p) #+ ggtitle(p)
    })

  ggpubr::ggarrange(plotlist = parameter_plots, nrow=1, ncol = length(parameters))
}


#' library(TOSCA)
#' library(dplyr)
#' library(ggplot2)
#' data("exampleData_CNA")
#' mutations = exampleData_CNA$Mutations
#' parameters = exampleData_CNA$Parameters
#' samples = exampleData_CNA$Samples
#' therapies = exampleData_CNA$Therapies
#'
#' x = init(mutations=mutations, samples=samples, therapies=therapies, parameters=parameters)
#' fit = TOSCA::fit(x, model_name='CNA', n_iterations = 1000, n_chains = 4, warm_up = 500)
#' plot_prior_vs_posterior(fit)
