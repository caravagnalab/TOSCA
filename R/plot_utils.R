# Plot Prior vs Posterior - all parameters

get_prior_distribution_type = function(par){
  prior_distributions = data.frame(
    'parameter'= c("t_eca","t_driver","t_mrca", "mrca", "omega", "mu_driver","s"),
    'prior_distribution' = c("uniform", "uniform", "uniform", "beta", "gamma", "gamma","beta")
  )

  prior_distributions %>% dplyr::filter(parameter == par) %>% dplyr::pull(prior_distribution)

}

eplot = function (){
  ggplot2::ggplot(data = data.frame(x = 0, y = 0, label = "X"),
                  ggplot2::aes(x = x, y = y, label = label)) + my_ggplot_theme() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black",
                                                        fill = NA, linetype = "dashed"), panel.background = ggplot2::element_rect(fill = "gainsboro"),
                   axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

my_ggplot_theme = function (cex = 1)
{
  cex_opt = getOption("CNAqc_cex", default = 1)
  ggplot2::theme_light(base_size = 10 * cex_opt) + ggplot2::theme(legend.position = "bottom",
                                                                  legend.key.size = ggplot2::unit(0.3 * cex_opt, "cm"),
                                                                  panel.background = ggplot2::element_rect(fill = "white"))
}



## Diagnostics
# plot_chains = function(x, pars=NA){
#   bayesplot::mcmc_trace(get_inferred_parameters(x))
#   if (!is.na(pars)){
#     bayesplot::mcmc_trace(get_inferred_parameters(x), pars = pars)
#   }
# }
# plot_density=function(x, pars=NA){
#   bayesplot::mcmc_dens(get_inferred_parameters(x))
#   if (!is.na(pars)){
#     bayesplot::mcmc_trace(get_inferred_parameters(x), pars = pars)
#   }
# }
#
# get_convergence =function(x){
#   r_hat<- as.data.frame(x$tosca_fit$summary)$rhat
#   if ( sum(r_hat < 1.1) == length(r_hat)){return(TRUE)}else{return(FALSE)}
# }
#get_diverget_transition = function(patient){}






