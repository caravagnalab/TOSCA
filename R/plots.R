# Plot Prior vs Posterior - all parameters

plot_prior_vs_posterior = function(x, parameter){

  posterior = get_inferred_parameters(x)

  if (!(parameter %in% colnames(posterior)) & parameter != "mrca") return(CNAqc:::eplot())

  prior = get_prior_distribution_type(parameter)
  alpha = get_prior_hyperparameters(x, parameter)[["alpha"]]
  beta = get_prior_hyperparameters(x, parameter)[["beta"]]

  if (prior == 'gamma') draws = rgamma(1000000, alpha, beta)
  if (prior == 'beta') draws = rbeta(1000000, alpha, beta)

  if (parameter == "mrca") { posterior = as.data.frame(posterior[["rho_mrca"]]) } else { posterior = as.data.frame(posterior[[parameter]]) }

  lb = min(min(draws), min(posterior[,1]))
  ub = max(max(draws), max(posterior[,1]))

  if (lb < 1e-10 & ub>2){
    density_post_vs_prior = ggplot2::ggplot() +
      CNAqc:::my_ggplot_theme() +
      geom_histogram(data= posterior, aes(x= posterior[,1], y = ..density..),
                     bins = 100,
                     alpha= .6)+
      geom_vline(
        xintercept = mean(posterior[,1]),
        color = 'indianred3',
        linetype = 'dashed',
        size = 1
      )
  }else{
    density_post_vs_prior = ggplot2::ggplot() +
      geom_density(aes(x = draws), alpha=.6, color= 'white', fill= 'grey') +
      # ggplot2::theme_bw() +
      CNAqc:::my_ggplot_theme() +
      geom_histogram(data= posterior, aes(x= posterior[,1], y = ..density..),
                     bins = 100,
                     alpha= .6)+
      geom_vline(
        xintercept = mean(posterior[,1]),
        color = 'indianred3',
        linetype = 'dashed',
        size = 1
      )
  }

  density_post_vs_prior = density_post_vs_prior +
    ggplot2::scale_x_continuous(name=bquote(.(rlang::sym(parameter))), breaks = scales::pretty_breaks(n=3) )

  L = ggplot2::ggplot_build(density_post_vs_prior)$layout$panel_params[[1]]

  density_post_vs_prior = density_post_vs_prior +
    ggrepel::geom_label_repel(
      ggplot2::aes(
        x = c(L$x.range[2] * .9, NA),
        y = c(L$y.range[2] * .9, NA),
        label = paste0("MAP ", round(mean(posterior[,1]), 3)),
      ),
      #ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      fill = 'indianred3',
      color = 'white',
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    )

  return(density_post_vs_prior)
}

get_prior_distribution_type = function(par){
  prior_distributions = data.frame(
    'parameter'= c("t_eca","t_driver","t_mrca", "mrca", "omega", "mu_driver"),
    'prior_distribution' = c("uniform", "uniform", "uniform", "beta", "gamma", "gamma")
  )

  prior_distributions %>% filter(parameter == par) %>% pull(prior_distribution)

}

get_parameters_of_model = function(model_name){

}

get_mutations_of_model = function(model_name){}

# Posterior Predictive checks
plot_posterior_predictive_checks = function(x, mut1= c("m_clock", "clock", NA), mut2 = c("m_driver","2",NA)){

  estimates = get_inferred_parameters(x)
  mut1_real = get_mutation(x, mut1[1], mut1[2], mut1[3])
  mut2_real = get_mutation(x, mut2[1], mut2[2], mut2[3])
  # to modify for alpha and beta!
  if (mut1[1] == "m_CNA" | mut1[1] == "m_CNA"){
    print("Error, you still need to write this!")
  }
  rep_name1 = paste0(mut1[1], "_rep")
  rep_name2 = paste0(mut2[1], "_rep")

  ppc = estimates %>% dplyr::select(rep_name1,rep_name2)

  ppc = ppc %>%
      dplyr::rename('m_rep_1' = rep_name1,
                    'm_rep_2' = rep_name2) %>%
      mutate(pp = (abs(m_rep_1 - mut1_real) + abs(m_rep_1 - mut2_real) )) %>%
      group_by(m_rep_1, m_rep_2, pp) %>%
      summarise(n = n())

  ggplot(ppc,aes(m_rep_1, m_rep_2, fill= pp, alpha = n)) +
      ggplot2::geom_tile() +
      CNAqc:::my_ggplot_theme() +
      scale_fill_gradient(
        high = "goldenrod",
        low = "purple4"
      )+
      geom_point(
        data =
          data.frame(x = mut1_real, y = mut2_real),
        aes(x = x, y = y),
        inherit.aes = FALSE,
        color = 'indianred3',
        size = 3,
        alpha = 1
      ) +
      geom_vline(
        xintercept = mut1_real,
        color = 'indianred3',
        linetype = 'dashed',
        size = .5
      ) +
      geom_hline(
        yintercept = mut2_real,
        color = 'indianred3',
        linetype = 'dashed',
        size = .5
      ) +
      ylab(bquote(.(rlang::sym('m'))['cl'])) +
      scale_alpha(range = c(.5, 1)) +
      guides(
        fill = 'none',
        alpha = 'none'
      ) + xlab(mut1[1]) + ylab(mut2[1])

  }
plot_expected_N = function(x){
  posterior = get_inferred_parameters(x) %>% as_tibble()
  N = exp(posterior$omega*(posterior$t_mrca-get_sample(x, sample ='2')))
  ggplot() + CNAqc:::my_ggplot_theme() +
    geom_histogram(aes(x=N))
}

# Plot clinical timeline + posterior times
plot_timing = function(x)
{
  clinical_timeline = x$clinical_records
  estimates = get_inferred_parameters(x)

  endpoints = clinical_timeline %>% filter(Clinical.name=="Sample")
  therapies = clinical_timeline %>% filter(Clinical.name!="Sample")

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) %>%
    apply(2, convert_date_real) %>%
    as_tibble()
  times = timing_estimates$variable %>% unique()

  for (i in 1:ncol(timing_estimates)){
    timing_estimates[[i]] = as.Date(timing_estimates[[i]])
  }

  timing_estimates = timing_estimates %>% reshape2::melt() %>% as_tibble()
  var_colors = ggsci::pal_npg()(ncol(timing_estimates))
  names(var_colors) = times

  posterior_plot = ggplot() +
    geom_histogram(
      data = timing_estimates %>% dplyr::rename(Date = value),
      aes(Date, fill = variable, ..density..),
      inherit.aes = FALSE,
      bins = 150
    ) +
    geom_point(
      data = endpoints,
      aes(x = as.Date(convert_date_real(Clinical.value.start)), y = 0),
      inherit.aes = FALSE,
      size = 3
    ) + CNAqc:::my_ggplot_theme()+
    theme(legend.position = 'bottom')+
    scale_fill_manual(values = get_inferred_times_colors())
  # +
  #   geom_label(
  #     data = endpoints,
  #     aes(x = as.Date(convert_date_real(Clinical.value.start)), y = .1, label= paste0(Clinical.name, ' ',Clinical.type))
  #   )
  # +
  #   geom_segment(
  #     data = endpoints,
  #     aes(x = as.Date(convert_date_real(Clinical.value.start)),
  #         xend = as.Date(convert_date_real(Clinical.value.start)),
  #         y=0,yend=.08),
  #     linetype = 'dashed'
  #   ) + CNAqc:::my_ggplot_theme()+
  #   theme(legend.position = 'bottom')+
  #   scale_fill_manual(values = get_inferred_times_colors())

  posterior_plot = posterior_plot + geom_rect(
    aes(xmin = as.Date(convert_date_real(therapies$Clinical.value.start[1])),
        xmax = as.Date(convert_date_real(therapies$Clinical.value.end[1]))),
    ymin = 0,
    ymax = Inf,
    fill = 'indianred',
    colour = "white",
    size = 0.5,
    alpha = .5
  )

  posterior_plot +
    geom_segment(data = therapies[2:nrow(therapies),],
               aes(x=as.Date(convert_date_real(Clinical.value.start)),
                   xend=as.Date(convert_date_real(Clinical.value.start)),
                   y=0, yend=Inf), color = "indianred", alpha=.5)

  # for (th in 2:nrow(therapies)){
  #   posterior_plot = posterior_plot +
  #     geom_vline(
  #     aes(
  #       xintercept = as.Date(convert_date_real(therapies$Clinical.value.start[th])),
  #       #xmax = as.Date(convert_date_real(therapies$Clinical.value.end[th]))),
  #     #ymin = 0,
  #     #ymax = Inf,
  #     #fill = 'indianred',
  #     colour = "indianred",
  #     size = 0.5,
  #     alpha = .5
  #   ))
  # }
  # posterior_plot
}

## Diagnostics
plot_chains = function(x, pars=NA){
  bayesplot::mcmc_trace(get_inferred_parameters(x))
  if (!is.na(pars)){
    bayesplot::mcmc_trace(get_inferred_parameters(x), pars = pars)
  }
}
plot_density=function(x, pars=NA){
  bayesplot::mcmc_dens(get_inferred_parameters(x))
  if (!is.na(pars)){
    bayesplot::mcmc_trace(get_inferred_parameters(x), pars = pars)
  }
}

get_convergence =function(x){
  r_hat<- as.data.frame(x$tosca_fit$summary)$rhat
  if ( sum(r_hat < 1.1) == length(r_hat)){return(TRUE)}else{return(FALSE)}
}
get_diverget_transition = function(patient){}






