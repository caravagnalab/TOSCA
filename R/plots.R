# Plot Prior vs Posterior - all parameters
get_prior_distribution_type = function(par){
  prior_distributions = data.frame(
    'parameter'= c("t_eca","t_driver","t_mrca", "mrca", "omega", "mu_driver"),
    'prior_distribution' = c("uniform", "uniform", "uniform", "beta", "gamma", "gamma")
  )

  prior_distributions %>% filter(parameter == par) %>% pull(prior_distribution)

}

# need to adjust this function for parameters with multiple entries!
plot_prior_vs_posterior_single_parameter = function(x, parameter){

  posterior = get_inferred_parameters(x)

  if (!(parameter %in% colnames(posterior)) & parameter != "mrca") return(CNAqc:::eplot())

  if (grepl("\\[", parameter)){
    prior = "gamma"
    par_name = strsplit(parameter, split = "\\[")[[1]][1]
    index = as.integer(strsplit(strsplit(parameter, split = "\\[")[[1]][2], split = "\\]")[1])
    alpha = get_prior_hyperparameters(x, par_name)[["alpha"]][index]
    beta = get_prior_hyperparameters(x, par_name)[["beta"]][index]
    } else {
      prior = get_prior_distribution_type(parameter)
      alpha = get_prior_hyperparameters(x, parameter)[["alpha"]]
      beta = get_prior_hyperparameters(x, parameter)[["beta"]]

      }

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

plot_prior_vs_posterior = function(x, model){
  parameters = c()
  posterior = get_inferred_parameters(x)

  is_omega = grepl("omega", colnames(posterior)) %>% sum()
  is_mu_driver = grepl("mu_driver", colnames(posterior)) %>% sum()
  if (is_omega == 1) parameters = c("omega", parameters)
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
    plot_prior_vs_posterior_single_parameter(x, p) + ggtitle(p)
  })

  ggpubr::ggarrange(plotlist = parameter_plots, nrow=1, ncol = length(parameters))
}

# Posterior Predictive checks
check_ppc = function(x){

  posterior = get_inferred_parameters(x)#x$tosca_fit$posterior
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = get_inference_data(x, model = x$tosca_fit$model_info$model_name,
                     fixed_omega = x$tosca_fit$model_info$fixed_omega,
                    fixed_mu = x$tosca_fit$model_info$fixed_mu)
  variables = names(training_data)[grepl("m_", names(training_data))]

  ppc_df = data.frame()

  for (v in variables){
    # print(v)

    true_value = training_data[[v]]

    if (length(true_value)>0){

      if (v %in% c("m_th_step","m_th_cauchy", "m_alpha", "m_beta",
                   "m_alpha_tetraploid_step","m_beta_tetraploid_step",
                   "m_alpha_tetraploid_cauchy","m_beta_tetraploid_cauchy")){

        for (i in 1:length(true_value)){
          rep_draws = posterior[[paste0(v,"_rep[",i,"]")]]
          coverage <- mean(true_value[i] > (rep_draws - 2*sd(rep_draws)) &
                             true_value[i] < (rep_draws + 2*sd(rep_draws)))
          # cat(paste0(v, '_', i,':'))
          # cat(sprintf("Approximate posterior predictive coverage: %.1f%%\n", coverage * 100))
          # cat('\n')
          ppc_df = rbind(ppc_df,
                         data.frame(
                           "variable" = paste0(v,"_",i),
                           "true_value"=true_value[i],
                           "ppc_coverage"=coverage,
                           "pass" = coverage > .6
                         ))
        }

          }else{

      rep_draws = posterior[[paste0(v, "_rep")]]
      coverage <- mean(true_value > (rep_draws - 2*sd(rep_draws)) &
                       true_value < (rep_draws + 2*sd(rep_draws)))
      # cat(paste0(v, ':'))
      # cat(sprintf("Approximate posterior predictive coverage: %.1f%%\n", coverage * 100))
      # cat('\n')
      ppc_df = rbind(ppc_df,
                     data.frame(
                       "variable" = v,
                       "true_value"=true_value,
                       "ppc_coverage"=coverage,
                       "pass" = coverage > .6
                     ))

      }
    }
  }

  return(ppc_df %>% select(variable, pass))

}

plot_ppc_single_mut = function(x, mut1_real, mut2_real, rep_name1, rep_name2){

  estimates = get_inferred_parameters(x)
  # mut1_real = get_mutation(x, mut1[1], mut1[2], mut1[3])
  # mut2_real = get_mutation(x, mut2[1], mut2[2], mut2[3])
  #
  # rep_name1 = paste0(mut1[1], "_rep")
  # rep_name2 = paste0(mut2[1], "_rep")

  ppc = estimates %>% dplyr::select(rep_name1,rep_name2)

  ppc = ppc %>%
      dplyr::rename('m_rep_1' = rep_name1,
                    'm_rep_2' = rep_name2) %>%
      mutate(pp = (abs(m_rep_1 - mut1_real) + abs(m_rep_1 - mut2_real) )) %>%
      group_by(m_rep_1, m_rep_2, pp) %>%
      summarise(n = n())

  ggplot(ppc,aes(m_rep_1, m_rep_2, color= pp, alpha = n), size=.05) +
      ggplot2::geom_point() +
      CNAqc:::my_ggplot_theme() +
      scale_color_gradient(
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
      ) + xlab(rep_name1) + ylab(rep_name2)+
    theme(legend.position = "none")

  }

plot_ppc = function(x){

  posterior = get_inferred_parameters(x)
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = get_inference_data(x, model = x$tosca_fit$model_info$model_name,
                                     fixed_omega = x$tosca_fit$model_info$fixed_omega,
                                     fixed_mu = x$tosca_fit$model_info$fixed_mu)
  variables = names(training_data)[grepl("m_", names(training_data))]

  ppc_plot_list = list()
  seen = c()

  for (v1 in variables){
    #print(v1)
    for (v2 in variables){
      #print(v2)

      if (v1 != v2){

        true_value1 = training_data[[v1]]
        true_value2 = training_data[[v2]]

        if (!(v1 %in% seen)){

          if (length(true_value1)==1 & length(true_value2)==1){

            rep_name1 = paste0(v1, "_rep")
            rep_name2 = paste0(v2, "_rep")
            ppc_plot = plot_ppc_single_mut(x, true_value1, true_value2, rep_name1, rep_name2)
            ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot


          }else{

            if (length(true_value1)>1 & length(true_value2)==1){
              for (i in 1:length(true_value1)){
                rep_name1 = paste0(v1,"_rep[",i,"]")
                rep_name2 = paste0(v2, "_rep")
                ppc_plot = plot_ppc_single_mut(x, true_value1[i], true_value2, rep_name1, rep_name2)
                ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot
              }
            }

            if (length(true_value1)==1 & length(true_value2)>1){
              for (i in 1:length(true_value1)){
                rep_name2 = paste0(v2,"_rep[",i,"]")
                rep_name1 = paste0(v1, "_rep")
                ppc_plot = plot_ppc_single_mut(x, true_value1, true_value2[i], rep_name1, rep_name2)
                ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot
              }
            }

          }

          seen = c(seen, v1, v2)

        }

        }


    }

  }

  ggpubr::ggarrange(plotlist = ppc_plot_list) + ggtitle("Posterior Predictive Checks")


  }

plot_expected_N = function(x){

  posterior = get_inferred_parameters(x) %>% as_tibble()
  N_rel = exp(posterior$omega*(get_sample(x, sample ='2')-posterior$t_mrca))
  N_pre = exp(posterior$omega*(get_sample(x, sample ='1')-posterior$t_mrca_primary))

  N_rel_plot = ggplot() + CNAqc:::my_ggplot_theme() +
    geom_histogram(aes(x=log(N_rel)))+
    scale_x_log10()
  N_pre_plot = ggplot() + CNAqc:::my_ggplot_theme() +
    geom_histogram(aes(x=log(N_pre)))+
    scale_x_log10()

  ggpubr::ggarrange(plotlist = list(N_rel_plot, N_pre_plot), nrow = 2)
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
  #times = timing_estimates$variable %>% unique()

  for (i in 1:ncol(timing_estimates)){
    timing_estimates[[i]] = as.Date(timing_estimates[[i]])
  }

  timing_estimates = timing_estimates %>% reshape2::melt() %>% as_tibble()
  times = timing_estimates$variable %>% unique()
  var_colors = ggsci::pal_npg()(length(times))
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
    scale_fill_manual(values = var_colors)
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






