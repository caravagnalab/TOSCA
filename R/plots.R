# Plot Prior vs Posterior - all parameters

#' Return prior distribution type for each parameter
#'
#' @param par Name of the parameter of interest.
#'
#' @return Name of the prior distribution used in the model
#'
#' @examples
get_prior_distribution_type = function(par){
  prior_distributions = data.frame(
    'parameter'= c("t_eca","t_driver","t_mrca", "mrca", "omega", "mu_driver","s"),
    'prior_distribution' = c("uniform", "uniform", "uniform", "beta", "gamma", "gamma","beta")
  )

  prior_distributions %>% dplyr::filter(parameter == par) %>% dplyr::pull(prior_distribution)

}

#' Empty plot
#'
#' @return Creates and empty plot
#'
#' @examples
eplot = function (){
  ggplot2::ggplot(data = data.frame(x = 0, y = 0, label = "X"),
                  ggplot2::aes(x = x, y = y, label = label)) + CNAqc:::my_ggplot_theme() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black",
                                                        fill = NA, linetype = "dashed"), panel.background = ggplot2::element_rect(fill = "gainsboro"),
                   axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

#' Customized ggplot theme
#'
#' @param cex
#'
#' @return Customized ggplot theme
#'
#' @examples
my_ggplot_theme = function (cex = 1)
{
  cex_opt = getOption("CNAqc_cex", default = 1)
  ggplot2::theme_light(base_size = 10 * cex_opt) + ggplot2::theme(legend.position = "bottom",
                                                                  legend.key.size = ggplot2::unit(0.3 * cex_opt, "cm"),
                                                                  panel.background = ggplot2::element_rect(fill = "white"))
}

#' Prior vs Posterior distribution of inferred parameter
#'
#' @param x TOSCA object
#' @param parameter parameter name
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_prior_vs_posterior_single_parameter = function(x, parameter){

  posterior = get_inferred_parameters(x)

  if (!(parameter %in% colnames(posterior)) & parameter != "mrca") return(eplot())

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

  if (prior == 'gamma') draws = stats::rgamma(1000000, alpha, beta)
  if (prior == 'beta') draws = stats::rbeta(1000000, alpha, beta)

  if (parameter == "mrca") { posterior = as.data.frame(posterior[["rho_mrca"]]) } else { posterior = as.data.frame(posterior[[parameter]]) }

  lb = min(min(draws), min(posterior[,1]))
  ub = max(max(draws), max(posterior[,1]))

  if (lb < 1e-10 & ub>2){
    density_post_vs_prior = ggplot2::ggplot() +
      my_ggplot_theme() +
      ggplot2::geom_histogram(data= posterior, ggplot2::aes(x= posterior[,1], y = ..density..),
                     bins = 100,
                     alpha= .6)+
      ggplot2::geom_vline(
        xintercept = mean(posterior[,1]),
        color = 'indianred3',
        linetype = 'dashed',
        size = 1
      )
  }else{
    density_post_vs_prior = ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = draws), alpha=.6, color= 'white', fill= 'grey') +
      # ggplot2::theme_bw() +
      my_ggplot_theme() +
      ggplot2::geom_histogram(data= posterior, ggplot2::aes(x= posterior[,1], y = ..density..),
                     bins = 100,
                     alpha= .6)+
      ggplot2::geom_vline(
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

#' Produce collective plot of prior vs posterior for all inferred parameters
#'
#' @param x TOSCA obj
#'
#' @return plot
#' @export
#'
#' @examples
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
plot_prior_vs_posterior = function(x){
  parameters = c()
  posterior = get_inferred_parameters(x)

  is_omega = grepl("omega", colnames(posterior)) %>% sum()
  is_s = grepl("s", colnames(posterior)) %>% sum()
  is_mu_driver = grepl("mu_driver", colnames(posterior)) %>% sum()
  if (is_omega >= 1) parameters = c("omega", parameters)
  if (is_s == 1) parameters = c("s", parameters)
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
    plot_prior_vs_posterior_single_parameter(x, p) #+ ggtitle(p)
  })

  ggpubr::ggarrange(plotlist = parameter_plots, nrow=1, ncol = length(parameters))
}

# Posterior Predictive checks
#' Posterior Predictive checks
#'
#' @param x TOSCA object
#'
#' @return dataframe with columns: name, the name of the mutation cluster, pass: whether the ppc test was passed, i.e. the real number of mutations falls within 1 sd from the mean of the posterior predictive distribution
#'
#' @examples
check_ppc = function(x){

  posterior = get_inferred_parameters(x)#x$tosca_fit$posterior
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = get_inference_data(x,
                                     model = x$tosca_fit$model_info$model_name,
                                     #fixed_omega = x$tosca_fit$model_info$fixed_omega,
                                     #fixed_mu = x$tosca_fit$model_info$fixed_mu,
                                     dormancy = x$tosca_fit$model_info$dormancy
                    )
  variables = names(training_data)[grepl("m_", names(training_data))]

  ppc_df = data.frame()

  for (v in variables){
    #print(v)

    true_value = training_data[[v]]

    if (length(true_value)>0){

      if (v %in% c("m_th_step","m_th_cauchy", "m_alpha", "m_beta") &
          x$tosca_fit$model_info$dormancy==F ){

        for (i in 1:length(true_value)){
          rep_draws = posterior[[paste0(v,"_rep[",i,"]")]]
          coverage <- mean(true_value[i] > (rep_draws - 2*sd(rep_draws)) &
                             true_value[i] < (rep_draws + 2*sd(rep_draws)))
          # cat(paste0(v, '_', i,':'))
          # cat(sprintf("Approximate posterior predictive coverage: %.1f%%\n", coverage * 100))
          # cat('\n')
          original_name = get_original_mutation_name(x, name = v, index = i)
          ppc_df = rbind(ppc_df,
                         data.frame(
                           "name" = original_name,
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

      original_name = get_original_mutation_name(x, name = v, index = NA)

      ppc_df = rbind(ppc_df,
                     data.frame(
                       "name" = original_name,
                       "variable" = v,
                       "true_value"=true_value,
                       "ppc_coverage"=coverage,
                       "pass" = coverage > .6
                     ))

      }
    }
  }

  return(ppc_df %>% dplyr::select(name, pass))

}

#' Plot the posterior predictive distribution and compares it to the real number of mutations
#'
#' @param x TOSCA object
#' @param mut1_real real value of the first mutation type of interest
#' @param mut2_real real value of the second mutation type of interest
#' @param rep_name1 name of the posterior predictive variable corresponding to the first mutation type of interest, as reported in the model
#' @param rep_name2 name of the posterior predictive variable corresponding to the second mutation type of interest, as reported in the model
#'
#' @return plot
#' @export
#'
#' @examples
plot_ppc_single_mut = function(x, mut1_real, mut2_real, rep_name1, rep_name2){

  estimates = get_inferred_parameters(x)
  # mut1_real = get_mutation(x, mut1[1], mut1[2], mut1[3])
  # mut2_real = get_mutation(x, mut2[1], mut2[2], mut2[3])
  #
  # rep_name1 = paste0(mut1[1], "_rep")
  # rep_name2 = paste0(mut2[1], "_rep")

  ppc = estimates %>% dplyr::select(rep_name1,rep_name2)

  n1 = strsplit(rep_name1, "_rep")[[1]][1]
  n2 = strsplit(rep_name2, "_rep")[[1]][1]
  if (grepl("\\[", rep_name1)) index1 = strsplit(strsplit(rep_name1, "\\[")[[1]][2], "\\]")[[1]][1] else index1 = NA
  if (grepl("\\[", rep_name2)) index2 = strsplit(strsplit(rep_name2, "\\[")[[1]][2], "\\]")[[1]][1] else index2 = NA
  original_name_1 = get_original_mutation_name(x, n1, index1)
  original_name_2 = get_original_mutation_name(x, n2, index2)

  ppc = ppc %>%
      dplyr::rename('m_rep_1' = rep_name1,
                    'm_rep_2' = rep_name2) %>%
    dplyr::mutate(pp = (abs(m_rep_1 - mut1_real) + abs(m_rep_1 - mut2_real) )) %>%
    dplyr::group_by(m_rep_1, m_rep_2, pp) %>%
    dplyr::summarise(n = n())

  ggplot2::ggplot(ppc,ggplot2::aes(m_rep_1, m_rep_2, color= pp, alpha = n), size=.05) +
      ggplot2::geom_point() +
      CNAqc:::my_ggplot_theme() +
    ggplot2::scale_color_gradient(
        high = "goldenrod",
        low = "purple4"
      )+
    ggplot2::geom_point(
        data =
          data.frame(x = mut1_real, y = mut2_real),
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = 'indianred3',
        size = 3,
        alpha = 1
      ) +
    ggplot2::geom_vline(
        xintercept = mut1_real,
        color = 'indianred3',
        linetype = 'dashed',
        size = .5
      ) +
    ggplot2::geom_hline(
        yintercept = mut2_real,
        color = 'indianred3',
        linetype = 'dashed',
        size = .5
      ) +
    ggplot2::ylab(bquote(.(rlang::sym('m'))['cl'])) +
    ggplot2::scale_alpha(range = c(.5, 1)) +
    ggplot2::guides(
        fill = 'none',
        alpha = 'none'
      ) + ggplot2::xlab(original_name_1) + ggplot2::ylab(original_name_2)+
    ggplot2::theme(legend.position = "none")

  }

#' Plot posterior predictive checks
#'
#' @param x TOSCA object
#'
#' @return plot of the posterior predictive distribution versus the real value of all mutation groups provided in input
#' @export
#'
#' @examples
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
#' plot_ppc(fit)
plot_ppc = function(x){

  posterior = get_inferred_parameters(x)
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = get_inference_data(x, model = x$tosca_fit$model_info$model_name,
                                     #fixed_omega = x$tosca_fit$model_info$fixed_omega,
                                     #fixed_mu = x$tosca_fit$model_info$fixed_mu,
                                     dormancy = x$tosca_fit$model_info$dormancy)
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

#' Posterior distribution of the number of cells collected in the first and second samples.
#'
#' @param x TOSCA object
#'
#' @return Posterior distribution of the number of cells collected in the first and second samples, considering the growth rate (omega) and the time of birth of the MRCA (t_mrca) of each sample.
#' @export
#'
#' @examples
plot_expected_N = function(x){

  posterior = get_inferred_parameters(x) %>% dplyr::as_tibble()
  N_rel = exp(posterior$omega*(get_sample(x, sample ='2')-posterior$t_mrca))
  N_pre = exp(posterior$omega*(get_sample(x, sample ='1')-posterior$t_mrca_primary))

  primary_name = x$clinical_records %>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "1") %>% dplyr::pull(Clinical.original.name)
  relapse_name = x$clinical_records %>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "2") %>% dplyr::pull(Clinical.original.name)

  N_rel_plot = ggplot2::ggplot() + my_ggplot_theme() +
    ggplot2::geom_histogram(ggplot2::aes(x=log(N_rel)))+
    ggplot2::scale_x_log10() +
    ggplot2::xlab(paste0("Tumor size at ", relapse_name, "\n(log-scale)"))
  N_pre_plot = ggplot2::ggplot() + my_ggplot_theme() +
    ggplot2::geom_histogram(ggplot2::aes(x=log(N_pre)))+
    ggplot2::xlab(paste0("Tumor size at ", primary_name, "\n(log-scale)"))

  ggpubr::ggarrange(plotlist = list(N_rel_plot, N_pre_plot), nrow = 2)
}

# Plot clinical timeline + posterior times
#' Posterior distribution of the inferred times
#'
#' @param x
#'
#' @return Posterior distributio plot of the inferred times, mapped on the clinical history
#' @export
#'
#' @examples
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
#' plot_timing(fit)
plot_timing = function(x)
{
  clinical_timeline = x$clinical_records
  estimates = get_inferred_parameters(x)

  endpoints = clinical_timeline %>% dplyr::filter(Clinical.name=="Sample")
  therapies = clinical_timeline %>% dplyr::filter(Clinical.name!="Sample")

  # 1. time posterior plot
  timing_estimates = estimates %>%
    dplyr::select(starts_with('t_')) %>%
    apply(2, convert_date_real) %>%
    dplyr::as_tibble()
  #times = timing_estimates$variable %>% unique()

  for (i in 1:ncol(timing_estimates)){
    timing_estimates[[i]] = as.Date(timing_estimates[[i]])
  }

  timing_estimates = timing_estimates %>% reshape2::melt() %>% dplyr::as_tibble()
  times = timing_estimates$variable %>% unique()
  therapy_names = therapies$Clinical.original.name %>% unique()
  var_colors = ggsci::pal_npg()(length(times)+length(therapy_names))
  times_colors = var_colors[1:length(times)]
  names(times_colors) = times
  clinical_colors = var_colors[length(times)+1:length(therapy_names)]
  names(clinical_colors) = therapy_names
  times_colors_df = data.frame("variable" = names(times_colors), "color"=times_colors)
  timing_estimates = dplyr::left_join(timing_estimates, times_colors_df, by = "variable")

  posterior_plot = ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = timing_estimates %>% dplyr::rename(Date = value),
      ggplot2::aes(Date, fill = color, ..density..),
      inherit.aes = FALSE,
      bins = 150
    ) +
    ggplot2::geom_point(
      data = endpoints,
      ggplot2::aes(x = as.Date(convert_date_real(Clinical.value.start)), y = 0),
      inherit.aes = FALSE,
      size = 3
    ) +
    my_ggplot_theme()+
    ggplot2::theme(legend.position = 'bottom')+
    ggplot2::scale_fill_identity()
    #scale_fill_manual(values = times_colors)

  hist_data <- ggplot2::ggplot_build(posterior_plot)$data[[1]]
  ymin <- 0
  ymax <- max(hist_data$y)

  # 5% above bottom
  ylab_pos <- ymin + 0.1 * (ymax - ymin)


  if ("t_dormancy_start" %in% timing_estimates$variable){
    MAP_dormancy_start = timing_estimates %>% dplyr::filter(variable == "t_dormancy_start") %>% dplyr::pull(value) %>% mean()
    MAP_dormancy_end = timing_estimates %>% dplyr::filter(variable == "t_dormancy_end") %>% dplyr::pull(value) %>% mean()
    timing_estimates = timing_estimates %>% dplyr::filter(!(variable %in% c("t_dormancy_start","t_dormancy_end")))

    posterior_plot = ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = timing_estimates %>% dplyr::rename(Date = value),
        ggplot2::aes(Date, fill = variable, ..density..),
        inherit.aes = FALSE,
        bins = 150
      ) +
      ggplot2::geom_point(
        data = endpoints,
        ggplot2::aes(x = as.Date(convert_date_real(Clinical.value.start)), y = 0),
        inherit.aes = FALSE,
        size = 3
      ) + my_ggplot_theme()+
      ggplot2::theme(legend.position = 'bottom')+
      ggplot2::scale_fill_manual(values = var_colors)+
      ggplot2::geom_rect(
        ggplot2::aes(xmin = MAP_dormancy_start,
            xmax = MAP_dormancy_end),
        ymin = 0,
        ymax = Inf,
        fill = 'grey',
        colour = "white",
        size = 0.5,
        alpha = .5
      )
  }

  therapies = therapies %>% dplyr::mutate(Duration = Clinical.value.end-Clinical.value.start) %>% dplyr::mutate(short=ifelse(Duration < 30/365, T, F))
  new_col = data.frame(Clinical.original.name = names(clinical_colors), colors = clinical_colors)
  therapies = dplyr::left_join(therapies,new_col, by="Clinical.original.name")

  posterior_plot = posterior_plot + ggplot2::geom_rect(
    data = therapies %>% dplyr::filter(short == F),
    ggplot2::aes(xmin = as.Date(convert_date_real(Clinical.value.start)),
        xmax = as.Date(convert_date_real(Clinical.value.end)),
        fill=colors),
    ymin = 0,
    ymax = Inf,
    #fill = c,
    colour = "white",
    size = 0.5,
    alpha = .5
  ) +
    ggplot2::geom_segment(
      data = therapies %>% dplyr::filter(short == T),
      ggplot2::aes(x=as.Date(convert_date_real(Clinical.value.start)),
          xend=as.Date(convert_date_real(Clinical.value.start)),
          y=0, yend=Inf), color = therapies %>% dplyr::filter(short == T) %>% dplyr::pull(colors), alpha=1)  +
    ggplot2::guides(color = "none")

  posterior_plot = posterior_plot +
    ggplot2::geom_segment(
      data = endpoints,
      ggplot2::aes(x = as.Date(convert_date_real(Clinical.value.start)),
          xend = as.Date(convert_date_real(Clinical.value.start)),
          y=0,
          yend = ylab_pos),
      size = .5, linetype="dashed"
    )+
    ggplot2::geom_label(
      data = endpoints,
      ggplot2::aes(x = as.Date(convert_date_real(Clinical.value.start)), y = ylab_pos, label=Clinical.original.name),
      size = 3
    )

  dummy_guide <- function(
    labels = NULL,
    ...,
    title = NULL,
    key   = draw_key_point,
    guide_args = list(),
    min_value_time
  ) {
    aesthetics <- list(...)
    n <- max(lengths(aesthetics), 0)
    labels <- labels %||% seq_len(n)

    aesthetics$alpha <- aesthetics$alpha %||% rep(1, n)

    guide_args$override.aes <- guide_args$override.aes %||% aesthetics
    guide <- do.call(guide_legend, guide_args)

    ggplot2::update_geom_defaults("point", list(dummy = "x"))

    dummy_geom <- ggplot2::geom_point(
      data = data.frame(
        x = rep(min_value_time, n),
        y = rep(Inf, n),
        dummy = factor(labels, levels = labels)   # preserve order!
      ),
      ggplot2::aes(x, y, dummy = dummy),
      alpha = 0,
      key_glyph = key
    )

    fills <- aesthetics$fill

    dummy_scale <- ggplot2::discrete_scale(
      "dummy", "dummy_scale",
      palette = function(x) {
        stats::setNames(fills, labels)[x]
      },
      name = title,
      guide = guide
    )

    list(dummy_geom, dummy_scale)
  }


  posterior_plot +
    dummy_guide(
    labels = c(names(times_colors), names(clinical_colors)),
    fill   = c(times_colors, clinical_colors),
    colour = NA,
    title  = "Inferred times and treatments",
    key = draw_key_polygon,
    min_value_time = min(timing_estimates$value, na.rm = TRUE)
  )

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






