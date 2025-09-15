# Posterior Predictive checks
#' Posterior Predictive checks
#'
#' @param x TOSCA object
#'
#' @return dataframe with columns: name, the name of the mutation cluster, pass: whether the ppc test was passed, i.e. the real number of mutations falls within 1 sd from the mean of the posterior predictive distribution
check_ppc = function(x){

  posterior = TOSCA:::get_inferred_parameters(x)#x$tosca_fit$posterior
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = TOSCA:::get_inference_data(x,
                                     model = x$Fit$model_info$model_name,
                                     #fixed_omega = x$tosca_fit$model_info$fixed_omega,
                                     #fixed_mu = x$tosca_fit$model_info$fixed_mu,
                                     dormancy = x$Fit$model_info$dormancy
  )
  variables = names(training_data)[grepl("m_", names(training_data))]

  ppc_df = data.frame()

  for (v in variables){
    # print(v)

    true_value = training_data[[v]]

    if (length(true_value)>0){

      if (v %in% c("m_th_step","m_th_cauchy", "m_alpha", "m_beta") &
          x$Fit$model_info$dormancy==F ){

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

        original_name = TOSCA:::get_original_mutation_name(x, name = v, index = NA)

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
plot_ppc_single_mut = function(x, mut1_real, mut2_real, rep_name1, rep_name2){

  estimates = TOSCA:::get_inferred_parameters(x)
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
  original_name_1 = TOSCA:::get_original_mutation_name(x, n1, as.integer(index1))
  original_name_2 = TOSCA:::get_original_mutation_name(x, n2, index2)

  ppc = ppc %>%
    dplyr::rename('m_rep_1' = rep_name1,
                  'm_rep_2' = rep_name2) %>%
    dplyr::mutate(pp = (abs(m_rep_1 - mut1_real) + abs(m_rep_1 - mut2_real) )) %>%
    dplyr::group_by(m_rep_1, m_rep_2, pp) %>%
    dplyr::summarise(n = n())

  ggplot2::ggplot(ppc,ggplot2::aes(m_rep_1, m_rep_2, color= pp, alpha = n), size=.05) +
    ggplot2::geom_point() +
    my_ggplot_theme() +
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

  posterior = TOSCA:::get_inferred_parameters(x)
  m_rep_draws = colnames(posterior)[grepl('_rep', colnames(posterior))]
  training_data = TOSCA:::get_inference_data(x, model = x$Fit$model_info$model_name,
                                     #fixed_omega = x$tosca_fit$model_info$fixed_omega,
                                     #fixed_mu = x$tosca_fit$model_info$fixed_mu,
                                     dormancy = x$Fit$model_info$dormancy)
  variables = names(training_data)[grepl("m_", names(training_data))]

  ppc_plot_list = list()
  seen = c()

  for (v1 in variables){
    #print(v1)
    for (v2 in variables){
     # print(v2)

      if (v1 != v2){

        true_value1 = training_data[[v1]]
        true_value2 = training_data[[v2]]

        if (!(v1 %in% seen)){

          if (length(true_value1)==1 & length(true_value2)==1){

            rep_name1 = paste0(v1, "_rep")
            rep_name2 = paste0(v2, "_rep")
            ppc_plot = TOSCA:::plot_ppc_single_mut(x, true_value1, true_value2, rep_name1, rep_name2)
            ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot


          }else{

            if (length(true_value1)>1 & length(true_value2)==1){
              for (i in 1:length(true_value1)){
                rep_name1 = paste0(v1,"_rep[",i,"]")
                rep_name2 = paste0(v2, "_rep")
                ppc_plot = TOSCA:::plot_ppc_single_mut(x, true_value1[i], true_value2, rep_name1, rep_name2)
                #mut1_real, mut2_real, rep_name1, rep_name2
                ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot
              }
            }

            if (length(true_value1)==1 & length(true_value2)>1){
              for (i in 1:length(true_value1)){
                rep_name2 = paste0(v2,"_rep[",i,"]")
                rep_name1 = paste0(v1, "_rep")
                ppc_plot = TOSCA:::plot_ppc_single_mut(x, true_value1, true_value2[i], rep_name1, rep_name2)
                ppc_plot_list[[length(ppc_plot_list)+1]] = ppc_plot
              }
            }

          }

          seen = c(seen, v1, v2)

        }

      }


    }

  }

  ggpubr::ggarrange(plotlist = ppc_plot_list, nrow = 2, ncol=round(length(ppc_plot_list)/2)) + ggplot2::ggtitle("Posterior Predictive Checks")


}
