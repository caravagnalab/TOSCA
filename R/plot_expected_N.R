#' Posterior distribution of the number of cells collected in the first and second samples.
#'
#' @param x TOSCA object
#'
#' @return Posterior distribution of the number of cells collected in the first and second samples, considering the growth rate (omega) and the time of birth of the MRCA (t_mrca) of each sample.
#' @export
plot_expected_N = function(x){

  posterior = TOSCA:::get_inferred_parameters(x) %>% dplyr::as_tibble()
  N_rel = exp(posterior$omega*(TOSCA:::get_sample(x, sample =2)-posterior$t_mrca))
  N_pre = exp(posterior$omega*(TOSCA:::get_sample(x, sample =1)-posterior$t_mrca_primary))

  if (nrow(x$Input$Samples)>2) i = 1 else i=0
  primary_name = (x$Input$Samples %>% pull(Name))[1+i] #%>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "1") %>% dplyr::pull(Clinical.original.name)
  relapse_name = (x$Input$Samples %>% pull(Name))[2+i] #x$clinical_records %>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "2") %>% dplyr::pull(Clinical.original.name)

  N_rel_plot = ggplot2::ggplot() + TOSCA:::my_ggplot_theme() +
    ggplot2::geom_histogram(ggplot2::aes(x=log(N_rel)))+
    ggplot2::scale_x_log10() +
    ggplot2::xlab(paste0("Tumor size at ", relapse_name, "\n(log-scale)"))
  N_pre_plot = ggplot2::ggplot() +  TOSCA:::my_ggplot_theme() +
    ggplot2::geom_histogram(ggplot2::aes(x=log(N_pre)))+
    ggplot2::xlab(paste0("Tumor size at ", primary_name, "\n(log-scale)"))

  ggpubr::ggarrange(plotlist = list(N_rel_plot, N_pre_plot), nrow = 1)
}
