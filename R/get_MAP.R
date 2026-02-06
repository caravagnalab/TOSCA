get_index_from_drug_name = function(x, drug_name){
  mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="alpha_th_step")
  mu_alpha = left_join(
    x$Input$Therapies %>% filter(Class=="Mutagenic"),
    mu_alpha %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
    dplyr::arrange(as.Date(Start))
  which(mu_alpha$Name == drug_name)
}

get_index_from_cna_length = function(x, len){
  index = x$Input$Mutations %>% dplyr::filter(Type == "alpha") %>%
    dplyr::arrange(Length)
  which(index$Length == len)
}

#' Get mean, mode and q5, q95
#'
#' @param x TOSCA obj
#' @param parameter string, optional. Parameter of interest. If nothing is provided the function returns the mean, mode and q5, q95, of all parameters. Possible values: t_eca, t_mrca, t_cna, t_driver, omega, mu_th_step, mu_driver.
#' @param cna_length if parameter == t_cna, length of the CNA of interest
#' @param drug_name if parameter == mu_th_step, name of the drug of interest
#' @param quantiles vector of quantiles, default (5, 95)
#'
#' @return a dataframe summarising the inference results
#' @export
#'
#' @examples
#' data("exampleFit")
#' get_fit_summary(exampleFit, parameter = "t_cna", cna_length = 4.5e+08, quantiles=c(5, 95))
#' get_fit_summary(exampleFit, parameter = "mu_th_step", drug_name = "Drug 2", quantiles=c(2, 98))
get_fit_summary = function(x, parameter=NULL,
                           cna_length=NULL, drug_name=NULL, quantiles=c(5, 95)){

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

  summary = summary %>% dplyr::filter(variable %in% vars_of_interest) %>% select(variable, median, mean, rhat, ess_bulk)

  quantiles_df = lapply(vars_of_interest, function(f){
      q = TOSCA:::get_inferred_parameters(x) %>% pull(f) %>% quantile(probs = c(quantiles[1]/100, quantiles[2]/100))
      q = data.frame("variable" = f,
                 "q1" = q[[paste0(quantiles[1],"%")]],
                 "q2" = q[[paste0(quantiles[2],"%")]]
                 )
      if (f %in% timing_vars){
        q = q %>% dplyr::mutate(q1 = TOSCA:::convert_date_real(x=x, date = q1),
                                q2 = TOSCA:::convert_date_real(x=x, date = q2))
      }
      colnames(q) = c("variable", paste0("q",quantiles[1]), paste0("q",quantiles[2]))
      q
      #names(q) = f
  }) %>% Reduce(rbind, .)


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

  summary = summary %>%
    dplyr::mutate(
      mean = ifelse(variable %in% timing_vars, TOSCA:::convert_date_real(x=x, date = mean), mean),
      median = ifelse(variable %in% timing_vars, TOSCA:::convert_date_real(x=x, date = median), median))
      #q5 = ifelse(variable %in% timing_vars, TOSCA:::convert_date_real(x=x, date = q5), q5),
      #q95 = ifelse(variable %in% timing_vars, TOSCA:::convert_date_real(x=x, date = q95), q95))
  quantiles_df$variable <- ifelse(quantiles_df$variable %in% names(new_names_par),
                             new_names_par[quantiles_df$variable],
                             quantiles_df$variable)
  quantiles_df$variable <- ifelse(quantiles_df$variable %in% names(new_names_ppc),
                             new_names_ppc[quantiles_df$variable],
                             quantiles_df$variable)

  summary$variable <- ifelse(summary$variable %in% names(new_names_par),
                             new_names_par[summary$variable],
                             summary$variable)
  summary$variable <- ifelse(summary$variable %in% names(new_names_ppc),
                             new_names_ppc[summary$variable],
                             summary$variable)

  summary = dplyr::left_join(summary, quantiles_df, by = "variable")
  # Select only variables of interest (timing, rate, ppc)
  # Convert the timings to real times
  # convert the names of rates and ppc

  if (!is.null(parameter)){
    if (!is.null(cna_length)) {
      index = TOSCA:::get_index_from_cna_length(x, cna_length)
      parameter = paste0(parameter, "[", index,"]")
    }
    if (!is.null(drug_name)) {
      parameter = paste0("mu_", drug_name)
      }
    #index = TOSCA:::get_index_from_drug_name(x, drug_name)
    #if (!is.null(cna_length) | !is.null(cna_length)) parameter = paste0(parameter, "[", index,"]")
    summary = summary %>% dplyr::filter(variable == parameter)
  }

  summary %>% dplyr::as_tibble()
}

#' Computes the days of the event of interest from a date of interest
#'
#' @param x TOSCA obj
#' @param parameter the inferred time of interest
#' @param from the date of interest (string, "YY-mm-dd")
#' @param cna_length if the inferred date refers to a CNA, the length of the CNA
#' @param quantiles
#'
#' @return the distance in days between an arbitrary date (parameter from) and a chosen inferred time
#' @export
#'
#' @examples
#' data("exampleFit")
#' days_from(exampleFit, parameter = "t_cna", from = "2010-01-04", cna_length = 4.5e+08, quantiles=c(5, 95))
days_from = function(x, parameter, from, cna_length=NULL, quantiles = c(5,95)){
  #TOSCA:::get_index_from_cna_length(x, cna_length)
  from_chosen_date = function(f) as.Date(f) - as.Date(from)

  get_fit_summary(x=x, parameter=parameter, cna_length=cna_length, quantiles=quantiles) %>%
    dplyr::select(!c("rhat", "ess_bulk", "variable")) %>%
    dplyr::mutate(across(everything(), from_chosen_date))

}

# get_MAP_t_eca = function(x, map = "mean"){
#   if(map=="mean") map = mean(x$Fit$posterior$t_eca) else map = median(x$Fit$posterior$t_eca)
#   TOSCA:::convert_date_real(date = map, x = x)
# }
# get_MAP_t_cna = function(x, index=NA, map = "mean"){
#   if (is.na(index)) var="t_cna" else var = paste0("t_cna[",index,"]")
#   if(map=="mean") map = mean(x$Fit$posterior[[var]]) else map = median(x$Fit$posterior[[var]])
#   TOSCA:::convert_date_real(date = map, x = x)
# }
# get_MAP_t_mrca = function(x, map = "mean"){
#   if(map=="mean") map = mean(x$Fit$posterior$t_mrca) else map = median(x$Fit$posterior$t_mrca)
#   TOSCA:::convert_date_real(date = map, x = x)
# }
# get_MAP_growth_rate = function(x, map = "mean"){
#   if(map=="mean") map = mean(x$Fit$posterior$omega) else map = median(x$Fit$posterior$omega)
#   map
# }
# get_MAP_mutation_rate = function(x, type = "th_step", name=NA, map = "mean"){
#   # type = "th_step", "driver
#   if (is.na(index)) var=paste0("mu_", type) else var = paste0("mu_",type, "[",index,"]")
#   # mu_th_step[1]
#   if(map=="mean") map = mean(x$Fit$posterior[[var]]) else map = median(x$Fit$posterior[[var]])
#   map
# }
