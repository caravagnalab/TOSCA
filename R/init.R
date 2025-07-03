#' Initialise a TOSCA object
#'
#' @param mutations dataframe containing input mutations
#' @param clinical_records
#' @param parameters
#'
#' @return
#' @export
#'
#' @examples
init = function(mutations, clinical_records, parameters, delta_omega = 0.02, growth_rate_variance=10){

  # check_input_mutations(mutations)
  # check_input_clinical(clinical_records)
  # parameters = check_parameters(mutations, clinical_records, parameters, delta_omega = delta_omega, growth_rate_variance=growth_rate_variance)

  x = list('mutations' = mutations, 'clinical_records' = clinical_records, 'parameters' = parameters)
  class(x)= "TOSCA"
  return(x)
}
