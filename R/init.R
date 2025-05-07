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
init = function(mutations, clinical_records, parameters){

  check_input_mutations(mutations)
  check_input_clinical(clinical_records)
  check_parameters(parameters)

  x = list('mutations' = mutations, 'clinical_records' = clinical_records, 'parameters' = parameters)
  class(x)= "TOSCA"
  return(x)
}
