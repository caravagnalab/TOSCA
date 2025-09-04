#' Initialise TOSCA object
#'
#' @param mutations data.frame with the information on mutation falling on the timed branches.
#'  1. Name (str), name of choosing for the cluster of mutations
#'  2. Length (int), length of the genomic segment they sit on
#'  3. Karyptype (str), karyotype of the genomic segment they sit as "Major:minor" (ex "1:1")
#'  4. Type (str), one between the following: "clock-like primary", "clock-like relapse", "alpha", "beta", "driver", the name of the mutagenic drug to which they are associated as reported in the therapies dataframe
#'  5. Value (int), number of mutations
#' @param parameters data.frame with the parameters required for inference.
#'  1. Name (str),
#'  2. Value (real),
#'  3. Index (int)
#' @param samples data.frame with the name and collection date of the longitudinal samples.
#'  1. Name (str), name of the sample
#'  2. Date (date) date of collection in the format "YYYY-mm-dd"
#' @param therapies data.frame with the name and exposure period to therapies affecting the mutation accumulation process.
#'  1. Name (str)
#'  2. Class (str), one of the following : "Mutagenic step", "Mutagenic cauchy", "Driver responsive", "Chemotherapy inducing dormancy"
#'  3. Start (date), start of the cicle of therapy in the format "YYYY-mm-dd"
#'  4. End (date), end of the cicle of therapy in the format "YYYY-mm-dd"
#'
#' @return TOSCA object
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
#'
init = function(mutations, parameters, samples, therapies){

  dfs = list(mutations, parameters, samples, therapies)
  dfs_names = c("mutations", "parameters", "samples", "therapies")

  for (df in dfs){
    if (!is.data.frame(df)) {
      stop("Input must be a data.frame.")
    }
  }

  for (t in 1:length(dfs)){
    check_required_cols(dfs[[t]], type=dfs_names[t])
  }

  TOSCA:::check_clinical_input(samples)
  TOSCA:::check_clinical_input(therapies)
  TOSCA:::check_genomic_input(mutations)
  TOSCA:::check_parameters_input(parameters)

  x = list(
    "Input" = list(
    "Samples" = samples,
    "Mutations" = mutations,
    "Parameters" = parameters,
    "Therapies" = therapies
    )

    )

  class(x)= "TOSCA"
  return(x)
}

