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
init = function(mutations, parameters, samples, therapies){

  dfs = c(mutations, parameters, samples, therapies)
  dfs_names = c("mutations", "parameters", "samples", "therapies")

  for (df in dfs){
    if (!is.data.frame(df)) {
      stop("Input must be a data.frame.")
    }
  }

  for (t in 1:length(dfs)){
    check_required_cols(dfs[t], type=dfs_names[t])
  }

  transformed_input = check_genomic_input(mutations, parameters)
  transformed_clinical_records = check_clinical_input(samples, therapies)
  transformed_mutations = transformed_input[[1]]
  transformed_parameters = transformed_input[[2]]

  x = list(
    'input_mutations' = mutations, 'mutations' = transformed_mutations,
    'input_clinical_records' = clinical_records, 'clinical_records' = transformed_clinical_records,
    'input_parameters' = parameters, 'parameters' = transformed_parameters
    )

  class(x)= "TOSCA"
  return(x)
}

check_clinical_input = function(samples, therapies){

  clinical_records = data.frame()
  samples = samples %>% dplyr::arrange(Date)
  for (c in 1:2){
    new_entry = data.frame(
      "Clinical.name"	= "Sample",
      "Clinical.type"	= as.character(c),
      "Clinical.value.start"= TOSCA:::convert_real_date(samples$Date[c]),
      "Clinical.value.end" = NA,
      "Clinical.index"=NA,
      "Clinical.original.name"=samples$Name[c]
    )
    clinical_records = rbind(clinical_records, new_entry)
  }

  therapies = therapies %>% dplyr::arrange(Start)

  add_therapy_type = function(clinical_records, therapies, type){

    therapy_names = therapies %>% filter(Type == type) %>% pull(Name) %>% unique()
    n_therapy_types = length(therapy_names)

    if (n_therapy_types > 0){

      if (type == "Mutagenic") type="Therapy step"
      if (type == "Driver responsive") type="Driver"
      if (type == "Chemotherapy inducing dormancy") type="Chemotherapy"

      for (t in 1:n_therapy_types){
        #therapy_df = data.frame()
        sub_df = therapies %>% filter(Name==therapy_names[t])
        for (tt in 1:nrow(sub_df)){
          new_entry = data.frame(
            "Clinical.name"	= "Therapy step",
            "Clinical.type"	= type,
            "Clinical.value.start"= TOSCA:::convert_real_date(sub_df$Start[tt]),
            "Clinical.value.end" = TOSCA:::convert_real_date(sub_df$End[tt]),
            "Clinical.index"=as.character(tt),
            "Clinical.original.name"=sub_df$Name[tt]
          )
          clinical_records = rbind(clinical_records, new_entry)
        }
      }
    }

    return(clinical_records)
  }

  clinical_records = add_therapy_type(clinical_records, therapies, type="Mutagenic")
  clinical_records = add_therapy_type(clinical_records, therapies, type="Driver responsive")
  clinical_records = add_therapy_type(clinical_records, therapies, type="Chemotherapy inducing dormancy")

  return(clinical_records)
}

check_genomic_input = function(mutations, parameters, therapies){

  mutations_new = data.frame()
  parameters_new = data.frame()

  l_diploid = mutations %>% filter(Karyptype=="1:1") %>% pull(Length) %>% max()
  l_CNA = mutations %>% filter(Karyptype!="1:1") %>% pull(Length) %>% unique()
  n_cna = length(l_CNA)
  drug_names = therapies %>% dplyr::arrange(Start) %>% pull(Name) %>% unique()

  for (m in 1:nrow(mutations)){

    # source
    if (grepl("Drug", mutations$Type[m])) source = "step" else source = NA
    # coeff
    if (mutations$Karyptype[m]=="1:1") coeff = NA
    if (mutations$Karyptype[m]=="2:0") coeff = "2"
    if (mutations$Karyptype[m]=="2:2") coeff = "4"
    # name
    if (grepl("clock", mutations$Type[m])) name = "m_clock"
    if (mutations$Type[m] %in% c("alpha", "beta")) name = "m_cna"
    if (mutations$Type[m] %in% drug_names) name = "m_th"
    if (mutations$Type[m] == "driver") name = "m_driver"
    # type
    if (grepl("primary", mutations$Type[m])) type = "primary"
    if (grepl("relapse", mutations$Type[m])) type = "relapse"
    if (mutations$Type[m] %in% c("alpha", "beta")) type = mutations$Type[m]
    if (mutations$Type[m] %in% drug_names) type = NA
    if (mutations$Type[m] == "driver") type = NA
    # index
    index = NA
    if (mutations$Type[m] %in% drug_names) {
      index = which(drug_names == mutations$Type[m])
      }
    if (mutations$Type[m] %in% c("alpha", "beta")){
      l = mutations$Length[m]
      index = which(l_CNA == l)
    }
    new_entry_mutations = data.frame(
      "Mutation.name"=name,
      "Mutation.type"=type,
      "Mutation.index"=index,
      "Mutation.value"=mutations$Value[m],
      "Mutation.coeff"=coeff,
      "Mutation.source" = source,
      "Mutation.original.name"= mutations$Name[m]
    )

    mutations_new = rbind(mutations_new, new_entry_mutations)

  }

  # add lengths
  parameters_lengths = data.frame("Parameter.name"= "l_diploid","Parameter.value"=l_diploid,"Parameter.index"=NA)
  if (n_cna > 0){
    for (n in 1:n_cna){
      parameters_lengths = data.frame("Parameter.name"= "l_CNA","Parameter.value"=l_CNA[n],"Parameter.index"=as.character(n))
    }
  }
  # add coefficients
  parameters_coeff = data.frame()
  if (n_cna > 0){
    for (n in 1:n_cna){
      karyo = mutations %>% filter(Length==l_CNA[n]) %>% pull(Karyptype) %>% unique()
      if (karyo == "2:0") coeff="2" else coeff="4"
      parameters_coeff = data.frame("Parameter.name"= "coeff","Parameter.value"=coeff,"Parameter.index"=as.character(n))
    }
  }
  # add compulsory parameters if not present
  if ("m_driver" %in% mutations_new$Mutation.name)

    }

check_required_cols = function(df, type){
  if (type == "samples") required_cols <- c("Name", "Date")
  if (type == "parameters") required_cols <- c("Name", "Value","Index")
  if (type == "mutations") required_cols <- c( "Name","Length","Karyptype","Type","Value")
  if (type == "therapies") required_cols <- c( "Name","Type","Start","End")
  required_cols <- c("Name", "Value")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }
}

check_mandatory_parameters = function(parameters){

  required_params <- c("mu_clock", )
  missing_params <- setdiff(required_params, df$Name)
  if (length(missing_params) > 0) {
    stop(paste("Missing required parameter(s):", paste(missing_params, collapse = ", ")))
  }

  # 4. Extract values
  a_val <- df$Value[df$Name == "a"]
  b_val <- df$Value[df$Name == "b"]

  # 5. Validate "a" is integer
  if (length(a_val) != 1 || is.na(a_val) || a_val %% 1 != 0) {
    stop("'a' must be a single integer value.")
  }

  # 6. Validate "b" is numeric
  if (length(b_val) != 1 || is.na(b_val) || !is.numeric(b_val)) {
    stop("'b' must be a single numeric value.")
  }

  # If all checks pass, return TRUE (or the validated values)
  return(TRUE)
}



