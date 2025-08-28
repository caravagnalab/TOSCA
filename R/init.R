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

  dfs = list(mutations, parameters, samples, therapies)
  dfs_names = c("mutations", "parameters", "samples", "therapies")

  for (df in dfs){
    #print(df)
    if (!is.data.frame(df)) {
      stop("Input must be a data.frame.")
    }
  }

  for (t in 1:length(dfs)){
    #print(dfs_names[t])
    #print(dfs[t])
    check_required_cols(dfs[[t]], type=dfs_names[t])
  }

  transformed_clinical_records = check_clinical_input(samples, therapies)
  transformed_input = check_genomic_input(mutations, parameters, transformed_clinical_records)
  transformed_mutations = transformed_input[[1]]
  transformed_parameters = transformed_input[[2]]

  x = list(
    'clinical_records' = transformed_clinical_records,
    'mutations' = transformed_mutations,
    'parameters' = transformed_parameters,
    "inputs" = list('input_mutations' = mutations,
                    'input_samples' = samples, "input_therapies"="therapies",
                    'input_parameters' = parameters)
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

  add_therapy_name = function(clinical_records, therapies, type){

    therapy_names = therapies %>% filter(Class == type) %>% pull(Name) %>% unique()
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
            "Clinical.name"	= type,
            "Clinical.type"	= t,
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

  clinical_records = add_therapy_name(clinical_records, therapies, type="Mutagenic")
  clinical_records = add_therapy_name(clinical_records, therapies, type="Driver responsive")
  clinical_records = add_therapy_name(clinical_records, therapies, type="Chemotherapy inducing dormancy")

  return(clinical_records)
}

check_genomic_input = function(mutations, parameters, transformed_clinical_records){

  mutations_new = data.frame()
  parameters_new = data.frame()

  l_diploid = mutations %>% filter(Karyptype=="1:1") %>% pull(Length) %>% max()
  l_CNA = mutations %>% filter(Karyptype!="1:1") %>% pull(Length) %>% unique()
  n_cna = length(l_CNA)
  drug_names = transformed_clinical_records %>% filter(Clinical.name == "Therapy step") %>%
    dplyr::arrange(Clinical.value.start) %>% pull(Clinical.original.name) %>% unique()

  ### Create mutations dataframe
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


  ### Create Parameters dataframe
  colnames(parameters) = c("Parameter.name","Parameter.value","Parameter.index")

  # add lengths
  parameters_lengths = data.frame("Parameter.name"= "l_diploid","Parameter.value"=l_diploid,"Parameter.index"=NA)
  if (n_cna > 0){
    for (n in 1:n_cna){
      parameters_lengths = rbind(
        parameters_lengths,
        data.frame("Parameter.name"= "l_CNA","Parameter.value"=l_CNA[n],"Parameter.index"=as.character(n)))
    }
  }

  # add coefficients
  parameters_coeff = data.frame()
  if (n_cna > 0){
    for (n in 1:n_cna){
      karyo = mutations %>% filter(Length==l_CNA[n]) %>% pull(Karyptype) %>% unique()
      if (karyo == "2:0") coeff="2" else coeff="4"
      parameters_coeff = rbind(
        parameters_coeff,
        data.frame("Parameter.name"= "coeff","Parameter.value"=coeff,"Parameter.index"=as.character(n)))
    }
  }

  # check compulsory parameters
  if (!("mu_clock" %in% parameters$Parameter.name)) stop(paste("Missing required parameter: mu_clock"))
  mu_clock = parameters %>% filter("Parameter.name" == "mu_clock") %>% pull("Parameter.value")

  # add compulsory parameters if not present
  # apart from mu_clock, which is compulsory, all mu, if not provided are approximated
  # based on the number of mutation and the exposure time
  get_omega_approximation = function(transformed_clinical_records){
    if ("omega_alpha" %in% parameters$Parameter.name){
      omega_alpha = parameters %>% filter(Parameter.name=="omega_alpha") %>% pull(Parameter.value)
      omega_beta = parameters %>% filter(Parameter.name=="omega_beta") %>% pull(Parameter.value)
      omega = omega_alpha/omega_beta
    }
    sample_1= transformed_clinical_records %>% filter(Clinical.name == "Sample", Clinical.type == "1") %>% pull(Clinical.value.start)
    sample_2= transformed_clinical_records %>% filter(Clinical.name == "Sample", Clinical.type == "2") %>% pull(Clinical.value.start)
    omega_lb = log(10e6) / (sample_2 - sample_1)
    omega_ub = log(10e10) / (sample_2-(sample_2 - 30/365))
    omega = mean(c(omega_lb, omega_ub))
    return(omega)
  }
  if (!("omega_alpha" %in% parameters$Parameter.name)){
    omega_alpha_approx = get_omega_approximation(transformed_clinical_records) / 10
    omega_beta_approx = .1
    parameters = rbind(parameters,
                       data.frame("Parameter.name"= c("omega_alpha", "omega_beta"),
                                  "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
  }
  if ("m_driver" %in% mutations_new$Mutation.name & !("mu_clock_driver" %in% parameters$Name)){
    parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_clock_driver"),
                                              "Parameter.value"=mu_clock,"Parameter.index"= NA))
  }
  if ("m_driver" %in% mutations_new$Mutation.name & !("mu_driver" %in% parameters$Name)){
    exposure_lb = transformed_clinical_records %>% filter(Clinical.type=="Driver") %>% pull(Clinical.value.start)
    exposure_lb = exposure_lb[1]
    exposure_ub = transformed_clinical_records %>% filter(Clinical.type=="Driver") %>% pull(Clinical.value.end)
    exposure_ub = exposure_ub[length(exposure_ub)]

    m_driver = mutations_new %>% filter(Mutation.name=="m_driver") %>% pull(Mutation.value)
    omega = get_omega_approximation(transformed_clinical_records)
    mu_driver = m_driver / (2*l_diploid*omega*(exposure_ub-exposure_lb))

    parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_driver", "omega_beta"),
                                  "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
  }
  for (th in drug_names){
    th_index = transformed_clinical_records %>% filter(Clinical.original.name==th) %>% pull(Clinical.type) %>% unique()
    if (parameters %>% filter(Parameter.name=="mu_th_step", Parameter.index==th_index) %>% nrow()==0){
      exposure_lb = transformed_clinical_records %>% filter(Clinical.original.name==th) %>% pull(Clinical.value.start)
      exposure_lb = exposure_lb[1]
      exposure_ub = transformed_clinical_records %>% filter(Clinical.original.name==th) %>% pull(Clinical.value.end)
      exposure_ub = exposure_ub[length(exposure_ub)]

      m_th = mutations_new %>% filter(Mutation.name=="m_th", Mutation.index==th_index) %>% pull(Mutation.value)
      omega = get_omega_approximation(transformed_clinical_records)
      mu_driver = m_driver / (2*l_diploid*omega*(exposure_ub-exposure_lb))

      parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_driver", "omega_beta"),
                                                "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
    }
  }
  if (!("k_step" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= "k_step",
                                              "Parameter.value"=10,"Parameter.index"=NA))
  }
  if (!("N_min" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= c("N_min","N_min","N_max","N_max"),
                                              "Parameter.value"=c(1e6,1e6,1e10,1e10),
                                              "Parameter.index"=c("1","2","1","2")))
  }
  if (!("exponential_growth" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= "exponential_growth",
                                              "Parameter.value"=1,"Parameter.index"=NA))
  }
  if (!("mrca_alpha" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= c("mrca_alpha","mrca_beta"),
                                              "Parameter.value"=c(1,1),"Parameter.index"=c(1,1)))
  }
  if (!("phi_driver" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= "phi_driver",
                                              "Parameter.value"=0.05,"Parameter.index"=NA))
  }
  for (th in drug_names){
    th_index = transformed_clinical_records %>% filter(Clinical.original.name==th) %>% pull(Clinical.type) %>% unique()
    if (parameters %>% filter(Parameter.name=="phi_th_step", Parameter.index==th_index) %>% nrow()==0){
      parameters = rbind(parameters, data.frame("Parameter.name"= "phi_th_step",
                                                "Parameter.value"=0.05,"Parameter.index"=th_index))
    }
  }
  n_phi_cna = sum(parameters$Parameter.name == "phi_cna")
  if (n_phi_cna!=n_cna) {
    for (c in 1:n_cna){
      index = mutations_new %>% filter(Mutation.name=="m_cna") %>% pull(Mutation.index) %>% unique()
      index = index[c]
      if (parameters %>% filter(Parameter.name=="phi_th_step", Parameter.index==th_index) %>% nrow()==0){
        parameters = rbind(parameters, data.frame("Parameter.name"= "phi_cna",
                                                  "Parameter.value"=0.05,"Parameter.index"=index))
      }
    }
  }
  if (!("phi_clock" %in% parameters$Parameter.name)) {
    parameters = rbind(parameters, data.frame("Parameter.name"= "phi_clock",
                                              "Parameter.value"=0.05,"Parameter.index"=NA))
  }

  parameters_new = rbind(
    parameters,
    parameters_lengths,
    parameters_coeff
  )

  return(list(mutations_new, parameters_new))
    }

check_required_cols = function(df, type){
  if (type == "samples") required_cols <- c("Name", "Date")
  if (type == "parameters") required_cols <- c("Name", "Value","Index")
  if (type == "mutations") required_cols <- c( "Name","Length","Karyptype","Type","Value")
  if (type == "therapies") required_cols <- c( "Name","Class","Start","End")
  #required_cols <- c("Name", "Value")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }
}

# check_mandatory_parameters = function(parameters){
#
#   required_params <- c("mu_clock", )
#   missing_params <- setdiff(required_params, df$Name)
#   if (length(missing_params) > 0) {
#     stop(paste("Missing required parameter(s):", paste(missing_params, collapse = ", ")))
#   }
#
#   # 4. Extract values
#   a_val <- df$Value[df$Name == "a"]
#   b_val <- df$Value[df$Name == "b"]
#
#   # 5. Validate "a" is integer
#   if (length(a_val) != 1 || is.na(a_val) || a_val %% 1 != 0) {
#     stop("'a' must be a single integer value.")
#   }
#
#   # 6. Validate "b" is numeric
#   if (length(b_val) != 1 || is.na(b_val) || !is.numeric(b_val)) {
#     stop("'b' must be a single numeric value.")
#   }
#
#   # If all checks pass, return TRUE (or the validated values)
#   return(TRUE)
# }



