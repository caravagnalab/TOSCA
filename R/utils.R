#' Check clinical input
#'
#' @param clinical_input
#'
#' @return The function returns an error of the format of the "samples" dataframe has errors
#'
#' @examples
check_clinical_input = function(clinical_input){

  if (nrow(clinical_input) > 0){

  if (class(clinical_input$Name) != "character"){
    stop("Sample names and Therapy names must be strings")
  }


  if ("Date" %in% colnames(clinical_input)){
    for (d in 1:nrow(clinical_input)){
      #print(d)
      if (!TOSCA:::check_valid_date(clinical_input$Date[d])){
        stop("Invalid date format. Date format must be: yyyy-mm-dd")
      }
    }
  }else{
    for (d in 1:nrow(clinical_input)){
      if (!check_valid_date(clinical_input$Start[d])){
        stop("Invalid date format. Date format must be: yyyy-mm-dd")
      }
    }

    for (d in 1:nrow(clinical_input)){
      if (!check_valid_date(clinical_input$End[d])){
        stop("Invalid date format. Date format must be: yyyy-mm-dd")
      }
    }
  }

  if ("Class" %in% colnames(clinical_input)){

    therapy_types = c("Mutagenic", "Mutagenic cauchy", "Driver responsive", "Chemotherapy inducing dormancy")

    for (name in clinical_input$Class){
      if (!(name %in% therapy_types)){
        stop(paste0("Unrecognised class for ",name, ". Select one of the following classes: ", paste(therapy_types, collapse= " , ")))
      }
    }
  }
  }

}

#' Check that the input mutation dataframe is valid
#'
#' @param mutations
#'
#' @return
#'
#' @examples
check_genomic_input = function(mutations){
  if (class(mutations$Name) != "character"){stop("Mutation names must be strings")}
  if (class(mutations$Length) != "numeric"){stop("Genome segments lengths must be integer numbers")}
  for (k in 1:nrow(mutations)){
    if (!TOSCA:::check_valid_karyotype(mutations$Karyotype[k])) stop(paste0("Invalid karyotype: ", k))
  }
  mutation_types = c("clock-like primary", "clock-like relapse", "alpha", "beta", "driver")
  non_standard_mutations = mutations$Type[!(mutations$Type %in% mutation_types)]
  if (length(non_standard_mutations > 0)) warning(paste0("The following mutations do not follow the conventional naming scheme: ",
                                                         paste(non_standard_mutations, collapse = " , "),
                                                         ". If these mutations are associated to a therapy-related mutational process, make sure that the name reported in 'Type' is the same as the one reported in the 'Therapy' dataframe for the corresponding mutational process. If that is not the case, please choose a valid mutation type between the following: ",
                                                         paste(mutation_types, collapse=" , ")))
  if (class(mutations$Value) != "numeric") stop("Mutations must be integer numbers")
}

#' Checks that the parameters dataframe is valid
#'
#' @param parameters
#'
#' @return
#'
#' @examples
check_parameters_input = function(parameters){
  valid_parameters_names = c(
    "mu_clock",
    "mu_driver_clock",
    "mu_driver",
    "mu_th_step",
    "location_cauchy",
    "scales_cauchy",
    "omega_alpha",
    "omega_beta",
    "exponential_growth",
    "N_min",
    "N_max",
    "mrca_alpha",
    "mrca_beta",
    "phi_clock",
    "phi_driver",
    "phi_cna",
    "phi_th_step",
    "phi_th_cauchy",
    "wgd",
    "k_step",
    "alpha_th_step",
    "beta_th_step",
    "phi_cna_beta",
    "phi_cna_alpha",
    "mu_driver_alpha",
    "mu_driver_beta"
  )
  invalid_parameters = parameters$Name[!(parameters$Name %in% valid_parameters_names)]
  if (length(invalid_parameters > 1)){
    stop(paste0("Invalid parameters: ", paste(invalid_parameters, collapse="\n"),
                                              ", please select a valid name in: ",
                                              paste(valid_parameters_names, collapse="\n")))
  }
}

#' Check that the input karytypes in the mutations dataframe are valid
#'
#' @param karyo
#'
#' @return
#'
#' @examples
check_valid_karyotype = function(karyo){

  k1= strsplit(karyo, ":")[[1]][1]
  k2= strsplit(karyo, ":")[[1]][2]

  class(karyo) == "character" & as.integer(k1) >= as.integer(k2) & as.integer(k1) <= 2 & as.integer(k2) <= 2
}

#' Check that the input date is valid
#'
#' @param str_date
#'
#' @return
#'
#' @examples
check_valid_date = function(str_date){
  y = strsplit(str_date, "-")[[1]][1]
  condition_y = nchar(y) == 4 & as.integer(y) >= 1900
  m = strsplit(str_date, "-")[[1]][2]
  condition_m = as.integer(m) <= 12
  d = strsplit(str_date, "-")[[1]][3]
  condition_d = as.integer(d) <= 31

  condition_y & condition_m & condition_d
}

#' Check that for the required columns for each dataframe are present
#'
#' @param df
#' @param type
#'
#' @return
#'
#' @examples
check_required_cols = function(df, type){
  if (type == "samples") required_cols <- c("Name", "Date")
  if (type == "parameters") required_cols <- c("Name", "Value","Index")
  if (type == "mutations") required_cols <- c( "Name","Length","Karyotype","Type","Value")
  if (type == "therapies") required_cols <- c( "Name","Class","Start","End")
  #required_cols <- c("Name", "Value")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }
}

#' Converts date into real number
#'
#' @param x
#' @param ref_year
#'
#' @return
#'
#' @examples
convert_date_real = function(date, x, ref_year = 1900) {

  if (x$Input$Samples %>% nrow() > 2){
    ref_year = as.integer(strsplit((x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Date))[1], "-")[[1]][1])
  }

  y = ref_year + date %>% floor()
  py = (date - (date %>% floor())) * 365
  m = (py / 30) %>% floor
  d = (py - (m * 30)) %>% floor
  date_string = paste(y, m +1, d + 1, sep = '-')

  date_string = ifelse (d>30, paste(y, m +1, 30, sep = '-'), date_string)
  date_string = ifelse ((m==1 & d>27), paste(y, m +1, 28, sep = '-'), date_string)
  date_string = ifelse (m>=12, paste(y+1, 1, 1, sep = '-'), date_string)

  return(date_string)
}

#' Converts real number into date
#'
#' @param date
#' @param ref_year
#'
#' @return
#'
#' @examples
convert_real_date = function(x, date = NULL, ref_year = 1900) {

  if (x$Input$Samples %>% nrow() > 2){
    ref_year = as.integer(strsplit((x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Date))[1], "-")[[1]][1])
  }

  ref_month = 1
  ref_day = 1

  year = as.integer(strsplit(date, '-')[[1]][1])
  month = as.integer(strsplit(date, '-')[[1]][2])
  day = as.integer(strsplit(date, '-')[[1]][3])

  return((year - ref_year) + (month / 12 - ref_month / 12) + (day / 365 - ref_day / 365))
}

# check_genomic_input = function(mutations, parameters, transformed_clinical_records){
#
#   mutations_new = data.frame()
#   parameters_new = data.frame()
#
#   l_diploid = mutations %>% dplyr::filter(Karyptype=="1:1") %>% dplyr::pull(Length) %>% max()
#   l_CNA = mutations %>% dplyr::filter(Karyptype!="1:1") %>% dplyr::pull(Length) %>% unique()
#   n_cna = length(l_CNA)
#   drug_names = transformed_clinical_records %>% dplyr::filter(Clinical.name == "Therapy step") %>%
#     dplyr::arrange(Clinical.value.start) %>% dplyr::pull(Clinical.original.name) %>% unique()
#
#   ### Create mutations dataframe
#   for (m in 1:nrow(mutations)){
#
#     # source
#     if (grepl("Drug", mutations$Type[m])) source = "step" else source = NA
#     # coeff
#     if (mutations$Karyptype[m]=="1:1") coeff = NA
#     if (mutations$Karyptype[m]=="2:0") coeff = "2"
#     if (mutations$Karyptype[m]=="2:2") coeff = "4"
#     # name
#     if (grepl("clock", mutations$Type[m])) name = "m_clock"
#     if (mutations$Type[m] %in% c("alpha", "beta")) name = "m_cna"
#     if (mutations$Type[m] %in% drug_names) name = "m_th"
#     if (mutations$Type[m] == "driver") name = "m_driver"
#     # type
#     if (grepl("primary", mutations$Type[m])) type = "primary"
#     if (grepl("relapse", mutations$Type[m])) type = "relapse"
#     if (mutations$Type[m] %in% c("alpha", "beta")) type = mutations$Type[m]
#     if (mutations$Type[m] %in% drug_names) type = NA
#     if (mutations$Type[m] == "driver") type = NA
#     # index
#     index = NA
#     if (mutations$Type[m] %in% drug_names) {
#       index = which(drug_names == mutations$Type[m])
#     }
#     if (mutations$Type[m] %in% c("alpha", "beta")){
#       l = mutations$Length[m]
#       index = which(l_CNA == l)
#     }
#     new_entry_mutations = data.frame(
#       "Mutation.name"=name,
#       "Mutation.type"=type,
#       "Mutation.index"=index,
#       "Mutation.value"=mutations$Value[m],
#       "Mutation.coeff"=coeff,
#       "Mutation.source" = source,
#       "Mutation.original.name"= mutations$Name[m]
#     )
#
#     mutations_new = rbind(mutations_new, new_entry_mutations)
#
#   }
#
#
#   ### Create Parameters dataframe
#   #colnames(parameters) = c("Parameter.name","Parameter.value","Parameter.index")
#   parameters = parameters %>% dplyr::rename(Parameter.name=Name, Parameter.value=Value,Parameter.index=Index)
#
#   # add lengths
#   parameters_lengths = data.frame("Parameter.name"= "l_diploid","Parameter.value"=l_diploid,"Parameter.index"=NA)
#   if (n_cna > 0){
#     for (n in 1:n_cna){
#       parameters_lengths = rbind(
#         parameters_lengths,
#         data.frame("Parameter.name"= "l_CNA","Parameter.value"=l_CNA[n],"Parameter.index"=as.character(n)))
#     }
#   }
#
#   # add coefficients
#   parameters_coeff = data.frame()
#   if (n_cna > 0){
#     for (n in 1:n_cna){
#       karyo = mutations %>% dplyr::filter(Length==l_CNA[n]) %>% dplyr::pull(Karyptype) %>% unique()
#       if (karyo == "2:0") coeff="2" else coeff="4"
#       parameters_coeff = rbind(
#         parameters_coeff,
#         data.frame("Parameter.name"= "coeff","Parameter.value"=coeff,"Parameter.index"=as.character(n)))
#     }
#   }
#
#   # check compulsory parameters
#   if (!("mu_clock" %in% parameters$Parameter.name)) stop(paste("Missing required parameter: mu_clock"))
#   mu_clock = parameters %>% dplyr::filter(Parameter.name == "mu_clock") %>% dplyr::pull(Parameter.value)
#
#   # add compulsory parameters if not present
#   # apart from mu_clock, which is compulsory, all mu, if not provided are approximated
#   # based on the number of mutation and the exposure time
#   get_omega_approximation = function(transformed_clinical_records){
#     if ("omega_alpha" %in% parameters$Parameter.name){
#       omega_alpha = parameters %>% dplyr::filter(Parameter.name=="omega_alpha") %>% dplyr::pull(Parameter.value)
#       omega_beta = parameters %>% dplyr::filter(Parameter.name=="omega_beta") %>% dplyr::pull(Parameter.value)
#       omega = omega_alpha/omega_beta
#     }
#     sample_1= transformed_clinical_records %>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "1") %>% dplyr::pull(Clinical.value.start)
#     sample_2= transformed_clinical_records %>% dplyr::filter(Clinical.name == "Sample", Clinical.type == "2") %>% dplyr::pull(Clinical.value.start)
#     omega_lb = log(10e6) / (sample_2 - sample_1)
#     omega_ub = log(10e10) / (sample_2-(sample_2 - 30/365))
#     omega = mean(c(omega_lb, omega_ub))
#     return(omega)
#   }
#   if (!("omega_alpha" %in% parameters$Parameter.name)){
#     omega_alpha_approx = get_omega_approximation(transformed_clinical_records) / 10
#     omega_beta_approx = .1
#     parameters = rbind(parameters,
#                        data.frame("Parameter.name"= c("omega_alpha", "omega_beta"),
#                                   "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
#   }
#   if ("m_driver" %in% mutations_new$Mutation.name & !("mu_clock_driver" %in% parameters$Parameter.name)){
#     parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_clock_driver"),
#                                               "Parameter.value"=mu_clock,"Parameter.index"= NA))
#   }
#   if ("m_driver" %in% mutations_new$Mutation.name & !("mu_driver" %in% parameters$Parameter.name)){
#     cycles_of_driver = transformed_clinical_records %>% dplyr::filter(Clinical.type=="Driver")
#     exposure = 0
#     for (c in 1:length()){
#       exposure_lb = transformed_clinical_records %>% dplyr::filter(Clinical.type=="Driver") %>% dplyr::pull(Clinical.value.start)
#       exposure_lb = exposure_lb[c]
#       exposure_ub = transformed_clinical_records %>% dplyr::filter(Clinical.type=="Driver") %>% dplyr::pull(Clinical.value.end)
#       exposure_ub = exposure_ub[c]
#       exposure = exposure + (exposure_ub-exposure_lb)
#     }
#
#     m_driver = mutations_new %>% dplyr::filter(Mutation.name=="m_driver") %>% dplyr::pull(Mutation.value)
#     omega = get_omega_approximation(transformed_clinical_records)
#     mu_driver = m_driver / (2*l_diploid*omega*(exposure))
#
#     parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_driver", "omega_beta"),
#                                               "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
#   }
#   for (th in drug_names){
#     th_index = transformed_clinical_records %>% dplyr::filter(Clinical.original.name==th) %>% dplyr::pull(Clinical.type) %>% unique()
#     if (parameters %>% filter(Parameter.name=="mu_th_step", Parameter.index==th_index) %>% nrow()==0){
#       exposure_lb = transformed_clinical_records %>% dplyr::filter(Clinical.original.name==th) %>% dplyr::pull(Clinical.value.start)
#       exposure_lb = exposure_lb[1]
#       exposure_ub = transformed_clinical_records %>% dplyr::filter(Clinical.original.name==th) %>% dplyr::pull(Clinical.value.end)
#       exposure_ub = exposure_ub[length(exposure_ub)]
#
#       m_th = mutations_new %>% dplyr::filter(Mutation.name=="m_th", Mutation.index==th_index) %>% dplyr::pull(Mutation.value)
#       omega = get_omega_approximation(transformed_clinical_records)
#       mu_driver = m_driver / (2*l_diploid*omega*(exposure_ub-exposure_lb))
#
#       parameters = rbind(parameters, data.frame("Parameter.name"= c("mu_driver", "omega_beta"),
#                                                 "Parameter.value"=c(omega_alpha_approx, omega_beta_approx),"Parameter.index"=c(NA,NA)))
#     }
#   }
#   if (!("k_step" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= "k_step",
#                                               "Parameter.value"=10,"Parameter.index"=NA))
#   }
#   if (!("N_min" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= c("N_min","N_min","N_max","N_max"),
#                                               "Parameter.value"=c(1e6,1e6,1e10,1e10),
#                                               "Parameter.index"=c("1","2","1","2")))
#   }
#   if (!("exponential_growth" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= "exponential_growth",
#                                               "Parameter.value"=1,"Parameter.index"=NA))
#   }
#   if (!("mrca_alpha" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= c("mrca_alpha","mrca_beta"),
#                                               "Parameter.value"=c(1,1),"Parameter.index"=c(1,1)))
#   }
#   if (!("phi_driver" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= "phi_driver",
#                                               "Parameter.value"=0.05,"Parameter.index"=NA))
#   }
#   for (th in drug_names){
#     th_index = transformed_clinical_records %>% dplyr::filter(Clinical.original.name==th) %>% dplyr::pull(Clinical.type) %>% unique()
#     if (parameters %>% dplyr::filter(Parameter.name=="phi_th_step", Parameter.index==th_index) %>% nrow()==0){
#       parameters = rbind(parameters, data.frame("Parameter.name"= "phi_th_step",
#                                                 "Parameter.value"=0.05,"Parameter.index"=th_index))
#     }
#   }
#   n_phi_cna = sum(parameters$Parameter.name == "phi_cna")
#   if (n_phi_cna!=n_cna) {
#     for (c in 1:n_cna){
#       index = mutations_new %>% dplyr::filter(Mutation.name=="m_cna") %>% dplyr::pull(Mutation.index) %>% unique()
#       index = index[c]
#       if (parameters %>% dplyr::filter(Parameter.name=="phi_th_step", Parameter.index==th_index) %>% nrow()==0){
#         parameters = rbind(parameters, data.frame("Parameter.name"= "phi_cna",
#                                                   "Parameter.value"=0.05,"Parameter.index"=index))
#       }
#     }
#   }
#   if (!("phi_clock" %in% parameters$Parameter.name)) {
#     parameters = rbind(parameters, data.frame("Parameter.name"= "phi_clock",
#                                               "Parameter.value"=0.05,"Parameter.index"=NA))
#   }
#
#   parameters_new = rbind(
#     parameters,
#     parameters_lengths,
#     parameters_coeff
#   )
#
#   return(list(mutations_new, parameters_new))
# }

# check_clinical_input = function(samples, therapies){
#
#   clinical_records = data.frame()
#   samples = samples %>% dplyr::arrange(Date)
#   for (c in 1:2){
#     new_entry = data.frame(
#       "Clinical.name"	= "Sample",
#       "Clinical.type"	= as.character(c),
#       "Clinical.value.start"= convert_real_date(samples$Date[c]),
#       "Clinical.value.end" = NA,
#       "Clinical.index"=NA,
#       "Clinical.original.name"=samples$Name[c]
#     )
#     clinical_records = rbind(clinical_records, new_entry)
#   }
#
#   therapies = therapies %>% dplyr::arrange(Start)
#
#   add_therapy_name = function(clinical_records, therapies, type){
#
#     therapy_names = therapies %>% dplyr::filter(Class == type) %>% dplyr::pull(Name) %>% unique()
#     n_therapy_types = length(therapy_names)
#
#     if (n_therapy_types > 0){
#
#       if (type == "Mutagenic") type="Therapy step"
#       if (type == "Driver responsive 1") type="Driver 1"
#       if (type == "Driver responsive 2") type="Driver 2"
#       if (type == "Chemotherapy inducing dormancy") type="Chemotherapy"
#
#       for (t in 1:n_therapy_types){
#         #therapy_df = data.frame()
#         sub_df = therapies %>% dplyr::filter(Name==therapy_names[t])
#         for (tt in 1:nrow(sub_df)){
#           new_entry = data.frame(
#             "Clinical.name"	= type,
#             "Clinical.type"	= t,
#             "Clinical.value.start"= convert_real_date(sub_df$Start[tt]),
#             "Clinical.value.end" = convert_real_date(sub_df$End[tt]),
#             "Clinical.index"=as.character(tt),
#             "Clinical.original.name"=sub_df$Name[tt]
#           )
#           clinical_records = rbind(clinical_records, new_entry)
#         }
#       }
#     }
#
#     return(clinical_records)
#   }
#
#   clinical_records = add_therapy_name(clinical_records, therapies, type="Mutagenic")
#   clinical_records = add_therapy_name(clinical_records, therapies, type="Driver responsive")
#   clinical_records = add_therapy_name(clinical_records, therapies, type="Chemotherapy inducing dormancy")
#
#   return(clinical_records)
# }

# Discussion on prior choosing

# Growth rate :
# The clonal expansion begins sometimes after the last therapy with clonal signatures.
# To be conservative and set a lower bound, we make the end of therapy coincide with the MRCA and compute omega_LB as :
# LOWER BOUND :
#   omega_LB = log(N) / delta_t
#   N = 10^6 (minimum)
#   delta_t = Sampling_2 - last_therapy_end
# UPPER BOUND :
#   omega_UB = log(N) / delta_t
#   N = 10^13 (maximum)
#   delta_t = Sampling_2 - small_delta (=1 month)
# We then compute the hyperparameters for the gamma distribution as:
# E = (lb_omega + ub_omega) / 2 # (expected value)
# V = 10 # (input variance)
# beta_omega = V / E
# alpha_omega = E*beta_omega

expected_growth_rate = function(N_min, N_max, last_th_end, Sampling_2, delta = 0.02){
  omega_LB = log(N_min) / (Sampling_2 - last_th_end)
  omega_UB = log(N_max) / (Sampling_2 - delta) # 0.02 ~ one week
  E = (omega_LB + omega_UB) / 2
  E
}

get_hyperparameters_growth_rate = function(E, V=10){
  beta_omega = V / E
  alpha_omega = E*beta_omega
  return(c('alpha_omega'=alpha_omega, 'beta_omega'=beta_omega))
}

# Mutation rate under GCV
# The lower bound mu_th[n]_LB is obtained assuming that the cell has cycled throughout therapy_n exposure:
# mu_th[n]_LB = muts_th[n] / (2*length_genome*omega*delta_therapy)
# The upper bound assuming that all therapy[n] mutations were accumulated during a single cell division:
# mu_th[n]_UB = muts_th[n] / (2*length_genome)

expected_mutation_rate_th = function(muts_th, omega, delta_therapy, length_genome=3e9){
  mu_th_LB = muts_th / (2*length_genome*omega*delta_therapy)
  mu_th_UB = muts_th / (2*length_genome)
  mu_th = mean(mu_th_LB, mu_th_UB)
  mu_th
}

get_hyperparameters_mutation_rate = function(mu_th, delta_therapy, var=.05){
  alpha_mu_th = mu_th*delta_therapy*var
  beta_mu_th = delta_therapy*var
  return(c('alpha'=alpha_mu_th, 'beta'=beta_mu_th))
}
