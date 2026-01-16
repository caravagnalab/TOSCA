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

check_valid_karyotype = function(karyo){

  k1= strsplit(karyo, ":")[[1]][1]
  k2= strsplit(karyo, ":")[[1]][2]

  class(karyo) == "character" & as.integer(k1) >= as.integer(k2) & as.integer(k1) <= 2 & as.integer(k2) <= 2
}

check_valid_date = function(str_date){
  y = strsplit(str_date, "-")[[1]][1]
  condition_y = nchar(y) == 4 & as.integer(y) >= 1900
  m = strsplit(str_date, "-")[[1]][2]
  condition_m = as.integer(m) <= 12
  d = strsplit(str_date, "-")[[1]][3]
  condition_d = as.integer(d) <= 31

  condition_y & condition_m & condition_d
}

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

#' Converts real into date
#'
#' @param x TOSCA obj
#' @param ref_year year of reference
#' @param date real number to transform to date
#'
#' @return Converts date into real number
#'
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

#' Converts date into real
#'
#' @param x TOSCA obj
#' @param date date to transform
#' @param ref_year year of reference
#'
#' @return Converts real number into date
#'
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

get_first_sample_name = function(x){
  samples = x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Name)
  if (nrow(x$Input$Samples) == 3) samples[2] else samples[1]
}

get_first_sample_date = function(x){
  samples = x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Date) %>% as.character()
  if (nrow(x$Input$Samples) == 3) samples[2] else samples[1]
}

get_second_sample_name = function(x){
  samples = x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Name)
  if (nrow(x$Input$Samples) == 3) samples[3] else samples[2]
}

get_second_sample_date = function(x){
  samples = x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Date) %>% as.character()
  if (nrow(x$Input$Samples) == 3) samples[3] else samples[2]
}

get_therapies_ordered = function(x){
  x$Input$Therapies %>% dplyr::arrange(as.Date(Start)) %>%
    dplyr::mutate(event = paste0(Name, " (", as.Date(End) - as.Date(Start), " days)")) %>% pull(event) %>% paste(collapse=" --> ")
}

get_clinical_history_length = function(x){
  days = as.Date(TOSCA:::get_second_sample_date(x)) - as.Date(TOSCA:::get_first_sample_date(x))
  years = as.integer((days / 365))
  months = as.integer((days - as.integer((days / 365))*365)/30)
  days =  as.integer(days - as.integer((days / 365))*365 - months*30)
  paste0("~ ",years, " years, ", months, " months and ", days, " days")
}

get_number_of_cna = function(x){
  cnas = x$Input$Mutations %>% dplyr::filter(Type == "alpha") %>% nrow()
  cnas
}

get_possible_dormancy = function(x){
  if ((x$Input$Therapies %>% dplyr::filter(Class == "Chemotherapy inducing dormancy") %>% nrow()) >= 1) "with" else "without"
}

get_driver = function(x){
  if ((x$Input$Mutations %>% dplyr::filter(Type == "driver") %>% nrow()) == 1) "with" else "without"
}

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
