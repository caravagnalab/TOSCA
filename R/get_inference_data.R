get_inference_data_cna = function(x){

  data = list()

  data[['m_clock_primary']] = TOSCA:::get_mutation(x, type = "clock-like primary")
  data[['m_clock']] = TOSCA:::get_mutation(x, type = "clock-like relapse")
  data[['l_diploid']] = TOSCA:::get_length(x)
  data[['mu_clock']] = TOSCA:::get_mutation_rate(x, type = "clock")

  data[['n_cna']] = TOSCA:::get_n_cna(x)
  data[['m_alpha']] = TOSCA:::get_mutation(x, type = "alpha")
  data[['m_beta']] = TOSCA:::get_mutation(x, type = "beta")
  data[['l_CNA']] = TOSCA:::get_length(x, diploid = F)
  # data[['coeff_alpha']] = TOSCA:::get_coeff(x)[["alpha"]]
  # data[['coeff_beta']] = TOSCA:::get_coeff(x)[["beta"]]
  data[['coeff']] = TOSCA:::get_coeff(x)[["beta"]]

  data[['n_th_step']]= TOSCA:::get_n_therapy_cycles(x, class = 'Mutagenic')
  data[['n_th_step_type']]= TOSCA:::get_n_therapy_classes(x, class = "Mutagenic")
  data[['start_th_step']] = TOSCA:::get_start_therapy(x, class= "Mutagenic")
  data[['end_th_step']] =  TOSCA:::get_end_therapy(x, class= "Mutagenic")
  data[['type_th_step']]= TOSCA:::get_therapy_class_index(x, class= "Mutagenic")
  # data[['mu_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")
  data[['alpha_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")[["alpha"]]
  data[['beta_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")[["beta"]]
  data[['m_th_step']]= TOSCA:::get_mutation(x, type = "Mutagenic") # get_m_th(x, type = 'step')

  # data[['n_th_cauchy']]= TOSCA:::get_n_therapy_cycles(x, class = 'Mutagenic cauchy')
  # data[['n_th_cauchy_type']]= TOSCA:::get_n_therapy_classes(x, class = "Mutagenic cauchy")
  # data[['location_th_cauchy']]= TOSCA:::get_start_therapy(x, class= "Mutagenic cauchy")
  # data[['type_th_cauchy']]= TOSCA:::get_therapy_class_index(x, class= "Mutagenic cauchy")
  # data[['scales_th_cauchy']] = TOSCA:::get_cauchy_scales(x)
  # data[['m_th_cauchy']]=  TOSCA:::get_mutation(x, type = "Mutagenic", cauchy=T)

  data[['omega_alpha']] = TOSCA:::get_parameter(x, "omega_alpha")
  data[['omega_beta']] = TOSCA:::get_parameter(x, "omega_beta")
  data[['k_step']] = TOSCA:::get_parameter(x, "k_step")

  data[['Sample_1']] = TOSCA:::get_sample(x, sample=1)
  data[['Sample_2']] = TOSCA:::get_sample(x, sample=2)
  # data[['max_therapy']] = TOSCA:::get_max_th(x)

  data[['exponential_growth']] = TOSCA:::get_parameter(x, "exponential_growth")
  data[['N_min']] = TOSCA:::get_N(x, which="min")
  data[['N_max']] = TOSCA:::get_N(x, which="max")
  #data[['alpha_mrca']] = TOSCA:::get_parameter(x, "mrca_alpha")
  #data[['beta_mrca']] = TOSCA:::get_parameter(x, "mrca_beta")

  data[['phi_clock']] = TOSCA:::get_phi(x, "clock")
  data[['phi_th_step']] = TOSCA:::get_phi(x, "phi_th_step")
  # data[['phi_th_cauchy']] = TOSCA:::get_phi(x, "phi_th_cauchy")
  data[['phi_cna_alpha']] = TOSCA:::get_phi(x, "phi_cna")[["alpha"]]
  data[['phi_cna_beta']] = TOSCA:::get_phi(x, "phi_cna")[["beta"]]

  data
}

get_inference_data_driver = function(x){

  data = list()

  data[['m_clock_primary']] = TOSCA:::get_mutation(x, type = "clock-like primary")
  data[['m_clock']] = TOSCA:::get_mutation(x, type = "clock-like relapse")
  data[['l_diploid']] = TOSCA:::get_length(x)
  data[['mu_clock']] = TOSCA:::get_mutation_rate(x, type = "clock")

  data[['driver_type']] = TOSCA:::get_driver_type(x)
  data[['cycles_driver']] = TOSCA:::get_n_therapy_cycles(x, "Driver responsive")
  data[['driver_start']] = TOSCA:::get_start_therapy(x, class= "Driver responsive")
  data[['driver_end']] = TOSCA:::get_end_therapy(x, class= "Driver responsive")
  data[['m_driver']] = TOSCA:::get_mutation(x, type = "driver")
  data[['mu_driver_clock']] = TOSCA:::get_mutation_rate(x, type = "driver_clock")
  #data[['mu_driver']] = TOSCA:::get_mutation_rate(x, type = "driver")

  data[['n_th_step']]= TOSCA:::get_n_therapy_cycles(x, class = 'Mutagenic')
  data[['n_th_step_type']]= TOSCA:::get_n_therapy_classes(x, class = "Mutagenic")
  data[['start_th_step']] = TOSCA:::get_start_therapy(x, class= "Mutagenic")
  data[['end_th_step']] =  TOSCA:::get_end_therapy(x, class= "Mutagenic")
  data[['type_th_step']]= TOSCA:::get_therapy_class_index(x, class= "Mutagenic")
  data[['alpha_th_step']] =  TOSCA:::get_mutation_rate(x, type = "th_step")[["alpha"]]
  data[['beta_th_step']] =  TOSCA:::get_mutation_rate(x, type = "th_step")[["beta"]]
  # data[['mu_th_step']] =  TOSCA:::get_mutation_rate(x, type = "th_step")
  data[['m_th_step']]= TOSCA:::get_mutation(x, type = "Mutagenic")

  data[['n_th_cauchy']]= TOSCA:::get_n_therapy_cycles(x, class = 'Mutagenic cauchy')
  data[['n_th_cauchy_type']]= TOSCA:::get_n_therapy_classes(x, class = "Mutagenic cauchy")
  data[['location_th_cauchy']]= TOSCA:::get_start_therapy(x, class= "Mutagenic cauchy")
  data[['type_th_cauchy']]= TOSCA:::get_therapy_class_index(x, class= "Mutagenic cauchy")
  # data[['scales_th_cauchy']] = TOSCA:::get_cauchy_scales(x)
  data[['alpha_th_cauchy']] = get_cauchy_scales(x)[["alpha"]]
  data[['beta_th_cauchy']] = get_cauchy_scales(x)[["beta"]]
  data[['m_th_cauchy']]=  TOSCA:::get_mutation(x, type = "Mutagenic", cauchy=T)

  data[['omega_alpha']] = TOSCA:::get_parameter(x, "omega_alpha")
  data[['omega_beta']] = TOSCA:::get_parameter(x, "omega_beta")
  data[['mu_driver_alpha']] =  TOSCA:::get_mutation_rate(x, type = "driver")[["alpha"]]
  data[['mu_driver_beta']] =  TOSCA:::get_mutation_rate(x, type = "driver")[["beta"]]
  data[['k_step']] = TOSCA:::get_parameter(x, "k_step")

  data[['Sample_1']] = TOSCA:::get_sample(x, sample=1)
  data[['Sample_2']] = TOSCA:::get_sample(x, sample=2)
  data[['max_therapy']] = TOSCA:::get_max_th(x)

  data[['exponential_growth']] = TOSCA:::get_parameter(x, "exponential_growth")
  data[['N_min']] = TOSCA:::get_N(x, which="min")
  data[['N_max']] = TOSCA:::get_N(x, which="max")
  data[['mrca_alpha']] = TOSCA:::get_parameter(x, "mrca_alpha")
  data[['mrca_beta']] = TOSCA:::get_parameter(x, "mrca_beta")

  data[['phi_clock']] = TOSCA:::get_phi(x, "clock")
  data[['phi_driver']] = TOSCA:::get_phi(x, "driver")
  data[['phi_th_step']] = TOSCA:::get_phi(x, "phi_th_step")
  data[['phi_th_cauchy']] = TOSCA:::get_phi(x, "phi_th_cauchy")

  data
}

get_inference_data_dormancy = function(x){

  data = list()
  data[["n_cna"]] = TOSCA:::get_n_cna(x)
  data[["m_clock_primary"]] = TOSCA:::get_mutation(x, type = "clock-like primary")
  data[["m_clock"]] = TOSCA:::get_mutation(x, type = "clock-like relapse")
  data[["l_diploid"]] =  TOSCA:::get_length(x)
  data[["mu_clock"]] = TOSCA:::get_mutation_rate(x, type = "clock")

  data[["m_alpha"]] = TOSCA:::get_mutation(x, type = "alpha")
  data[["m_beta"]] = TOSCA:::get_mutation(x, type = "beta")
  data[["l_CNA"]] = TOSCA:::get_length(x, diploid = F)
  data[["coeff"]] = TOSCA:::get_coeff(x)[["beta"]]

  data[['n_th_step']]= TOSCA:::get_n_therapy_cycles(x, class = 'Mutagenic')
  data[['n_th_step_type']]= TOSCA:::get_n_therapy_classes(x, class = "Mutagenic")
  data[['start_th_step']] = TOSCA:::get_start_therapy(x, class= "Mutagenic")
  data[['end_th_step']] =  TOSCA:::get_end_therapy(x, class= "Mutagenic")
  data[['type_th_step']]= TOSCA:::get_therapy_class_index(x, class= "Mutagenic")
  # data[['mu_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")
  data[['alpha_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")[["alpha"]]
  data[['beta_th_step']] = TOSCA:::get_mutation_rate(x, type = "th_step")[["beta"]]
  data[['m_th_step']]= TOSCA:::get_mutation(x, type = "Mutagenic") # get_m_th(x, type = 'step')

  data[['omega_alpha']] = TOSCA:::get_parameter(x, "omega_alpha")
  data[['omega_beta']] = TOSCA:::get_parameter(x, "omega_beta")
  data[['k_step']] = TOSCA:::get_parameter(x, "k_step")

  data[['Sample_1']] = TOSCA:::get_sample(x, sample=1)
  data[['Sample_2']] = TOSCA:::get_sample(x, sample=2)
  data[['chemo_start']]= TOSCA:::get_start_therapy(x, class= "Chemotherapy inducing dormancy")
 data[['chemo_end']]= TOSCA:::get_start_therapy(x, class= "Chemotherapy inducing dormancy")

  data[['exponential_growth']] = TOSCA:::get_parameter(x, "exponential_growth")
  data[['N_min']] = TOSCA:::get_N(x, which="min")
  data[['N_max']] = TOSCA:::get_N(x, which="max")


  data[['phi_clock']] = TOSCA:::get_phi(x, "clock")
  data[['phi_th_step']] = TOSCA:::get_phi(x, "phi_th_step")
  # data[['phi_th_cauchy']] = TOSCA:::get_phi(x, "phi_th_cauchy")
  data[['phi_cna_alpha']] = TOSCA:::get_phi(x, "phi_cna")[["alpha"]]
  data[['phi_cna_beta']] = TOSCA:::get_phi(x, "phi_cna")[["beta"]]

  data
}

get_inference_data = function(x, model, dormancy = F){

  if (model == "CNA" & !dormancy) data = get_inference_data_cna(x)
  if (model == "CNA" & dormancy) data = get_inference_data_dormancy(x)
  if (model == "Driver") data = get_inference_data_driver(x)

  data

}
