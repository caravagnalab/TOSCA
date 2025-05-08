check_input_mutations = function(mutations){
  condition1= class(mutations) == "data.frame"
  condition2= (colnames(mutations) == c('Mutation.type','Number.of.mutations')) %>% sum() == ncol(mutations)
  condition3= (c('m_clock', 'm_alpha', 'm_beta', 'm_th_1') %in% mutations$Mutation.type) %>% sum() == nrow(mutations)
  condition4= class(mutations$Number.of.mutations) == 'numeric'

  if(!condition1){
    stop("Mutations input must be a dataframe")
  }
  if(!condition2){
    stop("Missing columns: Mutation.type,Number.of.mutations")
  }
  if(!condition3){
    stop("Missing mutation types: 'm_clock', 'm_alpha', 'm_beta', 'm_th_1'")
  }
  if(!condition4){
    stop("Number.of.mutations must be numerical")
  }

}

check_input_clinical = function(clinical_records){
  condition1 = class(clinical_records) == "data.frame"
  condition2 = (colnames(clinical_data) == c('Timepoint','Start','End')) %>% sum() == ncol(clinical_data)
  condition3 = (c('Sample 1', 'Sample 2', 'Therapy 1') %in% clinical_data$Timepoint) %>% sum() == 3
  condition4 =
    (sapply(clinical_data$Start, function(x){
    year=strsplit(x, '-')[[1]][1]
    month=strsplit(x, '-')[[1]][2]
    day=strsplit(x, '-')[[1]][3]
    c1 = nchar(year) == 4
    c2 = as.integer(month) < 12
    c1 & c2
  }) %>% sum()) == ncol(clinical_data)

  if(!condition1){
    stop("Clinical records input must be a dataframe")
  }
  if(!condition2){
    stop("Missing columns: 'Timepoint','Start','End'")
  }
  if(!condition3){
    stop("Missing events: 'Sample 1', 'Sample 2', 'Therapy 1'")
  }
  if(!condition4){
    stop("Erroneous date format: YYYY-M-D")
  }
}

check_parameters = function(mutations, clinical_records, parameters, delta_omega = 0.02, growth_rate_variance=10){

  if (!('N_max' %in% parameters$Param.name)){parameters = rbind(parameters, data.frame('Param.name'='N_max', 'Value'=10**13))}
  if (!('N_min' %in% parameters$Param.name)){parameters = rbind(parameters, data.frame('Param.name'='N_min', 'Value'=10**6))}
  if (!('diploid_length' %in% parameters$Param.name)){parameters = rbind(parameters, data.frame('Param.name'='diploid_length', 'Value'=3e9))}
  if (!('k' %in% parameters$Param.name)){parameters = rbind(parameters, data.frame('Param.name'='k', 'Value'=1e4))}
  if (!('mu' %in% parameters$Param.name)){
    stop("Missing required parameter: mutation rate")
  }
  if (!('CNA_length' %in% parameters$Param.name)){
    stop("Missing required parameter: length of the CNA region")
  }
  # Compute best alpha/beta hyperparameters for growth rate
  nmin = parameters %>% filter(Param.name=='N_min') %>% pull(Value)
  nmax = parameters %>% filter(Param.name=='N_max') %>% pull(Value)
  last_th_end = clinical_records %>% filter(grepl('Therapy', Timepoint)) %>% arrange(End) %>% pull(End)
  end_last_therapy = convert_real_date(last_th_end[length(last_th_end)])
  last_sample = convert_real_date(clinical_records %>% filter(Timepoint=='Sample 2') %>% pull(Start))
  exp_omega = expected_growth_rate(nmin, nmax, end_last_therapy, last_sample, delta = delta_omega)
  if (!('omega_alpha' %in% parameters$Param.name)){
    alpha_beta = get_hyperparameters_growth_rate(exp_omega, V=growth_rate_variance)
    parameters = rbind(parameters,
                       data.frame('Param.name'=c('omega_alpha','omega_beta'), 'Value'=c(alpha_beta['alpha_omega'], alpha_beta['beta_omega'])))
  }
  # Compute best alpha/beta hyperparameters for mutation rate under therapy
  n_therapies = clinical_records %>% filter(grepl('Therapy', Timepoint)) %>% nrow()
  for (i in 1:n_therapies){
    par_name = paste0('alpha_mu_th_',i)
    if (!(par_name %in% parameters$Param.name)){

      th_start = convert_real_date(clinical_records %>% filter(Timepoint == paste0('Therapy ',i)) %>% pull(Start))
      th_end = convert_real_date(clinical_records %>% filter(Timepoint == paste0('Therapy ',i)) %>% pull(End))
      delta_therapy = th_end - th_start
      l_diploid = parameters %>% filter(Param.name == 'diploid_length') %>% pull(Value)
      muts_th = mutations %>% filter(Mutation.type == paste0('m_th_',i)) %>% pull(Number.of.mutations)
      expected_mu = expected_mutation_rate_th(muts_th, exp_omega, delta_therapy, l_diploid)
      alpha_beta = get_hyperparameters_mutation_rate(expected_mu, delta_therapy, var=.05)
      parameters = rbind(parameters,
                         data.frame('Param.name'=c(paste0('alpha_mu_th_',i),paste0('beta_mu_th_',i)),
                                    'Value'=c(alpha_beta['alpha'], alpha_beta['beta'])))
    }
  }

}

convert_date_real = function(x) {
  y = 2000 + x %>% floor()
  py = (x - (x %>% floor())) * 365
  m = (py / 30) %>% floor
  d = (py - (m * 30)) %>% floor
  date_string = paste(y, m +1, d + 1, sep = '-')

  date_string = ifelse (d>30, paste(y, m +1, 30, sep = '-'), date_string)
  date_string = ifelse ((m==1 & d>27), paste(y, m +1, 28, sep = '-'), date_string)
  date_string = ifelse (m>=12, paste(y+1, 1, 1, sep = '-'), date_string)

  return(date_string)
}

convert_real_date = function(date = NULL, ref_year = 2000) {
  ref_month = 1
  ref_day = 1

  year = as.integer(strsplit(date, '-')[[1]][1])
  month = as.integer(strsplit(date, '-')[[1]][2])
  day = as.integer(strsplit(date, '-')[[1]][3])

  return((year - ref_year) + (month / 12 - ref_month / 12) + (day / 365 - ref_day / 365))
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
