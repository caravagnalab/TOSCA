#### Getters for input data

## Get mutations
# get_m_clock = function(x, type = "relapse"){
#   as.integer(x$mutations %>% filter(Mutation.name=='m_clock', Mutation.type == type) %>% pull(Mutation.value))
# }
# get_m_driver = function(x){
#   as.integer(x$mutations %>% filter(Mutation.name=='m_driver') %>% pull(Mutation.value))
# }
# get_m_CNA = function(x, type = 'alpha'){
#   as.integer(x$mutations %>% filter(Mutation.type==type) %>% pull(Mutation.value))
# }
# get_m_th = function(x, type = 'step'){
#   as.integer(x$mutations %>% filter(Mutation.type==type) %>% pull(Mutation.value))
# }

get_n_cna = function(x){
  (x$mutations %>% dplyr::filter(Mutation.name=='m_cna') %>% nrow())/2
}

get_mutation = function(x, name, type=NA, index=NA, source=NA, coeff=NA){
  m = x$mutations %>% dplyr::filter(Mutation.name==name)
  if (!is.na(type)) m = m %>% dplyr::filter(Mutation.type == type)
  if (!is.na(index)) m = m %>% dplyr::filter(Mutation.index == index)
  if (!is.na(source)) m = m %>% dplyr::filter(Mutation.source == source)
  if (!is.na(coeff)) m = m %>% dplyr::filter(Mutation.coeff == coeff)
  m %>% dplyr::pull(Mutation.value) %>% as.integer()
}

## Get Clinical Data
get_sample = function(x, sample='1'){
  x$clinical_records %>% dplyr::filter(Clinical.name=='Sample', Clinical.type==sample) %>% dplyr::pull(Clinical.value.start)
}

get_therapy_driver = function(x){
  start = x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% dplyr::pull(Clinical.value.start)
  end = x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% dplyr::pull(Clinical.value.end)
  list('start'=start,'end'=end)
}

get_n_th = function(x, name = 'Therapy step'){
  x$clinical_records %>% dplyr::filter(Clinical.name==name) %>% nrow()
}

get_n_th_type = function(x, name = 'Therapy step'){
  n_th_type = length(x$clinical_records %>% filter(Clinical.name== name) %>% dplyr::pull(Clinical.type) %>% unique() )
  if (is.null(n_th_type)){n_th_type=0}
  n_th_type
}

start_th = function(x, type='Therapy step'){
  x$clinical_records %>% dplyr::filter(Clinical.name==type) %>% dplyr::pull(Clinical.value.start)
}

end_th_step = function(x){
  x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy step') %>% dplyr::pull(Clinical.value.end)
}

start_chemo = function(x){
  chemo_names = x$clinical_records$Clinical.name[grepl("Chemotherapy", x$clinical_records$Clinical.name)]
  start = x$clinical_records %>% dplyr::filter(Clinical.name %in% chemo_names) %>% dplyr::pull(Clinical.value.start) #%>% arrange()
  start[1]
}

end_chemo = function(x){
  chemo_names = x$clinical_records$Clinical.name[grepl("Chemotherapy", x$clinical_records$Clinical.name)]
  end = x$clinical_records %>% dplyr::filter(Clinical.name %in% chemo_names) %>% dplyr::pull(Clinical.value.end) #%>% arrange()
  end[length(end)]
}

get_type_th_step = function(x, name='Therapy step'){
  as.integer(x$clinical_records %>% dplyr::filter(Clinical.name==name) %>% dplyr::pull(Clinical.type) )#%>% unique()
}


## Get parameters
get_length = function(x, model=NA){
  if (is.na(model)){
    return(x$parameters %>% dplyr::filter(Parameter.name=="l_diploid") %>% dplyr::pull(Parameter.value) %>% unique() %>% as.double())
  }
  if (model=="CNA"){
    return(x$parameters %>% dplyr::filter(Parameter.name=="l_CNA") %>% dplyr::pull(Parameter.value) %>% unique() %>% as.double())
  }
  if (model=="WGD"){
    return(x$parameters %>% dplyr::filter(Parameter.name=="l_tetraploid") %>% dplyr::pull(Parameter.value) %>% as.double())
  }
  # if (karyotype=="2:0" & model=="WGD"){
  #   return(x$parameters %>% filter(Parameter.name=="l_cnloh") %>% pull(Parameter.value))
  # }
}

get_coeff_CNA = function(x){
  coeff_beta = x$parameters %>% dplyr::filter(Parameter.name=="coeff_CNA") %>% dplyr::pull(Parameter.value)
  coeff_alpha = c()
  for (c in coeff_beta){
    if (c=="2"){coeff_alpha=c(coeff_alpha, 1)}
    if (c=="4"){coeff_alpha=c(coeff_alpha, 2)}
  }
  list('alpha'=coeff_alpha, 'beta'=as.integer(coeff_beta))
}

# get_mu_clock = function(x){
#   x$parameters %>% filter(Parameter.name=="mu_clock") %>% pull(Parameter.value)
# }
#
# get_mu_driver_clock = function(x){
#   x$parameters %>% filter(Parameter.name=="mu_clock_driver") %>% pull(Parameter.value)
# }

get_mutation_rate = function(x, type, index=NA){
  # type = clock, clock_driver, driver, th
  # index = 1...n_th_step_type
  mu = x$parameters %>% dplyr::filter(Parameter.name==paste0("mu_", type))
  if (!is.na(index)) mu = mu %>% dplyr::filter(Parameter.index==index)
  mu %>% dplyr::pull(Parameter.value) %>% as.double()
}

get_driver_type = function(x){
  # (x$mutations %>% dplyr::filter(Mutation.name=='m_driver') %>% dplyr::pull(Mutation.type) %>% as.integer()) -1
  2
}

get_cycles_drivers = function(x){
  x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% nrow()
}

get_prior_hyperparameters = function(x, name){
  # name = mu_driver*, mu_th_step*, scale_th_cauchy*, omega*, mrca*, s*
  if (name %in% x$parameters$Parameter.name){
    x$parameters %>% dplyr::filter(Parameter.name==name) %>% dplyr::pull(Parameter.value) %>% as.double()
  }else{
    alpha= x$parameters %>% dplyr::filter(Parameter.name==paste0(name,'_alpha')) %>% dplyr::pull(Parameter.value) %>% as.double()
    beta= x$parameters %>% dplyr::filter(Parameter.name==paste0(name,'_beta')) %>% dplyr::pull(Parameter.value) %>% as.double()
    list('alpha'=alpha, 'beta'=beta)
  }
}

get_k_step = function(x){
  x$parameters %>% dplyr::filter(Parameter.name=='k_step') %>% dplyr::pull(Parameter.value) %>% as.double()
}

get_max_th = function(x){

  max_index = x$clinical_records %>% dplyr::filter(Clinical.name %in% c("Therapy step", "Therapy cauchy", "Therapy driver"))

  if (nrow(max_index) > 0){
    max_index = max_index %>% dplyr::pull(Clinical.type) %>% as.integer() %>% max() %>% as.character()
    return(min(x$clinical_records %>% dplyr::filter(Clinical.name!="Sample", Clinical.type == max_index) %>% dplyr::pull(Clinical.value.start)))
  }else{
    return(x$clinical_records %>% dplyr::filter(Clinical.name=="Sample", Clinical.type == 1) %>% dplyr::pull(Clinical.value.start))
  }
  #min(dates)
}

get_exponential_growth = function(x){
  x$parameters %>% dplyr::filter(Parameter.name=='exponential_growth') %>% dplyr::pull(Parameter.value) %>% as.integer()
}

get_N_min = function(x){
  N_min = x$parameters %>% dplyr::filter(Parameter.name=='N_min') %>% dplyr::pull(Parameter.value) %>% as.double()
  if (length(N_min) < 2){return(c(N_min, N_min))}else{N_min}
}

get_N_max= function(x){
  N_max = x$parameters %>% dplyr::filter(Parameter.name=='N_max') %>% dplyr::pull(Parameter.value) %>% as.double()
  if (length(N_max) < 2){return(c(N_max, N_max))}else{N_max}
}

get_cauchy_scales = function(x){
  x$parameters %>% dplyr::filter(Parameter.name == "cauchy_scales") %>% dplyr::pull(Parameter.value) %>% as.double()
}

get_phi = function(x, name){
  # name : th_cauchy, th_step, cna, clock, driver
  x$parameters %>% dplyr::filter(Parameter.name == paste0("phi_",name)) %>% dplyr::pull(Parameter.value) %>% as.double()
}

get_extra_therapy = function(x){
  if (nrow(x$clinical_records %>% dplyr::filter(!(Clinical.name %in% c("Sample", "Chemotherapy")))) > 0){
    1
  }else{
    0
  }
}

get_wgd = function(x){
  if ("WGD" %in% x$mutations$Mutation.source) 1
  else 0
}

get_n_clock = function(x){
  get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA) %>% length()
}

get_n_step = function(x){
  get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) %>% length()
}

get_inference_data_driver = function(x){

  data = list()

  data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
  data[['m_clock']] = get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA)
  data[['l_diploid']] = get_length(x, model=NA)
  data[['mu_clock']] = get_mutation_rate(x, type = "clock")

  data[['driver_type']] = get_driver_type(x)
  data[['cycles_driver']] = get_cycles_drivers(x)
  data[['driver_start']] = get_therapy_driver(x)[["start"]]
  data[['driver_end']] = get_therapy_driver(x)[["end"]]
  data[['m_driver']] = get_mutation(x, name="m_driver", type=NA, index=NA, source=NA)
  data[['mu_driver_clock']] = get_mutation_rate(x, type = "driver_clock")
  data[['mu_driver']] = get_mutation_rate(x, type = "driver")

  data[['n_th_step']]= get_n_th(x, name = 'Therapy step')
  data[['n_th_step_type']]= get_n_th_type(x, name = 'Therapy step')
  data[['start_th_step']] = start_th(x, type='Therapy step')
  data[['end_th_step']] =  end_th_step(x)
  data[['type_th_step']]= get_type_th_step(x, name='Therapy step')
  data[['mu_th_step']] = get_mutation_rate(x, type = "th_step")
  data[['m_th_step']]= get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) # get_m_th(x, type = 'step')

  data[['n_th_cauchy']]= get_n_th(x, name = 'Therapy cauchy')
  data[['n_th_cauchy_type']]= get_n_th_type(x, name = 'Therapy cauchy')
  data[['location_th_cauchy']]= start_th(x, type='Therapy cauchy')
  data[['type_th_cauchy']]= get_type_th_step(x, name='Therapy cauchy')
  data[['scales_th_cauchy']] = get_cauchy_scales(x)
  data[['m_th_cauchy']]= get_mutation(x, name="m_th", type=NA, index=NA, source="cauchy", coeff=NA) #get_m_th(x, type = 'Cauchy')



  data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
  data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
  data[['k_step']] = get_k_step(x)

  data[['Sample_1']] = get_sample(x, sample='1')
  data[['Sample_2']] = get_sample(x, sample='2')
  data[['max_therapy']] = get_max_th(x)

  data[['exponential_growth']] = get_exponential_growth(x)
  data[['N_min']] = get_N_min(x)
  data[['N_max']] = get_N_max(x)
  data[['alpha_mrca']] = get_prior_hyperparameters(x, name='mrca')[["alpha"]]
  data[['beta_mrca']] = get_prior_hyperparameters(x, name='mrca')[["beta"]]
  data[['phi_clock']] = get_phi(x, "clock")
  data[['phi_driver']] = get_phi(x, "driver")
  data[['phi_th_step']] = get_phi(x, "th_step")
  data[['phi_th_cauchy']] = get_phi(x, "th_cauchy")

  data
}

get_inference_data_cna = function(x){

  data = list()

  data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
  data[['m_clock']] = get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA)
  data[['l_diploid']] = get_length(x, model=NA)
  data[['mu_clock']] = get_mutation_rate(x, type = "clock")

  data[['n_cna']] = get_n_cna(x)
  data[['m_alpha']] = get_mutation(x, name="m_cna", type="alpha", index=NA, source=NA)
  data[['m_beta']] = get_mutation(x, name="m_cna", type="beta", index=NA, source=NA)
  data[['l_CNA']] = get_length(x, model="CNA")
  data[['coeff_alpha']] = get_coeff_CNA(x)[["alpha"]]
  data[['coeff_beta']] = get_coeff_CNA(x)[["beta"]]

  data[['n_th_step']]= get_n_th(x, name = 'Therapy step')
  data[['n_th_step_type']]= get_n_th_type(x, name = 'Therapy step')
  data[['start_th_step']] = start_th(x, type='Therapy step')
  data[['end_th_step']] =  end_th_step(x)
  data[['type_th_step']]= get_type_th_step(x, name='Therapy step')
  data[['mu_th_step']] = get_mutation_rate(x, type = "th_step")
  data[['m_th_step']]= get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) # get_m_th(x, type = 'step')

  data[['n_th_cauchy']]= get_n_th(x, name = 'Therapy cauchy')
  data[['n_th_cauchy_type']]= get_n_th_type(x, name = 'Therapy cauchy')
  data[['location_th_cauchy']]= start_th(x, type='Therapy cauchy')
  data[['type_th_cauchy']]= get_type_th_step(x, name='Therapy cauchy')
  data[['scales_th_cauchy']] = get_cauchy_scales(x)
  data[['m_th_cauchy']]= get_mutation(x, name="m_th", type=NA, index=NA, source="cauchy", coeff=NA) #get_m_th(x, type = 'Cauchy')



  data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
  data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
  data[['k_step']] = get_k_step(x)

  data[['Sample_1']] = get_sample(x, sample='1')
  data[['Sample_2']] = get_sample(x, sample='2')
  data[['max_therapy']] = get_max_th(x)

  data[['exponential_growth']] = get_exponential_growth(x)
  data[['N_min']] = get_N_min(x)
  data[['N_max']] = get_N_max(x)
  data[['alpha_mrca']] = get_prior_hyperparameters(x, name='mrca')[["alpha"]]
  data[['beta_mrca']] = get_prior_hyperparameters(x, name='mrca')[["beta"]]

  data[['phi_clock']] = get_phi(x, "clock")
  data[['phi_th_step']] = get_phi(x, "th_step")
  data[['phi_th_cauchy']] = get_phi(x, "th_cauchy")
  data[['phi_cna']] = get_phi(x, "cna")

  data
}

get_inference_data_dormancy = function(x){

  data = list()

  data[['wgd']] = get_wgd(x)
  data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
  data[["n_clock"]] = get_n_clock(x)
  data[['m_clock']] = get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA)
  data[['l_diploid']] = get_length(x, model=NA)
  data[['mu_clock']] = get_mutation_rate(x, type = "clock")

  data[['m_alpha']] = get_mutation(x, name="m_cna", type="alpha", index=NA, source=NA)
  data[['m_beta']] = get_mutation(x, name="m_cna", type="beta", index=NA, source=NA)
  data[['l_CNA']] = get_length(x, model="CNA")
  data[['coeff_alpha']] = get_coeff_CNA(x)[["alpha"]]
  data[['coeff_beta']] = get_coeff_CNA(x)[["beta"]]

  data[['extra_therapy']] = get_extra_therapy(x)
  data[["n_steps"]] = get_n_step(x)
  data[['start_th_step']] = start_th(x, type='Therapy step')
  data[['end_th_step']] =  end_th_step(x)
  data[['mu_th_step']] = get_mutation_rate(x, type = "th_step")
  data[['m_th_step']]= get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) # get_m_th(x, type = 'step')

  data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
  data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
  data[['k_step']] = get_k_step(x)

  data[['Sample_1']] = get_sample(x, sample='1')
  data[['Sample_2']] = get_sample(x, sample='2')
  data[['chemo_start']] = start_chemo(x)
  data[['chemo_end']] = end_chemo(x)

  data[['exponential_growth']] = get_exponential_growth(x)
  data[['N_min']] = get_N_min(x)
  data[['N_max']] = get_N_max(x)
  data[['alpha_mrca']] = get_prior_hyperparameters(x, name='mrca')[["alpha"]]
  data[['beta_mrca']] = get_prior_hyperparameters(x, name='mrca')[["beta"]]

  data[['phi_clock']] = get_phi(x, "clock")
  data[['phi_th_step']] = get_phi(x, "th_step")
  data[['phi_cna']] = get_phi(x, "cna")

  data
}

get_inference_data = function(x, model, dormancy = F){

  if (model == "CNA" & !dormancy) data = get_inference_data_cna(x)
  if (model == "CNA" & dormancy) data = get_inference_data_dormancy(x)
  if (model == "Driver") data = get_inference_data_driver(x)

  data

}

## Get data for inference
# get_inference_data = function(x, model, dormancy){
#
#   data = list()
#   data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
#   data[['l_diploid']] = get_length(x, model=NA)
#   data[['m_clock']] = get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA)
#   data[['mu_clock']] = get_mu_clock(x)
#
#   if (model == 'Driver'){
#   # mutations associated to driver
#       data[['driver_type']] = get_driver_type(x)
#       data[['cycles_driver']] = get_cycles_drivers(x)
#       data[['driver_start']] = get_therapy_driver(x)[["start"]]
#       data[['driver_end']] = get_therapy_driver(x)[["end"]]
#       data[['m_driver']] = get_mutation(x, name="m_driver", type=NA, index=NA, source=NA)
#       data[['mu_driver_clock']] = get_mu_driver_clock(x)
#       data[['mu_driver']] = get_prior_hyperparameters(x, name='mu_driver')
#
#   }
#
#   if (model == 'CNA'){
#     # mutations on CNA
#     data[['n_cna']] = get_n_cna(x)
#     data[['m_alpha']] = get_mutation(x, name="m_cna", type="alpha", index=NA, source=NA)
#     data[['m_beta']] = get_mutation(x, name="m_cna", type="beta", index=NA, source=NA)
#     data[['l_CNA']] = get_length(x, model="CNA")
#     data[['coeff_alpha']] = get_coeff_CNA(x)[["alpha"]]
#     data[['coeff_beta']] = get_coeff_CNA(x)[["beta"]]
#   }
#
#   if (model %in% c('Driver', 'CNA')){
#     data[['n_th_step']]= get_n_th(x, name = 'Therapy step')
#     data[['n_th_step_type']]= get_n_th_type(x, name = 'Therapy step')
#     data[['start_th_step']] = start_th(x, type='Therapy step')
#     data[['end_th_step']] =  end_th_step(x)
#     data[['type_th_step']]= get_type_th_step(x, name='Therapy step')
#     data[['alpha_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["alpha"]]
#     data[['beta_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["beta"]]
#   }
#
#   if (model %in% c('Driver', 'CNA')){
#       data[['m_th_step']]= get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) # get_m_th(x, type = 'step')
#   }
#
#   if (model %in% c('Driver', 'CNA')){
#     # mutations associated to cauchy
#     data[['n_th_cauchy']]= get_n_th(x, name = 'Therapy cauchy')
#     data[['n_th_cauchy_type']]= get_n_th_type(x, name = 'Therapy cauchy')
#     data[['location_th_cauchy']]= start_th(x, type='Therapy cauchy')
#     data[['type_th_cauchy']]= get_type_th_step(x, name='Therapy cauchy')
#     data[['alpha_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["alpha"]]
#     data[['beta_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["beta"]]
#     if (model %in% c('Driver', 'CNA')){
#       data[['m_th_cauchy']]= get_mutation(x, name="m_th", type=NA, index=NA, source="cauchy", coeff=NA) #get_m_th(x, type = 'Cauchy')
#     }
#   }
#
#   if (model == 'WGD'){
#     data[['m_alpha']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="clock", coeff="4")
#     data[['m_beta']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="clock", coeff="4")
#     data[['l_wgd']]= get_length(x, model="WGD")
#
#     # data[['m_alpha_tetraploid_step']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="step", coeff="4")
#     # data[['m_beta_tetraploid_step']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="step", coeff="4")
#     # data[['m_alpha_tetraploid_cauchy']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="cauchy", coeff="4")
#     # data[['m_beta_tetraploid_cauchy']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="cauchy", coeff="4")
#     # data[['m_alpha_tetraploid_clock']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="clock", coeff="4")
#     # data[['m_beta_tetraploid_clock']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="clock", coeff="4")
#
#     # data[['alpha_cnloh_step']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="step", coeff="2")
#     # data[['beta_cnloh_step']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="step", coeff="2")
#     # data[['alpha_cnloh_cauchy']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="cauchy", coeff="2")
#     # data[['beta_cnloh_cauchy']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="cauchy", coeff="2")
#     # data[['alpha_cnloh_clock']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="clock", coeff="2")
#     # data[['beta_cnloh_clock']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="clock", coeff="2")
#
#     #data[['l_tetraploid']]= get_length(x, model="WGD")
#     # data[['l_cnloh']]= get_length(x, karyotype="2:0", model="WGD")
#
#   }
#
#   # other parameters
#   if (fixed_omega==T){
#     data[['omega']] = get_prior_hyperparameters(x, name='omega')
#   }else{
#     data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
#     data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
#   }
#
#   data[['k_step']] = get_k_step(x)
#   data[['Sample_1']] = get_sample(x, sample='1')
#   data[['Sample_2']] = get_sample(x, sample='2')
#
#   if (model %in% c('Driver', 'CNA')){
#     data[['max_therapy']] = get_max_th(x)
#   }
#
#   if (model == "WGD" & dormancy){
#     data[['chemo_start']] = start_chemo(x)
#     data[['chemo_end']] = end_chemo(x)
#   }
#
#   data[['exponential_growth']] = get_exponential_growth(x)
#   data[['N_min']] = get_N_min(x)
#   data[['N_max']] = get_N_max(x)
#
#   if (model %in% c('CNA', 'Driver')){
#     data[['alpha_mrca']] = get_prior_hyperparameters(x, name='mrca')[["alpha"]]
#     data[['beta_mrca']] = get_prior_hyperparameters(x, name='mrca')[["beta"]]
#   }
#
#   data
# }

get_model <- function(model_name='Driver', dormancy=F) {

  if (model_name=='CNA' & dormancy) model_name = "CNA_dormancy"

  all_paths <- list(
    "Driver" = "Driver.stan",
    "CNA" = "CNA.stan",
    "CNA_dormancy" = "CNA_with_dormancy.stan"
  )

  if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")

  model_path <- system.file("stan", all_paths[[model_name]], package = "TOSCA", mustWork = T)
  tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
  model
}

#### Getters for inferred data
get_inferred_parameters = function(x){
  # posterior::as_draws_df(x$tosca_fit$draws())
  x$tosca_fit$posterior
}

get_inferred_times_colors = function(){
  c("t_eca"="#ae4532ff", "t_driver"="#d48b3eff", "t_mrca"="#7876adff",
    "t_cna"="#dc5895ff", "t_wgd"="#438455ff")
}

get_original_mutation_name = function(x, name, index){
  if (name == "m_clock_primary") original_name = x$mutations %>% filter(Mutation.name=="m_clock", Mutation.type=="primary") %>% pull(Mutation.original.name)
  if (name == "m_clock") original_name = x$mutations %>% filter(Mutation.name=="m_clock", Mutation.type=="relapse") %>% pull(Mutation.original.name)
  if (name == "m_alpha") original_name = x$mutations %>% filter(Mutation.name=="m_cna", Mutation.type=="alpha", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_beta") original_name = x$mutations %>% filter(Mutation.name=="m_cna", Mutation.type=="beta", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_th_step") original_name = x$mutations %>% filter(Mutation.name=="m_th", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_driver") original_name = x$mutations %>% filter(Mutation.name=="m_driver") %>% pull(Mutation.original.name)
  return(original_name)
}

## Get inferred times

## Get inferred parameters
