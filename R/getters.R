#### Getters for input data

## Get mutations
get_m_clock = function(x){
  return(x$mutations %>% filter(Mutation.type=='m_clock') %>% pull(Number.of.mutations))
}
get_m_chemo = function(x){
  return(x$mutations %>% filter(Mutation.type=='m_chemo') %>% pull(Number.of.mutations))
}
get_m_mutag = function(x){
  return(x$mutations %>% filter(Mutation.type=='m_mutag') %>% pull(Number.of.mutations))
}
get_m_alpha = function(x){
  return(x$mutations %>% filter(Mutation.type=='m_alpha') %>% pull(Number.of.mutations))
}
get_m_beta = function(x){
  return(x$mutations %>% filter(Mutation.type=='m_beta') %>% pull(Number.of.mutations))
}
get_m_th = function(x, th='1'){
  return(x$mutations %>% filter(Mutation.type==paste0('m_th_',th)) %>% pull(Number.of.mutations))
}

## Get Clinical Data
get_therapy_start = function(x, th='1'){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint==paste0('Therapy_',th)) %>% pull(Start)))
}
get_therapy_end = function(x, th='1'){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint==paste0('Therapy_',th)) %>% pull(End)))
}
get_chemo_start = function(x){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint=='Chemotherapy') %>% pull(Start)))
}
get_chemo_end = function(x){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint=='Chemotherapy') %>% pull(End)))
}
get_sampling_time = function(x, time='1'){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint==paste0('Sample_',time)) %>% pull(Start)))
}

## Get parameters
get_mutation_rate = function(x){
  return(x$parameters %>% filter(Param.name=='mu') %>% pull(Value))
}
get_N_max = function(x){
  return(x$parameters %>% filter(Param.name=='N_max') %>% pull(Value))
}
get_N_min = function(x){
  return(x$parameters %>% filter(Param.name=='N_min') %>% pull(Value))
}
get_k = function(x){
  return(x$parameters %>% filter(Param.name=='k') %>% pull(Value))
}
get_k_sm = function(x){
  return(x$parameters %>% filter(Param.name=='k_sm') %>% pull(Value))
}
get_omega_alpha = function(x){
  return(x$parameters %>% filter(Param.name=='omega_alpha') %>% pull(Value))
}
get_omega_beta = function(x){
  return(x$parameters %>% filter(Param.name=='omega_beta') %>% pull(Value))
}
get_mu_th_alpha = function(x, th='1'){
  return(x$parameters %>% filter(Param.name==paste0('alpha_mu_th_',th)) %>% pull(Value))
}
get_mu_th_beta = function(x, th='1'){
  return(x$parameters %>% filter(Param.name==paste0('beta_mu_th_',th)) %>% pull(Value))
}
get_diploid_length = function(x){
  x$parameters %>% filter(Param.name=="diploid_length") %>% pull(Value)
}
get_CNA_length = function(x){
  x$parameters %>% filter(Param.name=="CNA_length") %>% pull(Value)
}
get_Major_allele = function(x){
  x$parameters %>% filter(Param.name=="Major") %>% pull(Value)
}
get_Minor_allele = function(x){
  x$parameters %>% filter(Param.name=="Minor") %>% pull(Value)
}
get_cauchy_centre = function(x){
  return(convert_real_date(x$clinical_records %>% filter(Timepoint=='Mutagenic exposure') %>% pull(Start)))
}
get_cauchy_alpha = function(x){
  alpha=x$parameters %>% filter(Param.name=="Cauchy alpha") %>% pull(Value)
  if (is.null(alpha)){
    alpha=10
  }
  beta
}
get_cauchy_beta = function(x){
  beta=x$parameters %>% filter(Param.name=="Cauchy beta") %>% pull(Value)
  if (is.null(beta)){
    beta=10
  }
  beta
}

## Get data for inference
get_inference_data = function(x){
  data[["m_clock"]]
  data[["l_diploid"]]
  data[["mu_clock"]]

  data[["n_cna"]]
  data[["[n_cna] m_alpha"]] # vector
  data[["[n_cna] m_beta"]] # vector
  data[["[n_cna] l_CNA"]] # vector
  data[["[n_cna] coeff"]]

  data[["n_driver"]]
  data[["[n_driver] m_driver"]] # vector
  data[["[n_driver] mu_driver"]] # vector

  data[["n_th_step"]]
  data[["[n_th_step] start_th_step"]] # vector
  data[["[n_th_step] end_th_step"]] # vector
  data[["[n_th_step] alpha_th_step"]] # vector
  data[["[n_th_step] beta_th_step"]] # vector
  data[["[n_th_step] m_th_step"]] # vector

  data[["n_th_cauchy"]]
  data[["[n_th_cauchy] location_th_cauchy"]] # vector
  data[["[n_th_step] scales_th_cauchy"]] # vector
  data[["[n_th_step] m_th_cauchy"]] # vector

  data[["omega_alpha"]]
  data[["omega_beta"]]
  data[["N_min"]]
  data[["N_max"]]
  data[["k_step"]]
  data[["k_softmax"]]

  data[["Sample_1"]]
  data[["Sample_2"]]
}

#### Getters for inferred data

## Get inferred times

## Get inferred parameters
