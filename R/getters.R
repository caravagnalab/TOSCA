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
get_inference_data_standard = function(x){

  n_therapies = length(grep('Therapy', x$clinical_records$Timepoint))

  data = list(
    'm_clock'= as.integer(get_m_clock(x)),
    'm_alpha'= as.integer(get_m_alpha(x)),
    'm_beta'= as.integer(get_m_beta(x))
  )
  for (i in 1:n_therapies){
    data[[paste0('m_th_',i)]] = as.integer(get_m_th(x, th=i))
  }

  data[["mu"]]= get_mutation_rate(x)
  data[["diploid_length"]] = get_diploid_length(x)
  data[["CNA_length"]] = get_CNA_length(x)
  data[["Major"]] = get_Major_allele(x)
  data[["Minor"]] = get_Minor_allele(x)
  data[["N_max"]] = get_N_max(x)
  data[["N_min"]] = get_N_min(x)
  data[["k"]] = get_k(x)
  data[["omega_alpha"]] = get_omega_alpha(x)
  data[["omega_beta"]] = get_omega_beta(x)
  for (i in 1:n_therapies){
    data[[paste0('alpha_mu_th_',i)]] = get_mu_th_alpha(x, th=i)
    data[[paste0('beta_mu_th_',i)]] = get_mu_th_beta(x, th=i)
  }
  data[["Sample_1"]] = get_sampling_time(x,time='1')
  data[["Sample_2"]] = get_sampling_time(x,time='2')

  for (i in 1:n_therapies){
    data[[paste0('Therapy_',i,'_start')]] = get_therapy_start(x, th=i)
    data[[paste0('Therapy_',i,'_end')]] = get_therapy_end(x, th=i)
  }
  data
}
get_inference_data_dormancy = function(x){

  data = list()

  data[["m_clock"]] = as.integer(get_m_clock(x))
  data[["m_chemo"]] = as.integer(get_m_chemo(x))
  data[["m_th_1"]] = as.integer(get_m_th(x, th='1'))
  data[["m_alpha"]] = as.integer(get_m_alpha(x))
  data[["m_beta"]] = as.integer(get_m_beta(x))
  data[["Major"]]= get_Major_allele(x)
  data[["Minor"]]= get_Minor_allele(x)
  data[["N_min"]] = get_N_min(x)
  data[["N_max"]] = get_N_max(x)

  data[["mu"]] = get_mutation_rate(x)
  data[["omega_alpha"]] = get_omega_alpha(x)
  data[["omega_beta"]] = get_omega_beta(x)
  data[["alpha_mu_th_1"]] = get_mu_th_alpha(x, th='1')
  data[["beta_mu_th_1"]] = get_mu_th_beta(x, th='1')

  data[["k"]] = get_k(x)
  data[["k_sm"]] = get_k_sm(x)
  data[["CNA_length"]] = get_CNA_length(x)
  data[["diploid_length"]] = get_diploid_length(x)
  data[["Sample_2"]] = get_sampling_time(x, time='2')
  data[["Chemo_start"]] = get_chemo_start(x)
  data[["Therapy_1_start"]] = get_therapy_start(x, th='1')
  data[["Chemo_end"]] = get_chemo_end(x)
  data[["Therapy_1_end"]] = get_therapy_end(x, th='1')

  data
}
get_inference_data_cauchy = function(x){

  data = list()

  data[["m_clock"]]= as.integer(get_m_clock(x))
  data[["m_chemo"]]= as.integer(get_m_chemo(x))
  data[["m_th_1"]]= as.integer(get_m_th(x, th='1'))
  data[["m_mutag"]]= as.integer(get_m_th(x)) # cauchy
  data[["Major"]]= get_Major_allele(x)
  data[["Minor"]]= get_Minor_allele(x)
  data[["m_alpha"]]= as.integer(get_m_alpha(x))
  data[["m_beta"]]= as.integer(get_m_beta(x))
  data[["N_min"]]= get_N_min(x)
  data[["N_max"]]= get_N_max(x)
  data[["mu"]]= get_mutation_rate(x)
  data[["omega_alpha"]]= get_omega_alpha(x)
  data[["omega_beta"]]= get_omega_alpha(x)
  data[["alpha_mu_th_1"]]= get_mu_th_alpha(x, th='1')
  data[["beta_mu_th_1"]]= get_mu_th_beta(x, th='2')

  data[["k"]]= get_k(x)
  data[["k_sm"]]= get_k_sm(x)
  data[["CNA_length"]]= get_CNA_length(x)
  data[["diploid_length"]]= get_diploid_length(x)
  data[["Sample_2"]]= get_sampling_time(x, time='2')

  data[["Chemo_start"]]= get_chemo_start(x)
  data[["Therapy_1_start"]]= get_therapy_start(x, th='1')
  data[["Chemo_end"]]= get_chemo_end(x)
  data[["Therapy_1_end"]]= get_therapy_end(x, th='1')

  data[["cauchy_centre"]]= get_cauchy_centre(x)
  data[["cauchy_centre_alpha"]]= get_cauchy_alpha(x)
  data[["cauchy_centre_beta"]]= get_cauchy_beta(x)

  data
}

get_inference_data = function(x, model='standard'){

  if (model == 'standard'){
    data = get_inference_data_standard(x)
  }
  if (model == 'dormancy'){
    data = get_inference_data_dormancy(x)
  }
  if (model == 'cauchy'){
    data = get_inference_data_cauchy(x)
  }

  data
}

#### Getters for inferred data

## Get inferred times

## Get inferred parameters
