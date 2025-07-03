#### Getters for input data

## Get mutations
get_m_clock = function(x){
  as.integer(x$mutations %>% filter(Mutation.name=='m_clock') %>% pull(Mutation.value))
}
get_m_driver = function(x){
  as.integer(x$mutations %>% filter(Mutation.name=='m_driver') %>% pull(Mutation.value))
}
get_m_CNA = function(x, type = 'alpha'){
  as.integer(x$mutations %>% filter(Mutation.name==type) %>% pull(Mutation.value))
}
get_m_th = function(x, type = 'Step'){
  as.integer(x$mutations %>% filter(Mutation.name==type) %>% pull(Mutation.value))
}

## Get Clinical Data
get_sample = function(x, sample='1'){
  x$clinical_records %>% filter(Clinical.name=='Sample', Clinical.type==sample) %>% pull(Clinical.value.start)
}

get_therapy_driver = function(x){
  start = x$clinical_records %>% filter(Clinical.name=='Therapy driver') %>% pull(Clinical.value.start)
  end = x$clinical_records %>% filter(Clinical.name=='Therapy driver') %>% pull(Clinical.value.end)
  list('start'=start,'end'=end)
}

get_n_th = function(x, name = 'Therapy step'){
  x$clinical_records %>% filter(Clinical.name==name) %>% nrow()
}

get_n_th_type = function(x, name = 'Therapy step'){
  n_th_type = x$clinical_records %>% filter(Clinical.name== name) %>% pull(Clinical.type) %>% unique() %>% nrow()
  if (is.null(n_th_type)){n_th_type=0}
  n_th_type
}

start_th = function(x, type='Therapy step'){
  x$clinical_records %>% filter(Clinical.name==type) %>% pull(Clinical.value.start)
}

end_th_step = function(x){
  x$clinical_records %>% filter(Clinical.name=='Therapy step') %>% pull(Clinical.value.end)
}

get_type_th_step = function(x, name='Therapy step'){
  as.integer(x$clinical_records %>% filter(Clinical.name==name) %>% pull(Clinical.type))
}


## Get parameters
get_l_diploid = function(x){
  x$parameters %>% filter(Parameter.name=="l_diploid") %>% pull(Parameter.value)
}

get_mu_clock = function(x){
  x$parameters %>% filter(Parameter.name=="mu_clock") %>% pull(Parameter.value)
}

get_mu_driver_clock = function(x){
  x$parameters %>% filter(Parameter.name=="mu_clock_driver") %>% pull(Parameter.value)
}

get_driver_type = function(x){
  (x$mutations %>% filter(Mutation.name=='m_driver') %>% pull(Mutation.type) %>% as.integer()) -1
}

get_cycles_drivers = function(x){
  x$clinical_records %>% filter(Clinical.name=='Therapy driver') %>% nrow()
}

get_prior_hyperparameters = function(x, name){
  # name = mu_driver*, mu_th_step*, scale_th_cauchy*, omega*
  alpha= x$parameters %>% filter(Parameter.name==paste0(name,'_alpha')) %>% pull(Parameter.value)
  beta= x$parameters %>% filter(Parameter.name==paste0(name,'_beta')) %>% pull(Parameter.value)
  list('alpha'=alpha, 'beta'=beta)
}

get_k_step = function(x){
  x$parameters %>% filter(Parameter.name=='k_step') %>% pull(Parameter.value)
}

get_max_th = function(x){
  dates=c(x$clinical_records %>% filter(Clinical.name!="Sample") %>% pull(Clinical.value.start), x$clinical_records %>% filter(Clinical.name!="Sample") %>% pull(Clinical.value.end))
  max(dates)
}

get_exponential_growth = function(x){
  x$parameters %>% filter(Parameter.name=='exponential_growth') %>% pull(Parameter.value)
}

get_N_min = function(x){
  x$parameters %>% filter(Parameter.name=='N_min') %>% pull(Parameter.value)
}

get_N_max= function(x){
  x$parameters %>% filter(Parameter.name=='N_max') %>% pull(Parameter.value)
}


## Get data for inference
get_inference_data = function(x){
  data = list(
    # Clock-like mutations
    'm_clock' = get_m_clock(x),
    'l_diploid' = get_l_diploid(x),
    'mu_clock' = get_mu_clock(x),


    # mutations associated to driver
    'driver_type' = get_driver_type(x),
    'cycles_driver' = get_cycles_drivers(x),
    'driver_start' = get_therapy_driver(x)[["start"]],
    'driver_end' = get_therapy_driver(x)[["end"]],
    'm_driver' = get_m_driver(x),
    'mu_driver_alpha' = get_prior_hyperparameters(x, name='mu_driver')[["alpha"]],
    'mu_driver_beta' = get_prior_hyperparameters(x, name='mu_driver')[["beta"]],
    'mu_driver_clock' = get_mu_driver_clock(x),

    'n_th_step'= get_n_th(x, name = 'Therapy step'),
    'n_th_step_type'= get_n_th_type(x, name = 'Therapy step'),
    'start_th_step' = start_th(x, type='Therapy step'),
    'end_th_step' =  end_th_step(x),
    'type_th_step'= get_type_th_step(x, name='Therapy step'),
    'alpha_th_step'= get_prior_hyperparameters(x, name='mu_th_step')[["alpha"]],
    'beta_th_step'= get_prior_hyperparameters(x, name='mu_th_step')[["beta"]],
    'm_th_step'= get_m_th(x, type = 'Step'),

    # mutations associated to cauchy
    'n_th_cauchy'= get_n_th(x, name = 'Therapy cauchy'),
    'n_th_cauchy_type'= get_n_th_type(x, name = 'Therapy cauchy'),
    'location_th_cauchy'= start_th(x, type='Therapy cauchy'),
    'type_th_cauchy'= get_type_th_step(x, name='Therapy cauchy'),
    'alpha_th_cauchy'= get_prior_hyperparameters(x, name='scale_th_cauchy')[["alpha"]],
    'beta_th_cauchy'= get_prior_hyperparameters(x, name='scale_th_cauchy')[["beta"]],
    'm_th_cauchy'= get_m_th(x, type = 'Cauchy'),

    # other parameters
    'omega_alpha' = get_prior_hyperparameters(x, name='omega')[["alpha"]],
    'omega_beta' = get_prior_hyperparameters(x, name='omega')[["beta"]],
    'k_step' = get_k_step(x),

    'Sample_1' = get_sample(x, sample='1'),
    'Sample_2' = get_sample(x, sample='2'),
    'max_therapy' = get_max_th(x),
    'exponential_growth' = get_exponential_growth(x),
    'N_min' = get_N_min(x),
    'N_max' = get_N_max(x)
  )
  data
}

#### Getters for inferred data

## Get inferred times

## Get inferred parameters
