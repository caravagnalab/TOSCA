#### Getters for input data

## Get mutations
get_m_clock = function(x, type = "relapse"){
  as.integer(x$mutations %>% filter(Mutation.name=='m_clock', Mutation.type == type) %>% pull(Mutation.value))
}
get_m_driver = function(x){
  as.integer(x$mutations %>% filter(Mutation.name=='m_driver') %>% pull(Mutation.value))
}
get_m_CNA = function(x, type = 'alpha'){
  as.integer(x$mutations %>% filter(Mutation.type==type) %>% pull(Mutation.value))
}
get_n_cna = function(x){
  (x$mutations %>% filter(Mutation.name=='m_cna') %>% nrow())/2
}
get_m_th = function(x, type = 'step'){
  as.integer(x$mutations %>% filter(Mutation.type==type) %>% pull(Mutation.value))
}


get_phi = function(x,type = "clock"){

v = paste0("phi_",type)  

x$parameters %>% filter(Parameter.name == v) %>% pull(Parameter.value) 
  
  
}

get_mutation = function(x, name, type=NA, index=NA){
  m = x$mutations %>% filter(Mutation.name==name)
  if (!is.na(type)) m = m %>% filter(Mutation.type == type)
  if (!is.na(index)) m = m %>% filter(Mutation.index == index)
  m %>% pull(Mutation.value) %>% as.integer()
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
  n_th_type = length(x$clinical_records %>% filter(Clinical.name== name) %>% pull(Clinical.type) %>% unique() )
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
  as.integer(x$clinical_records %>% filter(Clinical.name==name) %>% pull(Clinical.type) )#%>% unique()
}


## Get parameters
get_l_diploid = function(x){
  x$parameters %>% filter(Parameter.name=="l_diploid") %>% pull(Parameter.value)
}

get_l_CNA = function(x){
  x$parameters %>% filter(Parameter.name=="l_CNA") %>% pull(Parameter.value)
}

get_coeff_CNA = function(x){
  x$parameters %>% filter(Parameter.name=="coeff_CNA") %>% pull(Parameter.value)
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
  # name = mu_driver*, mu_th_step*, scale_th_cauchy*, omega*, mrca*
  if (name %in% x$parameters$Parameter.name){
    x$parameters %>% filter(Parameter.name==name) %>% pull(Parameter.value)
  }else{
    alpha= x$parameters %>% filter(Parameter.name==paste0(name,'_alpha')) %>% pull(Parameter.value)
    beta= x$parameters %>% filter(Parameter.name==paste0(name,'_beta')) %>% pull(Parameter.value)
    list('alpha'=alpha, 'beta'=beta)
  }
}

get_k_step = function(x,what){

  x$parameters %>% filter(Parameter.name==paste0('k_step_',what)) %>% pull(Parameter.value)
   # x$parameters %>% filter(Parameter.name==paste0('k_step')) %>% pull(Parameter.value)
}

get_max_th = function(x){

  max_index = x$clinical_records %>% filter(Clinical.name!="Sample") %>% pull(Clinical.type) %>% as.integer() %>% 
    max() %>% as.character()

  x$clinical_records %>% filter(Clinical.name!="Sample", Clinical.type == max_index) %>% pull(Clinical.value.start) %>% min()
  #min(dates)
}

get_exponential_growth = function(x){
  x$parameters %>% filter(Parameter.name=='exponential_growth') %>% pull(Parameter.value)
}

get_N_min = function(x){
  N_min = x$parameters %>% filter(Parameter.name=='N_min') %>% pull(Parameter.value)
  if (length(N_min) < 2){return(c(N_min, N_min))}else{N_min}
}

get_N_max= function(x){
  N_max = x$parameters %>% filter(Parameter.name=='N_max') %>% pull(Parameter.value)
  if (length(N_max) < 2){return(c(N_max, N_max))}else{N_max}
}


## Get data for inference
get_inference_data = function(x, model='Driver', fixed_omega, fixed_mu){

  data = list()

  # Clock-like mutations
  data[['m_clock_primary']] = get_m_clock(x, type = "primary")
  data[['m_clock']] = get_m_clock(x, type = "relapse")
  data[['l_diploid']] = get_l_diploid(x)
  data[['mu_clock']] = get_mu_clock(x)

  if (model == 'Driver'){
  # mutations associated to driver
      data[['driver_type']] = get_driver_type(x)
      data[['cycles_driver']] = get_cycles_drivers(x)
      data[['driver_start']] = get_therapy_driver(x)[["start"]]
      data[['driver_end']] = get_therapy_driver(x)[["end"]]
      data[['m_driver']] = get_m_driver(x)

      if (fixed_mu==T){
        data[['mu_driver']] = get_prior_hyperparameters(x, name='mu_driver')
      }else{
        data[['mu_driver_alpha']] = get_prior_hyperparameters(x, name='mu_driver')[["alpha"]]
        data[['mu_driver_beta']] = get_prior_hyperparameters(x, name='mu_driver')[["beta"]]
      }

      data[['mu_driver_clock']] = get_mu_driver_clock(x)
      
      data[["phi_clock"]] = get_phi(x,type = "clock")
      data[["phi_driver"]] = get_phi(x,type = "driver")

      # if ('mu_driver_clock' %in% fixed_pars){
      #   data[['mu_driver_clock']] = get_mu_driver_clock(x)
      # }else{
      #   data[['mu_driver_clock_alpha']] = get_prior_hyperparameters(x, name='mu_driver_clock')[["alpha"]]
      #   data[['mu_driver_clock_beta']] = get_prior_hyperparameters(x, name='mu_driver_clock')[["beta"]]
      # }

  }
  
  if(model == "Driver_cumulative"){
    
    # data[["m_tail_driver"]] = get_m_driver_tail(x, type = "driver")
    # data[["m_tail_step"]] = get_m_driver_tail(x, type = "step")
    # data[["m_tail_cauchy"]] = get_m_driver_tail(x, type = "cauchy")
    # data[["f_min"]] = get_ccf(x)[["min"]]
    # data[["f_max"]] = get_ccf(x)[["max"]]
    data[["tau_driver"]] = get_tau(x,type = "driver")
    data[["tau_step"]] =  get_tau(x,type = "step")
    data[['driver_type']] = get_driver_type(x)
    data[['cycles_driver']] = get_cycles_drivers(x)
    data[['driver_start']] = get_therapy_driver(x)[["start"]]
    data[['driver_end']] = get_therapy_driver(x)[["end"]]
    data[['m_driver']] = get_m_driver(x)
    
    data[['mu_driver_alpha']] = get_prior_hyperparameters(x, name='mu_driver')[["alpha"]]
    data[['mu_driver_beta']] = get_prior_hyperparameters(x, name='mu_driver')[["beta"]]
    
    
    data[['mu_driver_clock']] = get_mu_driver_clock(x)
    
    data[["phi_clock"]] = get_phi(x,type = "clock")
    data[["phi_driver"]] = get_phi(x,type = "driver")
    
    
  }

  if (model == 'CNA'){
    # mutations on CNA
    data[['n_cna']] = get_n_cna(x)
    data[['m_alpha']] = get_m_CNA(x, type = 'alpha')
    data[['m_beta']] = get_m_CNA(x, type = 'beta')
    data[['l_CNA']] = get_l_CNA(x)
    data[['coeff']] = get_coeff_CNA(x)
  }

  data[['n_th_step']]= get_n_th(x, name = 'Therapy step')
  data[['n_th_step_type']]= get_n_th_type(x, name = 'Therapy step')
  data[['start_th_step']] = start_th(x, type='Therapy step')
  data[['end_th_step']] =  end_th_step(x)
  data[['type_th_step']]= get_type_th_step(x, name='Therapy step')
  data[['alpha_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["alpha"]]
  data[['beta_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["beta"]]
  data[['m_th_step']]= get_m_th(x, type = 'step')
  data[["phi_th_step"]] = get_phi(x, type = "step")

  # mutations associated to cauchy
  data[['n_th_cauchy']]= get_n_th(x, name = 'Therapy cauchy')
  data[['n_th_cauchy_type']]= get_n_th_type(x, name = 'Therapy cauchy')
  data[['location_th_cauchy']]= start_th(x, type='Therapy cauchy')
  data[['type_th_cauchy']]= get_type_th_step(x, name='Therapy cauchy')
  data[['alpha_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["alpha"]]
  data[['beta_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["beta"]]
  data[['m_th_cauchy']]= get_m_th(x, type = 'Cauchy')
  data[["phi_th_cauchy"]] = get_phi(x, type = 'Cauchy')

  # other parameters
  if (fixed_omega==T){
    data[['omega']] = get_prior_hyperparameters(x, name='omega')
  }else{
    data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
    data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
  }

  # data[['k_step']] = get_k_step(x)
  data[['k_step_in_th']] = get_k_step(x,what = "in_th")
  data[['k_step_out_th']] = get_k_step(x,what = "out_th")
  data[['k_step_in_driver']] = get_k_step(x,what = "in_driver")
  data[['k_step_out_driver']] = get_k_step(x,what = "out_driver")
  data[['Sample_1']] = get_sample(x, sample='1')
  data[['Sample_2']] = get_sample(x, sample='2')
  data[['max_therapy']] = get_max_th(x)
  data[['exponential_growth']] = get_exponential_growth(x)
  data[['N_min']] = get_N_min(x)
  data[['N_max']] = get_N_max(x)
  data[['alpha_mrca']] = get_prior_hyperparameters(x, name='mrca')[["alpha"]]
  data[['beta_mrca']] = get_prior_hyperparameters(x, name='mrca')[["beta"]]
  # data[['alpha_eca']] = get_prior_hyperparameters(x, name='eca')[["alpha"]]
  # data[['beta_eca']] = get_prior_hyperparameters(x, name='eca')[["beta"]]

  data
}

# get_m_driver_tail = function(x, type = "driver"){
#   
#   v = paste0("m_tail_",type)  
#   
#   x$mutations  %>% filter(Mutation.name == v) %>% pull(Mutation.value)
#   
#   
# }

# get_ccf = function(x){
#   
#   f = x$parameters  %>% filter(Parameter.name %in% c("f_min","f_max")) %>% pull(Parameter.value) 
#   
#   names(f) = c("min","max")
#   
#   f
#   
# }

get_tau = function(x,type = "driver"){
  
  v = paste0("tau_",type) 
  
  x$parameters  %>% filter(Parameter.name  == v) %>% pull(Parameter.value)  
  
  
}




get_model <- function(model_name='Driver', fixed_omega = F, fixed_mu = F) {

  if (model_name=='Driver' & fixed_omega & !fixed_mu) model_name = "Driver_fixed_omega"
  if (model_name=='Driver' & fixed_mu & !fixed_omega) model_name = "Driver_fixed_mu_driver"
  if (model_name=='Driver' & fixed_omega & fixed_mu) model_name = "Driver_fixed_mu_and_omega"
  if (model_name=="Driver_cumulative" & !fixed_omega & !fixed_mu) model_name = "Driver_cumulative"
  
  all_paths <- list(
    "Driver" = "Driver.stan",
    "Driver_fixed_mu_and_omega" = "Driver_fixed_mut_rate_and_growth.stan",
    "Driver_fixed_omega" = "Driver_fixed_growth.stan",
    "Driver_fixed_mu_driver" = "Driver_fixed_mu_driver.stan",
    "CNA" = "CNA.stan",
    "Driver_cumulative" = "Driver_cumulative.stan"
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
  c("t_eca"="#44998dff", "t_driver"="#cb3144ff", "t_mrca"="#64419bff","t_mrca_primary" = "grey")
}

## Get inferred times

## Get inferred parameters
