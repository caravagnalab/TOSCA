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
  (x$mutations %>% filter(Mutation.name=='m_cna') %>% nrow())/2
}

get_mutation = function(x, name, type=NA, index=NA, source=NA, coeff=NA){
  m = x$mutations %>% filter(Mutation.name==name)
  if (!is.na(type)) m = m %>% filter(Mutation.type == type)
  if (!is.na(index)) m = m %>% filter(Mutation.index == index)
  if (!is.na(source)) m = m %>% filter(Mutation.source == source)
  if (!is.na(coeff)) m = m %>% filter(Mutation.coeff == coeff)
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
get_length = function(x, model=NA){
  if (is.na(model)){
    return(x$parameters %>% filter(Parameter.name=="l_diploid") %>% pull(Parameter.value))
  }
  if (model=="CNA"){
    return(x$parameters %>% filter(Parameter.name=="l_CNA") %>% pull(Parameter.value))
  }
  if (model=="WGD"){
    return(x$parameters %>% filter(Parameter.name=="l_tetraploid") %>% pull(Parameter.value))
  }
  # if (karyotype=="2:0" & model=="WGD"){
  #   return(x$parameters %>% filter(Parameter.name=="l_cnloh") %>% pull(Parameter.value))
  # }
}

get_coeff_CNA = function(x){
  coeff_beta = x$parameters %>% filter(Parameter.name=="coeff_CNA") %>% pull(Parameter.value)
  coeff_alpha = c()
  for (c in coeff_beta){
    if (c=="2"){coeff_alpha=c(coeff_alpha, 1)}
    if (c=="4"){coeff_alpha=c(coeff_alpha, 2)}
  }
  list('alpha'=coeff_alpha, 'beta'=as.integer(coeff_beta))
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

get_k_step = function(x){
  x$parameters %>% filter(Parameter.name=='k_step') %>% pull(Parameter.value)
}

get_max_th = function(x){

  max_index = x$clinical_records %>% filter(Clinical.name!="Sample") %>% pull(Clinical.index) %>% as.integer() %>% max() %>% as.character()

  x$clinical_records %>% filter(Clinical.name!="Sample", Clinical.index == max_index) %>% pull(Clinical.value.start)
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
  data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
  data[['l_diploid']] = get_length(x, model=NA)
  if (model %in% c('Driver', 'CNA')){
    # Clock-like mutations
    #data[['m_clock_primary']] = get_mutation(x, name="m_clock", type="primary", index=NA, source=NA)
    data[['m_clock']] = get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA)
  }
  data[['mu_clock']] = get_mu_clock(x)

  if (model == 'Driver'){
  # mutations associated to driver
      data[['driver_type']] = get_driver_type(x)
      data[['cycles_driver']] = get_cycles_drivers(x)
      data[['driver_start']] = get_therapy_driver(x)[["start"]]
      data[['driver_end']] = get_therapy_driver(x)[["end"]]
      data[['m_driver']] = get_mutation(x, name="m_driver", type=NA, index=NA, source=NA)

      if (fixed_mu==T){
        data[['mu_driver']] = get_prior_hyperparameters(x, name='mu_driver')
      }else{
        data[['mu_driver_alpha']] = get_prior_hyperparameters(x, name='mu_driver')[["alpha"]]
        data[['mu_driver_beta']] = get_prior_hyperparameters(x, name='mu_driver')[["beta"]]
      }

      data[['mu_driver_clock']] = get_mu_driver_clock(x)

  }

  if (model == 'CNA'){
    # mutations on CNA
    data[['n_cna']] = get_n_cna(x)
    data[['m_alpha']] = get_mutation(x, name="m_cna", type="alpha", index=NA, source=NA)
    data[['m_beta']] = get_mutation(x, name="m_cna", type="beta", index=NA, source=NA)
    data[['l_CNA']] = get_length(x, model="CNA")
    data[['coeff_alpha']] = get_coeff_CNA(x)[["alpha"]]
    data[['coeff_beta']] = get_coeff_CNA(x)[["beta"]]
  }


    data[['n_th_step']]= get_n_th(x, name = 'Therapy step')
    data[['n_th_step_type']]= get_n_th_type(x, name = 'Therapy step')
    data[['start_th_step']] = start_th(x, type='Therapy step')
    data[['end_th_step']] =  end_th_step(x)
    data[['type_th_step']]= get_type_th_step(x, name='Therapy step')
    data[['alpha_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["alpha"]]
    data[['beta_th_step']]= get_prior_hyperparameters(x, name='mu_th_step')[["beta"]]

    if (model %in% c('Driver', 'CNA')){
      data[['m_th_step']]= get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) # get_m_th(x, type = 'step')
    }

    # mutations associated to cauchy
    data[['n_th_cauchy']]= get_n_th(x, name = 'Therapy cauchy')
    data[['n_th_cauchy_type']]= get_n_th_type(x, name = 'Therapy cauchy')
    data[['location_th_cauchy']]= start_th(x, type='Therapy cauchy')
    data[['type_th_cauchy']]= get_type_th_step(x, name='Therapy cauchy')
    data[['alpha_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["alpha"]]
    data[['beta_th_cauchy']]= get_prior_hyperparameters(x, name='scale_th_cauchy')[["beta"]]
    if (model %in% c('Driver', 'CNA')){
      data[['m_th_cauchy']]= get_mutation(x, name="m_th", type=NA, index=NA, source="cauchy", coeff=NA) #get_m_th(x, type = 'Cauchy')
    }

  if (model == 'WGD'){
    data[['alpha_tetraploid_step']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="step", coeff="4")
    data[['beta_tetraploid_step']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="step", coeff="4")
    data[['alpha_tetraploid_cauchy']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="cauchy", coeff="4")
    data[['beta_tetraploid_cauchy']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="cauchy", coeff="4")
    data[['alpha_tetraploid_clock']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="clock", coeff="4")
    data[['beta_tetraploid_clock']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="clock", coeff="4")

    # data[['alpha_cnloh_step']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="step", coeff="2")
    # data[['beta_cnloh_step']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="step", coeff="2")
    # data[['alpha_cnloh_cauchy']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="cauchy", coeff="2")
    # data[['beta_cnloh_cauchy']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="cauchy", coeff="2")
    # data[['alpha_cnloh_clock']]= get_mutation(x, name="m_wgd", type="alpha", index=NA, source="clock", coeff="2")
    # data[['beta_cnloh_clock']]= get_mutation(x, name="m_wgd", type="beta", index=NA, source="clock", coeff="2")

    data[['l_tetraploid']]= get_length(x, karyotype="2:2", model="WGD")
    # data[['l_cnloh']]= get_length(x, karyotype="2:0", model="WGD")

  }

  # other parameters
  if (fixed_omega==T){
    data[['omega']] = get_prior_hyperparameters(x, name='omega')
  }else{
    data[['omega_alpha']] = get_prior_hyperparameters(x, name='omega')[["alpha"]]
    data[['omega_beta']] = get_prior_hyperparameters(x, name='omega')[["beta"]]
  }

  data[['k_step']] = get_k_step(x)
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

get_model <- function(model_name='Driver', fixed_omega = F, fixed_mu = F) {

  if (model_name=='Driver' & fixed_omega & !fixed_mu) model_name = "Driver_fixed_omega"
  if (model_name=='Driver' & fixed_mu & !fixed_omega) model_name = "Driver_fixed_mu_driver"
  if (model_name=='Driver' & fixed_omega & fixed_mu) model_name = "Driver_fixed_mu_and_omega"

  all_paths <- list(
    "Driver" = "Driver.stan",
    "Driver_fixed_mu_and_omega" = "Driver_fixed_mut_rate_and_growth.stan",
    "Driver_fixed_omega" = "Driver_fixed_growth.stan",
    "Driver_fixed_mu_driver" = "Driver_fixed_mu_driver.stan",
    "CNA" = "CNA.stan",
    "WGD" = "WGD.stan"
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

## Get inferred times

## Get inferred parameters
