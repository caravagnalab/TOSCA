get_mutation_rate = function(x, type){
  if (!(type %in% c("clock", "driver_clock", "driver", "th_step"))) stop("Invalid mutation rate type")

  if (type %in% c("clock", "driver_clock")){
    mu = x$Input$Parameters %>% dplyr::filter(Name==paste0("mu_", type)) %>% pull(Value) %>% as.double()
    if (length(mu)==0 & type == "clock") stop("mu_clock missing!")
    if (length(mu)==0 & type == "driver_clock") mu = x$Input$Parameters %>% dplyr::filter(Name=="mu_clock") %>% pull(Value) %>% as.double()
  }

  if (type == "driver"){
    mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="mu_driver_alpha") %>% pull(Value) %>% as.double()
    mu_beta = x$Input$Parameters %>% dplyr::filter(Name=="mu_driver_beta") %>% pull(Value) %>% as.double()
    if (length(mu_alpha) == 0){
      mu = get_mu_prior_estimates(x, type="driver")
    }else{
      mu = list("alpha"=mu_alpha, "beta"=mu_beta)
    }
  }

  if (type == "th_step") {
    mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="alpha_th_step") #%>% pull(Value) %>% as.double()
    mu_beta = x$Input$Parameters %>% dplyr::filter(Name=="beta_th_step") #%>% pull(Value) %>% as.double()

    if (length(mu_alpha)==nrow(x$Input$Therapies %>% filter(Class=="Mutagenic"))){
      mu_alpha = left_join(
        x$Input$Therapies %>% filter(Class=="Mutagenic"),
        mu_alpha %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
        dplyr::arrange(as.Date(Start)) %>% pull(Value) %>% unique() %>% as.double()
      mu_beta = left_join(
        x$Input$Therapies %>% filter(Class=="Mutagenic"),
        mu_beta %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
        dplyr::arrange(as.Date(Start)) %>% pull(Value) %>% unique() %>% as.double()
    }else{
      therapies = x$Input$Therapies %>% filter(Class=="Mutagenic") %>% arrange(Start)
      mu_alpha = c()
      mu_beta = c()
      for (th in (therapies$Name %>% unique())){
        if (nrow(x$Input$Parameters %>% filter(Name == "alpha_th_step", Index == th))==1) {
          mu_alpha = c(mu_alpha, x$Input$Parameters %>% filter(Name == "alpha_th_step", Index == th) %>% pull(Value))
          mu_beta = c(mu_beta, x$Input$Parameters %>% filter(Name == "beta_th_step", Index == th) %>% pull(Value))
        }else{
          mu_alpha = c(mu_alpha, get_mu_prior_estimates(x, type="th_step", index = th)[["alpha"]])
          mu_beta = c(mu_beta, get_mu_prior_estimates(x, type="th_step", index = th)[["beta"]])
        }
      }

    }
    mu = list("alpha"=mu_alpha, "beta"=mu_beta)
  }
  return(mu)
}

get_growth_rate_estimate = function(x){
  # Second sample name
  sample_names=x$Input$Samples %>% dplyr::arrange(Date) %>% pull(Name)
  last_sample_name = sample_names[length(sample_names)]

  N = x$Input$Parameters %>% dplyr::filter(Name %in% c("N_min", "N_max"), Index == last_sample_name) %>% pull(Value)
  if (length(N)<2) N = c(1e6, 1e13)

  # Last therapy
  last_therapy = x$Input$Therapies %>% dplyr::arrange(Start) %>% pull(Start)
  last_therapy = last_therapy[length(last_therapy)]
  relapse = x$Input$Samples %>% filter(Name == last_sample_name) %>% pull(Date)

  # Omega
  lb_interval = as.integer(as.Date(relapse) - as.Date(last_therapy)) / 365
  Lower_bound = log(N[1]) / (lb_interval)
  ub_interval =  30/365
  Upper_bound = log(N[2]) / (ub_interval)

  omega = mean(c(Lower_bound, Upper_bound))
  omega
}

get_growth_rate_prior_estimate = function(x){
  omega = get_growth_rate_estimate(x)
  omega_alpha = omega*.1
  omega_beta = .1
  # hist(rgamma(1000,omega_alpha, omega_beta))
  list("alpha"=omega_alpha, "beta"=omega_beta)
}

get_mu_driver_estimate = function(x){
  # Exposure period
  starts = x$Input$Therapies %>% dplyr::filter(Class == "Driver responsive") %>% pull(Start)
  ends = x$Input$Therapies %>% dplyr::filter(Class == "Driver responsive") %>% pull(End)
  period = 0
  for (i in 1:length(starts)){
    period = period + as.integer((as.Date(ends[i])-as.Date(starts[i])))/365
  }
  m_driver = x$Input$Mutations %>% dplyr::filter(Type == "Driver") %>% pull(Value)
  len = x$Input$Mutations %>% dplyr::filter(Type == "Driver") %>% pull(Length)
  omega = get_growth_rate_estimate(x)

  mu_driver = m_driver / (2*len*omega*period)
  mu_driver
}

get_mu_th_estimate = function(x, index){
  # Exposure period
  starts = x$Input$Therapies %>% dplyr::filter(Name == index) %>% pull(Start)
  ends = x$Input$Therapies %>% dplyr::filter(Name == index) %>% pull(End)
  period = 0
  for (i in 1:length(starts)){
    period = period + as.integer((as.Date(ends[i])-as.Date(starts[i])))/365
  }
  m_th = x$Input$Mutations %>% dplyr::filter(Type == index) %>% pull(Value)
  len = x$Input$Mutations %>% dplyr::filter(Type == index) %>% pull(Length)
  omega = get_growth_rate_estimate(x)

  mu_th = m_th / (2*len*omega*period)
  mu_th
}

get_mu_prior_estimates = function(x, type, index = NA){
  if (type == "driver"){
    mu_driver = get_mu_driver_estimate(x)
    mu_alpha = mu_driver*(1/mu_driver)*2
    mu_beta = (1/mu_driver)*2
  }
  if (type == "th_step"){
    mu_th = get_mu_th_estimate(x, index = index)
    mu_alpha = mu_th*(1/mu_th)*2
    mu_beta = (1/mu_th)*2
  }

  list("alpha"=mu_alpha, "beta"=mu_beta)
}

get_cauchy_scales = function(x){
  alpha = left_join(
    x$Input$Therapies %>% filter(Class=="Mutagenic cauchy") %>% dplyr::arrange(as.Date(Start)),
    x$Input$Parameters %>% filter(Name == "alpha_th_cauchy") %>% dplyr::rename(Par_name = Name, Name = Index),
    by = "Name"
  ) %>% dplyr::pull(Value) %>% unique() %>% as.double()
  beta = left_join(
    x$Input$Therapies %>% filter(Class=="Mutagenic cauchy") %>% dplyr::arrange(as.Date(Start)),
    x$Input$Parameters %>% filter(Name == "beta_th_cauchy") %>% dplyr::rename(Par_name = Name, Name = Index),
    by = "Name"
  ) %>% dplyr::pull(Value) %>% unique() %>% as.double()
  return(list("alpha"=alpha, "beta"=beta))
}

get_parameter = function(x, name){
  par = x$Input$Parameters %>% dplyr::filter(Name == name) %>% pull(Value) %>% as.double()
  if (length(par)==0){
    if (name == "omega_alpha") par = TOSCA:::get_growth_rate_prior_estimate(x)[["alpha"]]
    if (name == "omega_beta") par = TOSCA:::get_growth_rate_prior_estimate(x)[["beta"]]
    if (name == "k_step") par = 100
    if (name == "exponential_growth") par = 1
    if (name == "mrca_alpha") par = 1
    if (name == "mrca_beta") par = 1
  }
  par
}

get_exponential_growth = function(x){
  left_join(
    x$Input$Samples,
    x$Input$Parameters %>% dplyr::filter(Name == "exponential_growth") %>%
      dplyr::rename(par_name = Name, Name= Index)) %>% pull(Value) %>% as.double()
}

get_N = function(x, which){
  samples = x$Input$Samples
  if (nrow(samples)==3) samples = (samples %>% dplyr::arrange(as.Date(Date)))[2:3,]
  left_join(
    samples,
    x$Input$Parameters %>% dplyr::filter(Name==paste0("N_", which)) %>% dplyr::rename(Par_name=Name, Name=Index),
    by = "Name"
  ) %>% dplyr::arrange(as.Date(Date)) %>% pull(Value) %>% as.double()
}

get_phi = function(x, name){
  if (name %in% c("clock", "driver")) return(x$Input$Parameters %>% dplyr::filter(Name==paste0("phi_", name)) %>% dplyr::pull(Value) %>% as.double())
  if (name == "phi_cna"){
    phi_cna_alpha = x$Input$Parameters %>% filter(Name == "phi_cna_alpha") %>% dplyr::arrange(Index) %>% dplyr::pull(Value) %>% as.double()
    phi_cna_beta = x$Input$Parameters %>% filter(Name == "phi_cna_beta") %>% dplyr::arrange(Index) %>% dplyr::pull(Value) %>% as.double()
    phi_cna = list("alpha"=phi_cna_alpha, "beta"=phi_cna_beta)
    return(phi_cna)
  }
  if (name == "phi_th_step"){
    return(left_join(
      x$Input$Therapies %>% filter(Class == "Mutagenic") %>% dplyr::arrange(as.Date(Start)),
      x$Input$Parameters %>% filter(Name == "phi_th_step") %>% dplyr::rename(Par_Name = Name, Name = Index),
      by = "Name"
    ) %>% dplyr::select(Name, Value) %>% unique() %>% dplyr::pull(Value) %>% as.double())
  }
  if (name == "phi_th_cauchy"){
    return(left_join(
      x$Input$Therapies %>% filter(Class == "Mutagenic cauchy") %>% dplyr::arrange(as.Date(Start)),
      x$Input$Parameters %>% filter(Name == "phi_th_cauchy") %>% dplyr::rename(Par_Name = Name, Name = Index),
      by = "Name"
    ) %>% dplyr::pull(Value) %>% as.double())
  }
}

get_driver_type = function(x){
  if (x$Input$Therapies %>% dplyr::filter(Class == "Driver responsive") %>% nrow() >= 1) 1 else 0
}
