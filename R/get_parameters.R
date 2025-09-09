get_mutation_rate = function(x, type){
  if (!(type %in% c("clock", "driver_clock", "driver", "th_step"))) stop("Invalid mutation rate type")

  if (type %in% c("clock", "driver_clock")){
    mu = x$Input$Parameters %>% dplyr::filter(Name==paste0("mu_", type)) %>% pull(Value) %>% as.double()
  }

  if (type == "driver"){
    mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="mu_driver_alpha") %>% pull(Value) %>% as.double()
    mu_beta = x$Input$Parameters %>% dplyr::filter(Name=="mu_driver_beta") %>% pull(Value) %>% as.double()
    mu = list("alpha"=mu_alpha, "beta"=mu_beta)
  }

  if (type == "th_step") {
    mu_alpha = x$Input$Parameters %>% dplyr::filter(Name=="alpha_th_step") #%>% pull(Value) %>% as.double()
    mu_beta = x$Input$Parameters %>% dplyr::filter(Name=="beta_th_step") #%>% pull(Value) %>% as.double()

    mu_alpha = left_join(
      x$Input$Therapies %>% filter(Class=="Mutagenic"),
      mu_alpha %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
      dplyr::arrange(as.Date(Start)) %>% pull(Value) %>% unique() %>% as.double()
    mu_beta = left_join(
      x$Input$Therapies %>% filter(Class=="Mutagenic"),
      mu_beta %>% dplyr::rename(par_name = Name, Name = Index), by = "Name") %>%
      dplyr::arrange(as.Date(Start)) %>% pull(Value) %>% unique() %>% as.double()
    mu = list("alpha"=mu_alpha, "beta"=mu_beta)

  }
  return(mu)
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
  x$Input$Parameters %>% dplyr::filter(Name == name) %>% pull(Value) %>% as.double()
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
