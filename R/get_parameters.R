get_mutation_rate = function(x, type){
  if (!(type %in% c("clock", "driver_clock", "driver", "th_step"))) stop("Invalid mutation rate type")
  mu = x$Input$Parameters %>% dplyr::filter(Name==paste0("mu_", type))
  if (type == "th_step") {
    mu = mu %>% dplyr::rename(Mutation_rate = Name, Name= Index)
    left_join(x$Input$Therapies %>% filter(Class=="Mutagenic"), mu, by = "Name") %>% dplyr::arrange(as.Date(Start)) %>% pull(Value) %>% unique() %>% as.double()
  }else{
    mu %>% pull(Value) %>% as.double()
  }
}

get_cauchy_scales = function(x){
  left_join(
    x$Input$Therapies %>% filter(Class=="Mutagenic cauchy") %>% dplyr::arrange(as.Date(Start)),
    x$Input$Parameters %>% filter(Name == "scales_cauchy") %>% dplyr::rename(Par_name = Name, Name = Index),
    by = "Name"
  ) %>% dplyr::pull(Value) %>% unique() %>% as.double()
}

get_parameter = function(x, name){
  x$Input$Parameters %>% dplyr::filter(Name == name) %>% pull(Value) %>% as.double()
}

get_N = function(x, which){
  left_join(
    x$Input$Samples,
    x$Input$Parameters %>% dplyr::filter(Name==paste0("N_", which)) %>% dplyr::rename(Par_name=Name, Name=Index),
    by = "Name"
  ) %>% dplyr::arrange(as.Date(Date)) %>% pull(Value) %>% as.double()
}

get_phi = function(x, name){
  if (name %in% c("clock", "driver")) return(x$Input$Parameters %>% dplyr::filter(Name==paste0("phi_", name)) %>% dplyr::pull(Value) %>% as.double())
  if (name == "phi_cna"){
    return(x$Input$Parameters %>% filter(Name == "phi_cna") %>% dplyr::arrange(Index) %>% dplyr::pull(Value) %>% as.double())
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
