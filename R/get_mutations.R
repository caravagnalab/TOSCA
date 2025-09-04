get_n_cna = function(x){
  (x$Input$Mutations %>% filter(Type %in% c("alpha", "beta")) %>% nrow())/2
}

get_mutation = function(x, type, cauchy=F){
  if (type == "Mutagenic"){
    if (!cauchy){
    drugs = x$Input$Therapies %>% dplyr::filter(Class == "Mutagenic")
    left_join(
      drugs %>% dplyr::arrange(as.Date(Start)),
      x$Input$Mutations %>% dplyr::filter(Type %in% drugs$Name) %>% dplyr::select(Type, Value) %>% dplyr::rename(Name=Type),
      by = "Name") %>% pull(Value) %>% unique() %>% as.integer()
    }else{
      drugs = x$Input$Therapies %>% dplyr::filter(Class == "Mutagenic cauchy")
      left_join(
        drugs %>% dplyr::arrange(as.Date(Start)),
        x$Input$Mutations %>% dplyr::filter(Type %in% drugs$Name) %>% dplyr::select(Type, Value) %>% dplyr::rename(Name=Type),
        by = "Name") %>% pull(Value) %>% unique() %>% as.integer()
    }
  }else{
  x$Input$Mutations %>% dplyr::filter(Type == type) %>% dplyr::arrange(Length) %>% pull(Value) %>% as.integer()
  }
}

get_length = function(x, diploid = T){
  if (diploid){
    l = x$Input$Mutations %>% dplyr::filter(Karyptype == "1:1") %>% dplyr::pull(Length) %>% unique()
    is_wgd_variable = x$Input$Parameters %>% dplyr::filter(Name=="wgd") %>% nrow() > 0
    if (is_wgd_variable)  is_wgd = x$Input$Parameters %>% dplyr::filter(Name=="wgd") %>% dplyr::pull(Value) else is_wgd=NA
    if (length(l)!=1 & is.na(is_wgd)) stop("The entry for the length of the diploid genome must be unique")
  }else{
    l = x$Input$Mutations %>% dplyr::filter(Karyptype != "1:1") %>% dplyr::arrange(Length) %>% pull(Length) %>% unique()
  }
  return(l)
}

get_coeff = function(x){
  karyotypes = x$Input$Mutations %>% dplyr::filter(Type %in% c("alpha", "beta")) %>% dplyr::arrange(Length) %>% dplyr::select(Karyptype) %>% unique()
  coeff_alpha = c()
  coeff_beta = c()
  for (k in karyotypes$Karyptype){
    if (k == "2:0") {
      coeff_alpha = c(coeff_alpha, 1)
      coeff_beta = c(coeff_beta, 2)
    }else{
      coeff_alpha = c(coeff_alpha, 2)
      coeff_beta = c(coeff_beta, 4)
    }
  }
  return(list("alpha"=coeff_alpha, "beta"=coeff_beta))
}
