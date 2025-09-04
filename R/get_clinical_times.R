get_drugs_order = function(x){
  x$Input$Therapies %>% dplyr::arrange(as.Date(Start)) %>% dplyr::pull(Name) %>% unique()
}

get_n_therapy_cycles = function(x, class = "Mutagenic"){
  x$Input$Therapies %>% dplyr::filter(Class==class) %>% nrow()
}

get_n_therapy_classes = function(x, class = "Mutagenic"){
  x$Input$Therapies %>% dplyr::filter(Class==class) %>% dplyr::pull(Name) %>% unique() %>% length()
}

get_start_therapy = function(x, class= "Mutagenic"){
  starts = x$Input$Therapies %>% dplyr::filter(Class==class)
  if (nrow(starts) >0){
    starts = starts %>% dplyr::arrange(as.Date(Start)) %>% dplyr::pull(Start)
    return(unname(sapply(starts, TOSCA:::convert_real_date)) %>% unlist())
  }else{
    return(starts = starts %>% dplyr::arrange(as.Date(Start)) %>% dplyr::pull(Start) %>% as.double())
  }
}

get_end_therapy = function(x, class= "Mutagenic"){
  starts = x$Input$Therapies %>% dplyr::filter(Class==class)
  if (nrow(starts) >0){
    starts = starts %>% dplyr::arrange(as.Date(Start)) %>% dplyr::pull(End)
    return(unname(sapply(starts, TOSCA:::convert_real_date)) %>% unlist())
  }else{
    return(starts = starts %>% dplyr::arrange(as.Date(Start)) %>% dplyr::pull(End)%>% as.double())
  }
}

get_therapy_class_index = function(x, class= "Mutagenic"){
  therapy_names = x$Input$Therapies %>% dplyr::filter(Class==class) %>% dplyr::arrange(as.Date(Start)) %>% pull(Name)
  match(therapy_names, unique(therapy_names))
}

get_sample = function(x, sample=1){
  samples = x$Input$Samples %>% dplyr::arrange(as.Date(Date)) %>% dplyr::pull(Date)
  TOSCA:::convert_real_date(samples[sample])
}

get_max_th = function(x){
  starts = x$Input$Therapies %>% dplyr::arrange(as.Date(Start)) %>%
    distinct(Name, .keep_all = TRUE) %>% dplyr::pull(Start)
  TOSCA:::convert_real_date(starts[length(starts)])
}

############## TODO

## Get Clinical Data


get_therapy_driver = function(x){
  start = x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% dplyr::pull(Clinical.value.start)
  end = x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% dplyr::pull(Clinical.value.end)
  list('start'=start,'end'=end)
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
