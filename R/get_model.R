get_model <- function(model_name='Driver', dormancy=F) {

  if (model_name=='CNA' & dormancy) model_name = "CNA_dormancy"

  all_paths <- list(
    "Driver" = "Driver.stan",
    "CNA" = "CNA.stan",
    "CNA_dormancy" = "CNA_dormancy.stan" # "CNA_with_dormancy.stan"
  )

  if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")

  model_path <- system.file("stan", all_paths[[model_name]], package = "TOSCA", mustWork = T)
  tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
  model
}

#### Getters for inferred data
get_inferred_parameters = function(x){
  x$Fit$posterior
}

get_inferred_times_colors = function(){
  c("t_eca"="#ae4532ff", "t_driver"="#d48b3eff", "t_mrca"="#7876adff",
    "t_cna"="#dc5895ff", "t_wgd"="#438455ff")
}

get_original_mutation_name = function(x, name, index){
  if (name == "m_clock_primary") original_name = x$Input$Mutations %>% dplyr::filter(Type=="clock-like primary") %>% pull(Name)
  if (name == "m_clock") original_name = x$Input$Mutations %>% dplyr::filter(Type=="clock-like relapse") %>% pull(Name)

  # if (name %in% c("m_alpha", "m_beta") & x$Fit$model_info$dormancy){
  #   if (name == "m_alpha") original_name = (x$Input$Mutations %>% dplyr::filter(Type=="alpha") %>% dplyr::arrange(Length) %>% dplyr::pull(Name))
  #   if (name == "m_beta") original_name = (x$Input$Mutations %>% dplyr::filter(Type=="beta") %>% dplyr::arrange(Length) %>% dplyr::pull(Name))
  # }else{
    if (name == "m_alpha") original_name = (x$Input$Mutations %>% dplyr::filter(Type=="alpha") %>% dplyr::arrange(Length) %>% dplyr::pull(Name))[index]
    if (name == "m_beta") original_name = (x$Input$Mutations %>% dplyr::filter(Type=="beta") %>% dplyr::arrange(Length) %>% dplyr::pull(Name))[index]
  #}

  if (name == "m_th_step") {
    drugs = x$Input$Therapies %>% dplyr::filter(Class == "Mutagenic") %>% dplyr::pull(Name) %>% unique()
    original_name = (left_join(
      x$Input$Mutations %>% dplyr::filter(Type %in% drugs) %>% dplyr::rename(OR_name = Name, Name = Type),
      x$Input$Therapies %>% filter(Class == "Mutagenic"), by="Name") %>% dplyr::arrange(as.Date(Start)) %>%
      dplyr::select(OR_name) %>% unique() %>% dplyr::pull(OR_name))[index]
  }

  if (name == "m_th_cauchy") {
    drugs = x$Input$Therapies %>% filter(Class == "Mutagenic cauchy") %>% dplyr::pull(Name) %>% unique()
    original_name = (left_join(
      x$Input$Mutations %>% filter(Type %in% drugs) %>% rename(OR_name = Name, Name = Type),
      x$Input$Therapies %>% filter(Class == "Mutagenic cauchy"), by = "Name") %>% dplyr::arrange(as.Date(Start)) %>%
        dplyr::select(OR_name) %>% unique() %>% dplyr::pull(OR_name))[index]
  }

  if (name == "m_driver") original_name = x$Input$Mutations %>% filter(Type=="driver") %>% pull(Name)
  return(original_name)
}
