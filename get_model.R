get_model <- function(model_name='Driver', dormancy=F) {

  if (model_name=='CNA' & dormancy) model_name = "CNA_dormancy"

  all_paths <- list(
    "Driver" = "Driver.stan",
    "CNA" = "CNA.stan",
    "CNA_dormancy" = "CNA_with_dormancy.stan"
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

get_original_mutation_name = function(x, name, index){
  if (name == "m_clock_primary") original_name = x$mutations %>% filter(Mutation.name=="m_clock", Mutation.type=="primary") %>% pull(Mutation.original.name)
  if (name == "m_clock") original_name = x$mutations %>% filter(Mutation.name=="m_clock", Mutation.type=="relapse") %>% pull(Mutation.original.name)
  if (name == "m_alpha") original_name = x$mutations %>% filter(Mutation.name=="m_cna", Mutation.type=="alpha", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_beta") original_name = x$mutations %>% filter(Mutation.name=="m_cna", Mutation.type=="beta", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_th_step") original_name = x$mutations %>% filter(Mutation.name=="m_th", Mutation.index == index) %>% pull(Mutation.original.name)
  if (name == "m_driver") original_name = x$mutations %>% filter(Mutation.name=="m_driver") %>% pull(Mutation.original.name)
  return(original_name)
}
