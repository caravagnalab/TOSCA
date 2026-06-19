#
#
# get_cycles_drivers = function(x){
#   x$clinical_records %>% dplyr::filter(Clinical.name=='Therapy driver') %>% nrow()
# }
#
# get_prior_hyperparameters = function(x, name){
#   # name = mu_driver*, mu_th_step*, scale_th_cauchy*, omega*, mrca*, s*
#   if (name %in% x$parameters$Parameter.name){
#     x$parameters %>% dplyr::filter(Parameter.name==name) %>% dplyr::pull(Parameter.value) %>% as.double()
#   }else{
#     alpha= x$parameters %>% dplyr::filter(Parameter.name==paste0(name,'_alpha')) %>% dplyr::pull(Parameter.value) %>% as.double()
#     beta= x$parameters %>% dplyr::filter(Parameter.name==paste0(name,'_beta')) %>% dplyr::pull(Parameter.value) %>% as.double()
#     list('alpha'=alpha, 'beta'=beta)
#   }
# }
#
# get_k_step = function(x){
#   x$parameters %>% dplyr::filter(Parameter.name=='k_step') %>% dplyr::pull(Parameter.value) %>% as.double()
# }
#
#
# get_exponential_growth = function(x){
#   x$parameters %>% dplyr::filter(Parameter.name=='exponential_growth') %>% dplyr::pull(Parameter.value) %>% as.integer()
# }
#
# get_N_min = function(x){
#   N_min = x$parameters %>% dplyr::filter(Parameter.name=='N_min') %>% dplyr::pull(Parameter.value) %>% as.double()
#   if (length(N_min) < 2){return(c(N_min, N_min))}else{N_min}
# }
#
# get_N_max= function(x){
#   N_max = x$parameters %>% dplyr::filter(Parameter.name=='N_max') %>% dplyr::pull(Parameter.value) %>% as.double()
#   if (length(N_max) < 2){return(c(N_max, N_max))}else{N_max}
# }
#
#
#
#
# get_extra_therapy = function(x){
#   if (nrow(x$clinical_records %>% dplyr::filter(!(Clinical.name %in% c("Sample", "Chemotherapy")))) > 0){
#     1
#   }else{
#     0
#   }
# }
#
# get_wgd = function(x){
#   if ("WGD" %in% x$mutations$Mutation.source) 1
#   else 0
# }
#
# get_n_clock = function(x){
#   get_mutation(x, name="m_clock", type="relapse", index=NA, source=NA) %>% length()
# }
#
# get_n_step = function(x){
#   get_mutation(x, name="m_th", type=NA, index=NA, source="step", coeff=NA) %>% length()
# }



## Get inferred times

## Get inferred parameters
