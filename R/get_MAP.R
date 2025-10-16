get_MAP_t_eca = function(x, map = "mean"){
  if(map=="mean") map = mean(x$Fit$posterior$t_eca) else map = median(x$Fit$posterior$t_eca)
  TOSCA:::convert_date_real(date = map, x = x)
}
get_MAP_t_cna = function(x, index=NA){
  if (is.na(index)) var="t_cna" else var = paste0("t_cna[",index,"]")
  if(map=="mean") map = mean(x$Fit$posterior[[var]]) else map = median(x$Fit$posterior[[var]])
  TOSCA:::convert_date_real(date = map, x = x)
}
get_MAP_t_mrca = function(x){
  if(map=="mean") map = mean(x$Fit$posterior$t_mrca) else map = median(x$Fit$posterior$t_mrca)
  TOSCA:::convert_date_real(date = map, x = x)
}
get_MAP_growth_rate = function(x){
  if(map=="mean") map = mean(x$Fit$posterior$omega) else map = median(x$Fit$posterior$omega)
  map
}
get_MAP_mutation_rate = function(x, type = "th_step", index=NA){
  # type = "th_step", "driver
  if (is.na(index)) var=paste0("mu_", type) else var = paste0("mu_",type, "[",index,"]")
  # mu_th_step[1]
  if(map=="mean") map = mean(x$Fit$posterior[[var]]) else map = median(x$Fit$posterior[[var]])
  map
}
