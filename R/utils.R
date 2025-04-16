check_input_mutations = function(mutations){
  condition1= class(mutations) == "data.frame"
  condition2= colnames(mutations) == c('Mutation.type','Number.of.mutations')
  condition3= c('m_clock', 'm_alpha', 'm_beta', 'm_th_1') %in% mutations$Mutation.type
}

check_input_clinical = function(clinical_records){}

create_parameters = function(mutations, clinical_records){}

# convert_date_real = function(x) {
#   y = 2000 + x %>% floor()
#   py = (x - (x %>% floor())) * 365
#   m = (py / 30) %>% floor
#   d = (py - (m * 30)) %>% floor
#   date_string = paste(y, m +1, d + 1, sep = '-')
#
#   date_string = ifelse (d>30, paste(y, m +1, 30, sep = '-'), date_string)
#   date_string = ifelse ((m==1 & d>27), paste(y, m +1, 28, sep = '-'), date_string)
#   date_string = ifelse (m>=12, paste(y+1, 1, 1, sep = '-'), date_string)
#
#   return(date_string)
# }
#
# convert_real_date = function(date = NULL, ref_year = 2000) {
#   ref_month = 1
#   ref_day = 1
#
#   year = as.integer(strsplit(date, '-')[[1]][1])
#   month = as.integer(strsplit(date, '-')[[1]][2])
#   day = as.integer(strsplit(date, '-')[[1]][3])
#
#   return((year - ref_year) + (month / 12 - ref_month / 12) + (day / 365 - ref_day / 365))
# }
