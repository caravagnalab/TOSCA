## Need to change these functions according to the model

generate_functions = function(x, model='standard'){
  'functions{
  real p_m2(real ti, real tf,real t_1,real t_2, real k){

    vector[2] a1;
    vector[2] a2;
    vector[2] a3;
    vector[2] a4;

    real f1;
    real f2;
    real f3;
    real f4;

    a1[1] = k*(tf-t_1);
    a1[2] = 1;
    a2[1] = k*(tf-t_2);
    a2[2] = 1;
    a3[1] = k*(ti-t_1);
    a3[2] = 1;
    a4[1] = k*(ti-t_2);
    a4[2] = 1;

    f1 = 1/k*(log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1));
    f2 = 1/k*(log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2));
    f3 = 1/k*(log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3));
    f4 = 1/k*(log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4));
    return f1 - f2 - f3 + f4;

  }
}'
}

generate_data = function(x, model='standard'){
  line1= 'data{
  '
  # Mutations
  mutations_lines= sapply(1:nrow(x$mutations), function(m){
    paste0('int <lower=0> ', x$mutations$Mutation.type[m], ';')
  }) %>% paste(collapse='\n')
  # Parametes
  parameters_lines= sapply(1:nrow(x$parameters), function(p){
    paste0('real <lower=0> ', x$parameters$Param.name[p], ';')
  }) %>% paste(collapse='\n')
  # Timepoints
  clinical_times_lines= sapply(1:nrow(x$clinical_records), function(t){
    if (grepl('Therapy', x$clinical_records$Timepoint[t])){
      paste0('real <lower=0> ', x$clinical_records$Timepoint[t], '_start;\nreal <lower= 0> ',x$clinical_records$Timepoint[t], '_end;')
    }else{
      paste0('real <lower=0> ', x$clinical_records$Timepoint[t], ';')
    }
  }) %>% paste(collapse='\n')
  linef= '
  }'

  paste(line1, mutations_lines, '\n',parameters_lines, '\n',clinical_times_lines, linef, collapse='\n')
}

generate_parameters = function(x, model='standard'){

  line1= 'parameters{
  '

  timeECA_line= 'real <lower=0, upper=Sample_1> t_eca;\n'
  last_th = clinical_records %>% filter(grepl('Therapy', Timepoint)) %>% arrange(desc(End)) %>% pull(Timepoint)
  last_th = last_th[length(last_th)]
  timeMRCA_line= paste0('real <lower=',last_th,'_end, upper=Sample_2> t_mrca;\n')
  timeCNA_line = 'real <lower=t_eca, upper=t_mrca> t_cna;\n'

  omega_line = 'real <lower=0> omega;\n'
  mu_th = x$parameters %>% filter(grepl('alpha_mu_th',Param.name)) %>% pull(Param.name) %>% strsplit('alpha_') %>% unlist()
  mu_th = mu_th[grepl('mu_th', mu_th)]
  mu_th_line = sapply(1:length(mu_th), function(m){
    paste0('real <lower=0> ',mu_th[m],';')
  }) %>% paste(collapse='\n')

  linef='
  }'

  paste(line1, timeECA_line, timeMRCA_line, timeCNA_line, omega_line, mu_th_line, linef, collapse='\n')
}

generate_tranfsormed_parameters = function(x, model='standard'){}

generate_model = function(x, model='standard'){

  line1 = 'model{
  '
  # Priors
  tECA_line = 't_eca ~ uniform(0, Sample_1);\n'
  last_th = clinical_records %>% filter(grepl('Therapy', Timepoint)) %>% arrange(desc(End)) %>% pull(Timepoint)
  last_th = last_th[length(last_th)]
  tMRCA_line= paste0('t_mrca ~ uniform(',last_th,'_end,Sample_2);\n')
  tCNA_line = 't_cna ~ uniform(t_eca, t_mrca);\n'

  omega_line = 'omega ~ gamma(omega_alpha, omega_beta);\n'
  mu_th = x$parameters %>%
    filter(grepl('alpha_mu_th',Param.name)) %>% pull(Param.name) %>% strsplit('alpha_') %>% unlist()
  mu_th = mu_th[grepl('mu_th', mu_th)]
  mu_th_line = sapply(1:length(mu_th), function(m){
    paste0(mu_th[m],' ~ gamma(alpha_mu_th_',m,', beta_mu_th_',m,');')
  }) %>% paste(collapse='\n')

  mu_th = x$parameters %>% filter(grepl('alpha_mu_th',Param.name)) %>% pull(Param.name) %>% strsplit('alpha_') %>% unlist()
  mu_th = mu_th[grepl('mu_th', mu_th)]

  # Likelihood
  m_clock_line = 'm_clock ~ poisson(2*omega*diploid_length*mu*(t_mrca - t_eca));\n'

  # alpha: always coeff. 1
  # beta: 2:0 -> coeff 2, 2:1 -> 3, 2:2 -> 4
  m_alpha_line_1 = 'm_alpha ~ poisson(omega*CNA_length*( mu*(t_cna-t_eca) + '
  m_alpha_line_2 = sapply(1:length(mu_th), function(a){
    paste0(mu_th[a], '*p_m2(t_eca, t_cna, Therapy_',a,'_start, Therapy_',a,'_end, k)')
  }) %>% paste(collapse=' + ')
  m_alpha_line = paste0(m_alpha_line_1, m_alpha_line_2, '));\n')

  coeff = x$parameters %>% filter(Param.name == 'Minor') %>% pull(Value)
  coeff = coeff +2
  m_beta_line_1 = paste0('m_beta ~ poisson(',coeff,'*omega*CNA_length*( mu*(t_mrca-t_cna) + ')
  m_beta_line_2 = sapply(1:length(mu_th), function(b){
    paste0(mu_th[b], '*p_m2(t_cna, t_mrca, Therapy_',b,'_start, Therapy_',b,'_end, k)')
  }) %>% paste(collapse=' + ')
  m_beta_line = paste0(m_beta_line_1, m_beta_line_2, '));\n')

  m_th_lines = sapply(1:length(mu_th), function(m){
    paste0('m_th_',m,' ~ poisson(2*omega*diploid_length*',mu_th[m],'*(Therapy_',m,'_end - Therapy_',m,'_start);')
  }) %>% paste(collapse='\n')

  # Penalty
  penalty_line='
  target += -N_min*exp(-omega*(Sample_1-t_eca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_1-t_eca))));
  target += -N_min*exp(-omega*(Sample_2-t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2-t_mrca))));\n'

  linef = '
  }'

  paste(line1, tECA_line, tMRCA_line, tCNA_line, omega_line, mu_th_line, '\n',
        m_clock_line, m_alpha_line, m_beta_line, penalty_line, linef, collapse = '\n')
}

generate_generated_quantities = function(x, model='standard'){
  line_1= 'generated quantities{\n'

  mu_th = x$parameters %>% filter(grepl('alpha_mu_th',Param.name)) %>% pull(Param.name) %>% strsplit('alpha_') %>% unlist()
  mu_th = mu_th[grepl('mu_th', mu_th)]

  m_clock_line = 'int m_clock_rep = poisson_rng(2*omega*diploid_length*mu*(t_mrca - t_eca));\n'
  m_alpha_line_1 = 'int m_alpha_rep = poisson_rng(omega*CNA_length*( mu*(t_cna-t_eca) + '
  m_alpha_line_2 = sapply(1:length(mu_th), function(a){
    paste0(mu_th[a], '*p_m2(t_eca, t_cna, Therapy_',a,'_start, Therapy_',a,'_end, k)')
  }) %>% paste(collapse=' + ')
  m_alpha_line = paste0(m_alpha_line_1, m_alpha_line_2, '));\n')

  coeff = x$parameters %>% filter(Param.name == 'Minor') %>% pull(Value)
  coeff = coeff +2
  m_beta_line_1 = paste0('int m_beta_rep = poisson_rng(',coeff,'*omega*CNA_length*( mu*(t_mrca-t_cna) + ')
  m_beta_line_2 = sapply(1:length(mu_th), function(b){
    paste0(mu_th[b], '*p_m2(t_cna, t_mrca, Therapy_',b,'_start, Therapy_',b,'_end, k)')
  }) %>% paste(collapse=' + ')
  m_beta_line = paste0(m_beta_line_1, m_beta_line_2, '));\n')

  # Penalty
  penalty_line='real N = exp(omega*(Sample_2 - t_mrca));\n'

  line_f= '\n}\n'

  paste(line_1, m_clock_line, m_alpha_line, m_beta_line, penalty_line, line_f, collapse = '\n')
}

generate_stan_code = function(x, model='standard', saving_dir = NULL){
  if (is.null(saving_dir)){
    dir.create(paste0(getwd(), '/TOSCA_fit'))
    saving_dir = paste0(getwd(), '/TOSCA_fit')
  }else{
    dir.create(saving_dir)
  }
  stan_code = paste0(
    generate_functions(x),'\n',
    generate_data(x),'\n',
    generate_parameters(x),'\n',
    generate_tranfsormed_parameters(x),'\n',
    generate_model(x),'\n',
    generate_generated_quantities(x)
  )
  writeLines(stan_code, paste0(saving_dir,"/stan_model.stan"))
  return(stan_code)
}
