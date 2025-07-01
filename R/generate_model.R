generate_functions = function(x){
  line_1= 'functions{'

  step_therapy = '
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

  }'

  cauchy_therapy = '
  real couchy_cdf(real location, real scale, real a,real b){
    real d= ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d);
	}'

  line_f='}'
}

number_of_cna = function(x){
  grepl('m_alpha',x$mutations$Mutation.type) %>% sum()
}

number_of_drivers = function(x){
  grepl('m_driver',x$mutations$Mutation.type) %>% sum()
}

number_of_exogenous_muts = function(x){
  grepl('m_th',x$mutations$Mutation.type) %>% sum()
}

number_of_step_therapies = function(x){
  grepl('Step',x$clinical_records$Type) %>% sum()
}

number_of_cauchy_therapies = function(x){
  grepl('Cauchy',x$clinical_records$Type) %>% sum()
}

generate_data = function(x){

  line_1 = 'data{'
  'int <lower=0> m_clock;
   real <lower=0> mu_clock;
   real <lower=0> l_diploid;
   '
  if (number_of_cna(x)>0){
    lines_cna = '
    int <lower=0> n_cnas;
    vector[n_cnas] int <lower=0> m_alpha;
    vector[n_cnas] int <lower=0> m_beta;
    vector[n_cnas] real <lower=0> l_CNA;
    '
  }

  if (number_of_drivers(x)>0){
    line_drivers = '
    int int <lower=0> n_drivers;
    vector[n_drivers] int <lower=0> m_driver;
    real <lower=0> mu_driver;
    '
  }

  if (number_of_exogenous_muts(x)>0){
    line_drivers = '
    int <lower=0> n_th;
    vector[n_th] real <lower=0> m_th;
    '
    line_step = '
    int <lower>
    vector
    '
  }


  line_f = '}'
}
