functions {

  real lambda_therapy_single(real ti, real tf, real t_therapy_i, real t_therapy_f, real k) {
  real a1 = k * (tf - t_therapy_i);
  real a2 = k * (tf - t_therapy_f);
  real a3 = k * (ti - t_therapy_i);
  real a4 = k * (ti - t_therapy_f);

  real f1 = log_sum_exp(0, a1) / k;
  real f2 = log_sum_exp(0, a2) / k;
  real f3 = log_sum_exp(0, a3) / k;
  real f4 = log_sum_exp(0, a4) / k;

  return f1 - f2 - f3 + f4;
  }


	real couchy_cdf_single(real location, real scale, real a,real b){
    real d= ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d);
	}

}

data{
  // Clock-like mutations
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;


  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  array[cycles_driver] real<lower=0> driver_start; // if driver have effect only associated to external therapy
  array[cycles_driver] real<lower=0> driver_end;
  int <lower=0> m_driver;
  real <lower=0> mu_driver_alpha;
  real <lower=0> mu_driver_beta;
  real <lower=0> mu_driver_clock; // if driver alters basal mutation rate of clock-like

  // mutations associated to step-like therapies
  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  array[n_th_step] real<lower=0> start_th_step;
  array[n_th_step] real<lower=0> end_th_step;
  array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)
  array[n_th_step_type] real<lower=0> alpha_th_step;
  array[n_th_step_type] real<lower=0> beta_th_step;
  array[n_th_step_type] int<lower=0> m_th_step;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  //vector<lower=0>[n_th_cauchy] cycles_th_cauchy;
  array[n_th_cauchy] real<lower=0> location_th_cauchy;
  array[n_th_cauchy] int<lower=0> type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy)
  array[n_th_cauchy_type] real<lower=0> alpha_th_cauchy;
  array[n_th_cauchy_type] real<lower=0> beta_th_cauchy;
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy;

  // other parameters
  // real <lower=0> omega_alpha;
  // real <lower=0> omega_beta;
  real <lower=0> omega;
  real <lower=0> k_step;
  // real <lower=0> k_softmax;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> max_therapy;
  int <lower=0, upper=1> exponential_growth;
  real N_min;
  real N_max;

  real <lower=0> alpha_mrca;
  real <lower=0> beta_mrca;
  // real <lower=0> alpha_eca;
  // real <lower=0> beta_eca;

}

parameters{
  real <lower=0, upper=Sample_1> t_eca;
  // real <lower=max_therapy, upper=Sample_2> t_mrca;
  real <lower=0, upper=1> rho_mrca;
  // real <lower=0, upper=1> rho_eca;
  real <lower=t_eca, upper=driver_end[cycles_driver]> t_driver;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  // real <lower=0> omega;
  real <lower=0> mu_driver;
}

transformed parameters{

  real <lower=max_therapy> t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);
  // real <lower=0, upper=Sample_1> t_eca = Sample_1 - rho_eca;
  // real <lower=t_eca, upper=driver_end[cycles_driver]> t_driver;

  array[n_th_step_type] real lambda_th_step;
  array[n_th_cauchy_type] real lambda_th_cauchy;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy[i] = 0;

  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    for (cycle in 1:n_th_step){
      if (type_th_step[cycle] == th_type){
        lambda_th_step[th_type] += lambda_therapy_single(t_eca, t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);
      }
    }
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    for (cycle in 1:n_th_cauchy){
      if (type_th_cauchy[cycle] == th_cauchy){
        lambda_th_cauchy[th_cauchy] += couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], t_eca, t_mrca);
      }
    }
  }
  }

  real lambda_driver=0;
    for (c in 1:cycles_driver){
    lambda_driver += lambda_therapy_single(t_driver,t_mrca, driver_start[c],driver_end[c],k_step);
    }

}

model{

  // Priors
  t_eca ~ uniform(0, Sample_1);
  // t_mrca ~ uniform(max_therapy, Sample_2);
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
  // rho_eca ~ beta(alpha_eca, beta_eca);
  t_driver ~ uniform(t_eca, t_mrca);

  for (m in 1:n_th_step_type){
      mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
    }

  for (ch in 1:n_th_cauchy_type){
      scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
    }


  // omega ~ gamma(omega_alpha, omega_beta);
  mu_driver ~ gamma(mu_driver_alpha,mu_driver_beta);

  // Likelihood
  m_clock ~ poisson(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    m_th_step[th_type] ~ poisson(2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type]);
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    m_th_cauchy[th_cauchy] ~ poisson(2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy]);
  }
  }

  if (driver_type==0){
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*lambda_driver);
  }

  if (exponential_growth==1){
    target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
  }

}

generated quantities{

  int<lower =0> m_clock_rep = poisson_rng(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  int <lower =0> m_driver_rep;

  if (driver_type==0){
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*lambda_driver);
  }


  if (exponential_growth==1){
    real N = exp(omega*(Sample_2-t_mrca));
  }

}
