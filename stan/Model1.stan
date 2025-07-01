functions {

  real lambda_therapy_single(real ti, real tf,real t_therapy_i,real t_therapy_f, real k){

    vector[2] a1;
    vector[2] a2;
    vector[2] a3;
    vector[2] a4;

    real f1;
    real f2;
    real f3;
    real f4;

    a1[1] = k*(tf-t_therapy_i);
    a1[2] = 1;
    a2[1] = k*(tf-t_therapy_f);
    a2[2] = 1;
    a3[1] = k*(ti-t_therapy_i);
    a3[2] = 1;
    a4[1] = k*(ti-t_therapy_f);
    a4[2] = 1;

    f1 = 1/k*(log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1));
    f2 = 1/k*(log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2));
    f3 = 1/k*(log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3));
    f4 = 1/k*(log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4));
    return f1 - f2 - f3 + f4;

  }

  real lambda_therapy(int n_cycles, real ti, real tf, vector[n_cycles] t_therapy_i, vector[n_cycles] t_therapy_f, real k){

    real f=0;

    for (c in 1:n_cycles){
      f += lambda_therapy_single(ti[c], tf[c], t_therapy_i, t_therapy_f, k);
    }
    return(f);

  }

	real couchy_cdf_single(real location, real scale, real a,real b){
    real d= ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d);
	}

	real couchy_cdf(int n_cycles, vector[n_cycles] location, vector[n_cycles] scale, real a,real b){
	  real c=0;
	  for (c in 1:n_cycles){
      c += couchy_cdf_single(location[c], scale[c], a, b);
	  }
    return(c);
	}

	real smooth_max(real a, real b, real k_softmax) {
    return log_sum_exp(k_softmax * a, k_softmax * b) / k_softmax;
  }

  real smooth_min(real x, real y, real k_softmax) {
  return -log_sum_exp(-k_softmax * x, -k_softmax * y) / k_softmax;
  }

}

data{
  // Clock-like mutations
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;


  // mutations on CNA
  int <lower=0> n_cna;
  vector<lower=0>[n_cna] m_alpha;
  vector<lower=0>[n_cna] m_beta;
  vector<lower=0>[n_cna] l_CNA;
  vector<lower=0>[n_cna] coeff;

  // mutations on WGD
  int <lower=0> n_cna_wgd;
  int <lower=0> n_sbs;
  matrix[n_cna_wgd, n_sbs] m_alpha_wgd;
  matrix[n_cna_wgd, n_sbs] m_beta_wgd;
  vector<lower=0>[n_cna_wgd] l_CNA_wgd;
  vector<lower=0>[n_cna_wgd] coeff_wgd;

  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  vector[cycles_driver] driver_start; // if driver have effect only associated to external therapy
  vector[cycles_driver] driver_end;
  int <lower=0> m_driver;
  real <lower=0> mu_driver_alpha;
  real <lower=0> mu_driver_beta;
  real <lower=0> mu_driver_clock; // if driver alters basal mutation rate of clock-like

  // mutations associated to step-like therapies
  int <lower=0> n_th_step;
  vector<lower=0>[n_th_step] cycles_th_step;
  vector<lower=0>[n_th_step*cycles_th_step] start_th_step;
  vector<lower=0>[n_th_step*cycles_th_step] end_th_step;
  vector<lower=0>[n_th_step*cycles_th_step] type_th_step; // vector with numbers identifying the therapy
  vector<lower=0>[n_th_step] alpha_th_step;
  vector<lower=0>[n_th_step] beta_th_step;
  vector<lower=0>[n_th_step] m_th_step;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  vector<lower=0>[n_th_cauchy] cycles_th_cauchy;
  vector<lower=0>[n_th_cauchy*cycles_th_cauchy] location_th_cauchy;
  vector<lower=0>[n_th_cauchy*cycles_th_cauchy] scales_th_cauchy;
  vector<lower=0>[n_th_cauchy*cycles_th_cauchy] type_th_cauchy;
  vector<lower=0>[n_th_cauchy] m_th_cauchy;

  // other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> N_min;
  real <lower=0> N_max;
  real <lower=0> k_step;
  real <lower=0> k_softmax;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;

}

parameters{
  real <lower=Sample_1-1, upper=Sample_1> t_eca;
  real <lower=max(location_th_cauchy[n_th_cauchy], start_th_step[n_th_step]), upper=Sample_2> t_mrca;
  vector<lower=t_eca, upper=t_mrca>[n_cna] t_cn;
  vector<lower=t_eca, upper=t_mrca>[n_driver] t_driver;
  vector<lower=0>[n_th_step] mu_th_step;
  real <lower=0> omega;
  real <lower=0> mu_driver;
}

model{

  // Priors
  t_eca ~ uniform(Sample_1-1, Sample_1);
  t_mrca ~ uniform(max(location_th_cauchy[n_th_cauchy], end_th_step[n_th_step]), Sample_2);

  for (cna in 1:n_cna){
    t_cn[cna] ~ uniform(t_eca, t_mrca);
  }
  for (driver in 1:n_driver){
    t_driver[driver] ~ uniform(t_eca, t_mrca);
  }
  for (m in 1:n_th_step){
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }
  omega ~ gamma(omega_alpha, omega_beta);
  mu_driver ~ gamma(mu_driver_alpha,mu_driver_beta);

  // Likelihood
  m_clock ~ poisson(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver));

  for (m in 1:n_th_step){
    m_th_step[m] ~ poisson(2*l_diploid*omega*mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_eca, t_mrca, start_th_step[m], end_th_step[m], k_step));
  }

  if (driver_type==0){
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*lambda_therapy(cycles_driver, t_eca, t_mrca, smooth_max(driver_start, t_driver, k_softmax), smooth_min(t_mrca, driver_end, k_softmax), k_step));
  }

  for (cauchy in 1:n_th_cauchy){
    m_th_cauchy[cauchy] ~ poisson(2*l_diploid*omega*mu_clock*couchy_cdf(cycles_th_cauchy[cauchy],cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_eca, t_mrca));
  }

  // add clock
  for (cna in 1:n_cna){
    real lambda_alpha_th_step=0;
    real lambda_beta_th_step=0;
    real lambda_alpha_driver=0;
    real lambda_beta_driver=0;
    real lambda_alpha_th_cauchy=0;
    real lambda_beta_th_cauchy=0;

    real lambda_alpha_clock = mu_clock*(t_cn[cna]-t_eca);
    real lambda_beta_clock = mu_clock*(t_mrca-t_cn[cna]);

    for (m in 1:n_th_step){
      lambda_alpha_th_step += mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_eca, t_cn[cna], start_th_step[m], end_th_step[m], k_step));
      lambda_beta_th_step += mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_cn[cna], t_mrca, start_th_step[m], end_th_step[m], k_step));
      }
    for (driver in 1:n_driver){
      lambda_alpha_driver += mu_driver[driver]*lambda_therapy(cycles_th_step[driver], t_eca, t_cn[cna], t_driver[driver], t_mrca, k_step));
      lambda_beta_driver += mu_driver[driver]*lambda_therapy(cycles_th_step[driver], t_cn[cna], t_mrca, t_driver[driver], t_mrca, k_step));
      }
    for (cauchy in 1:n_th_cauchy){
      lambda_alpha_th_cauchy += couchy_cdf(cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_eca, t_cn[cna]));
      lambda_beta_th_cauchy += couchy_cdf(cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_cn[cna], t_mrca));
      }

    m_alpha[cna] ~ poisson(l_CNA[cna]*omega*(lambda_alpha_clock + lambda_alpha_th_step + lambda_alpha_driver + lambda_alpha_th_cauchy));
    m_beta[cna] ~ poisson(l_CNA[cna]*omega*coeff*(lambda_beta_clock + lambda_beta_th_step + lambda_beta_driver + lambda_beta_th_cauchy));

  }

  target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
}

generated quantities{
  // Test the soft max!

  int m_clock_rep = poisson_rng(2*l_diploid*mu_clock*omega*(t_mrca-t_eca));
  vector[n_th_step] m_th_step_rep;
  for (m in 1:n_th_step){
    m_th_step_rep[m] = poisson_rng(2*l_diploid*omega*mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_eca, t_mrca, start_th_step[m], end_th_step[m], k_step));
  }
  vector[n_driver] m_driver_rep;
  for (driver in 1:n_driver){
    m_driver_rep[driver] = poisson_rng(2*l_diploid*omega*mu_driver[driver]*lambda_therapy(cycles_th_step[driver], t_eca, t_mrca, t_driver[driver], t_mrca, k_step));
  }
  vector[n_th_cauchy] m_th_cauchy_rep;
  for (cauchy in 1:n_th_cauchy){
    m_th_cauchy_rep[cauchy] = poisson_rng(2*l_diploid*omega*mu_clock*couchy_cdf(cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_eca, t_mrca));
  }

  vector[n_cna] m_alpha_rep;
  vector[n_cna] m_beta_rep;
  for (cna in 1:n_cna){
    real lambda_alpha_th_step_rep=0;
    real lambda_beta_th_step_rep=0;
    real lambda_alpha_driver_rep=0;
    real lambda_beta_driver_rep=0;
    real lambda_alpha_th_cauchy_rep=0;
    real lambda_beta_th_cauchy_rep=0;

    real lambda_alpha_clock_rep = mu_clock*(t_cn[cna]-t_eca);
    real lambda_beta_clock_rep = mu_clock*(t_mrca-t_cn[cna]);

    for (m in 1:n_th_step){
      lambda_alpha_th_step_rep += mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_eca, t_cn[cna], start_th_step[m], end_th_step[m], k_step));
      lambda_beta_th_step_rep += mu_th_step[m]*lambda_therapy(cycles_th_step[m], t_cn[cna], t_mrca, start_th_step[m], end_th_step[m], k_step));
      }
    for (driver in 1:n_driver){
      lambda_alpha_driver_rep += mu_driver[driver]*lambda_therapy(cycles_th_step[driver], t_eca, t_cn[cna], t_driver[driver], t_mrca, k_step));
      lambda_beta_driver_rep += mu_driver[driver]*lambda_therapy(cycles_th_step[driver], t_cn[cna], t_mrca, t_driver[driver], t_mrca, k_step));
      }
    for (cauchy in 1:n_th_cauchy){
      lambda_alpha_th_cauchy_rep += couchy_cdf(cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_eca, t_cn[cna]));
      lambda_beta_th_cauchy_rep += couchy_cdf(cycles_th_cauchy[cauchy], location_th_cauchy[cauchy], scales_th_cauchy[cauchy], t_cn[cna], t_mrca));
      }

    m_alpha_rep[cna] = poisson_rng(l_CNA[cna]*omega*(lambda_alpha_clock_rep + lambda_alpha_th_step_rep + lambda_alpha_driver_rep + lambda_alpha_th_cauchy_rep));
    m_beta_rep[cna] = poisson_rng(l_CNA[cna]*omega*coeff*(lambda_beta_clock_rep + lambda_beta_th_step_rep + lambda_beta_driver_rep + lambda_beta_th_cauchy_rep));

  }

  real N = exp(omega*(Sample_2-t_mrca));
}
