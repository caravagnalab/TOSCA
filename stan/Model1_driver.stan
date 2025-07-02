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


  real lambda_therapy(int n_cycles, real ti, real tf, vector[n_cycles] t_therapy_i, vector[n_cycles] t_therapy_f, real k){

    real f=0;

    for (c in 1:n_cycles){
      f += lambda_therapy_single(ti, tf, t_therapy_i[c], t_therapy_f[c], k);
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

}

data{
  // Clock-like mutations
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;


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
  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  // vector<lower=0>[n_th_step] cycles_th_step;
  vector<lower=0>[n_th_step] start_th_step;
  vector<lower=0>[n_th_step] end_th_step;
  vector<lower=0>[n_th_step] type_th_step; // vector with numbers identifying the therapy (1:n_th_type)
  vector<lower=0>[n_th_step_type] alpha_th_step;
  vector<lower=0>[n_th_step_type] beta_th_step;
  vector<lower=0>[n_th_step_type] m_th_step;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  //vector<lower=0>[n_th_cauchy] cycles_th_cauchy;
  vector<lower=0>[n_th_cauchy] location_th_cauchy;
  vector<lower=0>[n_th_cauchy] type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy_type)
  vector<lower=0>[n_th_cauchy_type] alpha_th_cauchy;
  vector<lower=0>[n_th_cauchy_type] beta_th_cauchy;
  vector<lower=0>[n_th_cauchy_type] m_th_cauchy;

  // other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;
  // real <lower=0> k_softmax;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> max_therapy;
  int <lower=0, upper=1> exponential_growth;
  real N_min;
  real N_max;

}

transformed data {
  // Precompute indices for step therapies
  array[n_th_step_type] int[] idx_step;
  for (m in 1:n_th_step_type) {
    int count = 0;
    for (i in 1:n_th_step) {
      if (type_th_step[i] == m) count += 1;
    }
    idx_step[m] = rep_array(0, count);
    int k = 1;
    for (i in 1:n_th_step) {
      if (type_th_step[i] == m) {
        idx_step[m][k] = i;
        k += 1;
      }
    }
  }

  // Precompute indices for cauchy therapies
  array[n_th_cauchy_type] int[] idx_cauchy;
  for (c in 1:n_th_cauchy_type) {
    int count = 0;
    for (i in 1:n_th_cauchy) {
      if (type_th_cauchy[i] == c) count += 1;
    }
    idx_cauchy[c] = rep_array(0, count);
    int k = 1;
    for (i in 1:n_th_cauchy) {
      if (type_th_cauchy[i] == c) {
        idx_cauchy[c][k] = i;
        k += 1;
      }
    }
  }
}


parameters{
  real <lower=0, upper=Sample_1> t_eca;
  real <lower=max_therapy, upper=Sample_2> t_mrca;
  vector<lower=t_eca, upper=t_mrca> t_driver;
  vector<lower=0>[n_th_step_type] mu_th_step;
   vector<lower=0>[n_th_cauchy_type] scales_th_cauchy;
  real <lower=0> omega;
  real <lower=0> mu_driver;
}

model{

  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca ~ uniform(max_therapy, Sample_2);

  t_driver ~ uniform(t_eca, t_mrca);


  for (m in 1:n_th_step_type){
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }

  for (ch in 1:n_th_cauchy_type){
     scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
  }


  omega ~ gamma(omega_alpha, omega_beta);
  mu_driver ~ gamma(mu_driver_alpha,mu_driver_beta);

  // Likelihood
  m_clock ~ poisson(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  for (m in 1:n_th_step_type) {
    vector[size(idx_step[m])] t_start_m;
    vector[size(idx_step[m])] t_end_m;

    for (j in 1:size(idx_step[m])) {
      t_start_m[j] = t_therapy_i[idx_step[m][j]];
      t_end_m[j] = t_therapy_f[idx_step[m][j]];
    }

    m_th_step[m] ~ poisson(2 * l_diploid * omega * mu_th_step[m] *
                         lambda_therapy(size(idx_step[m]), t_eca, t_mrca, t_start_m, t_end_m, k_step));
  }


  for (c in 1:n_th_cauchy_type) {
    vector[size(idx_cauchy[c])] location_c;

    for (j in 1:size(idx_cauchy[c])) {
      location_c[j] = location_th_cauchy[idx_cauchy[c][j]];
    }

    m_th_cauchy[c] ~ poisson(2 * l_diploid * omega * mu_clock *
                           couchy_cdf(size(idx_cauchy[c]), location_c, scales_th_cauchy[c], t_eca, t_mrca));
  }



  if (driver_type==0){
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*lambda_therapy(cycles_driver,t_driver,t_mrca,
    driver_start,driver_end,k_step));
  }

  if (exponential_growth==1){
    target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
  }

}

generated quantities{


  int<lower =0> m_clock_rep = poisson_rng(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  vector[n_th_step_type] m_th_step_rep;

  for (m in 1:n_th_step_type) {
    int n_m = size(idx_step[m]);
    vector[n_m] t_start_m;
    vector[n_m] t_end_m;

    for (j in 1:n_m) {
      t_start_m[j] = t_therapy_i[idx_step[m][j]];
      t_end_m[j] = t_therapy_f[idx_step[m][j]];
    }

    m_th_step_rep[m] = poisson_rng(2 * l_diploid * omega * mu_th_step[m] *
                                 lambda_therapy(n_m, t_eca, t_mrca, t_start_m, t_end_m, k_step));
  }

  vector[n_th_cauchy_type] m_th_cauchy_rep;

  for (c in 1:n_th_cauchy_type) {
    int n_c = size(idx_cauchy[c]);
    vector[n_c] location_c;

    for (j in 1:n_c) {
      location_c[j] = location_th_cauchy[idx_cauchy[c][j]];
    }

    m_th_cauchy_rep[c] = poisson_rng(2 * l_diploid * omega * mu_clock *
                                   couchy_cdf(n_c, location_c, scales_th_cauchy[c], t_eca, t_mrca));
  }



  int <lower =0> m_driver_rep;

  if (driver_type==0){
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*lambda_therapy(cycles_driver,t_driver,t_mrca,
    driver_start,driver_end,k_step));
  }



  if (exponential_growth==1){
    real N = exp(omega*(Sample_2-t_mrca));
  }

}
