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

  int <lower=0> m_clock_primary;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;

  // mutations associated to step-like therapies
  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  array[n_th_step] real<lower=0> start_th_step;
  array[n_th_step] real<lower=0> end_th_step;
  array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)
  array[n_th_step_type] real<lower=0> alpha_th_step;
  array[n_th_step_type] real<lower=0> beta_th_step;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  array[n_th_cauchy] real<lower=0> location_th_cauchy;
  array[n_th_cauchy] int<lower=0> type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy)
  array[n_th_cauchy_type] real<lower=0> alpha_th_cauchy;
  array[n_th_cauchy_type] real<lower=0> beta_th_cauchy;

  // WGD
  array[n_th_step_type] int <lower=0> alpha_tetraploid_step; // mutations in 2:2, for each therapy + clock
  array[n_th_step_type] int <lower=0> beta_tetraploid_step;
  array[n_th_cauchy_type] int <lower=0> alpha_tetraploid_cauchy;
  array[n_th_cauchy_type] int <lower=0> beta_tetraploid_cauchy;
  int <lower=0> alpha_tetraploid_clock;
  int <lower=0> beta_tetraploid_clock;

  // array[n_th_step_type] int <lower=0> alpha_cnloh_step; // mutations in 2:0, for each therapy + clock
  // array[n_th_step_type] int <lower=0> beta_cnloh_step;
  // array[n_th_cauchy_type] int <lower=0> alpha_cnloh_cauchy;
  // array[n_th_cauchy_type] int <lower=0> beta_cnloh_cauchy;
  // int <lower=0> alpha_cnloh_clock;
  // int <lower=0> beta_cnloh_clock;

  real <lower=0> l_tetraploid;
  // real <lower=0> l_cnloh;

  // other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> max_therapy;
  int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  real <lower=0> alpha_mrca;
  real <lower=0> beta_mrca;

}

parameters{
  real <lower=0, upper=Sample_1 - 1e-8 > t_eca;
  real <lower=t_eca, upper=Sample_1 > t_mrca_primary;
  real <lower=t_eca, upper=Sample_2> t_wgd;
  real <lower=0, upper=1> rho_mrca;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  real <lower=0> omega;
}

transformed parameters{

  real <lower=max_therapy> t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);

  array[n_th_step_type] real lambda_th_step_alpha;
  array[n_th_step_type] real lambda_th_step_beta;
  array[n_th_cauchy_type] real lambda_th_cauchy_alpha;
  array[n_th_cauchy_type] real lambda_th_cauchy_beta;

  for (i in 1:n_th_step_type) lambda_th_step_alpha[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy_alpha[i] = 0;
  for (i in 1:n_th_step_type) lambda_th_step_beta[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy_beta[i] = 0;

  // Step therapy mutations
  if (n_th_step_type > 0){
    for (th_type in 1:n_th_step_type){
      for (cycle in 1:n_th_step){

        if (type_th_step[cycle] == th_type){

          lambda_th_step_alpha[th_type] += lambda_therapy_single(t_eca, t_wgd, start_th_step[cycle], end_th_step[cycle], k_step);
          lambda_th_step_beta[th_type] += lambda_therapy_single(t_wgd, t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);

        }
      }
   }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    for (cycle in 1:n_th_cauchy){
      if (type_th_cauchy[cycle] == th_cauchy){

        lambda_th_cauchy_alpha[th_cauchy] += couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], t_eca, t_wgd);
        lambda_th_cauchy_beta[th_cauchy] += couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], t_wgd, t_mrca);

      }
    }
  }
  }


}

model{

  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
  t_wgd ~ uniform(t_eca, t_mrca);

  for (m in 1:n_th_step_type){
      mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
    }

  for (ch in 1:n_th_cauchy_type){
      scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
    }

  omega ~ gamma(omega_alpha, omega_beta);

  // Likelihood
  m_clock_primary ~ poisson(2*omega*l_diploid*mu_clock*(t_mrca_primary-t_eca));
  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    alpha_tetraploid_step[th_type] ~ poisson(2 * l_tetraploid * omega * mu_th_step[th_type] * lambda_th_step_alpha[th_type]);
    beta_tetraploid_step[th_type] ~ poisson(4 * l_tetraploid * omega * mu_th_step[th_type] * lambda_th_step_beta[th_type]);

    // alpha_cnloh_step[th_type] ~ poisson(l_cnloh * omega * mu_th_step[th_type] * lambda_th_step_alpha[th_type]);
    // beta_cnloh_step[th_type] ~ poisson(2 * l_cnloh * omega * mu_th_step[th_type] * lambda_th_step_beta[th_type]);
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    alpha_tetraploid_cauchy[th_cauchy] ~ poisson(2 * l_tetraploid * omega * mu_clock * lambda_th_cauchy_alpha[th_cauchy]);
    beta_tetraploid_cauchy[th_cauchy] ~ poisson(4 * l_tetraploid * omega * mu_clock * lambda_th_cauchy_beta[th_cauchy]);

    // alpha_cnloh_cauchy[th_cauchy] ~ poisson(l_cnloh * omega * mu_clock * lambda_th_cauchy_alpha[th_cauchy]);
    // beta_cnloh_cauchy[th_cauchy] ~ poisson(2 * l_cnloh * omega * mu_clock * lambda_th_cauchy_beta[th_cauchy]);
  }
  }

  alpha_tetraploid_clock ~ poisson(2 * l_tetraploid * omega * mu_clock * (t_wgd - t_eca));
  beta_tetraploid_clock ~ poisson(4 * l_tetraploid * omega *  mu_clock * (t_mrca - t_wgd));

  // if (l_cnloh > 0){
    // alpha_cnloh_clock ~ poisson(l_cnloh * omega * mu_clock * (t_wgd - t_eca));
    // beta_cnloh_clock ~ poisson(2 * l_cnloh * omega *  mu_clock * (t_mrca - t_wgd));
  // }



  if (exponential_growth==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }

}

generated quantities{

  // Step therapy mutations
  array[n_th_step_type] int <lower=0> alpha_tetraploid_step_rep; // mutations in 2:2, for each therapy + clock
  array[n_th_step_type] int <lower=0> beta_tetraploid_step_rep;
  array[n_th_cauchy_type] int <lower=0> alpha_tetraploid_cauchy_rep;
  array[n_th_cauchy_type] int <lower=0> beta_tetraploid_cauchy_rep;

  // array[n_th_step_type] int <lower=0> alpha_cnloh_step_rep; // mutations in 2:2, for each therapy + clock
  // array[n_th_step_type] int <lower=0> beta_cnloh_step_rep;
  // array[n_th_cauchy_type] int <lower=0> alpha_cnloh_cauchy_rep;
  // array[n_th_cauchy_type] int <lower=0> beta_cnloh_cauchy_rep;

  int <lower=0> alpha_tetraploid_clock_rep;
  int <lower=0> beta_tetraploid_clock_rep;
  // int <lower=0> alpha_cnloh_clock_rep;
  // int <lower=0> beta_cnloh_clock_rep;

  int m_clock_primary_rep = poisson_rng(2*omega*l_diploid*mu_clock*(t_mrca_primary-t_eca));

  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    alpha_tetraploid_step_rep[th_type] = poisson_rng(2 * l_tetraploid * omega * mu_th_step[th_type] * lambda_th_step_alpha[th_type]);
    beta_tetraploid_step_rep[th_type] = poisson_rng(4 * l_tetraploid * omega * mu_th_step[th_type] * lambda_th_step_beta[th_type]);

    // alpha_cnloh_step_rep[th_type] = poisson_rng(l_cnloh * omega * mu_th_step[th_type] * lambda_th_step_alpha[th_type]);
    // beta_cnloh_step_rep[th_type] = poisson_rng(2 * l_cnloh * omega * mu_th_step[th_type] * lambda_th_step_beta[th_type]);
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    alpha_tetraploid_cauchy_rep[th_cauchy] = poisson_rng(2 * l_tetraploid * omega * mu_clock * lambda_th_cauchy_alpha[th_cauchy]);
    beta_tetraploid_cauchy_rep[th_cauchy] = poisson_rng(4 * l_tetraploid * omega * mu_clock * lambda_th_cauchy_beta[th_cauchy]);

    // alpha_cnloh_cauchy_rep[th_cauchy] = poisson_rng(l_cnloh * omega * mu_clock * lambda_th_cauchy_alpha[th_cauchy]);
    // beta_cnloh_cauchy_rep[th_cauchy] = poisson_rng(2 * l_cnloh * omega * mu_clock * lambda_th_cauchy_beta[th_cauchy]);
  }
  }

  alpha_tetraploid_clock_rep = poisson_rng(2 * l_tetraploid * omega * mu_clock * (t_wgd - t_eca));
  beta_tetraploid_clock_rep = poisson_rng(4 * l_tetraploid * omega *  mu_clock * (t_mrca - t_wgd));

  // if (l_cnloh > 0){
  //   alpha_cnloh_clock_rep = poisson_rng(l_cnloh * omega * mu_clock * (t_wgd - t_eca));
  //   beta_cnloh_clock_rep = poisson_rng(2 * l_cnloh * omega *  mu_clock * (t_mrca - t_wgd));
  // }


  if (exponential_growth==1){
    real N_sample_2 = exp(omega*(Sample_2-t_mrca));
    real N_sample_1 = exp(omega*(Sample_1-t_eca));
  }

}
