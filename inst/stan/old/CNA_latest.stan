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
  int <lower=0> m_clock_primary;
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;

  // mutations on CNA
  int <lower=0> n_cna;
  array[n_cna] int<lower=0> m_alpha;
  array[n_cna] int<lower=0> m_beta;
  array[n_cna] real<lower=0> l_CNA;
  array[n_cna] int<lower=0> coeff_alpha; // 1 for 2:0, 2 for 2:2
  array[n_cna] int<lower=0> coeff_beta; // 2 for 2:0, 4 for 2:2

  // mutations associated to step-like therapies
  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  array[n_th_step] real<lower=0> start_th_step;
  array[n_th_step] real<lower=0> end_th_step;
  array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)
  // array[n_th_step_type] real<lower=0> alpha_th_step;
  // array[n_th_step_type] real<lower=0> beta_th_step;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_step_type] int<lower=0> m_th_step;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  array[n_th_cauchy] real<lower=0> location_th_cauchy;
  array[n_th_cauchy] int<lower=0> type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy)
  // array[n_th_cauchy_type] real<lower=0> alpha_th_cauchy;
  // array[n_th_cauchy_type] real<lower=0> beta_th_cauchy;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy;

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

  // Overdispersion parameters for all mutation types
  real<lower=0> phi_clock;
  array[n_th_step_type] real <lower=0> phi_th_step;
  array[n_th_cauchy_type] real <lower=0> phi_th_cauchy;
  array[n_cna] real <lower=0> phi_cna;

}

parameters{
  real <lower=0, upper=Sample_1 - 1e-8 > t_eca;
  real <lower=t_eca, upper=Sample_1 > t_mrca_primary;
  array[n_cna] real<lower=t_eca, upper=Sample_2 - 1e-8 > t_cna;
  real <lower=0, upper=1> rho_mrca;
  // array[n_th_step_type] real<lower=0> mu_th_step;
  // array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  real <lower=0> omega;
}

transformed parameters{

  real <lower=max_therapy> t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);

  array[n_th_step_type] real lambda_th_step;
  array[n_th_cauchy_type] real lambda_th_cauchy;

  array[n_cna] real lambda_alpha_clock;
  array[n_cna] real lambda_beta_clock;
  array[n_cna] real lambda_alpha_th_step;
  array[n_cna] real lambda_beta_th_step;
  array[n_cna] real lambda_alpha_th_cauchy;
  array[n_cna] real lambda_beta_th_cauchy;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy[i] = 0;

  for (i in 1:n_cna) {
    lambda_alpha_clock[i] = 0;
    lambda_beta_clock[i] = 0;

    lambda_alpha_th_step[i] = 0;
    lambda_beta_th_step[i] = 0;

    lambda_alpha_th_cauchy[i] = 0;
    lambda_beta_th_cauchy[i] = 0;
  }

  // clock-like rate
  for (c in 1:n_cna){
    lambda_alpha_clock[c] += coeff_alpha[c] * l_CNA[c] * omega * mu_clock * (t_cna[c] - t_eca);
    lambda_beta_clock[c] += coeff_beta[c] * l_CNA[c] * omega * mu_clock * (t_mrca - t_cna[c]);
  }

  // Step therapy mutations
  if (n_th_step_type > 0){
    for (th_type in 1:n_th_step_type){
      for (cycle in 1:n_th_step){

        if (type_th_step[cycle] == th_type){

          lambda_th_step[th_type] += lambda_therapy_single(t_eca, t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);

          // Update lambda CNA
          for (c in 1:n_cna){
            lambda_alpha_th_step[c] += coeff_alpha[c] * l_CNA[c] * omega * mu_th_step[th_type] * lambda_therapy_single(t_eca, t_cna[c], start_th_step[cycle], end_th_step[cycle], k_step);
            lambda_beta_th_step[c] += coeff_beta[c] * l_CNA[c] * omega * mu_th_step[th_type] * lambda_therapy_single(t_cna[c], t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);
          }

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
        for (c in 1:n_cna){
          lambda_alpha_th_cauchy[c] += coeff_alpha[c] * l_CNA[c] * omega * mu_clock * couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], t_eca, t_cna[c]);
          lambda_beta_th_cauchy[c] += coeff_beta[c] * l_CNA[c] * omega * mu_clock * couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], t_cna[c], t_mrca);
        }
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

  for (c in 1:n_cna){
    t_cna[c] ~ uniform(t_eca, t_mrca);
  }

  // Negative Binomial shape parameters (inverse overdispersion)
  real shape_clock = 1 / phi_clock;
  array[n_th_step_type] real shape_th_step;
  array[n_th_cauchy_type] real shape_th_cauchy;
  array[n_cna] real shape_cna;

  for (m in 1:n_th_step_type)
    shape_th_step[m] = 1 / phi_th_step[m];
  for (ch in 1:n_th_cauchy_type)
    shape_th_cauchy[ch] = 1 / phi_th_cauchy[ch];
  for (cna in 1:n_cna)
    shape_cna[cna] = 1 / phi_cna[cna];

  // for (m in 1:n_th_step_type){
  //     mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  //   }
  //
  // for (ch in 1:n_th_cauchy_type){
  //     scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
  //   }


  omega ~ gamma(omega_alpha, omega_beta);

  // Likelihood
  m_clock_primary ~ neg_binomial_2(
    2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1,
    shape_clock
    );
  m_clock ~ neg_binomial_2(
    2*l_diploid*omega*mu_clock*(t_mrca-t_eca) +.1,
    shape_clock
    );

  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    m_th_step[th_type] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] +.1,
      shape_th_step[th_type]
      );
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    m_th_cauchy[th_cauchy] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy] +.1,
      shape_th_cauchy[th_cauchy]
    );
  }
  }

  for (c in 1:n_cna){
    m_alpha[c] ~ neg_binomial_2(lambda_alpha_clock[c] + lambda_alpha_th_step[c] + lambda_alpha_th_cauchy[c] +.1,
                                shape_cna[c]);
    m_beta[c] ~ neg_binomial_2(lambda_beta_clock[c] + lambda_beta_th_step[c] + lambda_beta_th_cauchy[c] +.1,
                               shape_cna[c]);
  }

  if (exponential_growth==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }

}

generated quantities{

  // Negative Binomial shape parameters (inverse overdispersion)
  real shape_clock_rep = 1 / phi_clock;
  array[n_th_step_type] real shape_th_step_rep;
  array[n_th_cauchy_type] real shape_th_cauchy_rep;
  array[n_cna] real shape_cna_rep;

  for (m in 1:n_th_step_type)
    shape_th_step_rep[m] = 1 / phi_th_step[m];
  for (ch in 1:n_th_cauchy_type)
    shape_th_cauchy_rep[ch] = 1 / phi_th_cauchy[ch];
  for (cna in 1:n_cna)
    shape_cna_rep[cna] = 1 / phi_cna[cna];

  int<lower =0> m_clock_primary_rep = neg_binomial_2_rng(2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1, shape_clock_rep);
  int<lower =0> m_clock_rep = neg_binomial_2_rng(2*l_diploid*omega*mu_clock*(t_mrca-t_eca) +.1, shape_clock_rep);

  array[n_cna] int<lower=0> m_alpha_rep;
  array[n_cna] int<lower=0> m_beta_rep;
  for (c in 1:n_cna){
    m_alpha_rep[c] = neg_binomial_2_rng(lambda_alpha_clock[c] + lambda_alpha_th_step[c] + lambda_alpha_th_cauchy[c] +.1, shape_cna_rep[c]);
    m_beta_rep[c] = neg_binomial_2_rng(lambda_beta_clock[c] + lambda_beta_th_step[c] + lambda_beta_th_cauchy[c] +.1, shape_cna_rep[c]);
  }

  // Step therapy mutations
  array[n_th_step_type] int<lower=0> m_th_step_rep;
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    m_th_step_rep[th_type] = neg_binomial_2_rng(2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] +.1, shape_th_step_rep[th_type]);
  }
  }

  // Cauchy therapy mutations
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy_rep;
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    m_th_cauchy_rep[th_cauchy] = neg_binomial_2_rng(2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy] +.1, shape_th_cauchy_rep[th_cauchy]);
  }
  }


  if (exponential_growth==1){
    real N_sample_2 = exp(omega*(Sample_2-t_mrca));
    real N_sample_1 = exp(omega*(Sample_1-t_eca));
  }

}
