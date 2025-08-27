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
    real d = ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d * pi() * scale);
  
	}
	

  // Smooth max using log-sum-exp
  real smooth_max(real x, real y, real tau) {
    return tau * log_sum_exp(x / tau, y / tau);
  }

  // Smooth min via negated smooth max
  real smooth_min(real x, real y, real tau) {
    return -smooth_max(-x, -y, tau);
  }

  // Smooth positive part using softplus (≈ max(0, x))
  real smooth_positive(real x, real tau) {
    return tau * log1p_exp(x / tau);
  }

  real mutation_integral(
      real t_mrca,
      real f_min,
      real f_max,
      real omega,
      real mu,
      real L,
      real t1,
      real t2,
      real tau_smooth
  ) {
    real t_start = t_mrca + log1m(-f_max) / -omega; // = log(1/f_max)/omega
    real t_end   = t_mrca + log1m(-f_min) / -omega;

    real t_upper = smooth_min(t_end, t2, tau_smooth);
    real t_lower = smooth_max(t_start, t1, tau_smooth);

    real integral = exp(omega * t_upper) - exp(omega * t_lower);

    // Clip to non-negative using smooth approximation
    real clipped_integral = smooth_positive(integral, tau_smooth);

    return 2 * L * mu * clipped_integral;
  }


}

data{
  // Clock-like mutations
  int <lower=0> m_clock_primary;
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;
  real <lower=0> mu_driver_clock;
  real<lower=0, upper=1> f_min;
  real<lower=0, upper=1> f_max;



  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  array[cycles_driver] real<lower=0> driver_start; // if driver have effect only associated to external therapy
  array[cycles_driver] real<lower=0> driver_end;
  real <lower = 0> tau_driver;
  int <lower=0> m_driver;
  int <lower=0> m_tail_driver;
   // if driver alters basal mutation rate of clock-like

  // mutations associated to step-like therapies
  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  array[n_th_step] real<lower=0> start_th_step;
  array[n_th_step] real<lower=0> end_th_step;
  array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)
  array[n_th_step_type] real<lower=0> alpha_th_step;
  array[n_th_step_type] real<lower=0> beta_th_step;
  array[n_th_step_type] int<lower=0> m_th_step;
  array[n_th_step_type] real<lower=0> tau_step;
  array[n_th_step_type] int<lower=0> m_tail_step;
  array[n_th_step_type] real <lower=0> k_step_in_th;
  array[n_th_step_type] real <lower=0> k_step_out_th;

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  array[n_th_cauchy] real<lower=0> location_th_cauchy;
  array[n_th_cauchy] int<lower=0> type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy)
  array[n_th_cauchy_type] real<lower=0> alpha_th_cauchy;
  array[n_th_cauchy_type] real<lower=0> beta_th_cauchy;
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy;
  array[n_th_cauchy_type] int<lower=0> m_tail_cauchy;

  // other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step_in_driver;
  real <lower=0> k_step_out_driver;
  
  real <lower=0> tau_smooth;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> max_therapy;
  int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  real <lower=0> alpha_mrca;
  real <lower=0> beta_mrca;
  
  real <lower=0> mu_driver_alpha;
  real <lower=0> mu_driver_beta;
  
  // Overdispersion parameters for all mutation types
  real<lower=0> phi_clock;
  real<lower=0> phi_driver;
  array[n_th_step_type] real <lower=0> phi_th_step;
  array[n_th_cauchy_type] real <lower=0> phi_th_cauchy;
 
}


parameters{
  
  real <lower=0, upper=Sample_1> t_eca;
  real <lower=t_eca, upper=Sample_1> t_mrca_primary;
  real <lower=0, upper=1> rho_mrca;
   real <lower=0, upper=1> rho_driver;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  real <lower=0> omega;
  real <lower=0> mu_driver;
  
}


transformed parameters{

  real <lower=max_therapy,upper = Sample_2> t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);
  real <lower=t_eca,upper = t_mrca> t_driver = t_eca + rho_driver*(t_mrca-t_eca);
 

  array[n_th_step_type] real lambda_th_step;
  array[n_th_cauchy_type] real lambda_th_cauchy;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy[i] = 0;

// Step therapy mutations
if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {

         // Single block effect over entire therapy window
          lambda_th_step[th_type] += lambda_therapy_single(
            t_eca, t_mrca,
            start_th_step[cycle],
            end_th_step[cycle],
            k_step_in_th[th_type],
            k_step_out_th[th_type]
          );
        }

      }
    }
}



  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    for (cycle in 1:n_th_cauchy){
      if (type_th_cauchy[cycle] == th_cauchy){
        lambda_th_cauchy[th_cauchy] += couchy_cdf_single(location_th_cauchy[cycle], scales_th_cauchy[cycle], 
        t_eca, t_mrca);
      }
    }
  }
  }

 
   real lambda_driver=0;
    for (c in 1:cycles_driver){
    lambda_driver += lambda_therapy_single(t_driver,t_mrca, 
    driver_start[c],driver_end[c],k_step_in_driver,k_step_out_driver);
    }
  


// subclonal 

 array[n_th_step_type] real lambda_tail_th_step;
array[n_th_cauchy_type] real lambda_tail_th_cauchy;

// Initialize
for (i in 1:n_th_step_type)
  lambda_tail_th_step[i] = 0;


// Step therapy mutations
if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {
        // Add mutation integral for this therapy cycle
        lambda_tail_th_step[th_type] += mutation_integral(
          t_mrca,
          f_min,
          f_max,
          omega,
          mu_th_step[th_type],
          l_diploid,
          start_th_step[cycle],
          end_th_step[cycle],
          tau_smooth
        );
      }
    }
  }
}


   real lambda_tail_driver=0;
    for (c in 1:cycles_driver){
    lambda_tail_driver += mutation_integral(
          t_mrca,
          f_min,
          f_max,
          omega,
          mu_driver,
          l_diploid,
          driver_start[c],
          driver_end[c],
          tau_smooth
        );
    }

}


model {
  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
  rho_driver ~ uniform(0, 1);
  t_driver ~ uniform(t_eca, t_mrca);

  for (m in 1:n_th_step_type) {
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }
  for (ch in 1:n_th_cauchy_type) {
    scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
  }

  omega ~ gamma(omega_alpha, omega_beta);
  mu_driver ~ gamma(mu_driver_alpha, mu_driver_beta);

  // Negative Binomial shape parameters (inverse overdispersion)
  real shape_clock = 1 / phi_clock;
  real shape_driver = 1 / phi_driver;
  array[n_th_step_type] real shape_th_step;
  array[n_th_cauchy_type] real shape_th_cauchy;

  for (m in 1:n_th_step_type)
    shape_th_step[m] = 1 / phi_th_step[m];
  for (ch in 1:n_th_cauchy_type)
    shape_th_cauchy[ch] = 1 / phi_th_cauchy[ch];

  // Likelihood
  m_clock_primary ~ neg_binomial_2(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca),
    shape_clock
  );

  m_clock ~ neg_binomial_2(
    2 * l_diploid * omega * (mu_clock * (t_driver - t_eca) + mu_driver_clock * (t_mrca - t_driver)),
    shape_clock
  );

  for (th_type in 1:n_th_step_type) {
    m_th_step[th_type] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type],
      shape_th_step[th_type]
    );
     m_tail_step[th_type] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_tail_th_step[th_type] + 1,
      shape_th_step[th_type]
    );
  }

  for (th_cauchy in 1:n_th_cauchy_type) {
    m_th_cauchy[th_cauchy] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy],
      shape_th_cauchy[th_cauchy]
    );
    
    //  m_tail_cauchy[th_cauchy] ~ neg_binomial_2(
    //   2 * l_diploid * omega * mu_clock * lambda_tail_th_cauchy[th_cauchy] + 1,
    //   shape_th_cauchy[th_cauchy]
    // );
  }

  if (driver_type == 0) {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * (t_mrca - t_driver),
      shape_driver
    );
     m_tail_driver ~ neg_binomial_2(
      2 * l_diploid * mu_driver * (1/f_min - 1/f_max),
      shape_driver
    );
  } else {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * lambda_driver,
      shape_driver
    );
    m_tail_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * lambda_tail_driver + 1,
      shape_driver
    );
  
}

  if (exponential_growth == 1) {
   real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
   // real lambda2 = exp(-omega * (Sample_2 - t_mrca));
   real lambda2 = exp(-omega * (Sample_2 - t_driver));
  target += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
  // target += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
  // N_max[1] ~ exponential(lambda1);
  N_max[2] ~ exponential(lambda2);
}
  
}

generated quantities {
  
  int m_clock_primary_rep;
  int m_clock_rep;
  int m_driver_rep;
  int m_tail_driver_rep;
  array[n_th_step_type] int m_th_step_rep;
  array[n_th_step_type] int m_tail_step_rep;
  array[n_th_cauchy_type] int m_th_cauchy_rep;
  // array[n_th_cauchy_type] int m_tail_cauchy_rep;
  
  
  real shape_clock = 1 / phi_clock;
  real shape_driver = 1 / phi_driver;
  array[n_th_step_type] real shape_th_step;
  array[n_th_cauchy_type] real shape_th_cauchy;
  
   for (m in 1:n_th_step_type)
    shape_th_step[m] = 1 / phi_th_step[m];
  for (ch in 1:n_th_cauchy_type)
    shape_th_cauchy[ch] = 1 / phi_th_cauchy[ch];

  // Generate predictive draws using the same means as in the model block
  m_clock_primary_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca),
    shape_clock 
  );

  m_clock_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * (mu_clock * (t_driver - t_eca) + mu_driver_clock * (t_mrca - t_driver)),
   shape_clock 
  );

  for (th_type in 1:n_th_step_type) {
    m_th_step_rep[th_type] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type],
      shape_th_step[th_type]
    );

    m_tail_step_rep[th_type] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_tail_th_step[th_type],
      shape_th_step[th_type]
    );
  }

  for (th_cauchy in 1:n_th_cauchy_type) {
    m_th_cauchy_rep[th_cauchy] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy],
      shape_th_cauchy[th_cauchy]
    );

    // m_tail_cauchy_rep[th_cauchy] = neg_binomial_2_rng(
    //   2 * l_diploid * omega * mu_clock * lambda_tail_th_cauchy[th_cauchy],
    //   shape_th_cauchy[th_cauchy]
    // );
  }

  if (driver_type == 0) {
    m_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_driver * (t_mrca - t_driver),
      shape_driver 
    );
    m_tail_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * mu_driver * (1 / f_min - 1 / f_max),
      shape_driver 
    );
  } else {
    m_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_driver * lambda_driver,
    shape_driver 
    );
    m_tail_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_driver * lambda_tail_driver,
     shape_driver 
    );
  }
  
    real N_primary_rep = -1;
    real N_relapse_rep = -1;

   if (exponential_growth == 1) {
     
     
     // N_relapse_rep = exponential_rng(exp(-omega*(Sample_2 - t_mrca)));
     N_relapse_rep = exponential_rng(exp(-omega*(Sample_2 - t_driver)));
     N_primary_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
     
  }

}




