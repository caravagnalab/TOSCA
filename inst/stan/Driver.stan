functions {

//   real softplus(real x, real k) {
//     return log1p_exp(k * x) / k;
//   }
// 
//   real lambda_therapy_single(real ti, real tf, real t_therapy_i, real t_therapy_f, real k) {
//     real f1 = softplus(tf - t_therapy_i, k);
//     real f2 = softplus(tf - t_therapy_f, k);
//     real f3 = softplus(ti - t_therapy_i, k);
//     real f4 = softplus(ti - t_therapy_f, k);
// 
//     return f1 - f2 - f3 + f4;
//     
// }


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
  
  // // Smooth max using log-sum-exp
  // real smooth_max(real a, real b, real k) {
  //   return log_sum_exp(k * a, k * b) / k;
  // }
  // 
  // // Smooth min using negated log-sum-exp
  // real smooth_min(real a, real b, real k) {
  //   return -log_sum_exp(-k * a, -k * b) / k;
  // }
  // 
  // // Smooth version of max(x, 0)
  // real smooth_pos(real x, real k) {
  //   return smooth_max(x, 0, k);
  // }
  // 
  // // Smooth overlap between two time intervals: [th, tf] and [ts, te]
  // real lambda_therapy_single(real th, real tf, real ts, real te, real k) {
  //   real start = smooth_max(th, ts, k);
  //   real end_  = smooth_min(tf, te, k);
  //   return smooth_pos(end_ - start, k);
  // }


}

data{
  // Clock-like mutations
  int <lower=0> m_clock_primary;
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;


  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  array[cycles_driver] real<lower=0> driver_start; // if driver have effect only associated to external therapy
  array[cycles_driver] real<lower=0> driver_end;
  int <lower=0> m_driver;
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


  // other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> mu_driver_alpha;
  real <lower=0> mu_driver_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> min_mrca;
  real <lower=0> max_mrca;
  array[2] int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  real <lower=0> mrca_alpha;
  real <lower=0> mrca_beta;
  
  
  // Overdispersion parameters for all mutation types
  real<lower=0> phi_clock;
  real<lower=0> phi_driver;
  array[n_th_step_type] real <lower=0> phi_th_step;


}

transformed data{
  
 // Negative Binomial shape parameters (inverse overdispersion)
  real shape_clock = 1 / phi_clock;
  real shape_driver = 1 / phi_driver;
  array[n_th_step_type] real shape_th_step;

  for (m in 1:n_th_step_type)
    shape_th_step[m] = 1 / phi_th_step[m];
    
   
  
}


parameters{
  
  real <lower=0, upper=Sample_1> t_eca;
  real <lower=t_eca, upper=Sample_1> t_mrca_primary;
  real <lower=0, upper=1> rho_mrca;
  real <lower=0, upper=1> rho_driver;
  // real<lower=0, upper=1> mu_unit; 
  // real <lower=min_mrca, upper=Sample_2> t_mrca;
  // real <lower=t_eca, upper=t_mrca> t_driver;
  real <lower=0> mu_driver;
  array[n_th_step_type] real<lower=0> mu_th_step;
  real <lower=0> omega;
  
}

transformed parameters{

   real <lower= min_mrca ,upper = max_mrca> t_mrca = min_mrca + rho_mrca*(max_mrca-min_mrca);
   real <lower=t_eca,upper = t_mrca> t_driver = t_eca + rho_driver*(t_mrca -t_eca);
   
  
  array[n_th_step_type] real lambda_th_step;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;

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
            k_step
          );
        }

      }
    }
}



 real lambda_driver=0;
    for (c in 1:cycles_driver){
    lambda_driver += lambda_therapy_single(t_driver,t_mrca, 
    driver_start[c],driver_end[c],k_step);
    }
  


}


model {
  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  // t_mrca ~ uniform(min_mrca,Sample_2);
  // t_driver ~ uniform(t_eca,t_mrca);
  rho_mrca ~ beta(mrca_alpha, mrca_beta);
  rho_driver ~ uniform(0, 1);

  for (m in 1:n_th_step_type) {
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }

  omega ~ gamma(omega_alpha, omega_beta);
  
  mu_driver ~ gamma(mu_driver_alpha, mu_driver_beta);

 
  // Likelihood
  m_clock_primary ~ neg_binomial_2(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

  m_clock ~ neg_binomial_2(
    2 * l_diploid * omega * (mu_clock * (t_driver - t_eca) + mu_driver_clock * (t_mrca - t_driver)) + 0.1,
    shape_clock
  );
  
 
  for (th_type in 1:n_th_step_type) {
    m_th_step[th_type] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] + 0.1,
      shape_th_step[th_type]
    );
  }


  if (driver_type == 0) {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * (t_mrca - t_driver) + 0.1,
      shape_driver
    );
  } else {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * lambda_driver + 0.1,
      shape_driver
    );
}

  // if (exponential_growth == 1) {
  //   target += -N_min[1] * exp(-omega * (Sample_1 - t_mrca_primary)) +
  //     log(1 - exp(-(N_max[1] - N_min[1]) * exp(-omega * (Sample_1 - t_mrca_primary))));
  //   target += -N_min[2] * exp(-omega * (Sample_2 - t_mrca)) +
  //     log(1 - exp(-(N_max[2] - N_min[2]) * exp(-omega * (Sample_2 - t_mrca))));
  // }
  
  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    target += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
    // N_min[1] ~ exponential(lambda1);
    
  }
  
  if(exponential_growth[2] == 1){
    
    real lambda2 = exp(-omega * (Sample_2 - t_mrca));
    // N_min[2] ~ exponential(lambda2);
    target += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
    
  }
  
}

generated quantities {
  
  int<lower=0> m_clock_primary_rep;
  int<lower=0> m_clock_rep;
  int<lower=0> m_driver_rep;
  
  array[n_th_step_type] int<lower=0> m_th_step_rep;

 
 // Sample posterior predictive replicates from neg_binomial_2_rng
  m_clock_primary_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca),
    shape_clock
  );

  m_clock_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * (mu_clock * (t_driver - t_eca) + mu_driver_clock * (t_mrca - t_driver)),
    shape_clock
  );

  if (driver_type == 0) {
    m_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_driver * (t_mrca - t_driver),
      shape_driver
    );
  } else {
    m_driver_rep = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_driver * lambda_driver + 0.1,
      shape_driver
    );
  }


  for (m in 1:n_th_step_type) {
    m_th_step_rep[m] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_th_step[m] * lambda_th_step[m],
      shape_th_step[m]
    );
  }


  real N_primary_rep = -1;
  real N_relapse_rep = -1;

  //  if (exponential_growth == 1) {
  //    
  //    
  //    N_relapse_rep = exponential_rng(exp(-omega*(Sample_2 - t_mrca)));
  //    N_primary_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
  //    
  // }
  
   if (exponential_growth[2] == 1) {
     
     
     N_relapse_rep = exponential_rng(exp(-omega*(Sample_2 - t_mrca)));
   }
   
   if (exponential_growth[1] == 1) {
     
     N_primary_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
     
   }
     

}


