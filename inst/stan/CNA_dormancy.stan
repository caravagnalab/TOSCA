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

real delta(real t, real dorm_s, real dorm_e, real k) {
  real shift_weight = inv_logit(-k * (t - dorm_e));
  return t - (t - dorm_e) * shift_weight - dorm_s;
}

 real traslation(real t, real dorm_start, real dorm_end) {
    real t_cn = t;
    if (t > dorm_start) t_cn = t_cn + dorm_end - dorm_start;
    return t_cn;
  }
  

real smooth_mrca_weight(real t_mrca_tr,
real chemo_start,
real omega,
real k) {
  
  // Smooth transition weight using logistic
  real w = inv_logit(k * (t_mrca_tr - chemo_start));  // goes from 0 to 1
  
  real r = exp(-omega*(chemo_start - t_mrca_tr));
  
  // Smooth interpolation: r_mean when t_mrca_tr << chemo_start, 1 when t_mrca_tr >> chemo_start
  return (1 - w) * r + w * 1;
  
}


}

data{

 int <lower=1> n_cna;
 int <lower=0> m_clock_primary;
 int <lower=0> m_clock;
 real <lower=0> l_diploid;
 real <lower=0> mu_clock;

 array[n_cna] int <lower=0> m_alpha;
 array[n_cna] int <lower=0> m_beta;
 array[n_cna] real <lower=0> l_CNA;
 array[n_cna] int <lower=2,upper = 4> coeff;

  int <lower=0> n_th_step; // numero totale di terapie*cicli
  int <lower=0> n_th_step_type; // numero di tipi di terapia
  array[n_th_step] real<lower=0> start_th_step;
  array[n_th_step] real<lower=0> end_th_step;
  array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)
  array[n_th_step_type] real<lower=0> alpha_th_step;
  array[n_th_step_type] real<lower=0> beta_th_step;
  array[n_th_step_type] int<lower=0> m_th_step;


	real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> chemo_start;
  real <lower=0> chemo_end;

  array[2] int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  // Overdispersion parameters for all mutation types
   real<lower=0> phi_clock;
   array[n_th_step_type] real <lower=0> phi_th_step;
   array[n_cna] real<lower=0> phi_cna_alpha;
   array[n_cna] real<lower=0> phi_cna_beta;
  
}


parameters{

  real <lower=0, upper= Sample_1> t_eca;
  real <lower=t_eca, upper= Sample_1> t_mrca_primary;
  real <lower= chemo_end, upper= Sample_2> t_dormancy_end;
  // real <lower= t_dormancy_end, upper= Sample_2> t_mrca;
  real <lower= t_eca, upper= Sample_2> t_mrca_tr;
  // array[n_cna] real <lower= t_eca, upper= t_mrca - (t_dormancy_end - chemo_start)> t_cna_tr;
  array[n_cna] real <lower= t_eca, upper= t_mrca_tr> t_cna_tr;
  real <lower= 0> omega;
  array[n_th_step_type] real<lower=0> mu_th_step;
  
}

transformed parameters {
  
  array[n_th_step_type]      real lambda_th_step;

  array[n_cna] real lambda_alpha_clock;
  array[n_cna] real lambda_beta_clock;
  array[n_cna] real lambda_alpha_th_step;
  array[n_cna] real lambda_beta_th_step;

  // Initialize all lambda values to 0
  for (i in 1:n_th_step_type)     lambda_th_step[i] = 0;
  for (i in 1:n_cna) {
    lambda_alpha_clock[i]       = 0;
    lambda_beta_clock[i]        = 0;
    lambda_alpha_th_step[i]     = 0;
    lambda_beta_th_step[i]      = 0;
  }

  // CLOCK-LIKE mutations
  for (c in 1:n_cna) {
    real scale = l_CNA[c] * omega * mu_clock;
    real alpha = t_cna_tr[c] - t_eca;
    // real beta  = t_mrca-(t_dormancy_end - chemo_start) - t_cna_tr[c];
    real beta  = t_mrca_tr - t_cna_tr[c];

    if (coeff[c] == 2) {
      lambda_alpha_clock[c] += 1.0 * scale * alpha;
      lambda_beta_clock[c]  += 2.0 * scale * beta;
    } else if (coeff[c] == 3) {
      lambda_alpha_clock[c] += 1.0 * scale * alpha;
      lambda_beta_clock[c]  += 3.0 * scale * beta + 1.0 * scale * alpha;
    } else {
      lambda_alpha_clock[c] += 2.0 * scale * alpha;
      lambda_beta_clock[c]  += 4.0 * scale * beta;
    }
  }

  // STEP THERAPY mutations
  if (n_th_step_type > 0) {
    for (th_type in 1:n_th_step_type) {
      for (cycle in 1:n_th_step) {
        if (type_th_step[cycle] == th_type) {
          
          // lambda_th_step[th_type] += lambda_therapy_single(t_dormancy_end, t_mrca, start_th_step[cycle], 
          //     end_th_step[cycle], k_step);
          
          lambda_th_step[th_type] += lambda_therapy_single(t_eca, t_mrca_tr, 
                   start_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
                   end_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
                   k_step);
              
  
   for (c in 1:n_cna) {
     
            real alpha = lambda_therapy_single(t_eca, t_cna_tr[c], 
            start_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
            end_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
                           k_step);
            real beta  = lambda_therapy_single(t_cna_tr[c],
            // t_mrca-(t_dormancy_end - chemo_start), 
            t_mrca_tr,
            start_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
            end_th_step[cycle] - delta(start_th_step[cycle], chemo_start, t_dormancy_end, k_step), 
            k_step);
            
            real scale = l_CNA[c] * omega * mu_th_step[th_type];

            if (coeff[c] == 2) {
              lambda_alpha_th_step[c] += 1.0 * scale * alpha;
              lambda_beta_th_step[c]  += 2.0 * scale * beta;
            } else if (coeff[c] == 3) {
              lambda_alpha_th_step[c] += 1.0 * scale * alpha;
              lambda_beta_th_step[c]  += (3.0 * beta + 1.0 * alpha) * scale;
            } else {
              lambda_alpha_th_step[c] += 2.0 * scale * alpha;
              lambda_beta_th_step[c]  += 4.0 * scale * beta;
            }
          }
        }
      }
    }
  }


 }

model{
  
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  t_dormancy_end ~ uniform(chemo_end, Sample_2);
  // t_mrca ~ uniform(t_dormancy_end, Sample_2);
  t_mrca_tr ~ uniform(t_eca, Sample_2);

  omega ~ gamma(omega_alpha,omega_beta);
  
    for (m in 1:n_th_step_type){
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }
  

   // Global shapes
  real shape_clock    = 1 / phi_clock;
  
  // Therapy shapes
  array[n_th_step_type]       real shape_th_step;
  
  // CNA mutation shapes
  array[n_cna] real shape_cna_alpha;
  array[n_cna] real shape_cna_beta;
  
  // Compute therapy shapes
  for (m in 1:n_th_step_type)
  shape_th_step[m] = 1 / phi_th_step[m];
  
 
  // Compute CNA shapes
  for (c in 1:n_cna) {
    shape_cna_alpha[c] = 1 / phi_cna_alpha[c];
    shape_cna_beta[c]  = 1 / phi_cna_beta[c];
  }
  

m_clock_primary ~ neg_binomial_2(
  2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1,
  shape_clock
  );

if (m_clock > 0){
  m_clock ~  neg_binomial_2(
    // 2*mu_clock*omega*l_diploid*((chemo_start - t_eca) + (t_mrca-t_dormancy_end)) +.1,
     2*mu_clock*omega*l_diploid*(t_mrca_tr - t_eca) +.1,
    shape_clock
    );
}

    if (n_th_step_type > 0) {
        for (th_type in 1:n_th_step_type) {
          m_th_step[th_type] ~ neg_binomial_2(
            2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] + 0.1,
            shape_th_step[th_type]
            );
        }
      }


  for (c in 1:n_cna) {
    
   // t_cna_tr[c] ~ uniform(t_eca, t_mrca - (t_dormancy_end - chemo_start) );
   t_cna_tr[c] ~ uniform(t_eca, t_mrca_tr);
  
    m_alpha[c] ~ neg_binomial_2(
      lambda_alpha_clock[c] + lambda_alpha_th_step[c] + 0.1,
      shape_cna_alpha[c]
    );

    m_beta[c] ~ neg_binomial_2(
      lambda_beta_clock[c] + lambda_beta_th_step[c] + 0.1,
      shape_cna_beta[c]
    );
    
}
  
  
  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    target += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
    
  }
  
  if(exponential_growth[2] == 1){
    
    // real lambda2 = exp(-omega * (Sample_2 - t_mrca));
    // real lambda2 = exp(-omega * (Sample_2 - (t_dormancy_end - chemo_start) - t_mrca_tr));
    real kappa = smooth_mrca_weight(t_mrca_tr, chemo_start, omega,k_step);
    real lambda2 = exp(-omega * (Sample_2 - (t_dormancy_end - chemo_start) - t_mrca_tr))/kappa;
    target += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
    
  }
  
}

generated quantities {
  
  real <lower= t_eca, upper= Sample_2> t_mrca;
  
  t_mrca = traslation(t_mrca_tr, chemo_start, t_dormancy_end);
  
  array[n_cna] real <lower = t_eca,upper = t_mrca> t_cna;
  
  for (c in 1:n_cna){
    
   t_cna[c] = traslation(t_cna_tr[c], chemo_start, t_dormancy_end);

}

  
  // Posterior predictive distributions
  int m_clock_primary_rep;
  int m_clock_rep;
  array[n_th_step_type] int m_th_step_rep;
  array[n_cna] int m_alpha_rep;
  array[n_cna] int m_beta_rep;

  // Log-likelihoods for model evaluation (WAIC/LOO)
  real log_lik_m_clock_primary;
  real log_lik_m_clock;
  array[n_th_step_type] real log_lik_th_step;
  array[n_cna] real log_lik_alpha;
  array[n_cna] real log_lik_beta;

  // Derived quantities
  real dormancy_duration = t_dormancy_end - chemo_start;
 
  // Shape parameters
  real shape_clock = 1 / phi_clock;
  array[n_th_step_type] real shape_th_step;
  array[n_cna] real shape_cna_alpha;
  array[n_cna] real shape_cna_beta;

  for (m in 1:n_th_step_type)
    shape_th_step[m] = 1 / phi_th_step[m];

  for (c in 1:n_cna) {
    shape_cna_alpha[c] = 1 / phi_cna_alpha[c];
    shape_cna_beta[c]  = 1 / phi_cna_beta[c];
  }

  // Predict clock-like mutations
  m_clock_primary_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

  log_lik_m_clock_primary = neg_binomial_2_lpmf(
    m_clock_primary | 
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

  if (m_clock > 0) {
    m_clock_rep = neg_binomial_2_rng(
      // 2 * mu_clock * omega * l_diploid * ((chemo_start - t_eca) + (t_mrca - t_dormancy_end)) + 0.1,
      2*mu_clock*omega*l_diploid*(t_mrca_tr - t_eca) +.1,
      shape_clock
    );

    log_lik_m_clock = neg_binomial_2_lpmf(
      m_clock | 
      // 2 * mu_clock * omega * l_diploid * ((chemo_start - t_eca) + (t_mrca - t_dormancy_end)) + 0.1,
      2*mu_clock*omega*l_diploid*(t_mrca_tr - t_eca) +.1,
      shape_clock
    );
  } else {
    m_clock_rep = 0;
    log_lik_m_clock = 0;
  }

  // Predict therapy-related mutations and log-likelihoods
  for (th in 1:n_th_step_type) {
    real mean_th = 2 * l_diploid * omega * mu_th_step[th] * lambda_th_step[th] + 0.1;
    m_th_step_rep[th] = neg_binomial_2_rng(mean_th, shape_th_step[th]);
    log_lik_th_step[th] = neg_binomial_2_lpmf(m_th_step[th] | mean_th, shape_th_step[th]);
  }

  // Predict CNA-specific mutation counts and log-likelihoods
  for (c in 1:n_cna) {
    real mean_alpha = lambda_alpha_clock[c] + lambda_alpha_th_step[c] + 0.1;
    real mean_beta  = lambda_beta_clock[c] + lambda_beta_th_step[c] + 0.1;

    m_alpha_rep[c] = neg_binomial_2_rng(mean_alpha, shape_cna_alpha[c]);
    m_beta_rep[c]  = neg_binomial_2_rng(mean_beta,  shape_cna_beta[c]);

    log_lik_alpha[c] = neg_binomial_2_lpmf(m_alpha[c] | mean_alpha, shape_cna_alpha[c]);
    log_lik_beta[c]  = neg_binomial_2_lpmf(m_beta[c]  | mean_beta,  shape_cna_beta[c]);

  }
  
   real N_primary_rep = -1;
   real N_relapse_rep = -1;
  
  if (exponential_growth[2] == 1) {
    
     real kappa = smooth_mrca_weight(t_mrca_tr, chemo_start, omega, k_step);
     real lambda2 = exp(-omega * (Sample_2 - (t_dormancy_end - chemo_start) - t_mrca_tr))/kappa;
     
     N_relapse_rep = exponential_rng(lambda2);
  }
  
  if (exponential_growth[1] == 1) {
    
    N_primary_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
    
  }
  
  // Log-likelihood terms for branching process (exponential growth)
  array[2] real log_lik_branching;
  log_lik_branching[1] = 0;
  log_lik_branching[2] = 0;

  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    log_lik_branching[1] = -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
  }

  if (exponential_growth[2] == 1) {
     real kappa = smooth_mrca_weight(t_mrca_tr, chemo_start, omega, k_step);
     real lambda2 = exp(-omega * (Sample_2 - (t_dormancy_end - chemo_start) - t_mrca_tr))/kappa;
     log_lik_branching[2] = -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
  }
  
  
   
  real total_log_lik;
  
   // Total log-likelihood sum
  total_log_lik =
    log_lik_m_clock_primary +
    log_lik_m_clock +
    sum(log_lik_th_step) +
    sum(log_lik_alpha) +
    sum(log_lik_beta) +
    log_lik_branching[1] +
    log_lik_branching[2];


}



