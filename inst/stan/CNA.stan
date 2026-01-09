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
  array[n_cna] int<lower=2,upper = 4> coeff;

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
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  array[2] int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  
  real<lower=0> phi_clock;
  array[n_th_step_type] real <lower=0> phi_th_step;
  array[n_cna] real<lower=0> phi_cna_alpha;
  array[n_cna] real<lower=0> phi_cna_beta;

}

parameters{
  
  real <lower=0, upper= Sample_1 > t_eca;
  real <lower=t_eca, upper = Sample_1 > t_mrca_primary;
  real <lower=t_eca, upper = Sample_2 > t_mrca;
  array[n_cna] real <lower = t_eca, upper = t_mrca > t_cna;
  array[n_th_step_type] real<lower=0> mu_th_step;
  real <lower=0> omega;
  
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
    real alpha = t_cna[c] - t_eca;
    real beta  = t_mrca - t_cna[c];

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
          lambda_th_step[th_type] += lambda_therapy_single(t_eca, t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);

          for (c in 1:n_cna) {
            real alpha = lambda_therapy_single(t_eca, t_cna[c], start_th_step[cycle], end_th_step[cycle], k_step);
            real beta  = lambda_therapy_single(t_cna[c], t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);
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

model {
  
  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  t_mrca ~ uniform(t_eca, Sample_2);
  
  for (c in 1:n_cna){
    t_cna[c] ~ uniform(t_eca, t_mrca);
  }
  
  for (m in 1:n_th_step_type){
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }
  
  omega ~ gamma(omega_alpha, omega_beta);
  
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
  
  
  // Likelihood: clock-like mutations
  m_clock_primary ~ neg_binomial_2(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

if(m_clock > 0){
  m_clock ~ neg_binomial_2(
    2 * l_diploid * omega * (
      mu_clock * (t_mrca - t_eca) 
    ) + 0.1,
    shape_clock
  );
}

  // Step therapy mutations
  if (n_th_step_type > 0) {
    for (th_type in 1:n_th_step_type) {
      m_th_step[th_type] ~ neg_binomial_2(
        2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] + 0.1,
        shape_th_step[th_type]
      );
    }
  }

  // CNA-associated mutations (alpha and beta)
  for (c in 1:n_cna) {
    m_alpha[c] ~ neg_binomial_2(
      lambda_alpha_clock[c] + lambda_alpha_th_step[c] +  0.1,
      shape_cna_alpha[c]
    );

    m_beta[c] ~ neg_binomial_2(
      lambda_beta_clock[c] + lambda_beta_th_step[c] +  0.1,
      shape_cna_beta[c]
    );
  }
  
  
    
  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    target += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
     // N_min[1] ~ exponential(lambda1);
    
  }
  
  if(exponential_growth[2] == 1){
    
    real lambda2 = exp(-omega * (Sample_2 - t_mrca));
    target += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
    // N_min[2] ~ exponential(lambda2);
    
  }
  
}

generated quantities {
  
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
  
  int<lower=0> m_clock_primary_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

  int<lower=0> m_clock_relapse_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * (
      mu_clock * (t_mrca - t_eca)
    ) + 0.1,
    shape_clock
  );

  array[n_cna] int<lower=0> m_alpha_rep;
  array[n_cna] int<lower=0> m_beta_rep;

  for (c in 1:n_cna) {
    real mu_alpha = lambda_alpha_clock[c] + lambda_alpha_th_step[c] +  0.1;
    real mu_beta  = lambda_beta_clock[c] + lambda_beta_th_step[c] +  0.1;

    m_alpha_rep[c] = neg_binomial_2_rng(mu_alpha, shape_cna_alpha[c]);
    m_beta_rep[c]  = neg_binomial_2_rng(mu_beta,  shape_cna_beta[c]);
  }

  real N_primary_rep = -1;
  real N_relapse_rep = -1;

  if (exponential_growth[1] == 1) {
    N_primary_rep = exponential_rng(exp(-omega * (Sample_1 - t_mrca_primary)));
  }

  if (exponential_growth[2] == 1) {
    N_relapse_rep = exponential_rng(exp(-omega * (Sample_2 - t_mrca)));
  }
  
  real log_lik = 0;

  // Clock-like mutations
  log_lik += neg_binomial_2_lpmf(
    m_clock_primary | 
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1, 
    shape_clock
  );

  if (m_clock > 0) {
    log_lik += neg_binomial_2_lpmf(
      m_clock | 
      2 * l_diploid * omega * mu_clock * (t_mrca - t_eca) + 0.1, 
      shape_clock
    );
  }

  // Therapy mutations
  if (n_th_step_type > 0) {
    for (th_type in 1:n_th_step_type) {
      log_lik += neg_binomial_2_lpmf(
        m_th_step[th_type] |
        2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] + 0.1,
        shape_th_step[th_type]
      );
    }
  }

  // CNA mutations
  for (c in 1:n_cna) {
    log_lik += neg_binomial_2_lpmf(
      m_alpha[c] | 
      lambda_alpha_clock[c] + lambda_alpha_th_step[c] + 0.1,
      shape_cna_alpha[c]
    );

    log_lik += neg_binomial_2_lpmf(
      m_beta[c] | 
      lambda_beta_clock[c] + lambda_beta_th_step[c] + 0.1,
      shape_cna_beta[c]
    );
  }

  // Exponential growth contribution
  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    log_lik += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));
  }

  if (exponential_growth[2] == 1) {
    real lambda2 = exp(-omega * (Sample_2 - t_mrca));
    log_lik += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
  }

  
}


