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

  int generate_th_rng(int n_th_step_type, int index, real l_diploid, real omega,
  array[] real mu_th_step, array[] real lambda_th_step, array[] real shape_th_step){

    int m_th_step=1;
    if (n_th_step_type > 0) {
      //for (th_type in 1:n_th_step_type) {
        m_th_step = neg_binomial_2_rng(
          2 * l_diploid * omega * mu_th_step[index] * lambda_th_step[index] + 0.1,
          shape_th_step[index]
      );
    //}
  }
   return m_th_step;

  }

  int generate_fac_rng(int fac_n_th_step, real l_diploid, real omega, real fac_mu_th_step, real lambda_fac, real shape_fac){
    int m_fac=1;
    if (fac_n_th_step > 0){
      m_fac = neg_binomial_2_rng(2 * l_diploid * omega * fac_mu_th_step * lambda_fac, shape_fac);
    }
    return m_fac;
  }

}

data {
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

  // Therapies
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
  real<lower=0> fac_phi_th_step;
  array[n_th_step_type] real <lower=0> phi_th_step;
  array[n_cna] real<lower=0> phi_cna_alpha;
  array[n_cna] real<lower=0> phi_cna_beta;

  // Data specific to this model
  real <lower=0> first_clinical_event; // earlier date between the end of the last cycle of the first therapy (between mutagenic and chemo), and the first sample
  real <lower=0> chemo_start; // == dormancy start
  real <lower=0> fac; // First event After Chemo could be chemo end OR the end of the last cycle of the first therapy after dormancy IF the chemo finishes after that
}

parameters {
  real <lower=0, upper= first_clinical_event > t_eca; // by default before dormancy
  real <lower=t_eca, upper = Sample_1 > t_mrca_primary; // this branch doesn't go through dormancy

  real <lower= chemo_start, upper= fac> t_dormancy_end;
  real <lower= t_dormancy_end, upper=Sample_2> t_mrca;
  array[n_cna] real <lower = t_eca, upper = t_mrca - (t_dormancy_end-chemo_start) > t_cna_tr;

  array[n_th_step_type] real<lower=0> mu_th_step;
  real <lower=0> omega;
}

transformed parameters {
  real t_mrca_tr = t_mrca - (t_dormancy_end-chemo_start);

  array[n_cna] real lambda_alpha_clock;
  array[n_cna] real lambda_beta_clock;
  array[n_th_step_type] real lambda_th_step;
  array[n_cna] real lambda_alpha_th_step;
  array[n_cna] real lambda_beta_th_step;

  for (c in 1:n_cna) {
    lambda_alpha_clock[c] = 0;
    lambda_beta_clock[c] = 0;
    lambda_alpha_th_step[c] = 0;
    lambda_beta_th_step[c] = 0;
  }
  for (th_type in 1:n_th_step_type) {
    lambda_th_step[th_type] = 0;
  }


  array[n_cna] real t_cna;
  for (c in 1:n_cna){
    t_cna[c] = traslation(t_cna_tr[c], chemo_start, t_dormancy_end);
    }

  // Clock-like alpha and beta
  for (c in 1:n_cna) {
    real scale = l_CNA[c] * omega * mu_clock;
    real alpha = t_cna_tr[c] - t_eca;
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

  // Therapies Pre and Post Dormancy
  if (n_th_step_type > 0) {
    for (th_type in 1:n_th_step_type) {
      for (cycle in 1:n_th_step) {
        if (type_th_step[cycle] == th_type) {

          // the cycle happens whithin the dormancy inteval
          if (start_th_step[cycle] > chemo_start && end_th_step[cycle] < t_dormancy_end) lambda_th_step[th_type] += 0;
          // the cycle starts before dormancy and ends during
          if (start_th_step[cycle] < chemo_start && end_th_step[cycle] > chemo_start && end_th_step[cycle] < t_dormancy_end) lambda_th_step[th_type] += chemo_start - start_th_step[cycle];
          // the cycle starts during dormancy and ends after
          if (start_th_step[cycle] > chemo_start && start_th_step[cycle] < t_dormancy_end && end_th_step[cycle] > t_dormancy_end) lambda_th_step[th_type] += end_th_step[cycle] - t_dormancy_end;
          // the cycle starts and ends outside of dormancy
          if (start_th_step[cycle] > t_dormancy_end || end_th_step[cycle] < chemo_start) lambda_th_step[th_type] += end_th_step[cycle] - start_th_step[cycle];

          for (c in 1:n_cna) {
            real alpha_step;
            real beta_step;
            // the cycle happens whithin the dormancy inteval
            if (start_th_step[cycle] > chemo_start && end_th_step[cycle] < t_dormancy_end) {
              alpha_step = 0;
              beta_step  = 0;
              }
            // the cycle starts before dormancy and ends during
            if (start_th_step[cycle] < chemo_start && end_th_step[cycle] > chemo_start && end_th_step[cycle] < t_dormancy_end) {
              alpha_step = lambda_therapy_single(t_eca, t_cna[c], end_th_step[cycle], chemo_start, k_step);
              beta_step  = lambda_therapy_single(t_cna[c], t_mrca, end_th_step[cycle], chemo_start,k_step);
              }
            // the cycle starts during dormancy and ends after
            if (start_th_step[cycle] > chemo_start && start_th_step[cycle] < t_dormancy_end && end_th_step[cycle] > t_dormancy_end) {
              alpha_step = lambda_therapy_single(t_eca, t_cna[c], t_dormancy_end, end_th_step[cycle], k_step);
              beta_step  = lambda_therapy_single(t_cna[c], t_mrca, t_dormancy_end, end_th_step[cycle],k_step);
              }
            // the cycle starts and ends outside of dormancy
            if (start_th_step[cycle] > t_dormancy_end || end_th_step[cycle] < chemo_start) {
              alpha_step = lambda_therapy_single(t_eca, t_cna[c], start_th_step[cycle], end_th_step[cycle], k_step);
              beta_step  = lambda_therapy_single(t_cna[c], t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);
              }

            real scale_step = l_CNA[c] * omega * mu_th_step[th_type];

            if (coeff[c] == 2) {
              lambda_alpha_th_step[c] += 1.0 * scale_step * alpha_step;
              lambda_beta_th_step[c]  += 2.0 * scale_step * beta_step;
            } else if (coeff[c] == 3) {
              lambda_alpha_th_step[c] += 1.0 * scale_step * alpha_step;
              lambda_beta_th_step[c]  += (3.0 * beta_step + 1.0 * alpha_step) * scale_step;
            } else {
              lambda_alpha_th_step[c] += 2.0 * scale_step * alpha_step;
              lambda_beta_th_step[c]  += 4.0 * scale_step * beta_step;
            }
          }
          }

        }
      }
    }

}

model {
  // *** Prior ***
  t_eca ~ uniform(0, first_clinical_event); // by default before dormancy
  t_mrca_primary ~ uniform(t_eca, Sample_1); // this branch doesn't go through dormancy
  t_dormancy_end ~ uniform(chemo_start, fac);
  t_mrca ~ uniform(t_dormancy_end, Sample_2);
  for (c in 1:n_cna){
    t_cna_tr[c] ~ uniform(t_eca, t_mrca_tr);
  }
  omega ~ gamma(omega_alpha, omega_beta);
  for (m in 1:n_th_step_type){
    mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  }

  // *** Likelihood ***

  // Clock-like
  real shape_clock = 1 / phi_clock;
  m_clock_primary ~ neg_binomial_2(2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1,shape_clock);
  m_clock ~  neg_binomial_2(2*mu_clock*omega*l_diploid*(t_mrca_tr - t_eca) +.1, shape_clock);

  // Therapies Pre and Post Dormancy
  array[n_th_step_type] real shape_th_step;
  for (m in 1:n_th_step_type)
  shape_th_step[m] = 1 / phi_th_step[m];
  if (n_th_step_type > 0) {
    for (th_type in 1:n_th_step_type) {
      m_th_step[th_type] ~ neg_binomial_2(
        2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type] + 0.1,
        shape_th_step[th_type]
      );
    }
  }


  // Copy Numbers
  array[n_cna] real shape_cna_alpha;
  array[n_cna] real shape_cna_beta;
  for (c in 1:n_cna) {
    shape_cna_alpha[c] = 1 / phi_cna_alpha[c];
    shape_cna_beta[c]  = 1 / phi_cna_beta[c];

    m_alpha[c] ~ neg_binomial_2(lambda_alpha_clock[c] + lambda_alpha_th_step[c] + 0.1, shape_cna_alpha[c]);
    m_beta[c] ~ neg_binomial_2(lambda_beta_clock[c] + lambda_beta_th_step[c] + 0.1, shape_cna_beta[c]);
  }

  // *** Penalty ***
  if (exponential_growth[1] == 1) {
    real lambda1 = exp(-omega * (Sample_1 - t_mrca_primary));
    target += -lambda1 * N_min[1] + log1m_exp(-lambda1 * (N_max[1] - N_min[1]));

  }
  if(exponential_growth[2] == 1){
    real lambda2 = exp(-omega * (Sample_2 - t_mrca));
    target += -lambda2 * N_min[2] + log1m_exp(-lambda2 * (N_max[2] - N_min[2]));
  }

}

generated quantities {

  real shape_clock_rep = 1 / phi_clock;
  array[n_th_step_type] real shape_th_step_rep;
  for (m in 1:n_th_step_type)
  shape_th_step_rep[m] = 1 / phi_th_step[m];
  real shape_fac_rep = 1 / fac_phi_th_step;
  array[n_cna] real shape_cna_alpha_rep;
  array[n_cna] real shape_cna_beta_rep;
  for (c in 1:n_cna){
    shape_cna_alpha_rep[c] = 1 / phi_cna_alpha[c];
    shape_cna_beta_rep[c]  = 1 / phi_cna_beta[c];
  }

  // array[n_cna] real <lower = t_eca,upper = t_mrca> t_cna;

  int m_clock_primary_rep = neg_binomial_2_rng(2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1,shape_clock_rep);
  int m_clock_rep =  neg_binomial_2_rng(2*mu_clock*omega*l_diploid*(t_mrca_tr - t_eca) +.1, shape_clock_rep);
  array[n_th_step_type] int m_th_step_rep;
  for (t in 1:n_th_step_type){
    m_th_step_rep[t] = generate_th_rng(n_th_step_type, t, l_diploid, omega, mu_th_step, lambda_th_step, shape_th_step_rep);
  }
  array[n_cna] int m_alpha_rep;
  array[n_cna] int m_beta_rep;
  for (c in 1:n_cna){
    m_alpha_rep[c] = neg_binomial_2_rng(lambda_alpha_clock[c] + lambda_alpha_th_step[c] + 0.1, shape_cna_alpha_rep[c]);
    m_beta_rep[c] = neg_binomial_2_rng(lambda_beta_clock[c] + lambda_beta_th_step[c] + 0.1, shape_cna_beta_rep[c]);
  }

}
