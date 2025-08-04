functions {
  real lambda(real ti, real tf, real t_change, real sample, real omega, real lambda_pre, real mu, real k) {
    real a1 = k * (tf);
    real a2 = k * (tf - t_change);
    real a3 = k * (ti);
    real a4 = k * (ti - t_change);

    real f1 = log_sum_exp(0, a1) / k;
    real f2 = log_sum_exp(0, a2) / k;
    real f3 = log_sum_exp(0, a3) / k;
    real f4 = log_sum_exp(0, a4) / k;

    real pre = (f1 - f2 - f3 + f4) * lambda_pre;

    real b1 = k * (tf - t_change);
    real b2 = k * (tf - sample);
    real b3 = k * (ti - t_change);
    real b4 = k * (ti - sample);

    real g1 = log_sum_exp(0, b1) / k;
    real g2 = log_sum_exp(0, b2) / k;
    real g3 = log_sum_exp(0, b3) / k;
    real g4 = log_sum_exp(0, b4) / k;

    real post = (g1 - g2 - g3 + g4) * omega * mu;

    return pre + post;
  }
}

data {
  // Clock-like mutations
  int<lower=0> m_clock;
  real<lower=0> l_diploid;
  real<lower=0> mu_post;
  real<lower=0> lambda_pre;
  real<lower=0> alpha_omega;
  real<lower=0> beta_omega;
  real<lower=0> alpha_mrca;
  real<lower=0> beta_mrca;
  real<lower=0> k_step;

  int<lower=1> n_events;
  array[n_events] int<lower=0> copy_number_type;
  array[n_events] real<lower=0> alpha;
  array[n_events] real<lower=0> beta;
  array[n_events] real<lower=0> l_copy_number;

  real<lower=0> Sample;
  real<lower=0> phi;
  int<lower=0,upper=1> exponential_growth;

  // Added for exponential growth model
  real<lower=0> N_min;
  real<lower=0> N_max;
}

parameters {
  real<lower=0, upper=1> rho_mrca;
  real<lower=0, upper=Sample> t_driver;
  array[n_events] real<lower=0> t_copy_number;
  real<lower=0> omega;
}

transformed parameters {
  real<lower=0, upper=Sample> t_mrca;
  t_mrca = rho_mrca * Sample;
}

model {
  // Priors
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
  omega ~ gamma(alpha_omega, beta_omega);

  real shape = 1 / phi;

  // Likelihood: clock-like mutations
  m_clock ~ neg_binomial_2(
    2 * l_diploid * (omega * mu_post * (t_mrca - t_driver) + lambda_pre * t_driver),
    shape
  );

  // Likelihood: mutations in copy number segments
  for (i in 1:n_events) {
    if (copy_number_type[i] == 0) {
      alpha[i] ~ neg_binomial_2(
        l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );
      beta[i] ~ neg_binomial_2(
        2 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );

    } else if (copy_number_type[i] == 1) {
      alpha[i] ~ neg_binomial_2(
        l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );
      beta[i] ~ neg_binomial_2(
        l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step) +
          3 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );

    } else {
      alpha[i] ~ neg_binomial_2(
        2 * l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );
      beta[i] ~ neg_binomial_2(
        4 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step),
        shape
      );
    }
  }

  if (exponential_growth == 1) {
    target += -N_min * exp(-omega * (Sample - t_mrca)) +
      log(1 - exp(-(N_max - N_min) * exp(-omega * (Sample - t_mrca))));
  }
}

generated quantities {
  int m_clock_rep;
  array[n_events] int alpha_rep;
  array[n_events] int beta_rep;
  
  real shape = 1 / phi;

  // Replicate m_clock from same neg_binomial_2 model
  m_clock_rep = neg_binomial_2_rng(
    2 * l_diploid * (omega * mu_post * (t_mrca - t_driver) + lambda_pre * t_driver),
    shape
  );

  // Replicate alpha and beta counts for each event
  for (i in 1:n_events) {
    real lambda_alpha;
    real lambda_beta;

    if (copy_number_type[i] == 0) {
      lambda_alpha = l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step);
      lambda_beta  = 2 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step);
    } else if (copy_number_type[i] == 1) {
      lambda_alpha = l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step);
      lambda_beta  = l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step) +
                     3 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step);
    } else {
      lambda_alpha = 2 * l_copy_number[i] * lambda(0, t_copy_number[i], t_driver, Sample, omega, lambda_pre, mu_post, k_step);
      lambda_beta  = 4 * l_copy_number[i] * lambda(t_copy_number[i], t_mrca, t_driver, Sample, omega, lambda_pre, mu_post, k_step);
    }

    alpha_rep[i] = neg_binomial_2_rng(lambda_alpha, shape);
    beta_rep[i]  = neg_binomial_2_rng(lambda_beta, shape);
  }
  
   if (exponential_growth == 1) {
    real N_rep;
    N_rep = exponential_rng(exp(-omega*(Sample - t_mrca)));
  }
  
}



