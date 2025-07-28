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
	
	real cauchy_exp_integrand(real t, real xc, real[] theta, real[] x_r, int[] x_i) {
    real location = theta[1];
    real scale    = theta[2];
    real mu_eff   = theta[3];
    real omega    = theta[4];
    real t_mrca   = theta[5];
    
    real pdf = (1 / pi()) * (scale / (square(t - location) + square(scale)));
    
    return mu_eff * pdf * exp(omega * (t - t_mrca));
  }
	

  real step_integrand(real t, real xc, real[] theta, real[] x_r, int[] x_i) {
    
    real mu_eff = theta[1];
    real t1     = theta[2];
    real t2     = theta[3];
    real k_step = theta[4];
    real omega  = theta[5];
    real t_mrca = theta[6];

    real H1 = inv_logit(k_step * (t - t1));
    real H2 = inv_logit(k_step * (t - t2));
    real mu_t = mu_eff * (H1 - H2);

    return mu_t * exp(omega * (t - t_mrca));
    
  }

}

data{
  // Clock-like mutations
  int <lower=0> m_clock_primary;
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;
  real <lower=0> mu_driver_clock;


  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  int <lower=0,upper = 1> cumulative_driver;
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
  array[n_th_step_type] int<lower=0,upper = 1> cumulative_step;
  array[n_th_step_type] real<lower=0> tau_step;
  array[n_th_step_type] real<lower=0> m_tail_step;

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
  real <lower=0> k_step;

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

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters{
  
  real <lower=0, upper=Sample_1> t_eca;
  real <lower=t_eca, upper=Sample_1> t_mrca_primary;
  real <lower=0, upper=1> rho_mrca;
  real <lower=t_eca, upper=driver_end[cycles_driver]> t_driver;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  real <lower=0> omega;
  real <lower=0> mu_driver;
  
}

transformed parameters{

  real <lower=max_therapy> t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);

// clonal muts

  array[n_th_step_type] real lambda_th_step;
  array[n_th_cauchy_type] real lambda_th_cauchy;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy[i] = 0;

  array[n_th_step_type] real lambda_th_step;

// Initialize to zero
for (i in 1:n_th_step_type) {
  lambda_th_step[i] = 0;
}

// Step therapy clonal mutations
if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {

        if (cumulative_step[th_type] == 1) {
          // Cumulative dosing: apply effect for each daily dose
          real dose_time = start_th_step[cycle];
          while (dose_time <= end_th_step[cycle]) {
            lambda_th_step[th_type] += lambda_therapy_single(
              t_eca, t_mrca,
              dose_time,
              dose_time + tau_step[th_type],
              k_step
            );
            dose_time += 1.0 / 365.0;
          }

        } else {
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


// driver induced

 real lambda_driver = 0;

for (c in 1:cycles_driver) {
  
  if (cumulative_driver == 1) {
    real dose_time = driver_start[c];
    while(dose_time <= driver_end[c]){
      
      lambda_driver += lambda_therapy_single(t_driver, t_mrca, dose_time, dose_time + tau_driver, k_step);
      dose_time = dose_time + 1.0/365.0; 
    }
    
  } else {
    
    lambda_driver += lambda_therapy_single(
      t_driver, t_mrca,
      driver_start[c],
      driver_end[c],
      k_step
      );
      
  }
  
}


// subclonal 

// therapy step 

array[n_th_step_type] real lambda_tail_th_step = rep_array(0.0, n_th_step_type);
array[n_th_cauchy_type] real lambda_tail_th_cauchy = rep_array(0.0, n_th_cauchy_type);

if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {

        if (cumulative_step[th_type] == 1) {
          real dose_time = start_th_step[cycle];
          while (dose_time <= end_th_step[cycle]) {
            real theta[6];
            theta[1] = mu_th_step[th_type];              // <-- your per-type mu
            theta[2] = dose_time;
            theta[3] = dose_time + tau_step[th_type];
            theta[4] = k_step;
            theta[5] = omega;
            theta[6] = t_mrca;

            lambda_tail_th_step[th_type] += integrate_1d(
              step_integrand,
              t_mrca + log(1 / f_max) / omega,
              t_mrca + log(1 / f_min) / omega,
              0.0,
              theta,
              x_r, x_i,
              1e-8
            );

            dose_time += 1.0 / 365.0;
          }

        } else {
          real theta[6];
          theta[1] = mu_th_step[th_type];
          theta[2] = start_th_step[cycle];
          theta[3] = end_th_step[cycle];
          theta[4] = k_step;
          theta[5] = omega;
          theta[6] = t_mrca;

          lambda_tail_th_step[th_type] += integrate_1d(
            step_integrand,
            t_mrca + log(1 / f_max) / omega,
            t_mrca + log(1 / f_min) / omega,
            0.0,
            theta,
            x_r, x_i,
            1e-8
          );
        }
      }
    }
  }
}


array[n_th_cauchy_type] real lambda_tail_th_cauchy = rep_array(0.0, n_th_cauchy_type);

if (n_th_cauchy_type > 0) {
  for (th_cauchy in 1:n_th_cauchy_type) {
    for (cycle in 1:n_th_cauchy) {
      if (type_th_cauchy[cycle] == th_cauchy) {
        real theta[5];
        theta[1] = location_th_cauchy[cycle];
        theta[2] = scales_th_cauchy[cycle];
        theta[3] = mu_th_cauchy[th_cauchy];  // scaling factor per type
        theta[4] = omega;
        theta[5] = t_mrca;
        
        lambda_tail_th_cauchy[th_cauchy] += integrate_1d(
          cauchy_exp_integrand,
          t_mrca + log(1 / f_max) / omega,
          t_mrca + log(1 / f_min) / omega,
          0.0,
          theta,
          x_r, x_i,
          1e-8
        );
      }
    }
  }
}





// driver induced

real lambda_tail_driver = 0;

for (c in 1:cycles_driver) {

  if (cumulative_driver == 1) {

    real dose_time = driver_start[c];

    while (dose_time <= driver_end[c]) {

      real theta[6];
      theta[1] = mu_driver;
      theta[2] = dose_time;
      theta[3] = dose_time + tau_driver;
      theta[4] = k_step;
      theta[5] = omega;
      theta[6] = t_mrca;

      lambda_tail_driver += integrate_1d(
        step_integrand,
        t_mrca + log(1 / f_max) / omega,
        t_mrca + log(1 / f_min) / omega,
        0.0,
        theta,
        x_r, x_i,
        1e-8
      );

      dose_time += 1.0 / 365.0;  // advance by one day
    }

  } else {

    real theta[6];
    theta[1] = mu_driver;
    theta[2] = driver_start[c];
    theta[3] = driver_end[c];
    theta[4] = k_step;
    theta[5] = omega;
    theta[6] = t_mrca;

    lambda_tail_driver += integrate_1d(
      step_integrand,
      t_mrca + log(1 / f_max) / omega,
      t_mrca + log(1 / f_min) / omega,
      0.0,
      theta,
      x_r, x_i,
      1e-8
    );
  }
}



}


model {
  // Priors
  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
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
  }

  for (th_cauchy in 1:n_th_cauchy_type) {
    m_th_cauchy[th_cauchy] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy],
      shape_th_cauchy[th_cauchy]
    );
  }

  if (driver_type == 0) {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * (t_mrca - t_driver),
      shape_driver
    );
  } else {
    m_driver ~ neg_binomial_2(
      2 * l_diploid * omega * mu_driver * lambda_driver,
      shape_driver
    );
  
}

  if (exponential_growth == 1) {
    target += -N_min[1] * exp(-omega * (Sample_1 - t_mrca_primary)) +
      log(1 - exp(-(N_max[1] - N_min[1]) * exp(-omega * (Sample_1 - t_mrca_primary))));
    target += -N_min[2] * exp(-omega * (Sample_2 - t_mrca)) +
      log(1 - exp(-(N_max[2] - N_min[2]) * exp(-omega * (Sample_2 - t_mrca))));
  }
  
}

generated quantities {
  
  int<lower=0> m_clock_primary_rep;
  int<lower=0> m_clock_rep;
  int<lower=0> m_driver_rep;
  
  array[n_th_step_type] int<lower=0> m_th_step_rep;
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy_rep;

  // Calculate shape parameters as in model block
  real shape_clock = 1 / phi_clock;
  real shape_driver = 1 / phi_driver;
  array[n_th_step_type] real shape_th_step;
  array[n_th_cauchy_type] real shape_th_cauchy;

  for (m in 1:n_th_step_type) {
    shape_th_step[m] = 1 / phi_th_step[m];
  }
  for (ch in 1:n_th_cauchy_type) {
    shape_th_cauchy[ch] = 1 / phi_th_cauchy[ch];
  }

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
      2 * l_diploid * omega * mu_driver * lambda_driver,
      shape_driver
    );
  }

  for (m in 1:n_th_step_type) {
    m_th_step_rep[m] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_th_step[m] * lambda_th_step[m],
      shape_th_step[m]
    );
  }

  for (ch in 1:n_th_cauchy_type) {
    m_th_cauchy_rep[ch] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_clock * lambda_th_cauchy[ch],
      shape_th_cauchy[ch]
    );
  }
  
  
   if (exponential_growth == 1) {
     real N_primary_rep;
     real N_relapse_rep;
     
     N_primary_rep = exponential_rng(exp(-omega*(Sample_2 - t_mrca)));
     N_relapse_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
     
  }
  
}


