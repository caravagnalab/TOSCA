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
	
real cauchy_integrand(
  real t,
  array[] real locations,
  int n_locations,
  real scale,
  real omega,
  real t_mrca
) {
  real sum_peaks = 0;
  for (i in 1:n_locations) {
    sum_peaks += (1 / pi()) * (scale / (square(t - locations[i]) + square(scale)));
  }
  return sum_peaks * exp(omega * (t - t_mrca)) * pi() * scale;
}

real trapezoidal_rule_cauchy(
  real t0,
  real t1,
  int n_points,
  array[] real locations,
  int n_locations,
  real scale,
  real omega,
  real t_mrca
) {
  real h = (t1 - t0) / n_points;
  real sum = 0;

  for (i in 1:(n_points - 1)) {
    real t = t0 + i * h;
    sum += cauchy_integrand(t, locations, n_locations, scale, omega, t_mrca);
  }

  real f0 = cauchy_integrand(t0, locations, n_locations, scale, omega, t_mrca);
  real f1 = cauchy_integrand(t1, locations, n_locations, scale, omega, t_mrca);

  return h * (0.5 * f0 + sum + 0.5 * f1);
}

	

real step_integrand(real t, array[] real t1s, array[] real t2s, int n_bumps, real k_step, real omega, real t_mrca){
  real sum_bumps = 0;
  for (i in 1:n_bumps) {
    sum_bumps += inv_logit(k_step * (t - t1s[i])) - inv_logit(k_step * (t - t2s[i]));
  }
  return sum_bumps * exp(omega * (t - t_mrca));
}


real trapezoidal_rule_step(
  real t0,
  real t1,
  int n_points,
  array[] real t1s,
  array[] real t2s,
  int n_bumps,
  real k_step,
  real omega,
  real t_mrca
) {
  real h = (t1 - t0) / n_points;
  real sum = 0;
  
  for (i in 1:(n_points - 1)) {
    real t = t0 + i * h;
    sum += step_integrand(t, t1s, t2s, n_bumps, k_step, omega, t_mrca);
  }
  
  real f0 = step_integrand(t0, t1s, t2s, n_bumps, k_step, omega, t_mrca);
  real f1 = step_integrand(t1, t1s, t2s, n_bumps, k_step, omega, t_mrca);
  
  return h * (0.5 * f0 + sum + 0.5 * f1);
  
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
  
  int n_points = 1000;   

  // Pre-allocate arrays with maximum size, since Stan requires sizes fixed or known
  int n_bumps_driver = 0;
  for (c in 1:cycles_driver) {
    real dose_time = driver_start[c];
    while (dose_time <= driver_end[c]) {
      n_bumps_driver  += 1;
      dose_time += 1.0 / 365.0;
    }
  }

  array[n_bumps_driver] real driver_t1s;
  array[n_bumps_driver] real driver_t2s;
 
 int k = 0;
 
  for (c in 1:cycles_driver) {
    real dose_time = driver_start[c];
    while (dose_time <= driver_end[c]) {
       k += 1;
      driver_t1s[k] = dose_time;
      driver_t2s[k] = dose_time + tau_driver;
      dose_time += 1.0 / 365.0;
      
    }
  }
  
  array[n_th_cauchy_type] int n_peaks_cauchy;

  for (th_cauchy in 1:n_th_cauchy_type) {
    n_peaks_cauchy[th_cauchy] = 0;
    for (cycle in 1:n_th_cauchy) {
      if (type_th_cauchy[cycle] == th_cauchy) {
        n_peaks_cauchy[th_cauchy] += 1;
      }
    }
  }
 
 array[n_th_step_type] int n_bumps_th;
 
if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    // Temporary arrays to hold the t1s and t2s for this therapy type
    
    n_bumps_th[th_type] = 0;
    
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {
        real dose_time = start_th_step[cycle];

        while (dose_time <= end_th_step[cycle]) {
          n_bumps_th[th_type] += 1;
          dose_time += 1.0 / 365.0;
        }
      }
    }  
  
  
}

}

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


// Step therapy clonal mutations
if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {

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
  
     real dose_time = driver_start[c];
     
    while(dose_time <= driver_end[c]){
      
      lambda_driver += lambda_therapy_single(t_driver, t_mrca, dose_time, dose_time + tau_driver, k_step);
      dose_time = dose_time + 1.0/365.0; 
    }
 
  
}


// subclonal 

// therapy step 

array[n_th_step_type] real lambda_tail_th_step = rep_array(0.0, n_th_step_type);

if (n_th_step_type > 0) {
  for (th_type in 1:n_th_step_type) {
    // Temporary arrays to hold the t1s and t2s for this therapy type
    array[n_bumps_th[th_type]] real t1s_tmp;
    array[n_bumps_th[th_type]] real t2s_tmp;
    int n_bumps = 0;
    
    for (cycle in 1:n_th_step) {
      if (type_th_step[cycle] == th_type) {
        real dose_time = start_th_step[cycle];

        while (dose_time <= end_th_step[cycle]) {
          n_bumps += 1;
          t1s_tmp[n_bumps] = dose_time;
          t2s_tmp[n_bumps] = dose_time + tau_step[th_type];
          dose_time += 1.0 / 365.0;
        }
      }
    }

    // Only integrate if at least one bump
   
      real t0 = t_mrca + log(1 / f_max) / omega;
      real t1 = t_mrca + log(1 / f_min) / omega;

      lambda_tail_th_step[th_type] = trapezoidal_rule_step(
        t0,
        t1,
        n_points,
        t1s_tmp,
        t2s_tmp,
        n_bumps_th[th_type],
        k_step,
        omega,
        t_mrca
      );
    }
}


// cauchy

array[n_th_cauchy_type] real lambda_tail_th_cauchy = rep_array(0.0, n_th_cauchy_type);

if (n_th_cauchy_type > 0) {
  for (th_cauchy in 1:n_th_cauchy_type) {
    
    array[n_peaks_cauchy[th_cauchy]] real locs_tmp;
    int n_peaks = 0;
    
    for (cycle in 1:n_th_cauchy) {
      if (type_th_cauchy[cycle] == th_cauchy) {
        n_peaks += 1;
        locs_tmp[n_peaks] = location_th_cauchy[cycle];
      }
    }

      real t0 = t_mrca + log(1 / f_max) / omega;
      real t1 = t_mrca + log(1 / f_min) / omega;

      lambda_tail_th_cauchy[th_cauchy] = trapezoidal_rule_cauchy(
        t0,
        t1,
        n_points,
        locs_tmp,
        n_peaks_cauchy[th_cauchy],
        scales_th_cauchy[th_cauchy],
        omega,
        t_mrca
      );
    }
  }



// driver induced


// Perform a single integration using trapezoidal rule

real lambda_tail_driver = 0;

real t0 = t_mrca + log(1 / f_max) / omega;
real t1 = t_mrca + log(1 / f_min) / omega;

lambda_tail_driver = trapezoidal_rule_step(
  t0,
  t1,
  n_points,
  driver_t1s,
  driver_t2s,
  n_bumps_driver,
  k_step,
  omega,
  t_mrca
);



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
    
     m_tail_cauchy[th_cauchy] ~ neg_binomial_2(
      2 * l_diploid * omega * mu_clock * lambda_tail_th_cauchy[th_cauchy] + 1,
      shape_th_cauchy[th_cauchy]
    );
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
    target += -N_min[1] * exp(-omega * (Sample_1 - t_mrca_primary)) +
      log(1 - exp(-(N_max[1] - N_min[1]) * exp(-omega * (Sample_1 - t_mrca_primary))));
    target += -N_min[2] * exp(-omega * (Sample_2 - t_mrca)) +
      log(1 - exp(-(N_max[2] - N_min[2]) * exp(-omega * (Sample_2 - t_mrca))));
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
  array[n_th_cauchy_type] int m_tail_cauchy_rep;
  
  
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

    m_tail_cauchy_rep[th_cauchy] = neg_binomial_2_rng(
      2 * l_diploid * omega * mu_clock * lambda_tail_th_cauchy[th_cauchy],
      shape_th_cauchy[th_cauchy]
    );
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
     
     N_relapse_rep = exponential_rng(exp(-omega*(Sample_2 - t_mrca)));
     N_primary_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca_primary)));
     
  }

}




