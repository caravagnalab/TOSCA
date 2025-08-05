// functions {
// 
// 
// 
// 	real step_integrand(real t, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
//     
//     real t1     = theta[1];
//     real t2     = theta[2];
//     real k_step = theta[3];
//     real omega  = theta[4];
//     real t_mrca = theta[5];
// 
//     real H1 = inv_logit(k_step * (t - t1));
//     real H2 = inv_logit(k_step * (t - t2));
//     real mu_t = (H1 - H2);
// 
//     return mu_t * exp(omega * (t - t_mrca));
//     
//   }
// 
// }
// 
// data{
//   
//   // Clock-like mutations
//   real <lower=0> l_diploid;
//   real<lower=0, upper=1> f_min;
//   real<lower=0, upper=1> f_max;
//   real<lower=0, upper=1> Sample_1;
// 
//   // mutations associated to driver
//   int <lower=0> cycles_driver;
//   int <lower=0,upper = 1> cumulative_driver;
//   array[cycles_driver] real<lower=0> driver_start; // if driver have effect only associated to external therapy
//   array[cycles_driver] real<lower=0> driver_end;
//   real <lower = 0> tau_driver;
//   int <lower=0> m_tail_driver;
//   real <lower = 0> N;
//    // if driver alters basal mutation rate of clock-like
// 
//   // other parameters
//   real <lower=0> mu_driver_alpha;
//   real <lower=0> mu_driver_beta;
//   real <lower=0> k_step;
//   real <lower=0> phi_driver;
//  
// }
// 
// transformed data {
//   array[0] real x_r;
//   array[0] int x_i;
//   real <lower = 0> omega = 35;
// }
// 
// parameters{
//   
//   real <lower=0> mu_driver;
//   real <lower = 0> t_mrca;
//   
// }
// 
// transformed parameters{
// 
//   
// // driver induced
// 
// real lambda_tail_driver = 0;
// 
// for (c in 1:cycles_driver) {
// 
//   if (cumulative_driver == 1) {
// 
//     real dose_time = driver_start[c];
// 
//     while (dose_time <= driver_end[c]) {
// 
//       array[5] real theta;
//       theta[1] = dose_time;
//       theta[2] = dose_time + tau_driver;
//       theta[3] = k_step;
//       theta[4] = omega;
//       theta[5] = t_mrca;
// 
//       lambda_tail_driver += integrate_1d(
//         step_integrand,
//         t_mrca + log(1 / f_max) / omega,
//         t_mrca + log(1 / f_min) / omega,
//         theta,
//         x_r, x_i
//       );
// 
//       dose_time += 1.0 / 365.0;  // advance by one day
//     }
// 
//   } else {
// 
//     array[5] real theta;
//     theta[1] = driver_start[c];
//     theta[2] = driver_end[c];
//     theta[3] = k_step;
//     theta[4] = omega;
//     theta[5] = t_mrca;
// 
//     lambda_tail_driver += integrate_1d(
//       step_integrand,
//       t_mrca + log(1 / f_max) / omega,
//       t_mrca + log(1 / f_min) / omega,
//       theta,
//       x_r, x_i
//     );
//   }
//   
// }
// 
// 
// 
// }
// 
// 
// model {
//   // Priors
//   
//   mu_driver ~ gamma(mu_driver_alpha, mu_driver_beta);
//   
//   t_mrca ~  uniform(0,Sample_1);
// 
//   // Negative Binomial shape parameters (inverse overdispersion)
//   
//  real shape_driver = 1 / phi_driver;
// 
//  m_tail_driver ~ neg_binomial_2(
//       2 * l_diploid * omega * mu_driver * lambda_tail_driver,
//       shape_driver
//     );
//     
//     N ~ exponential(exp(-omega*(Sample_1 - t_mrca)));
//   
// }
// 
// generated quantities {
//   
//   int m_tail_driver_rep;
//   real N_rep;
// 
//    m_tail_driver_rep = neg_binomial_2_rng(
//       2 * l_diploid * omega * mu_driver * lambda_tail_driver,
//       phi_driver
//     );
//     
//     N_rep = exponential_rng(exp(-omega*(Sample_1 - t_mrca)));
// 
// }


functions {
  // Logistic bump for one driver interval
  real bump(real t, real t1, real t2, real k_step) {
    real H1 = inv_logit(k_step * (t - t1));
    real H2 = inv_logit(k_step * (t - t2));
    return H1 - H2;
  }

  // Combined integrand summing all driver bumps at time t
  real combined_integrand(real t, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    int cycles_driver = x_i[1];
    int cumulative_driver = x_i[2];
    real k_step = theta[1];
    real omega = theta[2];
    real t_mrca = theta[3];
    real tau_driver = theta[4];

    real sum_bumps = 0;

    if (cumulative_driver == 1) {
      for (c in 1:cycles_driver) {
        real start_c = x_r[2 * c - 1];
        real end_c = x_r[2 * c];
        real dose_time = start_c;

        while (dose_time <= end_c) {
          sum_bumps += bump(t, dose_time, dose_time + tau_driver, k_step);
          dose_time += 1.0 / 365.0;
        }
      }
    } else {
      for (c in 1:cycles_driver) {
        real start_c = x_r[2 * c - 1];
        real end_c = x_r[2 * c];
        sum_bumps += bump(t, start_c, end_c, k_step);
      }
    }

    return sum_bumps * exp(omega * (t - t_mrca));
  }
}

data {
  real<lower=0> l_diploid;
  real<lower=0, upper=1> f_min;
  real<lower=0, upper=1> f_max;
  real<lower=0> Sample_1;

  int<lower=0> cycles_driver;
  int<lower=0, upper=1> cumulative_driver;
  array[cycles_driver] real<lower=0> driver_start;
  array[cycles_driver] real<lower=0> driver_end;
  real<lower=0> tau_driver;
  int<lower=0> m_tail_driver;
  real<lower=0> N;

  real<lower=0> mu_driver_alpha;
  real<lower=0> mu_driver_beta;
  real<lower=0> k_step;
  real<lower=0> phi_driver;
}

transformed data {
  array[2 * cycles_driver] real x_r;
  array[2] int x_i;
  real omega = 35;

  for (i in 1:cycles_driver) {
    x_r[2 * i - 1] = driver_start[i];
    x_r[2 * i] = driver_end[i];
  }

  x_i[1] = cycles_driver;
  x_i[2] = cumulative_driver;
}

parameters {
  real<lower=0> mu_driver;
  real<lower=0,upper = Sample_1> t_mrca;
}

transformed parameters {
  real lambda_tail_driver;
  array[4] real theta;

  theta[1] = k_step;
  theta[2] = omega;
  theta[3] = t_mrca;
  theta[4] = tau_driver;

  real t_lower = t_mrca + log(1 / f_max) / omega;
  real t_upper = t_mrca + log(1 / f_min) / omega;

  lambda_tail_driver = integrate_1d(
    combined_integrand,
    t_lower,
    t_upper,
    theta,
    x_r,
    x_i
  );
}

model {
  
  mu_driver ~ gamma(mu_driver_alpha, mu_driver_beta);
  t_mrca ~ uniform(0, Sample_1);

  real shape_driver = 1 / phi_driver;

  m_tail_driver ~ neg_binomial_2(
    2 * l_diploid * omega * mu_driver * lambda_tail_driver,
    shape_driver
  );

  N ~ exponential(exp(-omega * (Sample_1 - t_mrca)));
  
}

generated quantities {
  int m_tail_driver_rep;
  real N_rep;
  real shape_driver = 1 / phi_driver;

  m_tail_driver_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_driver * lambda_tail_driver,
    shape_driver
  );

  N_rep = exponential_rng(exp(-omega * (Sample_1 - t_mrca)));
}




