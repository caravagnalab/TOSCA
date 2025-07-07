functions {

  real lambda_therapy_single(real ti, real tf,real t_therapy_i,real t_therapy_f, real k){

    vector[2] a1;
    vector[2] a2;
    vector[2] a3;
    vector[2] a4;

    real f1;
    real f2;
    real f3;
    real f4;

    a1[1] = k*(tf-t_therapy_i);
    a1[2] = 1;
    a2[1] = k*(tf-t_therapy_f);
    a2[2] = 1;
    a3[1] = k*(ti-t_therapy_i);
    a3[2] = 1;
    a4[1] = k*(ti-t_therapy_f);
    a4[2] = 1;

    f1 = 1/k*(log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1));
    f2 = 1/k*(log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2));
    f3 = 1/k*(log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3));
    f4 = 1/k*(log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4));
    return f1 - f2 - f3 + f4;

  }


  // real lambda_therapy(int n_cycles, real ti, real tf, array[n_cycles] t_therapy_i, array[n_cycles] t_therapy_f, real k){
  //
  //   real f=0;
  //
  //   for (c in 1:n_cycles){
  //     f += lambda_therapy_single(ti, tf, t_therapy_i[c], t_therapy_f[c], k);
  //   }
  //   return(f);
  //
  // }

	real couchy_cdf_single(real location, real scale, real a,real b){
    real d= ((1/pi()) * atan((b-location)/scale) + .5) - ((1/pi()) * atan((a-location)/scale) + .5);
  return(d);
	}

// 	real couchy_cdf(int n_cycles, array[n_cycles] location, array[n_cycles] scale, real a,real b){
// 	  real c=0;
// 	  for (i in 1:n_cycles){
//       c += couchy_cdf_single(location[i], scale[i], a, b);
// 	  }
//     return(c);
// 	}

}

data{
  // Clock-like mutations
  int <lower=0> m_clock;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;


  // mutations associated to driver
  int <lower=0, upper=1> driver_type; // 0 = endogeno, 1 = dipendente da esogeno
  int <lower=0> cycles_driver;
  array[cycles_driver] real<lower=0> driver_start; // if driver have effect only associated to external therapy
  array[cycles_driver] real<lower=0> driver_end;
  int <lower=0> m_driver;
  // real <lower=0> mu_driver_alpha;
  // real <lower=0> mu_driver_beta;
  real mu_driver;
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

  // mutations associated to cauchy
  int <lower=0> n_th_cauchy;
  int <lower=0> n_th_cauchy_type;
  //vector<lower=0>[n_th_cauchy] cycles_th_cauchy;
  array[n_th_cauchy] real<lower=0> location_th_cauchy;
  array[n_th_cauchy] int<lower=0> type_th_cauchy; // vector with numbers identifying the therapy (1:n_th_cauchy)
  array[n_th_cauchy_type] real<lower=0> alpha_th_cauchy;
  array[n_th_cauchy_type] real<lower=0> beta_th_cauchy;
  array[n_th_cauchy_type] int<lower=0> m_th_cauchy;

  // other parameters
  //real <lower=0> omega_alpha;
  //real <lower=0> omega_beta;
  real omega;
  real <lower=0> k_step;
  // real <lower=0> k_softmax;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> max_therapy;
  // int <lower=0, upper=1> exponential_growth;
  real N_min;
  real N_max;

  real alpha_mrca;
  real beta_mrca;

}

parameters{
  real <lower=0, upper=Sample_1> t_eca;
  // real <lower=max_therapy, upper=Sample_2> t_mrca;
  real rho_mrca;
  real <lower=t_eca, upper=t_mrca> t_driver;
  array[n_th_step_type] real<lower=0> mu_th_step;
  array[n_th_cauchy_type] real<lower=0> scales_th_cauchy;
  // real <lower=0> omega;
  //real <lower=0> mu_driver;
}

transformed parameters{
  real t_mrca = max_therapy + rho_mrca*(Sample_2-max_therapy);
  array[n_th_step_type] real lambda_th_step;
  array[n_th_cauchy_type] real lambda_th_cauchy;

  for (i in 1:n_th_step_type) lambda_th_step[i] = 0;
  for (i in 1:n_th_cauchy_type) lambda_th_cauchy[i] = 0;

  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    for (cycle in 1:n_th_step){
      if (type_th_step[cycle] == th_type){
        lambda_th_step[th_type] += lambda_therapy_single(t_eca, t_mrca, start_th_step[cycle], end_th_step[cycle], k_step);
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
      }
    }
  }
  }

}

model{

  // Priors
  t_eca ~ uniform(0, Sample_1);
  //t_mrca ~ uniform(max_therapy, Sample_2);
  rho_mrca ~ beta(alpha_mrca, beta_mrca);
  t_driver ~ uniform(t_eca, t_mrca);

  for (m in 1:n_th_step_type){
      mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
    }

  for (ch in 1:n_th_cauchy_type){
      scales_th_cauchy[ch] ~ gamma(alpha_th_cauchy[ch], beta_th_cauchy[ch]);
    }


  // omega ~ gamma(omega_alpha, omega_beta);
  //mu_driver ~ gamma(mu_driver_alpha,mu_driver_beta);

  // Likelihood
  m_clock ~ poisson(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  // Step therapy mutations
  if (n_th_step_type > 0){
  for (th_type in 1:n_th_step_type){
    m_th_step[th_type] ~ poisson(2 * l_diploid * omega * mu_th_step[th_type] * lambda_th_step[th_type]);
  }
  }

  // Cauchy therapy mutations
  if (n_th_cauchy_type > 0){
  for (th_cauchy in 1:n_th_cauchy_type){
    m_th_cauchy[th_cauchy] ~ poisson(2 * l_diploid * omega * mu_clock * lambda_th_cauchy[th_cauchy]);
  }
  }

  if (driver_type==0){
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    real lambda_driver=0;
    for (c in 1:cycles_driver){
    lambda_driver += lambda_therapy_single(t_driver,t_mrca, driver_start[c],driver_end[c],k_step);
    }
    m_driver ~ poisson(2*l_diploid*omega*mu_driver*lambda_driver);
  }

  // if (exponential_growth==1){
  //   target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
  // }

}

generated quantities{

  int<lower =0> m_clock_rep = poisson_rng(2*l_diploid*omega*(mu_clock*(t_driver-t_eca) + mu_driver_clock*(t_mrca-t_driver)));

  int <lower =0> m_driver_rep;

  if (driver_type==0){
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*(t_mrca-t_driver));
  }else{
    real lambda_driver_rep=0;
    for (c in 1:cycles_driver){
      lambda_driver_rep += lambda_therapy_single(t_driver,t_mrca, driver_start[c],driver_end[c],k_step);
    }
    m_driver_rep = poisson_rng(2*l_diploid*omega*mu_driver*lambda_driver_rep);
  }


  // if (exponential_growth==1){
  //   real N = exp(omega*(Sample_2-t_mrca));
  // }

}
