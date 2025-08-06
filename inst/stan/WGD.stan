data{

  // Clock-like mutations
  int <lower=0> m_clock_primary;
  real <lower=0> l_diploid;
  real <lower=0> mu_clock;

  int <lower=0> m_alpha;
  int <lower=0> m_beta;
  int <lower=0> l_wgd;

	// other parameters
  real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;

  // real <lower=0> chemo_start;
  // real <lower=0> chemo_end;

  int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;
}

parameters{

  real <lower=0, upper= Sample_1> t_eca;
  real <lower=t_eca, upper= Sample_1> t_mrca_primary;
  real <lower= t_eca, upper= Sample_2> t_mrca;
  real <lower=t_eca, upper=t_mrca> t_wgd;

  real <lower= 0> omega;

}

model{

  t_eca ~ uniform(0, Sample_1);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  t_mrca ~ uniform(t_eca, Sample_2);
  t_wgd ~ uniform(t_eca, t_mrca);

  omega ~ gamma(omega_alpha,omega_beta);

  m_clock_primary ~ poisson(2*mu_clock*omega*l_diploid*(t_mrca_primary-t_eca));

  m_alpha ~  poisson(2*omega*l_wgd*mu_clock*(t_wgd-t_eca));
  m_beta ~  poisson(4*omega*l_wgd*mu_clock*( t_mrca - t_wgd));

  if (exponential_growth==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }
}

generated quantities{

  int m_clock_primary_rep = poisson_rng(2*mu_clock*omega*l_diploid*(t_mrca_primary-t_eca));
  int m_alpha_rep =  poisson_rng(2*omega*l_wgd*mu_clock*(t_wgd-t_eca));
  int m_beta_rep =  poisson_rng(4*omega*l_wgd*mu_clock*(t_mrca - t_wgd));

  if (exponential_growth==1){
    real N_sample_2 = exp(omega*(Sample_2-t_mrca));
    real N_sample_1 = exp(omega*(Sample_1-t_eca));
  }

}
