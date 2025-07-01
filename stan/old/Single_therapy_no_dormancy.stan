functions{
  real p_m2(real ti, real tf,real t_1,real t_2, real k){

    vector[2] a1;
    vector[2] a2;
    vector[2] a3;
    vector[2] a4;

    real f1;
    real f2;
    real f3;
    real f4;

    a1[1] = k*(tf-t_1);
    a1[2] = 1;
    a2[1] = k*(tf-t_2);
    a2[2] = 1;
    a3[1] = k*(ti-t_1);
    a3[2] = 1;
    a4[1] = k*(ti-t_2);
    a4[2] = 1;

    f1 = 1/k*(log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1));
    f2 = 1/k*(log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2));
    f3 = 1/k*(log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3));
    f4 = 1/k*(log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4));
    return f1 - f2 - f3 + f4;

  }
}
data{
  real <lower=0> t_1;
  real <lower=0> t_2;
  real <lower=0> t_obs;

  real <lower=0> alpha_lam;
  real <lower=0> beta_lam;
  real <lower=0> alpha_omega;
  real <lower=0> beta_omega;

  real <lower=0> lambda;

  int <lower=0> m;
  int <lower=0> m12;
  int <lower=0> ma;
  int <lower=0> mb;

  real <lower=0> N_min;
  real <lower=0> N_max;
  real <lower=0> k;
  real <lower=0> n;
}
parameters{
  real <lower=0, upper=t_1> t_0;
  real <lower=t_0, upper=t_obs> t_n;
  real <lower=t_2, upper=t_obs> t_f;

  real <lower=0> lam_new;
  real <lower=0> omega;
}
model{
  // Priors
  t_0 ~ uniform(0, t_1);
  t_f ~ uniform(t_2, t_obs);
  t_n ~ uniform(t_0, t_f);

  lam_new ~ gamma(alpha_lam, beta_lam);
  omega ~ gamma(alpha_omega, beta_omega);

  // Likelihood
  m ~ poisson(omega*lambda*((t_1-t_0) + (t_f-t_2)));
  m12 ~ poisson(omega*lam_new*(t_2-t_1));
  ma ~ poisson(omega*(lambda*((t_n-t_0) - p_m2(t_0, t_n, t_1,t_2, k)) + lam_new*p_m2(t_0,t_n,t_1,t_2, k)));
  mb ~ poisson(n*omega*(lambda*((t_f-t_n) - p_m2(t_n, t_f, t_1,t_2, k)) + lam_new*p_m2(t_n,t_f,t_1,t_2, k)));

  // Penalty
  target += -N_min*exp(-omega*(t_1-t_0)) + log(1-exp(-(N_max-N_min)*exp(-omega*(t_1-t_0))));
  target += -N_min*exp(-omega*(t_obs-t_f)) + log(1-exp(-(N_max-N_min)*exp(-omega*(t_obs-t_f))));
}
generated quantities{
  int <lower=0> m_rep = poisson_rng(omega*lambda*((t_1-t_0) + (t_f-t_2)));
  int <lower=0> m12_rep = poisson_rng(omega*lam_new*(t_2-t_1));
  int <lower=0> ma_rep = poisson_rng(omega*(lambda*((t_n-t_0) - p_m2(t_0, t_n, t_1,t_2, k)) + lam_new*p_m2(t_0,t_n,t_1,t_2, k)));
  int <lower=0> mb_rep = poisson_rng(omega*n*(lambda*((t_f-t_n) - p_m2(t_n, t_f, t_1,t_2, k)) + lam_new*p_m2(t_n,t_f,t_1,t_2, k)));

  real <lower = 0> N_rep = exp(omega*(t_obs-t_f));
}
