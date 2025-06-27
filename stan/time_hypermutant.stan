functions{
  real lambda_mutagenesis(real ti, real tf,real t_therapy_i,real t_therapy_f, real k){

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
}

data{
  int <lower=0> m_clock;
  int <lower=0> m_mutag;

  real <lower=0> mu;
  real <lower=0> mu_mutag;
  real <lower=0> k;
  real <lower=0> omega;
  real <lower=0> N_min;
  real <lower=0> N_max;
  real <lower=0> diploid_length;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> Therapy_1_start;
  real <lower=0> Therapy_1_end;
}

parameters{
  real <lower= Sample_1 - 7, upper= Sample_1> t_eca;
  real <lower= t_eca, upper= Therapy_1_end> t_mutag;
  real <lower= Therapy_1_end, upper= Sample_2> t_mrca;
}

model{
  t_eca ~ uniform(Sample_1 - 7, Sample_1);
  t_mutag ~ uniform(t_eca, Therapy_1_end);
  t_mrca ~ uniform(Therapy_1_end, Sample_2);

  m_clock ~ poisson(2*omega*diploid_length*(mu*(t_mutag-t_eca) + mu_mutag*(t_mrca-t_mutag)));
  m_mutag ~ poisson(2*omega*diploid_length*mu_mutag*lambda_mutagenesis(t_mutag, t_mrca, Therapy_1_start, Therapy_1_end, k));

  target += -N_min*exp(-omega*(Sample_1 - t_eca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_1 - t_eca))));
  target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
}

generated quantities{
  int m_clock_rep = poisson_rng(2*omega*diploid_length*(mu*(t_mutag-t_eca) + mu_mutag*(t_mrca-t_mutag)));
  int m_mutag_rep = poisson_rng(2*omega*diploid_length*mu_mutag*lambda_mutagenesis(t_mutag, t_mrca, Therapy_1_start, Therapy_1_end, k));

  real N1 = -N_min*exp(-omega*(Sample_1 - t_eca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_1 - t_eca))));
  real N2 = -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
}
