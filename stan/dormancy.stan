functions {

  real lambda_therapy(real ti, real tf,real t_therapy_i,real t_therapy_f, real k){

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

  real delta(real t,real dorm_s,real dorm_e,real k){
    return  dorm_e + (t-dorm_e)/(1+exp(k*(t-dorm_e))) - dorm_s;
    }

  real traslation(real t, real dorm_start, real dorm_end ) {
    real t_cn =t;
    if (t> dorm_start){
      t_cn= t_cn + dorm_end - dorm_start;
      }
      return t_cn;
	}
}

data{

  int <lower = 0> m_clock;
	int <lower = 0> m_chemo;
	int <lower = 0> m_th_1;
	int <lower = 0> m_alpha;
	int <lower = 0> m_beta;
	int <lower = 0> Major;
	int <lower = 0> Minor;
	real <lower = 0> N_min;
	real <lower = 0> N_max;

	real <lower=0> mu;
	real <lower = 0> omega_alpha;
	real <lower = 0> omega_beta;
	real <lower = 0> alpha_mu_th_1;
	real <lower = 0> beta_mu_th_1;

	real <lower = 0> k;
	real <lower = 0> k_sm;
	real <lower = 0> CNA_length;
	real <lower = 0> diploid_length;
	real <lower = 0> Sample_2;
	//int <lower= 0> cycles;
	real <lower = 0> Chemo_start;
	real <lower = 0> Therapy_1_start;
	real <lower = 0> Chemo_end;
	real <lower = 0> Therapy_1_end;
	//real <lower =0> transplant;
}

parameters{

  real <lower=Chemo_start - 1, upper= Chemo_start> t_eca;
  real <lower= Chemo_start, upper= Chemo_end> t_dorm_start;
  real <lower= Chemo_end, upper= Therapy_1_end> t_dorm_end;
  real <lower= t_dorm_end, upper= Sample_2> t_mrca;
  real <lower= t_eca, upper= t_mrca - t_dorm_end+ t_dorm_start> t_cn_tr;

  real <lower= 0> omega;
  real <lower=0> mu_th_1;
}

transformed parameters{
  int coeff;
  if (Minor==0){
    coeff = 2;
  }
  if (Minor==1){
    coeff = 3;
  }
  if (Minor==1){
    coeff = 4;
  }
}

model{
  t_eca ~ uniform(Chemo_start - 1, Chemo_start);
  t_dorm_start~ uniform(Chemo_start, Chemo_end);
  t_dorm_end~ uniform(Chemo_end, Therapy_1_end);
  t_mrca ~ uniform(t_dorm_end,Sample_2);
  t_cn_tr ~ uniform(t_eca, t_mrca-t_dorm_end+t_dorm_start);

  omega ~ gamma(omega_alpha,omega_beta);
  mu_th_1 ~ gamma(alpha_mu_th_1,beta_mu_th_1);

  m_clock ~  poisson(2*mu*omega*diploid_length*((Chemo_start - t_eca) + (t_mrca-t_dorm_end)
  - lambda_therapy(t_dorm_end, t_mrca, Therapy_1_start, Therapy_1_end,k)));
  m_chemo ~  poisson(2*mu*omega*diploid_length*(t_dorm_start - Chemo_start));
  m_th_1  ~  poisson(2*mu_th_1*omega*diploid_length*lambda_therapy(t_dorm_end,t_mrca, Therapy_1_start, Therapy_1_end,k));
  m_alpha ~  poisson(omega*CNA_length*(mu*((t_cn_tr - t_eca) - lambda_therapy(t_eca,t_cn_tr,Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm),Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)) + mu_th_1*lambda_therapy(t_eca,t_cn_tr,Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm), k)));
  m_beta ~ poisson(coeff*omega*CNA_length*(mu*( t_mrca - (t_dorm_end - t_dorm_start) - t_cn_tr - lambda_therapy(t_cn_tr,t_mrca-(t_dorm_end - t_dorm_start), Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)) + mu_th_1*lambda_therapy(t_cn_tr,t_mrca-(t_dorm_end - t_dorm_start), Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)));

  target += -N_min*exp(-omega*(Sample_2 - (t_dorm_end- t_dorm_start) - t_cn_tr)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - (t_dorm_end- t_dorm_start) - t_cn_tr))));
  //target += -N_min*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max-N_min)*exp(-omega*(Sample_2 - t_mrca))));
  target += log(1-exp(-N_max*exp(-omega*(Chemo_start - t_eca ))));
}

generated quantities{

    real t_copynumber;
    real N;

    int m_alpha_rep =  poisson_rng(omega*CNA_length*(mu*((t_cn_tr - t_eca) - lambda_therapy(t_eca,t_cn_tr,Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm),Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)) + mu_th_1*lambda_therapy(t_eca,t_cn_tr,Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm), k)));
    int m_beta_rep = poisson_rng(2*omega*CNA_length*(mu*( t_mrca - (t_dorm_end - t_dorm_start) - t_cn_tr - lambda_therapy(t_cn_tr,t_mrca-(t_dorm_end - t_dorm_start), Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)) + mu_th_1*lambda_therapy(t_cn_tr,t_mrca-(t_dorm_end - t_dorm_start), Therapy_1_start - delta(Therapy_1_start,t_dorm_start,t_dorm_end,k_sm), Therapy_1_end - delta(Therapy_1_end,t_dorm_start,t_dorm_end,k_sm),k)));

    int m_clock_rep =  poisson_rng(2*mu*omega*diploid_length*((Chemo_start - t_eca) + (t_mrca-t_dorm_end) - lambda_therapy(t_dorm_end, t_mrca, Therapy_1_start, Therapy_1_end,k)));
    int m_th_1_rep  =  poisson_rng(2*mu_th_1*omega*diploid_length*lambda_therapy(t_dorm_end,t_mrca, Therapy_1_start, Therapy_1_end,k));

    t_cn = traslation(t_cn_tr,  t_dorm_start,  t_dorm_end );
    N = exp(omega*(Sample_2  - (t_dorm_end- t_dorm_start) - t_cn_tr));
}
