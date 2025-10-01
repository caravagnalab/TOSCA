functions {
  real lambda_therapy(real ti, real tf, real t_therapy_i, real t_therapy_f, real k) {
    vector[2] a1; vector[2] a2; vector[2] a3; vector[2] a4;
    real f1; real f2; real f3; real f4;
    a1[1] = k*(tf - t_therapy_i); a1[2] = 1;
    a2[1] = k*(tf - t_therapy_f); a2[2] = 1;
    a3[1] = k*(ti - t_therapy_i); a3[2] = 1;
    a4[1] = k*(ti - t_therapy_f); a4[2] = 1;
    f1 = (log(exp(a1[1]-max(a1)) + exp(-max(a1))) + max(a1)) / k;
    f2 = (log(exp(a2[1]-max(a2)) + exp(-max(a2))) + max(a2)) / k;
    f3 = (log(exp(a3[1]-max(a3)) + exp(-max(a3))) + max(a3)) / k;
    f4 = (log(exp(a4[1]-max(a4)) + exp(-max(a4))) + max(a4)) / k;
    return f1 - f2 - f3 + f4;
  }

  real delta(real t, real dorm_s, real dorm_e, real k) {
    return dorm_e + (t - dorm_e) / (1 + exp(k * (t - dorm_e))) - dorm_s;
  }

  real traslation(real t, real dorm_start, real dorm_end) {
    real t_cn = t;
    if (t > dorm_start) t_cn = t_cn + dorm_end - dorm_start;
    return t_cn;
  }

  // ---- RNG helper functions (only call from generated quantities) ----

  int generate_clock_rng(int wgd, real mu_clock, real omega, real l_diploid,
                     real chemo_start, real t_eca, real t_mrca,
                     real t_dormancy_end, real shape_clock) {
    int m_clock_rep;
    if (!wgd) {
      m_clock_rep = neg_binomial_2_rng(
        2 * mu_clock * omega * l_diploid *
        ((chemo_start - t_eca) + (t_mrca - t_dormancy_end)) + 0.1,
        shape_clock
      );
    } else {
      m_clock_rep = -1; // placeholder
    }
    return m_clock_rep;
  }

  int generate_alpha_rng(int extra_therapy, int coeff_alpha, real omega, real l_CNA,
                     real mu_clock, real t_cna_tr, real t_eca,
                     array[] real mu_th_step, array[] real start_th_step,
                     array[] real end_th_step,
                     real chemo_start, real t_dormancy_end, real k_step,
                     real shape_cna, array[] real lambda_alpha, int n_step_th) {
    int m_alpha_rep;
    if (extra_therapy) {

      //for (th in 1:n_step_th){
      // m_alpha_rep = neg_binomial_2_rng(lambda_alpha[1]+lambda_alpha[2],shape_cna);
      // *** ORIGINAL CLL MODEL
      // m_alpha_rep = neg_binomial_2_rng(lambda_alpha[1]+lambda_alpha[2],shape_cna);
      m_alpha_rep = neg_binomial_2_rng(lambda_alpha[1]+lambda_alpha[2]+lambda_alpha[3],shape_cna);
      // *** TO HERE
      //}

    } else {
      m_alpha_rep = neg_binomial_2_rng(
        coeff_alpha * omega * l_CNA * (mu_clock * (t_cna_tr - t_eca)),
        shape_cna
      );
    }
    return m_alpha_rep;
  }

  int generate_beta_rng(int extra_therapy, int coeff_beta, real omega, real l_CNA,
                    real mu_clock, real t_cna_tr, real t_mrca,
                    array[] real mu_th_step, array[] real start_th_step,
                    array[] real end_th_step,
                    real chemo_start, real t_dormancy_end, real k_step,
                    real shape_cna, array[] real lambda_beta, int n_step_th) {
    int m_beta_rep;
    if (extra_therapy) {
      //for (th in 1:n_step_th){
      // *** ORIGINAL CLL MODEL
      // m_beta_rep = neg_binomial_2_rng(lambda_beta[1]+lambda_beta[2],shape_cna);
      m_beta_rep = neg_binomial_2_rng(lambda_beta[1]+lambda_beta[2]+lambda_beta[3],shape_cna);
      // *** TO HERE
      //}
    } else {
      m_beta_rep = neg_binomial_2_rng(
        coeff_beta * omega * l_CNA *
        (mu_clock * (t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr)),
        shape_cna
      );
    }
    return m_beta_rep;
  }


  int generate_th_rng(int extra_therapy,
                  real mu_th_step, real omega, real l_diploid,
                  real t_eca, real t_mrca,
                  real start_th_step, real end_th_step,
                  real k_step, real shape_step_1,  int n_step_th) {
    int m_th_step_rep;
    if (extra_therapy) {
      //for (th in 1:n_step_th){
      m_th_step_rep = neg_binomial_2_rng(
        2 * mu_th_step * omega * l_diploid *
        lambda_therapy(t_eca, t_mrca, start_th_step, end_th_step, k_step) + 0.1,
        shape_step_1
      );
      //}
    } else {
      m_th_step_rep = -1; // placeholder
    }
    return m_th_step_rep;
  }

  // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***
  int generate2_th_rng(int extra_therapy,
                     real mu_th_step, real omega, real l_diploid,
                     real t_eca, real t_mrca,
                     real start_th_step_1, real start_th_step_2,
                     real end_th_step_1, real end_th_step_2,
                     real k_step, real shape_step_1, int n_step_th) {
  int m_th_step_rep;
  real lambda_sum;
  if (extra_therapy) {
    lambda_sum = lambda_therapy(t_eca, t_mrca, start_th_step_1, end_th_step_1, k_step)
               + lambda_therapy(t_eca, t_mrca, start_th_step_2, end_th_step_2, k_step);

    m_th_step_rep = neg_binomial_2_rng(
      2 * mu_th_step * omega * l_diploid * lambda_sum + 0.1,
      shape_step_1
    );
  } else {
    m_th_step_rep = -1;
  }
  return m_th_step_rep;
}


  // ** UNTIL HERE
}

data{

 int <lower=0,upper=1> wgd;
 int <lower=0> m_clock_primary;
 int <lower=0,upper=1> n_clock;
 array[n_clock] int <lower=0> m_clock;
 real <lower=0> l_diploid;
 real <lower=0> mu_clock;

 int <lower=0> m_alpha;
 int <lower=0> m_beta;
 real <lower=0> l_CNA;
 int <lower=0> coeff_alpha; // 1 for 2:0, 2 for 2:2
 int <lower=0> coeff_beta; // 2 for 2:0, 4 for 2:2

 // FOR SIMPLICITY AT THE MOMENT ONLY A SINGLE CYCLE FOR EACH THERAPY IS ALLOWED
 int <lower=0,upper=1> extra_therapy;
 int <lower=0> n_th_step; // numero totale di terapie*cicli
 //int <lower=0> n_th_step_type; // numero di tipi di terapia
 array[n_th_step] real<lower=0> start_th_step;
 array[n_th_step] real<lower=0> end_th_step;
 //array[n_th_step] int<lower=0> type_th_step; // vector with numbers identifying the therapy (1:n_th_step)

 // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***
 //array[n_th_step] real<lower=0> alpha_th_step;
 //array[n_th_step] real<lower=0> beta_th_step;
 //array[n_th_step] int<lower=0> m_th_step;

 array[2] real<lower=0> alpha_th_step;
 array[2] real<lower=0> beta_th_step;
 array[2] int<lower=0> m_th_step;

 // ** UNTIL HERE **

	real <lower=0> omega_alpha;
  real <lower=0> omega_beta;
  real <lower=0> k_step;

  real <lower=0> Sample_1;
  real <lower=0> Sample_2;
  real <lower=0> chemo_start;
  real <lower=0> chemo_end;
  real <lower=0> last_therapy;

  array[2] int <lower=0, upper=1> exponential_growth;
  array[2] real<lower=0> N_min;
  array[2] real<lower=0> N_max;

  real <lower=0> alpha_mrca;
  real <lower=0> beta_mrca;

  // Overdispersion parameters for all mutation types
  real <lower=0> phi_clock;
  // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***
  // array[n_th_step] real <lower=0> phi_th_step;
  array[2] real <lower=0> phi_th_step;
  // ** UNTIL HERE
  real <lower=0> phi_cna;
}

parameters{

  real <lower=0, upper= start_th_step[1]> t_eca;
  real <lower=t_eca, upper= Sample_1> t_mrca_primary;
  // real <lower= chemo_start, upper= chemo_end> chemo_start;
  // *** CLL ORIGINAL MODEL
  //real <lower= chemo_end, upper= last_therapy> t_dormancy_end;
  real <lower= chemo_end+30/365, upper= last_therapy> t_dormancy_end;
  // UNTIL HERE
  real <lower= last_therapy, upper= Sample_2> t_mrca;
  real <lower= t_eca, upper= t_mrca - (t_dormancy_end - chemo_start)> t_cna_tr;

  real <lower= 0> omega;

  // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***
  // array[n_th_step] real<lower=0> mu_th_step;
  array[2] real<lower=0> mu_th_step;
  // *** UNTIL HERE
}

transformed parameters{
  // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***

  // array[2] real lambda_alpha;
  // array[2] real lambda_beta;
  //
  // lambda_alpha[1] =  coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca) + mu_th_step[1]*lambda_therapy(t_eca, t_cna_tr, start_th_step[1],end_th_step[1], k_step)) +.1;
  // lambda_alpha[2] =  coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca) + mu_th_step[2]*lambda_therapy(t_eca, t_cna_tr, start_th_step[2] - delta(start_th_step[2], chemo_start, t_dormancy_end, k_step),end_th_step[2] - delta(end_th_step[2], chemo_start, t_dormancy_end, k_step), k_step)) +.1;
  //
  // lambda_beta[1] = coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr ) + mu_th_step[1]*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - chemo_start), start_th_step[1], end_th_step[1],k_step)) +.1;
  // lambda_beta[2] = coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr ) + mu_th_step[2]*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - chemo_start), start_th_step[2] - delta(start_th_step[2],chemo_start,t_dormancy_end,k_step), end_th_step[2] - delta(end_th_step[2],chemo_start,t_dormancy_end,k_step),k_step)) +.1;
  //

  array[3] real lambda_alpha;
  array[3] real lambda_beta;

  lambda_alpha[1] =  coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca) +
                     mu_th_step[1]*lambda_therapy(t_eca, t_cna_tr, start_th_step[1],end_th_step[1], k_step)) +.1;
  lambda_alpha[2] =  coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca) +
                     mu_th_step[1]*lambda_therapy(t_eca, t_cna_tr, start_th_step[2],end_th_step[2], k_step)) +.1;
  lambda_alpha[3] =  coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca) +
                     mu_th_step[2]*lambda_therapy(t_eca, t_cna_tr, start_th_step[3] - delta(start_th_step[3], chemo_start, t_dormancy_end, k_step),end_th_step[3] - delta(end_th_step[3], chemo_start, t_dormancy_end, k_step), k_step)) +.1;

  lambda_beta[1] = coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr ) + mu_th_step[1]*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - chemo_start), start_th_step[1], end_th_step[1],k_step)) +.1;
  lambda_beta[2] = coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr ) + mu_th_step[1]*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - chemo_start), start_th_step[2], end_th_step[2],k_step)) +.1;
  lambda_beta[3] = coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr ) + mu_th_step[2]*lambda_therapy(t_cna_tr,t_mrca-(t_dormancy_end - chemo_start), start_th_step[3] - delta(start_th_step[3],chemo_start,t_dormancy_end,k_step), end_th_step[3] - delta(end_th_step[3],chemo_start,t_dormancy_end,k_step),k_step)) +.1;

   // *** UNTIL HERE
}

model{
  t_eca ~ uniform(0, start_th_step[1]);
  t_mrca_primary ~ uniform(t_eca, Sample_1);
  // chemo_start~ uniform(chemo_start, chemo_end);
  // *** CLL  ORIGINAL MODEL
  // t_dormancy_end~ uniform(chemo_end, last_therapy);
  // ***  UNTIL HERE
  t_dormancy_end~ uniform(chemo_end+30/365, last_therapy);
  t_mrca ~ uniform(last_therapy, Sample_2);
  t_cna_tr ~ uniform(t_eca, t_mrca - (t_dormancy_end - chemo_start) );

  omega ~ gamma(omega_alpha,omega_beta);
  // *** CLL ORIGINAL MODEL
  //for (m in 1:n_th_step){
    //mu_th_step[m] ~ gamma(alpha_th_step[m], beta_th_step[m]);
  //}
  //for (m in 1:n_th_step){
    mu_th_step[1] ~ gamma(alpha_th_step[1], beta_th_step[1]);
    mu_th_step[2] ~ gamma(alpha_th_step[2], beta_th_step[2]);
  //}
  // *** TO HERE

  print("t_eca bounds: [0, ", start_th_step[1], "]");
  print("t_mrca_primary bounds: [", t_eca, ", ", Sample_1, "]");
  print("t_dormancy_end bounds: [", chemo_end, ", ", last_therapy, "]");
  print("t_mrca bounds: [", t_dormancy_end, ", ", Sample_2, "]");
  print("t_cna_tr bounds: [", t_eca, ", ", t_mrca - (t_dormancy_end - chemo_start), "]");
  print("sampled omega: ", omega);
  print("alpha, beta 1: ", alpha_th_step[1], ",", beta_th_step[1]);
  print("alpha, beta 2: ", alpha_th_step[2], ",", beta_th_step[2]);
  print("sampled mu_sbsB: ", mu_th_step[1]);
  print("sampled mu_sbsGan: ", mu_th_step[2]);

  real shape_clock = 1 / phi_clock;
  // if (extra_therapy){
  // real shape_step = 1 / phi_th_step[1];
  // }
  real shape_cna = 1 / phi_cna;

  m_clock_primary ~ neg_binomial_2(
    2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1,
    shape_clock
    );
  print("sampled m_clock_primary: ",2*l_diploid*omega*mu_clock*(t_mrca_primary-t_eca) +.1);

  if (!wgd){
  m_clock[1] ~  neg_binomial_2(
    2*mu_clock*omega*l_diploid*((chemo_start - t_eca) + (t_mrca-t_dormancy_end)) +.1,
    shape_clock
  );
  print("sampled m_clock: ", 2*mu_clock*omega*l_diploid*((chemo_start - t_eca) + (t_mrca-t_dormancy_end)) +.1);
  }

  if (extra_therapy){

    // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***

    // for (th in 1:n_th_step){
    // m_th_step[th]  ~  neg_binomial_2(
    //   2*mu_th_step[th]*omega*l_diploid*lambda_therapy(t_eca,t_mrca, start_th_step[th], end_th_step[th],k_step) +.1,
    //   1/phi_th_step[th]
    //   );
    // }

    m_th_step[1]  ~  neg_binomial_2(
      2*mu_th_step[1]*omega*l_diploid*(
        lambda_therapy(t_eca,t_mrca, start_th_step[1], end_th_step[1],k_step)+
        lambda_therapy(t_eca,t_mrca, start_th_step[2], end_th_step[2],k_step)
        ) +.1,
      1/phi_th_step[1]
      );
      print("sampled m_th_step[1]: ",2*mu_th_step[1]*omega*l_diploid*(lambda_therapy(t_eca,t_mrca, start_th_step[1], end_th_step[1],k_step)+lambda_therapy(t_eca,t_mrca, start_th_step[2], end_th_step[2],k_step)) +.1);

      m_th_step[2]  ~  neg_binomial_2(
      2*mu_th_step[2]*omega*l_diploid*lambda_therapy(t_eca,t_mrca, start_th_step[3], end_th_step[3],k_step) +.1,
      1/phi_th_step[2]
      );
      print("sampled m_th_step[2]: ",2*mu_th_step[2]*omega*l_diploid*lambda_therapy(t_eca,t_mrca, start_th_step[3], end_th_step[3],k_step));



    // m_alpha ~  neg_binomial_2(lambda_alpha[1]+lambda_alpha[2],shape_cna);
    // m_beta ~ neg_binomial_2(lambda_beta[1]+lambda_beta[2],shape_cna);

    m_alpha ~  neg_binomial_2(lambda_alpha[1]+lambda_alpha[2]+lambda_alpha[3],shape_cna);
    m_beta ~ neg_binomial_2(lambda_beta[1]+lambda_beta[2]+lambda_beta[3],shape_cna);
    print("alpha: ",lambda_alpha[1]+lambda_alpha[2]+lambda_alpha[3]);
    print("beta: ",lambda_beta[1]+lambda_beta[2]+lambda_beta[3]);

    // ** UNTIL HERE


  }else{
    m_alpha ~  neg_binomial_2(coeff_alpha*omega*l_CNA*(mu_clock*(t_cna_tr - t_eca)) +.1+.1,
                                  shape_cna);
    m_beta ~ neg_binomial_2(coeff_beta*omega*l_CNA*(mu_clock*( t_mrca - (t_dormancy_end - chemo_start) - t_cna_tr )) +.1,
                                shape_cna);

  }

  // if (exponential_growth==1){
  //   target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
  //   target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  // }

  if (exponential_growth[1]==1){
    target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    // target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }

  if (exponential_growth[2]==1){
    // target += -N_min[1]*exp(-omega*(Sample_1 - t_mrca_primary)) + log(1-exp(-(N_max[1]-N_min[1])*exp(-omega*(Sample_1 - t_mrca_primary))));
    target += -N_min[2]*exp(-omega*(Sample_2 - t_mrca)) + log(1-exp(-(N_max[2]-N_min[2])*exp(-omega*(Sample_2 - t_mrca))));
  }
}

generated quantities {
  real t_cna;
  real shape_clock = 1 / phi_clock;
  real shape_cna   = 1 / phi_cna;

  // safe shape for step[1]
  real shape_step_1;
  if (extra_therapy == 1 && n_th_step > 0) {
    shape_step_1 = 1 / phi_th_step[1];
  } else {
    shape_step_1 = 1; // dummy
  }

  int m_clock_primary_rep = neg_binomial_2_rng(
    2 * l_diploid * omega * mu_clock * (t_mrca_primary - t_eca) + 0.1,
    shape_clock
  );

  int m_clock_rep = generate_clock_rng(
    wgd, mu_clock, omega, l_diploid, chemo_start, t_eca, t_mrca, t_dormancy_end, shape_clock
  );

  int m_alpha_rep = generate_alpha_rng(
    extra_therapy, coeff_alpha, omega, l_CNA, mu_clock, t_cna_tr, t_eca,
    mu_th_step, start_th_step, end_th_step,
    chemo_start, t_dormancy_end, k_step, shape_cna, lambda_alpha, n_th_step
  );

  int m_beta_rep = generate_beta_rng(
    extra_therapy, coeff_beta, omega, l_CNA, mu_clock, t_cna_tr, t_mrca,
    mu_th_step, start_th_step, end_th_step,
    chemo_start, t_dormancy_end, k_step, shape_cna, lambda_beta, n_th_step
  );

  // *** UNCOMMENT HERE TO HAVE THE CLL MODEL WITH ONLY BENDAMUSTINE ***

  // array[n_th_step] int m_th_step_rep;
  array[2] int m_th_step_rep;

  // for (th in 1:n_th_step){
  // m_th_step_rep[th] = generate_th_rng(
  //   extra_therapy, mu_th_step[th], omega, l_diploid, t_eca, t_mrca,
  //   start_th_step[th], end_th_step[th], k_step, shape_step_1, n_th_step
  // );
  // }

   m_th_step_rep[1] = generate2_th_rng(
     extra_therapy, mu_th_step[1], omega, l_diploid, t_eca, t_mrca,
     start_th_step[1], start_th_step[2], end_th_step[1], end_th_step[2], k_step, shape_step_1, n_th_step
     );

   m_th_step_rep[2] = generate_th_rng(
     extra_therapy, mu_th_step[2], omega, l_diploid, t_eca, t_mrca,
     start_th_step[3], end_th_step[3], k_step, shape_step_1, n_th_step
     );

  // *** UNTIL HERE

  t_cna = traslation(t_cna_tr, chemo_start, t_dormancy_end);

  // if you want these always present
  real N_sample_1 = exp(omega * (Sample_1 - t_eca));
  real N_sample_2 = exp(omega * (Sample_2 - t_mrca));
}



