functions {
  vector skewnormBDFunction(real mu, real sigma, real alpha, int nsteps, vector days) {
    vector[nsteps] xvec;
    for (i in 1:nsteps) {
      xvec[i] = exp(skew_normal_lpdf(days[i] | mu, sigma, alpha));
    }
    return xvec / sum(xvec);
  }

  vector hazFunction(real mu, real std, int nsteps, int dstep) {
    real prob;
    vector[nsteps] hazprob;
    hazprob[1] = normal_cdf((1 * dstep) | mu, std);
    for (s in 2:nsteps) {
      prob = normal_cdf((s * dstep) | mu, std) - normal_cdf(((s - 1) * dstep) | mu, std);
      hazprob[s] = (1 - normal_cdf(((s - 1) * dstep) | mu, std)) < 0.000001
                     ? 1.0
                     : prob / (1 - normal_cdf(((s - 1) * dstep) | mu, std));
    }
    return hazprob;
  }

  matrix transMatrix(real xphi, real xgamma, real xzeta, real xchi, int stages, int dstep) {
    matrix[stages, stages] matarray = rep_matrix(0, stages, stages);
    matarray[1, 1] = pow(xphi, dstep) * (1 - xgamma) * (1 - xzeta);
    matarray[2, 1] = pow(xphi, dstep) * xgamma * (1 - xzeta);
   // matarray[3, 1] = pow(xphi, dstep) * xzeta;
    matarray[4, 1] = 1 - pow(xphi, dstep);

    matarray[2, 2] = pow(xphi, dstep) * (1 - xzeta);
    matarray[3, 2] = pow(xphi, dstep) * xzeta;
    matarray[4, 2] = 1 - pow(xphi, dstep);

    matarray[4, 4] = pow(xchi, dstep);
    matarray[5, 4] = 1 - pow(xchi, dstep);

    matarray[3, 3] = 1;
    matarray[5, 5] = 1;

    return matarray;
  }
}

data {
  int<lower=0> nsurveys;
  int<lower=0> nstates;
  int<lower=0> ndays;
  int<lower=0> nobsstates;
  int<lower=0> nvisits;
  int<lower=0> ncols;

  int<lower=1> K;              // spline basis functions
  matrix[ncols, K] X_spline;   // spline basis matrix

  array[nsurveys] int<lower=1, upper=ncols> colid_svy;
  array[nsurveys] int<lower=0> nvisits_svy;

  array[nvisits] int obs_vector_w;
  array[nvisits] int obs_vector_m;
  array[nvisits] int obs_vector_d;
  array[nvisits] int<lower=0> obsdays_svy;

  int<lower=0, upper=1> est_pobs_d;
  int<lower=0, upper=1> est_pobs_w;
  int<lower=0, upper=1> est_pobs_m;
  int<lower=0, upper=1> est_pcm;
}

transformed data {
  int daysPerBin = 5;
  int nIntervals_s = ndays %/% daysPerBin;

  array[nvisits] int obsdays_bin;
  for (v in 1:nvisits) {
    obsdays_bin[v] = obsdays_svy[v] %/% daysPerBin + 1; // add 1 for indexing
  }

  int nIntervals_c = min(20, nIntervals_s);
  vector[nIntervals_c] cohort_days;
  for (z in 1:nIntervals_c) {
    cohort_days[z] = z * daysPerBin - 0.5 * daysPerBin;
  }

  vector<lower=0>[ndays] daysvec;
  for (i in 1:ndays) daysvec[i] = i;

  real xchi = inv_logit(3.5);

  real muMoult = 23;
  real sdMoult = 5;
  vector[nIntervals_s] xgamma = hazFunction(muMoult, sdMoult, nIntervals_s, daysPerBin);
}

parameters {
  vector[K] beta_spline;       // spline coefficients for muBday_col
   real<lower=0> tau_spline;
  real muBday_global;          // Global mean birth date

  vector[ncols] nborn_col_raw;

  real<lower=0> mu_sdBday;       // Global mean for birth date SD
  real<lower=0> sigma_sdBday;    // SD of the random effects (Cauchy scale)
  vector[ncols] sdBday_col_raw;  // Raw colony-level effects

  real<lower=0> skewBday;

  real muLeave_raw;
  real<lower=0> sdLeave_raw;

  real muSurvInt;
  real<lower=0> sdSurvInt;
  vector[ncols] survInt_col_raw;
  real<lower=0> survSlope;

  real<lower=0.7, upper=1> pobs_d_est;
  real<lower=0.7, upper=1> pobs_w_est;
  real<lower=0.7, upper=1> pobs_m_est;
  real<lower=0.7, upper=1> pcm_est;
}

transformed parameters {
  vector[ncols] nborn_col = exp(5.0 + 2.0 * nborn_col_raw);
  vector[ncols] muBday_col = muBday_global + X_spline * beta_spline;  // Global mean + spline 
  //vector[ncols] muBday_col = X_spline * beta_spline;  // spline for birth date
  vector[ncols] sdBday_col = mu_sdBday * exp(sigma_sdBday * sdBday_col_raw);  // Hierarchical birth date SD
  real muLeave = exp(3.5 + 0.1 * muLeave_raw);
  real sdLeave = 1 + 5.0 * sdLeave_raw;
  vector[nIntervals_s] xzeta = hazFunction(muLeave, sdLeave, nIntervals_s, daysPerBin);

  vector[ncols] survInt_col = muSurvInt + sdSurvInt * survInt_col_raw;

  real pobs_d = est_pobs_d == 1 ? pobs_d_est : 0.85;
  real pobs_w = est_pobs_w == 1 ? pobs_w_est : 0.98;
  real pobs_m = est_pobs_m == 1 ? pobs_m_est : 0.91;
  real pcm    = est_pcm == 1 ? pcm_est : 0.84;
}

model {
  nborn_col_raw ~ std_normal();
   muBday_global ~ normal(45, 10);  // Global mean birth date around day 45
  // beta_spline ~ normal(0, 1);  // tighter prior for spline coefficients
  // 
 // Prior on penalty strength (based on brms default)
tau_spline ~ student_t(3, 0, 35.6);

// First-order difference penalty
for (k in 2:K)
  target += normal_lpdf(beta_spline[k] - beta_spline[k-1] | 0, tau_spline);
  // Hierarchical priors for birth date SD
  mu_sdBday ~ normal(2, 1);        // Global mean birth date SD
  sigma_sdBday ~ cauchy(0, 2.5);   // Cauchy prior for random effect scale
  sdBday_col_raw ~ std_normal();   // Standard normal for raw effects
  
  skewBday ~ normal(0, 5);

  muLeave_raw ~ std_normal();
  sdLeave_raw ~ gamma(1, 1);

  muSurvInt ~ normal(2, 1);
  sdSurvInt ~ gamma(1, 1);
  survInt_col_raw ~ std_normal();
  survSlope ~ normal(0, 0.5);

  if (est_pobs_d == 1) pobs_d_est ~ beta(35, 5);
  if (est_pobs_w == 1) pobs_w_est ~ beta(45, 5);
  if (est_pobs_m == 1) pobs_m_est ~ beta(45, 5);
  if (est_pcm == 1)    pcm_est ~ beta(35, 5);

  array[ncols] vector[ndays] xphi;
  for (c in 1:ncols)
    for (d in 1:ndays)
      xphi[c][d] = inv_logit(survInt_col[c] + survSlope * d);

  matrix[nstates, nobsstates + 1] obsMatrix = rep_matrix(0, nstates, nobsstates + 1);
  obsMatrix[1, 1] = pobs_w;
  obsMatrix[1, 4] = 1 - pobs_w;
  obsMatrix[2, 1] = pobs_m * (1 - pcm);
  obsMatrix[2, 2] = pobs_m * pcm;
  obsMatrix[2, 4] = 1 - pobs_m;
  obsMatrix[3, 4] = 1;
  obsMatrix[4, 3] = pobs_d;
  obsMatrix[4, 4] = 1 - pobs_d;
  obsMatrix[5, 4] = 1;

  int pos = 1;
  for (s in 1:nsurveys) {
    vector[nIntervals_c] Pb = skewnormBDFunction(
      muBday_col[colid_svy[s]], sdBday_col[colid_svy[s]], skewBday, nIntervals_c, cohort_days);

    vector[nIntervals_c] cohortSizeBin = nborn_col[colid_svy[s]] * Pb;

    array[nIntervals_s] matrix[nstates, nstates] matarray;
    for (i in 1:nIntervals_s) {
      matarray[i] = transMatrix(
        xphi[colid_svy[s]][(i - 1) * daysPerBin + 1], xgamma[i], xzeta[i], xchi, nstates, daysPerBin);
    }

    array[nIntervals_s] vector[nstates] stageN;
    for (a in 1:nIntervals_s) stageN[a] = rep_vector(0, nstates);

    for (h in 1:nIntervals_c) {
      vector[nstates] tempCohortVec = rep_vector(0, nstates);
      tempCohortVec[1] = cohortSizeBin[h];
      stageN[h] += tempCohortVec;

      for (a in 1:(nIntervals_s - h)) {
        tempCohortVec = matarray[a] * tempCohortVec;
        stageN[h + a] += tempCohortVec;
      }
    }

    for (v in 1:nvisits_svy[s]) {
      row_vector[nobsstates + 1] lambda_obs_vector = to_row_vector(stageN[obsdays_bin[pos]]) * obsMatrix;

      obs_vector_w[pos] ~ poisson(lambda_obs_vector[1]);
      obs_vector_m[pos] ~ poisson(lambda_obs_vector[2]);
      obs_vector_d[pos] ~ poisson(lambda_obs_vector[3]);
      pos += 1;
    }
  }
}

// generated quantities {
//   array[nvisits] int yrep_w;
//   array[nvisits] int yrep_m;
//   array[nvisits] int yrep_d;
// 
//   int pos = 1;
// 
//   array[ncols] vector[ndays] xphi;
//   for (c in 1:ncols)
//     for (d in 1:ndays)
//       xphi[c][d] = inv_logit(muSurvInt + sdSurvInt * 0 + survSlope * d);  // Use mean surv for yrep
// 
//   for (s in 1:nsurveys) {
//     vector[nIntervals_c] Pb = skewnormBDFunction(
//       muBday_col[colid_svy[s]], sdBday_col[colid_svy[s]], skewBday, nIntervals_c, cohort_days);
// 
//     vector[nIntervals_c] cohortSizeBin = nborn_col[colid_svy[s]] * Pb;
// 
//     array[nIntervals_s] matrix[nstates, nstates] matarray;
//     for (i in 1:nIntervals_s) {
//       matarray[i] = transMatrix(
//         xphi[colid_svy[s]][(i - 1) * daysPerBin + 1], xgamma[i], xzeta[i], xchi, nstates, daysPerBin);
//     }
// 
//     array[nIntervals_s] vector[nstates] stageN;
//     for (a in 1:nIntervals_s) stageN[a] = rep_vector(0, nstates);
// 
//     for (h in 1:nIntervals_c) {
//       vector[nstates] tempCohortVec = rep_vector(0, nstates);
//       tempCohortVec[1] = cohortSizeBin[h];
//       stageN[h] += tempCohortVec;
// 
//       for (a in 1:(nIntervals_s - h)) {
//         tempCohortVec = matarray[a] * tempCohortVec;
//         stageN[h + a] += tempCohortVec;
//       }
//     }
// 
//     for (v in 1:nvisits_svy[s]) {
//       row_vector[nobsstates + 1] lambda_obs_vector = to_row_vector(stageN[obsdays_bin[pos]]) * obsMatrix;
// 
//       yrep_w[pos] = poisson_rng(lambda_obs_vector[1]);
//       yrep_m[pos] = poisson_rng(lambda_obs_vector[2]);
//       yrep_d[pos] = poisson_rng(lambda_obs_vector[3]);
// 
//       pos += 1;
//     }
//   }
// }
