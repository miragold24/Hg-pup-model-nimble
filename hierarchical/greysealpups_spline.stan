functions {
  vector skewnormBDFunction(real mu, real sigma, real alpha, int nsteps, vector days){
    vector[nsteps] xvec;
    // Map mean/sd/shape -> (xi, omega, alpha) used by Stan's skew_normal
    real a2 = square(alpha);
    real c  = sqrt(2 / pi());
    real denom = sqrt(1 - (2/pi()) * (a2 / (1 + a2)));
    real omega_skn = sigma / denom;
    real xi_skn    = mu - omega_skn * (alpha / sqrt(1 + a2)) * c;

    for (i in 1:nsteps)
      xvec[i] = exp(skew_normal_lpdf(days[i] | xi_skn, omega_skn, alpha));
    return xvec / sum(xvec);
  }

  vector hazFunction(real mu, real std, int nsteps, int dstep) {
    real prob;
    vector[nsteps] hazprob;
    hazprob[1] = normal_cdf(1 * dstep | mu, std);
    for (s in 2:nsteps) {
      prob = normal_cdf(s * dstep | mu, std) - normal_cdf((s - 1) * dstep | mu, std);
      hazprob[s] = (1 - normal_cdf((s - 1) * dstep | mu, std)) < 1e-6
                     ? 1.0
                     : prob / (1 - normal_cdf((s - 1) * dstep | mu, std));
    }
    return hazprob;
  }

  matrix transMatrix(real xphi, real xgamma, real xzeta, real xchi, int stages, int dstep) {
    matrix[stages, stages] matarray = rep_matrix(0, stages, stages);
    // from state 1
    matarray[1, 1] = pow(xphi, dstep) * (1 - xgamma) * (1 - xzeta);
    matarray[2, 1] = pow(xphi, dstep) * xgamma * (1 - xzeta);
    matarray[4, 1] = 1 - pow(xphi, dstep);
    // from state 2
    matarray[2, 2] = pow(xphi, dstep) * (1 - xzeta);
    matarray[3, 2] = pow(xphi, dstep) * xzeta;
    matarray[4, 2] = 1 - pow(xphi, dstep);
    // from state 4
    matarray[4, 4] = pow(xchi, dstep);
    matarray[5, 4] = 1 - pow(xchi, dstep);
    // absorbing
    matarray[3, 3] = 1;
    matarray[5, 5] = 1;
    return matarray;
  }
}

data {
  int<lower=1> nsurveys;
  int<lower=1> nstates;
  int<lower=1> ndays;
  int<lower=1> nobsstates;
  int<lower=1> nvisits;
  int<lower=1> ncols;
  int<lower=1> nyears;

  int<lower=1> K;                   // smooth basis dimension
  matrix[ncols, K] X_spline;
  matrix[K, K] S1_block;
  matrix[K, K] S2_block;

  array[nsurveys] int<lower=1, upper=ncols> colid_svy;
  array[nsurveys] int<lower=1, upper=nyears> year_svy;
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
  for (v in 1:nvisits)
    obsdays_bin[v] = min(obsdays_svy[v] %/% daysPerBin + 1, nIntervals_s);

  int nIntervals_c = min(20, nIntervals_s);
  vector[nIntervals_c] cohort_days;
  for (z in 1:nIntervals_c)
    cohort_days[z] = z * daysPerBin - 0.5 * daysPerBin;

  real xchi = inv_logit(2.5);
  real muMoult = 23;
  real sdMoult = 5;
  vector[nIntervals_s] xgamma = hazFunction(muMoult, sdMoult, nIntervals_s, daysPerBin);
}

parameters {
  vector[ncols] nborn_col_raw;

  real muBday_global;
  real<lower=0> sigma_col;
  vector[ncols] z_col;


  vector[K] beta_spline_raw;
  vector<lower=1e-8>[2] lambda_spline;   // force positive

  // Birth schedule 
  real<lower=0> sdBday;     // true SD of birth dates within season
  real skewBday;            // shape (unconstrained)

  // Leave distribution (centered)
  real muLeave_raw;         // lognormal mean on log scale
  real<lower=0> sdLeave;    // centered

  // Survival
  real muSurvInt;
  real<lower=0> sdSurvInt;
  vector[ncols] survInt_col_raw;

  // Observation probabilities
  real<lower=0.7, upper=1> pobs_d_est;
  real<lower=0.7, upper=1> pobs_w_est;
  real<lower=0.7, upper=1> pobs_m_est;
  real<lower=0.7, upper=1> pcm_est;
}

transformed parameters {
  vector[ncols] nborn_col = exp(5.0 + 2.0 * nborn_col_raw);

  vector[ncols] col_eff = sigma_col * z_col;

  // Leave params + hazard
  real muLeave = exp(3.5 + 0.1 * muLeave_raw);
  vector[nIntervals_s] xzeta = hazFunction(muLeave, sdLeave, nIntervals_s, daysPerBin);

  // Survival curve by day
  vector[ncols] survInt_col = muSurvInt + sdSurvInt * survInt_col_raw;
  array[ncols] vector[ndays] xphi;
  for (c in 1:ncols)
    for (i in 1:ndays)
      xphi[c, i] = inv_logit(survInt_col[c] + 0.1 * (i - 1));

  // Observation probabilities
  real pobs_d = est_pobs_d == 1 ? pobs_d_est : 0.85;
  real pobs_w = est_pobs_w == 1 ? pobs_w_est : 0.98;
  real pobs_m = est_pobs_m == 1 ? pobs_m_est : 0.91;
  real pcm    = est_pcm    == 1 ? pcm_est    : 0.84;
}

model {
  // Priors
  muBday_global ~ normal(45, 10);

  // smoothing penalties (half-Cauchy via lower=0)
    lambda_spline[1] ~ cauchy(0,1);     #DM suggested smoothing penalty prior
    lambda_spline[2] ~ cauchy(0,1);  


  nborn_col_raw ~ std_normal();
  z_col ~ std_normal();
  sigma_col ~ normal(0, 1);

  sdBday ~ gamma(0.1, 0.1);
  skewBday ~ normal(1, 1);

  muLeave_raw ~ std_normal();
  sdLeave ~ gamma(0.1, 0.1);

  muSurvInt ~ normal(2, 1);
  sdSurvInt ~ gamma(1, 1);
  survInt_col_raw ~ std_normal();

  if (est_pobs_d == 1) pobs_d_est ~ beta(35, 5);
  if (est_pobs_w == 1) pobs_w_est ~ beta(45, 5);
  if (est_pobs_m == 1) pobs_m_est ~ beta(45, 5);
  if (est_pcm == 1)    pcm_est    ~ beta(35, 5);

  // Spline prior
  
  matrix[K, K] K1 = lambda_spline[1] * S1_block + lambda_spline[2] * S2_block;
  vector[K] zero = rep_vector(0,K);
  beta_spline_raw ~ multi_normal_prec(zero,K1);
  vector[ncols] beta_spline_mu = X_spline*beta_spline_raw ;  

  // Likelihood
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

  vector[nstates] tempCohortVec;
  array[nIntervals_s] vector[nstates] stageN;
  row_vector[nobsstates + 1] lambda_obs_vector;

  int pos = 1;
  for (s in 1:nsurveys) {
    vector[nIntervals_c] Pb = skewnormBDFunction(
      muBday_global + col_eff[colid_svy[s]] + beta_spline_mu[colid_svy[s]],
      sdBday, skewBday, nIntervals_c, cohort_days);

    vector[nIntervals_c] cohortSizeBin = nborn_col[colid_svy[s]] * Pb;

    array[nIntervals_s] matrix[nstates, nstates] matarray;
    for (i in 1:nIntervals_s)
      matarray[i] = transMatrix(xphi[colid_svy[s]][(i - 1) * daysPerBin + 1],
                                xgamma[i], xzeta[i], xchi, nstates, daysPerBin);

    for (a in 1:nIntervals_s) stageN[a] = rep_vector(0, nstates);

    for (h in 1:nIntervals_c) {
      tempCohortVec = rep_vector(0, nstates);
      tempCohortVec[1] = cohortSizeBin[h];
      stageN[h] += tempCohortVec;
      for (a in 1:(nIntervals_s - h)) {
        tempCohortVec = matarray[a] * tempCohortVec;
        stageN[h + a] += tempCohortVec;
      }
    }

    for (v in 1:nvisits_svy[s]) {
      lambda_obs_vector = to_row_vector(stageN[obsdays_bin[pos]]) * obsMatrix;
      obs_vector_w[pos] ~ poisson(lambda_obs_vector[1]);
      obs_vector_m[pos] ~ poisson(lambda_obs_vector[2]);
      obs_vector_d[pos] ~ poisson(lambda_obs_vector[3]);
      pos += 1;
    }
  }
}
