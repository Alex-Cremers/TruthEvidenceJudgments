data {
  int N;
  vector[N] evidence_zero;
  vector[N] evidence_one;
  vector[N] scaled_evidence;
  vector[N] truth;
  int<lower=0,upper=1> y[N];
}

parameters {
  real alpha; // intercept
  real beta_truth; // truth coef
  real beta_ev; // evidence coef
  real beta_inter; // truth x evidence coef
  real offset_z; // extra negative offset when evidence is 0
  real<lower=0> offset_o; // extra positive offset when evidence is 1
}

model {
  alpha ~ normal(0,5);
  beta_truth ~ normal(0,.5);
  beta_ev ~ normal(0,.1);
  beta_inter ~ normal(0,.1);
  offset_z ~ normal(0,.5);
  offset_o ~ normal(4,5);
  y ~ bernoulli_logit(alpha + beta_truth * truth + beta_ev * scaled_evidence + beta_inter * truth .* scaled_evidence  - offset_z*evidence_zero  + offset_o*evidence_one);
}




