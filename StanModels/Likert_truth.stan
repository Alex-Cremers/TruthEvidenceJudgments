data {
  int N;
  vector[N] evidence_zero;
  vector[N] evidence_one;
  vector[N] scaled_evidence;
  vector[N] truth;
  int<lower=1,upper=7> y[N]; 
}

parameters {
  real<lower=0> beta_truth; // truth coef
  real beta_ev; // evidence coef
  real beta_inter; // truth x evidence coef
  real offset_z; // extra negative offset when evidence is 0
  real offset_o; // extra positive offset when evidence is 1
  ordered[6] c; // cutpoints for the ordinal regression
}

transformed parameters {
  vector[N] predictor;
  predictor = beta_truth * truth + beta_ev * scaled_evidence + beta_inter * truth .* scaled_evidence  - evidence_zero*offset_z  + evidence_one*offset_o;
}

model {
  beta_truth ~ normal(5,5);
  beta_ev ~ normal(0,.1);
  beta_inter ~ normal(0,.1);
  offset_z ~ normal(0,.5);
  offset_o ~ normal(0,.5); 
  for (i in 1:N) {
    y[i] ~ ordered_logistic(predictor[i], c);
  }
}






