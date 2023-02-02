data {
    int d;
    int n;
    int tp;
    int tn;
}
parameters {
   real<lower=0, upper=1>p;
   real<lower=0, upper=1>se;
   real<lower=0, upper=1>sp;
}
model {
    p ~ beta(1, 1);
    se ~ beta(1, 1);
    sp ~ beta(1, 1);
    d ~ binomial(n, p);
    tp ~ binomial(d, se);
    tn ~ binomial(n-d, sp);  
}
