param N>= 0; 
param SO>= 0;
param PI {i in 1..N};
param lambda {i in 1..N};
param d {i in 1..N, k in 1..PI[i]};
param c {i in 1..N, j in 1..SO, k in 1..PI[i]};
param coef {i in 1..N, j in 1..SO, k in 1..PI[i]};
param ch {j in 1..SO};

var w {j in 1..SO} >= 0;

minimize Barycenter: sum{j in 1..SO} ch[j]*w[j];

subject to Offer {i in 1..N,k in 1..PI[i]}:
sum {j in 1..SO} coef[i,j,k]*w[j] = d[i,k];
