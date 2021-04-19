function x = prox_l1(b,lambda)
x = max(0,b-lambda)+min(0,b+lambda);