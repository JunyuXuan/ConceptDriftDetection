function out=snrnd(mu,sigma,lambda,n)

% sample n variables from a skew normal distribution

delta=lambda/sqrt(lambda^2+1);
u0=randn(n,1);
u1=delta*u0+sqrt(1-delta^2)*randn(n,1);

out=mu+sigma*sign(u0).*u1;
