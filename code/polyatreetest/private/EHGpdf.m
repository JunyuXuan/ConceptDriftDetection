function pdf = EHGpdf(x, N, m, n, theta1, theta2)

% pdf of Fisher non central hypergeometric distribution

n_theta = length(theta1);

n1 = m;
n2 = N - m;
n0 = repmat(max(0, n+m-N):min(m, n), n_theta, 1);
n_n0 = size(n0, 2);
theta1 = repmat(theta1, 1, n_n0);
theta2 = repmat(theta2, 1, n_n0);


logout = logbinopdf(n0, n1, theta1) + logbinopdf(n - n0, n2, theta2);

bias = max(logout, [], 2);
out = exp(logout - repmat(bias, 1, n_n0));
out = out./repmat(sum(out,2), 1, n_n0);
pdf = out(n0 == x);
pdf(isnan(pdf)) = 0;
end


function out = logbinopdf(n0, n1, theta)
    out = gammaln(n1+1) - gammaln(n0 +1) - gammaln(n1-n0+1) + n0.*log(theta) + (n1-n0).*log(1-theta);
    out(theta==1) = log(double(n0(theta==1)==n1));
    out(theta==0) = log(double(n0(theta==0)==0));
end