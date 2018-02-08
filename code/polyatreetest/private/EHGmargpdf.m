function pdf = EHGmargpdf(x, N, m, n, a, b)

% Marginal pdf for extended hypergeometric distribution
% Evaluated using Monte Carlo approximation

nsamples = 2000;
theta = betarnd(a, b, nsamples, 2);
% betarnd produces NaN sometimes with low a and b, so the next step
% sets Nan to 0 or 1 w.p. a/(a+b)
ind_nan = isnan(theta);
if sum(ind_nan)>0
    theta(ind_nan) = (rand(sum(ind_nan), 1)<(a/(a+b)));
end
pdf = mean(EHGpdf(x, N, m, n, theta(:, 1), theta(:, 2)));

end