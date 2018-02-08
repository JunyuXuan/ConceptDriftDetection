% experiments on synthetic data
clear

addpath('polyatreetest');

rng(1);

mu1 = 0.5;
sigma1 = 0.04;
mu2 = 0.5;
sigma2 = 0.04;
% mu2 = 0.55;
% sigma2 = 0.045;

allresults = zeros(4, 20);

trailnum = 1000;

for ii = 1 : trailnum

results = zeros(4, 20);

t = 1;
for datnum = 5 : 5: 100
    
    data1 = normrnd(mu1,sigma1, 1, datnum);
    data2 = normrnd(mu2,sigma2, 1, datnum);
    [driftKS, p, KSstatistic] = kstest2(data1, data2, 'alpha', 0.05);
    [drift, post, stats]      = PTtest(data1', data2', 'normalize', true);
    results(:, t) = [driftKS KSstatistic drift post];
    t = t + 1;    
end

results(4, :) = 1 - results(4, :) ;

allresults = allresults + results;

end

allresults = allresults/trailnum;

% plot 
xlabel = 5 : 5 : 100;
%xlabel = 50 : 50: 1000;

plot(xlabel, allresults(2, 1:length(xlabel))')
hold on
plot(xlabel, allresults(4, 1:length(xlabel))')




