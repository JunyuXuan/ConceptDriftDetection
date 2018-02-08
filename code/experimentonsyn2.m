% experiments on synthetic data
clear

rng(1);

mu = 0.5;
step = 0.04;
datnum = 50;

allresults = zeros(4, 21);

for ii = 1:1000
    
results = zeros(4, 10);

for t = 1 : 21
   
    if mod(t, 2) == 0
       i = t/2;
       data1 = normrnd(mu,0.1 + step*(i-1), 1, datnum);
       data2 = normrnd(mu,0.1 + step*(i+1-1), 1, datnum);
       [driftKS, p, KSstatistic] = kstest2(data1, data2, 'alpha', 0.05);
       [drift, post, stats]      = PTtest(data1', data2', 'normalize', true);
       %[drift, post, stats]      = PTtest(data1, data2, 'partition', 'empirical');
       results(:, t) = [driftKS KSstatistic drift post];
    else
       i = (t+1)/2;
       data1 = normrnd(mu,0.1 + step*(i-1), 1, datnum);
       data2 = normrnd(mu,0.1 + step*(i-1), 1, datnum);
       [driftKS, p, KSstatistic] = kstest2(data1, data2, 'alpha', 0.05);
       [drift, post, stats]      = PTtest(data1', data2', 'normalize', true);
       %[drift, post, stats]      = PTtest(data1, data2, 'partition', 'empirical');
       results(:, t) = [driftKS KSstatistic drift post];
    end
    
end

results(4, :) = 1 - results(4, :) ;

allresults = allresults + results;

end

allresults = allresults/1000;

xlabel = 1:21;

plot(xlabel, allresults(2, :)')
hold on
plot(xlabel, allresults(4, :)')




