function [theta_trial, polyatree, kolmo_smir, wilcoxon, LOR, h, h2] = powercurve(n, type, varargin)

%
% POWERCURVE computes power curves for the Polya tree, KS and Wilcoxon tests
%      [theta_trial, polyatree, kolmo_smir, wilcoxon, LOR, h, h2] = powercurve(n, type, varargin)
%
% Requires the statistics toolbox
%--------------------------------------------------------------------------
% INPUTS 
%   - n:    Number of samples from sample y1 and sample y2
%   - type: Type of power curve. Possible values are:
%           'mean': N(0,1) vs N(theta,1)
%           'var':  N(0,1) vs N(0, theta^2)
%           'mixture': N(0,1) vs 1/2 N(theta,1) + 1/2 N(-theta,1)
%           'skew': N(0,1) vs skewN(0,1,theta)
%           'tail': N(0,1) vs t(location=0,scale=1,df=1/theta)
%           'student_mean': t(0,1,4) vs t(theta,1,4)
%           'student_var':  t(0,1,4) vs t(0,theta,4)
%           'gamma':    Gamma(1,1) vs Gamma(shape=1/theta,inversescale=1/theta) 
%           'exponential':  Exp(1) vs Exp(theta)
%           'lognorm_mean': lognorm(0,1) vs lognorm(theta,1)
%           'lognorm_var': lognorm(0,1) vs lognorm(0,theta^2)
% Optional inputs:
%   - c:    Parameter c>0 of the Polya Tree (Default=1)
%   - threshold: Decision threshold between 0 and 1 (Default=0.5)
%   - normalize: binary variable to normalize or not the data (Default=true)
%   - estimate_c: logical if equal to 1, c is estimated with empirical
%               Bayes (Default=false)
%   - partition: partition of the Polya tree
%                either an object of the Matlab class ProbDist that defines the centering distribution of the partition
%                or 'empirical' to use the conditional test
%                (Default=normal(0,1))
%   - nsamples: Number of samples to obtain the power curves (default=100)
% 
% OUTPUTS
%   - post: P(H0 | y1,y2)
%   - h:    1 if H0 is rejected 0 otherwise
%   - stats: structure with the following fields:
%           - LOR:  Log ratio of marginal likelihoods
%           - LOR_level: Log ratio of marginal likelihoods at each level of the
%                        tree
%           - LM_H0: log-marginal likelihood under the null
%           - LM_H1: log-marginal likelihood under the alternative
%
% See also POWERCURVE
%--------------------------------------------------------------------------
% EXAMPLE
% y1 = randn(50, 1);
% y2 = randn(50, 1);
% [h, post, stats] = PTtest(y1, y2, 'estimate_c', true);
% [h, post, stats] = PTtest(y1, y2, 'estimate_c', true, 'partition', 'empirical');


%--------------------------------------------------------------------------
% Reference: C.C. Holmes, F. Caron, J. Griffin and D.A. Stephens. 
% Two-sample Bayesian nonparametric hypothesis testing. To appear in
% Bayesian Analysis, 2014. (First version: arXiv:0910.5060, 2009).
%
% Copyright INRIA, 2014
% Author: François Caron, University of Oxford
% caron@stats.ox.ac.uk
% June 2008; Last update July 2014
%--------------------------------------------------------------------------



param_names = {'c', 'threshold', 'normalize', 'estimate_c', 'partition', 'nsamples', 'doplot'};
param_default = {1, .5, true, false, [], 100, false};
[c, threshold, normalize, estimate_c, partition, nsamples, doplot] =...
    parse_options(varargin, param_names, param_default);

% Some checks on the inputs
if ~islogical(estimate_c)
    error('Option ''estimate_c'' should be of type logical');
end
if ~islogical(doplot)
    error('Option ''doplot'' should be of type logical');
end
if ~isnumeric(c) || c<0
    error('Option c should be a strictly positive scalar');
end

if nargin<4
    nsamples = 100;
end

isoctave = (exist('OCTAVE_VERSION')~=0); % Checks if octave running

switch(type)
    case {'mean', 'lognorm_mean', 'student_mean'}
        theta_trial = linspace(0, 2, 20);        
    case {'var', 'lognorm_var', 'student_var'}
        theta_trial= linspace(1, 5, 20);
    case 'mixture'
        theta_trial = linspace(0, 2, 20);
    case 'skew'
        theta_trial= linspace(0, 5, 20);
    case 'tail'
        theta_trial= linspace(10^-6, 10, 20);
    case 'gamma'
        theta_trial= linspace(1, 3, 20);
    case 'exponential'
        theta_trial= linspace(0, 3, 20);
    otherwise
        error('Unknown type')
end
n_trial = length(theta_trial);
polyatree = zeros(n_trial, 1);
LOR = zeros(n_trial, nsamples);
kolmo_smir = zeros(n_trial, 1);
wilcoxon = zeros(n_trial, 1);

hbar = waitbar(0,['Calculating the power curve for n=' num2str(n) ' and type=' type]);
for i=1:n_trial
    waitbar(i/n_trial,hbar,['Calculating the power curve for n=' num2str(n) ' and type=' type]);
    fprintf('Trial %d over %d\n',i, n_trial);
    for j=1:nsamples
        % Simulate data
        y1=randn(n,1);
        theta = theta_trial(i);
        switch(type)
            case 'mean'                
                y2 = theta + randn(n, 1);
            case 'var'
                y2 = theta * randn(n, 1);
            case 'mixture'
                u = rand(n, 1)>.5;
                y2 = u*(theta) - (1-u)*theta  + randn(n,1);
            case 'skew'
                y2 = snrnd(0, 1, theta,n);   
            case 'tail'
                y2 = trnd(1/theta, n,1);
            case 'student_mean'
                y1 = trnd(4, n,1);
                y2 = theta + trnd(4, n,1);
            case 'student_var'
                y1 = trnd(4, n,1);
                y2 = theta * trnd(4, n,1);
            case 'gamma'
                y1 = gamrnd(1,1,n,1);
                y2 = gamrnd(1/theta,theta,n,1);
            case 'exponential'
                y1 = exprnd(1,n,1);
                y2 = exprnd(1 + theta,n,1);
            case 'lognorm_mean'
                y1 = lognrnd(5, 1, n, 1);
                y2 = lognrnd(5 + theta, 1, n, 1);
            case 'lognorm_var'
                y1 = lognrnd(5, 1, n, 1);
                y2 = lognrnd(5, theta, n, 1);
        % Normalize data
        end
        
        % Kolmogorov-Smirnov test
        if ~isoctave
            % Matlab
            kolmo_smir(i) = kolmo_smir(i) + kstest2(y1,y2);
            [junk, temp] =  ranksum(y1,y2);
            % Wilcoxon test
            wilcoxon(i) = temp + wilcoxon(i);
        else
            % Octave
            pval = kolmogorov_smirnov_test_2(y1, y2);
            kolmo_smir(i) = kolmo_smir(i) + pval<.05;
            pval = u_test(y1, y2);
            wilcoxon(i) = wilcoxon(i) + pval<.05;
        end
        
        [~, ~, stats] = PTtest(y1, y2, 'c', c, 'threshold', threshold,...
            'normalize', normalize, 'estimate_c', estimate_c,...
            'partition', partition);
        LOR(i,j) = stats.LOR;        
    end
    
    % Find the threshold
    if i==1
        thres = quantile(LOR(1,:), .05, 2);
    end
    for j=1:nsamples
        polyatree(i)= polyatree(i) + (LOR(i,j)<=thres);
    end

end

polyatree = polyatree/nsamples;
kolmo_smir = kolmo_smir/nsamples;
wilcoxon = wilcoxon/nsamples;

%% Plots
if doplot
    % Power
    h = figure('name',['Power ' type]);
    plot(theta_trial, [polyatree, kolmo_smir, wilcoxon])
    xlabel('\theta')
    ylabel('Power')
    ylim([0, 1]) 
    legend('PT', 'KS', 'Wilcoxon')

    % Posterior probability of H1
    post = 1 - 1./(1+exp(-LOR));
    h2 = figure('name',['Posterior proba ' type]);
    plot(theta_trial, mean(post,2), 'b')
    quant = quantile(post, [.05, .95], 2);
    hold on
    plot(theta_trial, quant, 'r--')
    xlabel('\theta')
    ylabel('Pr(H_1 | Y)')
    ylim([0, 1]) 
else
    h=[];h2=[];
end

close(hbar)