function [h, post, stats] = PTtest(y1, y2, varargin)

%
% PTTEST performs a two-sample test based on a Polya tree prior
%       y1 ~ F1 ; y2 ~F2
%       H0: F1 = F2 vs H1: F1 ~=F2
%   [h, post, stats] = PTTEST(y1, y2, 'PARAM1', val1, 'PARAM2', val2,...)
%
% Requires the statistics toolbox
%--------------------------------------------------------------------------
% INPUTS 
%   - y1:   Vector of length n1. First sample.
%   - y2:   Vector of length n2. Second sample.
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

param_names = {'c', 'threshold', 'normalize', 'estimate_c', 'partition'};
param_default = {1, .5, true, false, []};
[c, threshold, normalize, estimate_c, partition] = parse_options(varargin, param_names, param_default);

% Some checks on the inputs
if ~islogical(estimate_c)
    error('Option ''estimate'' should be of type logical');
end
if ~islogical(normalize)
    error('Option ''normalize'' should be of type logical');
end
if ~isnumeric(c) || length(c)>2 || min(c)<0
    error('Option c should be a vector of length 1 or 2 of strictly positive scalars');
end
if (~isnumeric(threshold) || threshold<=0 || threshold>=1) 
    error('Parameter threshold should be a scalar in (0,1)');
end
if ~isempty(partition) && ~isa(partition, 'ProbDist') && ~strcmp(partition, 'empirical')
    error('Parameter partition must be of Matlab class ProbDist or equal to ''empirical''')
end

y = [y1; y2];
if strcmp(partition, 'empirical')
    isempirical = true;
    n1 = length(y1);
    n2 = length(y2);
    vect_a = logical([ones(1,n1), zeros(1, n2)]);% Locations of data from sample 1
    [~, ind12] = sort(y);
    vect_a = vect_a(ind12); % Locations of data from sample 1 in ordered data
else
    isempirical = false;
    partition_PD = partition;
end


if normalize
    % Normalize data
    y1 = (y1-median(y))/iqr(y);
    y2 = (y2-median(y))/iqr(y);
    
%     y1 = y1-median(y);
%     y2 = y2-median(y);
end



if estimate_c
    % Estimate the value of the parameter c by empirical Bayes
    c_all = [.01, .1, 1, 10, 100, 1000];
    nb_c_all = length(c_all);

    LM_H0 = zeros(nb_c_all, 1);
    LM_H1 = zeros(nb_c_all, 1);
    for i=1:length(c_all)
        const = c_all(i);
        if isempirical
            [~, ~, ~, LM_H1(i)] = empPTtestHelper(vect_a, n1+n2, 0, 0, 0, 0, const);
        else
           [~, ~, LM_H0(i), LM_H1(i)] = PTtesthelper(y1, y2, partition_PD, '', 0, const); 
        end
    end
    if isempirical
        [~, ind] = max(LM_H1);
        c = c_all(ind);
    else
        [~, ind_num] = max(LM_H0);
        [~, ind_den] = max(LM_H1);
        c = c_all([ind_num, ind_den]);
    end
end

if isempirical
    % Conditional test
    [LOR, LOR_level, LM_H0, LM_H1] = empPTtestHelper(vect_a, n1+n2, 0, 0, 0, 0, c);
else
    % Subjective test
    [LOR, LOR_level, LM_H0, LM_H1] = PTtesthelper(y1, y2, partition_PD, '', 0, c);
end
post = 1/(1+exp(-LOR));
h = post < threshold;
stats.LOR = LOR;
stats.LOR_level = LOR_level;
stats.LOR_num = LM_H0;
stats.LOR_den = LM_H1;
stats.c = c;