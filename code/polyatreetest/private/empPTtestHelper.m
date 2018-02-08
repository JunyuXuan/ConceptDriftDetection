function [LOR, LOR_level, LORnum, LORden] = empPTtestHelper(vect_a, n, m, dec_epsilon, LOR_level, printInfo, c)

%
% function  empPTtesthelper.m
% Recursive function
%
% INPUT
% vect_a            Binary vector of length the sample size. provides
%                   positions of data from the first sample when samples are sorted
% n                 Number of samples
% m                 Level of the tree
% dec_epsilon       level of the tree (decimal)
% c                 constant c for the coefficient of the Polya tree
% LOR_level         log-ratio at each level of the tree (only for debugging purposes)
% printInfo         tag
%
% OUTPUT
% LOR               Sum of log ratios
% LOR_level         Sum of log ratios at each level of the tree
%--------------------------------------------------------------------------

% Reference: C.C. Holmes, F. Caron, J. Griffin and D.A. Stephens. 
% Two-sample Bayesian nonparametric hypothesis testing. To appear in
% Bayesian Analysis, 2014. (First version: arXiv:0910.5060, 2009).
%
% Copyright INRIA, 2014
% Author: François Caron, University of Oxford
% caron@stats.ox.ac.uk
% June 2008; Last update May 2014
%--------------------------------------------------------------------------



if nargin<6
    printInfo=0;
end
if nargin<7
    c = 0.1;
end

if m==0
    LOR_level=zeros(200,1);
end

m = m +1; % level of the tree
alpha = c * m^2;

dec_epsilon0 = dec_epsilon*2;
dec_epsilon1 = (dec_epsilon0+1);

ind0 = floor(dec_epsilon0/2^m*n)+1:floor(dec_epsilon1/2^m*n);
ind1 = floor(dec_epsilon1/2^m*n)+1:floor((dec_epsilon1+1)/2^m*n);

n0_12 = length(ind0);
n1_12 = length(ind1);
n0_1 = sum(vect_a(ind0));
n1_1 = sum(vect_a(ind1));
n0_2 = n0_12 - n0_1;
n1_2 = n1_12 - n1_1;

% [n0_12, n1_12, n0_1, n1_1, n0_2, n1_2] 
   
% contrib = log(hygepdf(n0_1, n0_12+n1_12, n0_1+n1_1, n0_12)) ...
%     - log (walleniusmargpdf(n0_1, n0_12+n1_12, n0_1+n1_1, n0_12, alpha, alpha)); 
contrib_num = log(hygepdf(n0_1, n0_12+n1_12, n0_1+n1_1, n0_12));
contrib_den = log (EHGmargpdf(n0_1, n0_12+n1_12, n0_1+n1_1, n0_12, alpha, alpha));

contrib = contrib_num - contrib_den;  

if (n0_1>1 && n0_2>1)
    [contrib_left, LOR_level, contrib_left_num, contrib_left_den] = empPTtestHelper(vect_a, n, m, dec_epsilon0, LOR_level, printInfo, c);
else
    contrib_left = 0;
    contrib_left_num = 0;
    contrib_left_den = 0;
end
if (n1_1>1 && n1_2>1)
    [contrib_right, LOR_level, contrib_right_num, contrib_right_den] = empPTtestHelper(vect_a, n, m, dec_epsilon1, LOR_level, printInfo, c);
else
    contrib_right = 0;
    contrib_right_num = 0;
    contrib_right_den = 0;
end

LOR = contrib + contrib_left + contrib_right ;
LOR_level(m) = contrib + LOR_level(m) ;
LORnum = contrib_num + contrib_left_num + contrib_right_num;
LORden = contrib_den + contrib_left_den + contrib_right_den;
if printInfo
    fprintf('level=%d, contrib=%.1f %.1f %.1f\n', m, contrib)
    fprintf('%d, %d, %d, %d, %d, %d\n',n0_1,n1_1,n0_2,n1_2,n0_12,n1_12)
end
end