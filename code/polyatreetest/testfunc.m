% runs the functions PTtest and powercurve with different parameters

close all
clear all

n = 50;
y1 = randn(n, 1);
y2 = randn(n, 1);

% Subjective Polya tree test with empirical estimation of c
[h, post, stats] = PTtest(y1, y2, 'estimate_c', true);

% Conditional Polya tree test with empirical estimation of c
[h, post, stats] = PTtest(y1, y2, 'estimate_c', true, 'partition', 'empirical');


% Power curves for a variance shift
type = 'var';
nsamples = 20;

% Subjective Polya Tree test
[theta_trial, polyatree, kolmo_smir, wilcoxon, LOR, h, h2] = powercurve(n, type, ...
    'nsamples', nsamples, 'doplot', true);

% Conditional Polya tree test
[theta_trial, polyatree, kolmo_smir, wilcoxon, LOR, h, h2] = powercurve(n, type, ...
    'nsamples', nsamples, 'partition', 'empirical', 'doplot', true);

