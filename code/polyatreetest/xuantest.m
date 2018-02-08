mfiles_path = './'; % Indicate the path to the folder
addpath(mfiles_path); % Add to the Matlab/Octave path

% nsamples = 200; % Number of samples - Note: 1000 were used to create the figures in the paper
% 
% n = 50;
% %types = {'mean', 'var', 'mixture', 'tail', 'lognorm_mean', 'lognorm_var'};
% types = {'var'};
% for k=1:length(types)
%     type = types{k};
%     [theta_trial, polyatree, kolmo_smir, wilcoxon, LOR, h, h2] = powercurve(n, type, ...
%     'nsamples', nsamples, 'doplot', true);
%     figure(h)
%     title(type);
%     figure(h2)
%     title(type);
% end

y1 = randn(50, 1);
y2 = randn(50, 1);
[h, post, stats] = PTtest(y1, y2, 'partition', 'empirical')


