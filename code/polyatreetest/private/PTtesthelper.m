function [LOR, LOR_level, LOR_num, LOR_den] = PTtesthelper(y1, y2, partition_PD, epsilon, dec_epsilon, const, LOR_level, printInfo)

%
% function  PTtesthelper.m
% Recursive function
%
% INPUT
% y1, y2            two samples
% PD_centering      Probability distribution defining the centering
%                   partition (Default= normal(0,1))
% epsilon           level of the tree (binary)
% dec_epsilon       level of the tree (decimal)
% const             constant c for the coefficient of the Polya tree
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

% dec_epsilon=uint64(dec_epsilon);
if numel(const) ==2
    const_num = const(1);
    const_den = const(2);
else
    const_num = const;
    const_den = const;
end
    
    
level_max = floor(log2(length([y1;y2])));

if nargin<13
    printInfo=0;
end


if isempty(epsilon)
    LOR_level = zeros(level_max+1,1);
    LOR_all = zeros(level_max+1, 3);
end
% epsilon
epsilon0 = [epsilon '0'];
epsilon1 = [epsilon '1'];
dec_epsilon0 = dec_epsilon*2;
dec_epsilon1 = (dec_epsilon0+1);
inter0 = getInter(partition_PD, epsilon0, dec_epsilon0);
inter1 = getInter(partition_PD, epsilon1, dec_epsilon1);

n0_1 = num_y(y1,inter0);
n1_1 = num_y(y1,inter1);
n0_2 = num_y(y2,inter0);
n1_2 = num_y(y2,inter1);
n0_12 = n0_1+n0_2;
n1_12 = n1_1+n1_2;

alpha0_1 = coeff(epsilon0, const_den);
alpha1_1 = coeff(epsilon1, const_den);
alpha0_2 = coeff(epsilon0, const_den);
alpha1_2 = coeff(epsilon1, const_den);
alpha0_12 = coeff(epsilon0, const_num);
alpha1_12 = coeff(epsilon1, const_num);

if ((n0_1+n1_1)==0 || (n0_2+n1_2)==0 || length(epsilon)>=level_max) %((n0_12+n1_12)<=1) || length(epsilon)>=100 % safeguard for too many recursions
    % End of the tree
    contrib=0; 
    contrib_num = 0;
    contrib_den = 0;
    LOR=contrib;
    LOR_num = contrib_num;
    LOR_den = contrib_den;
    LOR_level(length(epsilon)+1) = contrib + LOR_level(length(epsilon)+1);
else
    contrib_num = - betaln(alpha0_12,alpha1_12)...
        +betaln(alpha0_12 + n0_12, alpha1_12 + n1_12);
    contrib_den = - betaln(alpha0_1,alpha1_1) - betaln(alpha0_2,alpha1_2)...
        + betaln(alpha0_1 + n0_1, alpha1_1 + n1_1)...
        + betaln(alpha0_2 + n0_2, alpha1_2 + n1_2);
    contrib = contrib_num - contrib_den;

    LOR_level(length(epsilon)+1) = contrib + LOR_level(length(epsilon)+1);    
    
    [contrib_left, LOR_level, contrib_left_num, contrib_left_den]= PTtesthelper(y1, y2, partition_PD, epsilon0, dec_epsilon0, const, LOR_level);
    [contrib_right, LOR_level, contrib_right_num, contrib_right_den]= PTtesthelper(y1, y2, partition_PD,epsilon1, dec_epsilon1, const, LOR_level);
    LOR = contrib + contrib_left + contrib_right ;
    LOR_num = contrib_num + contrib_left_num + contrib_right_num;
    LOR_den = contrib_den + contrib_left_den + contrib_right_den;
    if printInfo
        fprintf('epsilon=%s, contrib=%.5f\n',epsilon,contrib)
        fprintf('%d, %d, %d, %d, %d, %d\n',n0_1,n1_1,n0_2,n1_2,n0_12,n1_12)
    end
end
end



%% Subfunctions

function alpha=coeff(epsilon, const)

% Return coefficient alpha associated to interval inter and level epsilon
% under some prior G0

m=length(epsilon);
power_m = 2;
alpha = const*m^power_m;

% %     alpha = const; % ATTENTION MODIF!!!
% elseif inter(1)==inter(2)
%     alpha=0;
% else
%     log_alpha=log(const)+power_m*log(m);
%     log_alpha=log_alpha+log(G0(inter(2))-G0(inter(1)));
%     % alpha=const*m^2*2^m*(G0(inter(2))-G0(inter(1)));
%     alpha=exp(log_alpha);
% end
end
%%
function interval=getInter(PD, epsilon, epsilon_dec)

% Get the interval from the cdf (x, f), the level (epsilon, dec_epsilon)

m=length(epsilon);

if nargin<4
    epsilon_dec=bin2dec(epsilon);
end

f0_min=(epsilon_dec)/2^m;
f0_max=(epsilon_dec+1)/2^m;

if isempty(PD) ||  ~isa(PD,'ProbDistUnivParam')
    % partition based on gaussian 0,1 
    interval = norminv([f0_min, f0_max]);
    % partition based on student 0,1,4
%     interval = tinv([f0_min, f0_max], 4);    
else
    interval = icdf(PD, [f0_min, f0_max]);
%     
%     ind_min=find((f-f0_min)>=0,1);
%     ind_max=find((f-f0_max)>=0,1);
%     interval=[x(ind_min), x(ind_max)];

%     % partition based on gaussian 0,1
%     interval = norminv([f0_min, f0_max]);
    
end

if epsilon_dec==0
    interval(1)=-Inf;
elseif epsilon_dec==2^m-1
    interval(2)=Inf;
end
end

%%
function n=num_y(y,inter)

% Number of data in interval inter

n=sum(y<=inter(2) & y>inter(1));
end
