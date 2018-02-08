function varargout = parse_options(opt_param, param_names, param_default)

%
% PARSE_OPTIONS parses the set of optional parameters and returns their values
%   varargout = PARSE_OPTIONS(opt_param, param_names, param_default)
%
%   INPUTS
%   - opt_param:    cell of length 2*n {name_var1, value_var1, namevar2, value_var2, ...}
%   - param_names:  cell of length m containing the names of allowed optional parameters
%                   First value is the default value
%   - param_default:cell of length m of default values for optional parameters
%
%   OUTPUTS
%   - varargout:    Values of the optional parameters
%--------------------------------------------------------------------------
% EXAMPLE
% param_names = {'var1', 'var2', 'var3'};
% param_default = {1, true, 'b'};
% params = parsevar({'var1', 5, 'var3', 'c'}, param_names, param_default);
%--------------------------------------------------------------------------

%
% Reference: C.C. Holmes, F. Caron, J. Griffin and D.A. Stephens. 
% Two-sample Bayesian nonparametric hypothesis testing. To appear in
% Bayesian Analysis, 2014. (First version: arXiv:0910.5060, 2009).
%
% Copyright INRIA, 2014
% Author: François Caron, University of Oxford
% caron@stats.ox.ac.uk
% June 2008; Last update May 2014
%--------------------------------------------------------------------------

n_params = length(param_names);

if (nargout ~=n_params)
    error('The number of output arguments should be the same as the number of optional parameters')
end

% set default values
varargout = cell(n_params,1);
for i=1:n_params
    varargout{i} = param_default{i};
end

found = false(n_params, 1);
for i=1:2:length(opt_param)
    [tag, ind] = ismember(opt_param{i}, param_names);
  if ~ischar(opt_param{i}) || ~tag % Check if known optional parameter name
      warning(['Unknown optional parameter ' opt_param{i}]);      
  elseif found(ind)==true
      warning(['Value of the optional parameter ' opt_param{i} ' already specified'])
  else      
      varargout{ind} = opt_param{i+1};
      found(ind) = true;            
  end 
end
