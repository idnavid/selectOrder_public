function h = helper_functions()
% returns function handle h pointing to functions in this m-file. 
    h.estimate_sigma = @estimate_sigma;
end

function sigma_hat2 = estimate_sigma(V,M)
% estimate noise variance from eigenvalues of 
% Sigma + \sigma^2I.
% 
% Inputs:
%       V: eigenvalues of sample covariance matrix
%       M: number of bottom eigenvalues to be used. (default: 10)
%           NOTE: If data is too noisy, the method can be improved by
%           increasing M. 
% Outputs:
%       sigma_hat: estimated sigma^2. 
if nargin<2
    M = 10;
end
sigma_hat2 = sum(V(end-M:end).^2)/(sum(V(end-M:end)));
end