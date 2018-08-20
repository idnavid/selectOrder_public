% Small matlab version algorithm. 
% TODO:
%  - Implementation of other criteria (MDL, AIC, BIC, Laplace)

function h = estimate_dim()
% returns function handle h pointing to sum_squared function. 
    h.sum_squared = @sum_squared;
end

function [obj,dim] = sum_squared(V,N,mode)
% Estimate data dimension using sum-of-squared eigenvalues
% Based on method proposed in
% Seghouane, Shokouhi, "Consistent Esitmaiton of Dimensionality for Data 
% Driven Methods in fMRI Analysis", IEEE Transactions on Medical Imaging, 2018.
%
% Input: 
%       V:      squared eigenvalues
%       N:      number of samples
%       mode:   which penalty to use (default log)
%               options are: 2, log(N), log(log(N))
% Outputs:
%       dim:    estimated dimension 
% 
if nargin<3
    mode='log';
end
    
n = length(V);
obj = zeros(3,n);
for k = 1:n
    V_hat = abs(V - mean(V(k:end)));
    obj(1,k) = (N-1.)*sum(V_hat(k:end)) - 2*(n-k)*(n-k+1.)/2.;
    obj(2,k) = (N-1.)*sum(V_hat(k:end)) - log(N)*(n-k)*(n-k+1.)/2.;
    obj(3,k) = (N-1.)*sum(V_hat(k:end)) - log(log(N))*(n-k)*(n-k+1.)/2.;
end
switch mode
    case '2'
        obj = obj(1,:);
        [~,dim] = min(obj);
    case 'log'
        obj = obj(2,:);
        [~,dim] = min(obj);
    case 'loglog'
        obj = obj(1,:);
        [~,dim] = min(obj);
end
end

