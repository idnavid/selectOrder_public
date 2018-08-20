% Order estimation function for CCA 
function h = estimate_order()
% collection of function for order estimation
h.coherence_matrix = @coherence_matrix;
h.cca_order = @cca_order;
end
function [C,F,K,G,Rxx,Rxy,Ryy] = coherence_matrix(X,Y)
% Calculate sample correlation coefficients/coords using svd of coherence 
% matrix. Works for Real and Complex matrices. 
%
% Inputs: 
%       X,Y:  data matrices of size (dim x N), where N is number of
%             samples.
%
% Outputs:
%       C:              coherence matrix
%       F,G:            left and right eigenvectors of coherence matrix, 
%                       aka canonical coordinates of X,Y (respectively).
%       K:              diagonal matrix containing canonical correlatoin
%                       coefficients. 
%       Rxx,Rxy,Ryy:    sample covariance of X, cross-covariance of X and Y, 
%                       covariance of Y, respectively. 

assert(size(X,2)==size(Y,2),'Number of columns must be the same');
assert(size(X,1)+size(Y,1)<size(X,2),'Not enough sampls!');
N = size(X,2);
X = bsxfun(@minus,X,mean(X,2));
Y = bsxfun(@minus,Y,mean(Y,2));

Rxx = (X*X')/N;
Ryy = (Y*Y')/N;
Rxy = (X*Y')/N;
C = (Rxx^(-0.5))*Rxy*(Ryy^(-0.5));

[F,K,G] = svd(C);
end

function k_hat = cca_order(canon_corr,m,n,M,mode,plot_ic)
% Use information criteria to estimate the 
% optimum number of non-zero canonical correlations. 
% 
% Inputs:
%       canon_corr: list of canonical correlations in decreasing order. 
%       m:          dimension of first dataset. 
%       n:          dimension of second dataset. 
%                   Note that length(can_corr)=min(m,n)
%       M:          total number of samples used to estimation can_corr. 
%       mode:       formula used for estimation.
%                   options: aic, biv, cp, prop (default: bic)
%       plot_ic:  boolean value to plot information criterion. (default: false)
%                 this is for debugging purposes.
%
% Outputs:
%       k_hat: number of non-zero canonical correlations. 
%              1 =< k_hat =< length(can_corr)
%
% My primary reference used to for the information criteria is 
% Gunderson and Muirhead, "On estimating the dimensionality in canonical 
% correlation analysis." J. of Multiv. Analysis, vol 62,pp.121-136, 1997

assert(numel(canon_corr)==length(canon_corr)); 
assert(numel(canon_corr)==min(m,n)); 
canon_corr = sort(canon_corr,'descend');

if nargin<5
    mode = 'bic'; % default mode
end

if nargin<6
    plot_ic = false; % default plot option
end

information_criterion = zeros(length(canon_corr),1);
for k = 1:min(m,n)
    switch mode
        case 'aic'
            information_criterion(k,1) = cca_aic(m,n,M,k,canon_corr);
        case 'bic'
            information_criterion(k,1) = cca_bic(m,n,M,k,canon_corr);
        case 'cp'
            information_criterion(k,1) = cca_mcp(m,n,M,k,canon_corr);
        case 'prop'
            information_criterion(k,1) = cca_pro(m,n,M,k,canon_corr);
    end
end

if plot_ic
    figure
    plot(information_criterion)
end

[~,k_hat] = min(information_criterion);
end


function aic = cca_aic(m,n,M,k,canon_corr)
% calculate Akaike's information criterion for candidate order k. 
% This formula is due to Gunderson and Muirhead, 1997. 
% 
% Inputs:
%       m,n,M,canon_corr: as defined above in cca_order
%       k:                candidate order.
% Outputs: 
%    aic:                 scalar value of AIC(k). 


deg_of_freedom = (m - k)*(n - k); 
aic = -M*sum(log(1 - canon_corr(k+1:min(m,n)).^2)) - 2*deg_of_freedom;
end


function bic = cca_bic(m,n,M,k,canon_corr)
% calculate Bayesian information criterion for candidate order k. 
% 
% Inputs:
%       m,n,M,canon_corr: as defined above in cca_order
%       k:                candidate order.
% Outputs: 
%    bic:                 scalar value of BIC(k). 
%
% NOTE:
%  There's a slight ambiguity with the choice of name here. According to
%  Gunderson and Muirhead, 1997, which is my main reference for this code,
%  the following criterion is based on Schwartz's 1978 paper, famously
%  known as the Bayesian information criterion. However, I was originally
%  under the impression that this was the MDL, due to Zhang and Wong 1993. 

deg_of_freedom = (m - k)*(n - k); 
bic = -M*sum(log(1 - canon_corr(k+1:min(m,n)).^2)) - log(M-1)*deg_of_freedom;
end

function cp = cca_mcp(m,n,M,k,canon_corr)
% calculate Mallow's Cp for candidate order k. 
% 
% Inputs:
%       m,n,M,canon_corr: as defined above in cca_order
%       k:                candidate order.
% Outputs: 
%    cp:                 scalar value of Cp(k). 

deg_of_freedom = (m - k)*(n - k); 
cp = M*sum((canon_corr(k+1:min(m,n)).^2)./(1-canon_corr(k+1:min(m,n)).^2)) + ...
      -2*deg_of_freedom;
if false
    % NOTE: 
    %  This is an adjusted version of Mallow's C_p. I need to remember
    %  where I found this reference. 
    cp = M*sum((canon_corr(k+1:min(m,n)).^2)./(1-canon_corr(k+1:min(m,n)).^2)) + ...
        -2*(m-k)*(n-k)-2*(1/M)*(m-k)*(n-k)*(m+n+1+k-2*sum(canon_corr(1:k).^(-2)));
end
end

function pic = cca_pro(m,n,M,k,canon_corr)
% calculate Proposed criterion candidate order k. 
%
% Inputs:
%       m,n,M,canon_corr: as defined above in cca_order
%       k:                candidate order.
% Outputs: 
%    pic:                 scalar value. stands for propsed IC. 
pic = (M/2)*sum(canon_corr(k+1:min(m,n)).^2) - (m-k)*(n-k)*log(log(M));
end


