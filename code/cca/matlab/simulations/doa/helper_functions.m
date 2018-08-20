% DOA experiment Helper functions
function h = helper_functions()
h.plot_spectrum = @plot_spectrum;
h.draw_doa = @draw_doa;
h.gen_ccsg = @gen_ccsg;
end

function [J,sp] = plot_spectrum(Rxx,Ryy,C,rank_C,design,wavelength,precision,...
                                doas,figure_title,xticks)
% Calculates and plots DOA objective function for different theta. 
% Based on objective number 3 in 
% [1] Wang, Ge, Kirsteins. "Direction-of-arrival estimation using distributed arrays:
% A canonical coordinates perspective with limited array size and sample support." 
% ICASSP, 2010. 
%
% NOTE:
%  This function uses default_doa_grid, steering_matrix,
%  find_doa_from_specturum_1d, and plot_sp from doa-tools. 
%  See ../../../../tools/ or https://github.com/morriswmz/doa-tools.git
%  
% NOTE: 
%  Make sure to cite M. Wang, Z. Zhang, and A. Nehorai's papers from 
%  https://research.wmz.ninja/research.html
%  
% Inputs:
%       Rxx,Ryy:        sample covariance matrices for two arrays
%       C:              coherence matrix (Rxx^-0.5)*Rxy*(Ryy^-0.5)
%       rank_C:         estimated rank. i.e, number of non-zero canonical correlations
%       design:         doa-tools structure, containing positions/parameters of
%                       arrays. 
%       wavelength:     used to calculate steering vectors. 
%       precision:      precision of DOA objective function in radians. 
%       doas:           true direction of arrivals in radians. 
%       figure_title:   title to use in 3D plot
%       xticks:         included ticks for polar axis.
%       
%
% Outputs:
%       J:              DOA objective function (peaks correspond to DOA)
%       sp:             spectrum structure, as defined in doa-tools. 
[F,K,G] = svds(C,rank_C);
[doa_grid_rad, doa_grid_display, ~] = default_doa_grid(floor(1/precision), ...
                                                       'radian', 1);
A = steering_matrix(design, wavelength, doa_grid_rad);
nx = size(C,1);
ny = size(C,2);
Ax = A(1:nx,:);
Ay = A(1+nx:nx+ny,:);

% NOTE: J is calculated using objective 3 as described in [1]
numerator = abs(Ax'*(Rxx'^(-1/2))*F*G'*(Ryy^(-1/2))*Ay);
denominator = abs(Ax'*(Rxx'^(-1))*Ax*Ay'*(Ryy^(-1))*Ay);
J = diag(numerator)./diag(denominator);
[x_est, ~, resolved] = find_doa_from_spectrum_1d(doa_grid_display, ...
                               J, rank_C);

sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est;
sp.x_unit = 'radian';
sp.y = J(:)';
sp.resolved = resolved;
sp.discrete = false;
sp.true_positions = doas;

plot_sp(sp, 'title', figure_title, 'PlotType', 'Polar',...
            'ReuseFigure',true,'xticks',xticks);

% NOTE: Here I'm setting the view-point for the 
%       3D plot. The azimuth and elevation values
%       were chosen based on my personal preference
%       from trial and error. 
az = 80;
el = 26;
view(az, el);
end

function draw_doa(doas)
% plot lines in a figure corresponding to DOAs in x-y plane. 
% 
% Inputs:
%   doas: direction of arrival angles in radians. 
%
% Outputs:
%   none. 

hold on
x = 0;
y = 0;
L = 30; % Assume source is at distance L
for alpha = doas(:)'
    x2=x+(L*cos(alpha));
    y2=y+(L*sin(alpha));
    plot(x2,y2,'kx')
end
end


function X = gen_ccsg(m, n, cov)
% generate random m random source signals of dimension n. 
% The output of this is used as the latent (hidden) source signals
% received by the sensor arrays. So, 
% Inputs:
%   m:      corresponds to the number of sensors 
%   n:      corresponds to the number of sources. 
%   cov:    n x n covariance matrix determining the paper of each source
%           signal and the correlation between each source pair. 
% Outputs:
%   X:      source signals. If A is the steering matrix the arrays receive
%           A*X + noise; 
%
% NOTE: 
%  This function was adopted from doa-tools. 

X0 = randn(m, n) + 1j * randn(m, n);
if isscalar(cov)
    X = sqrt(cov/2) * X0;
elseif isvector(cov)
    X = bsxfun(@times, X0, cov(:));
else
    C = sqrtm(cov/2);
    X = C*X0;
end
end
