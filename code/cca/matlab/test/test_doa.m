% test_doa: a script to examine the effectiveness of accurate order
% selection for direction of arrival (DOA) estimation in a reduced-rank CCA
% framework.
for n_sources = 2
    wavelength = 1; % normalized wavelength
    d = wavelength / 2;
    n_array = 16;
    n_sensor = 2*n_array;
    doas = linspace(-pi/3, pi/4, n_sources);
    %doas(1) = deg2rad(30);
    %doas(2) = deg2rad(35);
    snapshot_count = 200;
    power_source = [1 0.7;0.7 1];
    power_noise = 0.1;
    doa_precision = 0.01;
    iter = 1;
    n_aic_cca = zeros(iter,1);
    n_mdl_cca = zeros(iter,1);
    n_cp_cca = zeros(iter,1);
    n_prop_cca = zeros(iter,1);
    ic_mean = 0;
    figure(1)
    for i = 1:iter
        %%
        positions_x = linspace(0,4,n_array);
        positions_y = 10+positions_x;
        positions = [positions_x,positions_y];
        design_nominal = design_array_1d('custom', positions, d);
        design_perturbed = design_nominal;
        design_perturbed.gain_errors = 0.5 + sqrt(0.01)*randn(n_sensor, 1);
        % add some phase errors
        design_perturbed.phase_errors = exp(1j*sqrt(0.1)*randn(n_sensor, 1));
        % add some position errors
        pos_err = randn(2, n_sensor, 1) * sqrt(0.5*d);
        pos_err(1,2:end) = pos_err(1,2:end) - pos_err(1,1);
        pos_err(2,2:end) = pos_err(2,2:end) - pos_err(2,1);
        pos_err(:,1) = 0;
        design_perturbed.position_errors = pos_err;
        design_ula = design_perturbed;
        A = steering_matrix(design_ula, wavelength, doas);
        [m, k] = size(A);
        S_internal = gen_ccsg(k, snapshot_count, power_source);
        X = A * S_internal + gen_ccsg(m, snapshot_count, power_noise);
        Y = X(n_array+1:2*n_array,:);
        X = X(1:n_array,:);
        
        Rxx = (X*X')/snapshot_count;
        Ryy = (Y*Y')/snapshot_count;
        Rxy = (X*Y')/snapshot_count;
        
        C = (Rxx^(-0.5))*Rxy*(Ryy^(-0.5));
        
        ic = ic_cca(C,snapshot_count);
        [~,n_aic_cca(i)] = min(ic(:,1));
        [~,n_mdl_cca(i)] = min(ic(:,2));
        [~,n_cp_cca(i)] = min(ic(:,3));
        [~,n_prop_cca(i)] = min(ic(:,4));
        ic_mean = ic_mean + ic/iter;
        subplot(131)
        plot_spectrum(Rxx,Ryy,C,n_array,design_perturbed,...
                      wavelength,doa_precision,n_array,n_array,doas,'full-rank');
        subplot(132)
        plot_spectrum(Rxx,Ryy,C,n_aic_cca(i),design_perturbed,...
                      wavelength,doa_precision,n_array,n_array,doas,'reduced-rank using AIC');
        subplot(133)
        plot_spectrum(Rxx,Ryy,C,n_prop_cca(i),design_perturbed,...
                      wavelength,doa_precision,n_array,n_array,doas,'reduced-rank using Proposed');
    end
    disp(['AIC ',num2str(prob_correct(n_aic_cca,n_sources))])
    disp(['MDL ',num2str(prob_correct(n_mdl_cca,n_sources))])
    disp(['C_p ',num2str(prob_correct(n_cp_cca,n_sources))])
    disp(['Prop ',num2str(prob_correct(n_prop_cca,n_sources))])
    disp('----------------')
    % NOTE: Cp is a bit larger than the rest, so I'm normalizing them.
    figure(3)
    plot(bsxfun(@rdivide,ic_mean,1+max(ic_mean)))
    legend('AIC','MDL','C_p','Proposed')
end

%% Helper functions
function prob=prob_correct(n_estimated,n_true)
prob=(100*sum(abs(n_estimated-n_true)<1)/length(n_estimated));
end

function [J,sp] = plot_spectrum(Rxx,Ryy,C,rank_C,design,wavelength,...
                                precision,nx,ny,doas,figure_title)

[F,~,G] = svds(C,rank_C);

[doa_grid_rad, doa_grid_display, ~] = default_doa_grid(floor(1/precision), ...
                                                       'radian', 1);
A = steering_matrix(design, wavelength, doa_grid_rad);
Ax = A(1:nx,:);
Ay = A(1+nx:nx+ny,:);
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
if nargin > 9
    sp.true_positions = doas;
end
fprintf('Estimated DOAs:\n');
disp(sp.x_est);
if nargin > 10
    plot_sp(sp, 'title', figure_title, 'PlotType', 'Polar',...
            'ReuseFigure',true);
    ax=gca;             %Get Axis Handle
else
    plot_sp(sp, 'PlotType', 'Polar');
end

% NOTE: Here I'm setting the view-point for the 
%       3D plot. The azimuth and elevation values
%       were chosen based on my personal preference
%       from trial and error. 
az = 89;
el = 24;
view(az, el);
end

function draw_doa(doas)
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

function ic = ic_cca(C,M)
[~,K,~] = svd(C);
[m,n] = size(C);
lambda = diag(K);
ic = zeros(min(m,n),2);
for k = 1:min(m,n)
    d = k*(m + n - k);
    % AIC
    ic(k,1) = -M*sum(log(1 - lambda(k+1:min(m,n)).^2)) + 2*d;
    % MDL
    ic(k,2) = -M*sum(log(1 - lambda(k+1:min(m,n)).^2)) + log(M)*d;
    % Cp - Mallow's Cp with correction from Fujikoshi and Veitch 1979
    ic(k,3) = M*sum((lambda(k+1:min(m,n)).^2)./(1-lambda(k+1:min(m,n)).^2)) + ...
        -2*(m-k)*(n-k) + ...
        -2*(1/M)*(m-k)*(n-k)*(m+n+1+k-2*sum(lambda(1:k).^(-2)));
    % Proposed
    ic(k,4) = M*sum(lambda(k+1:min(m,n)).^2) - 2*(m-k)*(n-k)*log(log(M));
end
end

function X = gen_ccsg(m, n, cov)
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

