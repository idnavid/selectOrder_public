% Experiment 2-1:
% Compare the detection accuracy across different SNR values for four
% sources using AIC, BIC, C_p, and Proposed. 
% This simulation uses arrays with randomly placed sensors. 

addpath(genpath('../../../../tools/'))
addpath('../../../code')

order = estimate_order;
helper = helper_functions;
%% 
snr_levels = 10:-2:-10;
noise_levels = 10.^(-snr_levels/10);
n_noises = length(noise_levels);
n_aic_cca = zeros(iter,n_noises);
n_bic_cca = zeros(iter,n_noises);
n_cp_cca = zeros(iter,n_noises);
n_prop_cca = zeros(iter,n_noises);

for j = 1:n_noises
    power_noise = noise_levels(j);
    wavelength = 1; % normalized wavelength
    d = wavelength / 2;
    n_array = 15;
    n_sensor = 2*n_array;
    array_distance = 10;
    n_sources = 8;
    snapshot_count = 300;
    power_source = eye(n_sources);
    doas = linspace(-pi/3, pi/4, n_sources); 
    iter = 1000;
    for i = 1:iter
        positions_x = linspace(0,4,n_array);
        positions_y = array_distance+positions_x;
        positions = [positions_x,positions_y];
        design_nominal = design_array_1d('custom', positions, d);
        design_perturbed = design_nominal;
        % Adding gain, phase, and position errors:
        design_perturbed.gain_errors = 0.5 + sqrt(0.01)*randn(n_sensor, 1);
        design_perturbed.phase_errors = exp(1j*sqrt(0.1)*randn(n_sensor, 1));
        pos_err = randn(2, n_sensor, 1) * sqrt(0.5*d);
        pos_err(1,2:end) = pos_err(1,2:end) - pos_err(1,1);
        pos_err(2,2:end) = pos_err(2,2:end) - pos_err(2,1);
        pos_err(:,1) = 0;
        
        design_perturbed.position_errors = pos_err;
        design_ula = design_perturbed;
        A = steering_matrix(design_ula, wavelength, doas);
        
        [m, k] = size(A);
        S_internal = helper.gen_ccsg(k, snapshot_count, power_source);
        X = A * S_internal + helper.gen_ccsg(m, snapshot_count, power_noise);
        Y = X(n_array+1:2*n_array,:);
        X = X(1:n_array,:);
        
        [C,F,K,G,Rxx,Rxy,Ryy] = order.coherence_matrix(X,Y);
        n_aic_cca(i,j)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'aic');
        n_bic_cca(i,j)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'bic');
        n_cp_cca(i,j)   = order.cca_order(diag(K),n_array,n_array,snapshot_count,'cp');
        n_prop_cca(i,j) = order.cca_order(diag(K),n_array,n_array,snapshot_count,'prop');
    end
end
%% 
figure
hold on 
results = zeros(4,length(snr_levels));
results(1,:) = prob_correct(n_aic_cca,n_sources);
results(2,:) = prob_correct(n_bic_cca,n_sources);
results(3,:) = prob_correct(n_cp_cca,n_sources);
results(4,:) = prob_correct(n_prop_cca,n_sources);
semilogx(snr_levels,results(1,:),'k:','Marker','o')
semilogx(snr_levels,results(2,:),'k:','Marker','s')
semilogx(snr_levels,results(3,:),'k:','Marker','*')
semilogx(snr_levels,results(4,:),'k','Linewidth',2)
ylim([0,101])
legend({'AIC','BIC','C_p','Proposed'})
set(gca, 'XDir','reverse')

