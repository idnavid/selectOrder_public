% Experiment 3-1:
% Compare the detection accuracy across for different number of sources 
% using AIC, BIC, C_p, and Proposed. 
% This simulation uses arrays with randomly placed sensors. 
clear
addpath(genpath('../../../../tools/'))
addpath('../../../code')

order = estimate_order;
helper = helper_functions;
%%
iter = 5000;
sources = 1:12;
acc_aic_cca  = zeros(1,length(sources));
acc_bic_cca  = zeros(1,length(sources));
acc_cp_cca   = zeros(1,length(sources));
acc_prop_cca = zeros(1,length(sources));

n_aic_cca = zeros(iter,length(sources));
n_bic_cca = zeros(iter,length(sources));
n_cp_cca = zeros(iter,length(sources));
n_prop_cca = zeros(iter,length(sources));

for j = 1:length(sources)
    n_sources = sources(j);
    wavelength = 1; % normalized wavelength
    d = wavelength / 2;
    n_array = 16;
    n_sensor = 2*n_array;
    array_distance = 10;
    snapshot_count = 600;
    snr_level = 3;
    power_noise = 10^(-snr_level/10);
    power_source = eye(n_sources);
    doas = linspace(-pi/3, pi/4, n_sources);
    doa_precision = 0.01;
    parfor i = 1:iter
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
        
        [~,~,K,~,~,~,~] = order.coherence_matrix(X,Y);
        n_aic_cca(i,j)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'aic');
        n_bic_cca(i,j)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'bic');
        n_cp_cca(i,j)   = order.cca_order(diag(K),n_array,n_array,snapshot_count,'cp');
        n_prop_cca(i,j) = order.cca_order(diag(K),n_array,n_array,snapshot_count,'prop');
    end
    acc_aic_cca(j)=prob_correct(n_aic_cca(:,j),n_sources);
    acc_bic_cca(j)=prob_correct(n_bic_cca(:,j),n_sources);
    acc_cp_cca(j)=prob_correct(n_cp_cca(:,j),n_sources);
    acc_prop_cca(j)=prob_correct(n_prop_cca(:,j),n_sources);    
end
%% 
figure
hold on 
results = zeros(4,length(sources));
results(1,:) = acc_aic_cca;
results(2,:) = acc_bic_cca;
results(3,:) = acc_cp_cca;
results(4,:) = acc_prop_cca;
semilogx(sources,results(1,:),'r:','Marker','o','MarkerFaceColor','r')
semilogx(sources,results(2,:),'m:','Marker','s','MarkerFaceColor','m')
semilogx(sources,results(3,:),'b:','Marker','^','MarkerFaceColor','b')
semilogx(sources,results(4,:),'k:','Marker','*','MarkerFaceColor','k')
ylim([0,101])
xlim([min(sources),max(sources)])
legend({'AIC','BIC','C_p','Proposed'})




