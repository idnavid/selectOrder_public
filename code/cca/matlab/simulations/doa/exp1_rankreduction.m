% Exp1: 
% calculate the DOA for 100 runs using three different approaches:
% - without rank-reduction (full-rank)
% - with rank-reduction using AIC
% - with rank-reduction using proposed

% I calculate order for bic and C_p as well, but don't use them here. 

addpath(genpath('../../../../../../ext_tools/'))
addpath('../../../../cca/matlab/')

order = estimate_order;
helper = helper_functions;
%% 
for n_sources = 3
    wavelength = 1; % normalized wavelength
    d = wavelength / 2;
    n_array = 16;
    n_sensor = 2*n_array;
    array_distance = 10;
    doas = linspace(-pi/3, pi/4, n_sources);
    snapshot_count = 200;
    power_source = eye(n_sources);%[1 0.5;0.5 1];
    power_noise = 0.9;
    doa_precision = 0.01;
    iter = 10;
    n_aic_cca = zeros(iter,1);
    n_bic_cca = zeros(iter,1);
    n_cp_cca = zeros(iter,1);
    n_prop_cca = zeros(iter,1);
    figure(1)
    % NOTE: 
    %  the boolean xtick tells helper_functions.plot_spectrum whether it
    %  should print xaxis labels for the figure. This is by default false
    %  when reusing the same figure, except for the last plot. Otherwise
    %  the plot_spectrum keeps overlaying the same ticks on top of
    %  eachother and it ends up looking smudged. 
    xtick = false;
    for i = 1:iter
        % NOTE: 
        %  In general the sensors in each array are clusterred together. 
        %  Here we start with a linear array (*_nominal). 
        %  Then some random purturbation is added. 
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
        n_aic_cca(i)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'aic');
        n_bic_cca(i)  = order.cca_order(diag(K),n_array,n_array,snapshot_count,'bic');
        n_cp_cca(i)   = order.cca_order(diag(K),n_array,n_array,snapshot_count,'cp');
        n_prop_cca(i) = order.cca_order(diag(K),n_array,n_array,snapshot_count,'prop');
        if i == iter
            xtick = true;
        end
        subplot(311)
        helper.plot_spectrum(Rxx,Ryy,C,n_array,design_perturbed,...
                             wavelength,doa_precision,doas,'full-rank',xtick);
        subplot(312)
        helper.plot_spectrum(Rxx,Ryy,C,n_aic_cca(i),design_perturbed,...
                             wavelength,doa_precision,doas,'reduced-rank using AIC',xtick);
        subplot(313)
        helper.plot_spectrum(Rxx,Ryy,C,n_prop_cca(i),design_perturbed,...
                             wavelength,doa_precision,doas,'reduced-rank using Proposed method',xtick);
    end
    disp(['AIC ',num2str(prob_correct(n_aic_cca,n_sources))])
    disp(['BIC ',num2str(prob_correct(n_bic_cca,n_sources))])
    disp(['C_p ',num2str(prob_correct(n_cp_cca,n_sources))])
    disp(['Prop ',num2str(prob_correct(n_prop_cca,n_sources))])
    disp('----------------')
    print('-bestfit','figures/exp1_results','-dpdf')
end


