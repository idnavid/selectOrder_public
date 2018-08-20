% testing MVDR using reduced rank CCA
clear
% close all
n = 1024; % samples
L = 100; % window length
h = 1; % hop length
k_fft = 1024; % number of frequency bins

x_true = 0;
y_true = 0;
n_comp = 3; % number of embeded sinusoids
% NOTE: figure out a better way to select frequencies. 
freqs = linspace(0.1,0.4,n_comp);
for i = 1:n_comp
    x_true = x_true + cos(2*pi*freqs(i)*(1:n)');
    phase_i = 2*pi*rand(1);
    y_true = y_true + cos(2*pi*freqs(i)*(1:n)' + phase_i);
end
S = 0;
S_reduced = 0; 
noise_var = 10;
p_orders = 2:1:4;
legends = {};
j = 0;
for p = p_orders
for iter = 1:500
x = x_true + noise_var*randn(n,1); 
y = y_true + noise_var*randn(n,1); 

[C,C_reduced]=coherence_matrix(x,y,n,L,h,p);

freq_bins = linspace(0,1,k_fft); 
S = S+mean(abs(fft(C,k_fft)),2);
S_reduced = S_reduced + mean(abs(fft(C_reduced,k_fft)),2);

end

S = (S/max(S));
S_reduced = (S_reduced/max(S_reduced));
hold on 
plot(freq_bins(1:k_fft/2),S_reduced(1:k_fft/2))
j = j+1;
legends{j} = num2str(p);
legend(legends);
pause(0.1)
end
plot(freq_bins(1:k_fft/2),S(1:k_fft/2))
legends{j+1} = 'full';
legend(legends);

for i = 1:n_comp
    text(freqs(i),max(S),num2str(freqs(i),'%.2f'))
end


function [C,C_reduced]=coherence_matrix(x,y,n,L,h,p)
% Calculate coherence matrix

X = buffer(x,L,L-h); 
Y = buffer(y,L,L-h); 

X = hamming_window(X);
Y = hamming_window(Y);

Rxx = (1/n)*(X*X');
Rxy = (1/n)*(X*Y');
Ryy = (1/n)*(Y*Y');
C = (Rxx^(-0.5))*Rxy*(Ryy^(-0.5));
[F,K,G] = svd(C);
K_reduced = 0*K;
K_reduced(1:p,1:p) = K(1:p,1:p); 
C_reduced = F*K_reduced*G';
end

function X = hamming_window(X)
[L,~] = size(X); 
X = bsxfun(@times,X,hamming(L));
end