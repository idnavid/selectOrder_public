clear
addpath('../')

h = estimate_dim();
helpers = helper_functions();
%%
q = 5; % number of sources
p = 50; % dimension
n = 1000; % number of samples
sigma2 = 0.5;
source_signal = randl(q,n);
s = source_signal;
A = randn(p,q);
information_criteria=0;
n_iter = 1;
for i = 1:n_iter
X = A*s + sqrt(sigma2)*randn(p,n);
[~,V,~] = svd((1/n)*(X*X'));
V = diag(V);
sigma2_hat = helpers.estimate_sigma(V,5);
[~,V,~] = svd(X*X'/n - sigma2_hat*eye(p));
[information_criteria_i,~]= h.sum_squared(diag(V),n);
information_criteria =  information_criteria + information_criteria_i/n_iter;
end
[~,q_hat]= min(information_criteria);
plot(information_criteria)
hold on 
plot(q_hat,information_criteria((q_hat)),'r*')
text(q_hat,information_criteria((q_hat))*-10,...
     ['estimated: ',num2str(q_hat);
      'true     : ',num2str(q)])
legend('Proposed IC')