## Blue-print for deriving AIC from KL divergence

minimize I(f,g(.|\theta)) = minimize \int f(x)\log(f(x)/g(x|\theta))dx

We don'te have I(f,g(.|\theta)), but I(f,g(.|\hat{\theta}_{ML})) is an upper bound. 

==> So, minimize I(f,g(.|\hat{\theta}_{ML})) = minimize E_t [E_x [\log(g(x|hat{\theta}_{ML}))] ]

Step 1: Taylor series expansion of \log(g(x|hat{\theta}_{ML})) around \theta_0

Step 2 : Taylor series expansion of \log(g(x|{\theta}_{0})) around \hat{\theta}_{ML}

Conclusion: if f \in {g}, then trace(I\Sigma) = k
