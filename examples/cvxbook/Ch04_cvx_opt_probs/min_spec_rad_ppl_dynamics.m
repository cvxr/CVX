% Sec. 4.5.4 example: Minimum spectral radius via Perron-Frobenius theory
% Boyd & Vandenberghe "Convex Optimization"
% Joelle Skaf - 01/29/06
% 
% Check pp. 165 - 167 for the detailed explanation of the population
% dynamics model p(t+1) = Ap(t), the underlying assumptions and the model 
% of the birth and survival rates as monomial functions of the 
% concentrations of 2 chemicals c1 and c2 in the environments. 
% The goal is minimize, over c1 and c2, the spectral radius of A which is 
% the Perron-Frobenius eigenvalue lambda of A and which corresponds to finding the
% fastest decay rate, or slowest growth rate, for the population. 
% The problem can be posed as a geometric program 
%           minimize    lambda
%               s.t.    b1v1 + b2v2 + b3v3 + b4v4 <= lambda*v1
%                       s1v1 <= lambda*v2
%                       s2v2 <= lambda*v3
%                       s3v3 <= lambda*v4
%                       1/2 <= ci <= 2 ,                        i = 1,2
%                       bi = bi^{nom}(c1/c1^{nom})^alpha_i*...
%                               (c2/c2^{nom})^beta_i ,          i = 1,...,4 
%                       si = si^{nom}(c1/c1^{nom})^gamma_i*...
%                               (c2/c2^{nom})^delta_i ,         i = 1,...,3
% with variables bi, si, ci, vi, lambda

% input data
n = 4;
alphas = [-0.3 -0.4 -0.2 -0.1];
betas = [0.9 0.8 -0.2 0.1];
gammas = [0.1 0.15 0.2];
deltas = [-0.5 -0.2 0.1];
b_nom = [0.5 1.8 1.2 0.2];
s_nom = [0.5 0.3 0.2];
c_nom = [1 1];

% GP formulation
cvx_begin 
    variables x(n) y(n-1) z(2) t(n) u
    minimize (u)
    logsumexp_sdp(x+t)<= u+t(1);
    y(1:3) + t(1:3) <= u + t(2:4);
    z <= log(2);
    z >= log(1/2);
    for i=1:n
        x(i) == log(b_nom(i))+alphas(i)*(z(1)-log(c_nom(1)))+ ...
                        betas(i)*(z(2) - log(c_nom(2)));
    end
    for j=1:n-1
        y(j) == log(s_nom(j))+gammas(j)*(z(1)-log(c_nom(1)))+ ...
                        deltas(j)*(z(2) - log(c_nom(2)));
    end
cvx_end

% recovering original variables
b = exp(x); 
s = exp(y);
c = exp(z);
v = exp(t);
lambda = exp(u);

% displaying results
if lambda < 1 
    disp('The fastest decay rate of the bacteria population is: ');
    disp(lambda);
else 
    disp('The slowest growth rate of the bacteria population is: ');
    disp(lambda);
end
disp('The concentration of chemical 1 achieving this result is: ');
disp(c(1));
disp('The concentration of chemical 2 achieving this result is: ');
disp(c(2));
