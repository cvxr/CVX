echo on

% simple equality-constrained norm minimization problem

n = 20;
A = randn(2*n,n);
b = randn(2*n,1);
C = randn(0.5*n,n);
d = randn(0.5*n,1);
p = 1;
cvx_begin
   variable x(n)
   dual variable y
   minimize( norm( A * x - b, p ) )
   subject to 
        y : C * x == d;
cvx_end

% You can change p to 2 or Inf as well, and it will still work.
% Also try norm_lagest( A * x - b, k ) for 1 <= k <= 2 * n

echo off

