% Equality constrained norm minimization.
%
% This script constructs a random equality-constrained norm minimization
% problem and solves it using CVX. You can also change p to +2 or +Inf
% to produce different results. Alternatively, you an replace
%     norm( A * x - b, p )
% with
%     norm_largest( A * x - b, 'largest', p )
% for 1 <= p <= 2 * n.

echo on

p = 1;
n = 20;
A = randn(2*n,n);
b = randn(2*n,1);
C = randn(0.5*n,n);
d = randn(0.5*n,1);
cvx_begin
   variable x(n)
   dual variable y
   minimize( norm( A * x - b, p ) )
   subject to
        y : C * x == d;
cvx_end

echo off

