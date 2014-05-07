% Builds and solves a simple least-squares problem using cvx

echo on

n = 100;
A = randn(2*n,n);
b = randn(2*n,1);

cvx_begin
   variable x1(n)
   minimize( norm( A*x1-b ) )
cvx_end
v1 = cvx_optval;

cvx_begin
   variable x2(n)
   minimize( sum_square( A*x2-b ) )
cvx_end
v2 = cvx_optval;

% The differences between the unsquared and squared versions should be small
norm( x1 - x2 ) / norm( x1 ) %#ok
abs( v1 - sqrt(v2) ) / v1 %#ok

echo off




