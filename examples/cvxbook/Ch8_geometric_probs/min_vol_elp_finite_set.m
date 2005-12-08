% MIN_VOL_ELP_FINITE_SET Minimum volume ellipsoid covering a finite set
%                        (a figure is generated)
% Sec. 8.4.1, Boyd & Vandenberghe "Convex Optimization"
% Almir Mutapcic - 10/05
%
% Given a finite set of points x_i in R^2, we find the minimum volume 
% ellipsoid (described by matrix A and vector b) that covers all of
% the points by solving the optimization problem:
%
%           minimize     log det A^{-1}
%           subject to   || A x_i + b || <= 1   for all i
%
% Note: CVX still needs to handle log det problems

x = [ 0.55  0.0;
      0.25  0.35
     -0.2   0.2
     -0.25 -0.1
     -0.0  -0.3
      0.4  -0.2 ]';
[n, m] = size(x);

% Solving with CVX
cvx_begin
    variable A(n,n) symmetric
    variable b(n)
    minimize ( log det inv(A) )
    for k = 1:m
      norm ( A*x(:,k) + b ) <= 1
    end
cvx_end

% Plotting
