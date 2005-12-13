function cvx_optval = lp_opt_val( b )

% LP_OPT_VAL  An example of an incomplete problem specification.
%
% This function uses cvx to compute the value of a linear program.
% The value of the linear program depends upon the value of the
% function parameter 'b'. If b is a constant, then the problem is
% solved and a numeric value is returned. But in fact, this function
% is convex and nonincreasing in b. Therefore, it can be used in
% cvx constraints and objective functions, as long as b is affine
% or concave.

% b is passed in as argument
A = [1 2; 3 4; -3 1; -1 -2];
c = [-1 3]';

% A and c are numerically defined here
cvx_begin
    [m,n]=size(A);
    variable x(n)
    minimize (c'*x)
    subject to
        A*x <= b;
cvx_end
