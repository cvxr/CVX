function cvx_optval = lp_opt_val( b )

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
