% CHEBYSHEV_CENTER.M    computes the Chebyshev center of a polyhedron 
% Sec. 4.3.1, Boyd & Vandenberghe "Convex Optimization"
% Joëlle Skaf - 08/16/05
%
% The goal is to find the largest Euclidian ball (i.e. its center and 
% radius) that lies in a polyhedron described by linear inequalites in this
% fashion: P = {x : a_i'*x <= b_i, i=1,...,m}

% Input data
randn('state',0);
n = 10;
m = 2*n;
A = randn(m,n);
b = A*rand(n,1) + 2*rand(m,1);
norm_ai = sum(A.^2,2).^(.5);

fprintf(1,'Computing Chebyshev center...');
cvx_begin
    variable r(1)
    variable x_c(n) 
    dual variable y
    maximize ( r )
    y: A*x_c + r*norm_ai <= b
cvx_end

% Displaying results
fprintf(1,'Done! \n');
fprintf(1,'The Chebyshev center coordinates are: \n');
disp(x_c);
fprintf(1,'The radius of the largest Euclidian ball is: \n');
disp(r);
