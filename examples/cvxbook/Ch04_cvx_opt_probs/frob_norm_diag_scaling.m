% Section 4.5.4 example: Frobenius norm diagonal scaling
% Boyd & Vandenberghe "Convex Optimization"
% Joelle Skaf - 01/29/06
%
% Given a square matrix M, the goal is to find a vector (with dii > 0)
% such that ||DMD^{-1}||_F is minimized, where D = diag(d). 
% The problem can be cast into a geometric progam: 
%           minimize        sum_{i,j=1}^{n} Mij^2*di^2/dj^2 

% data generation 
randn('state',0);
n = 3;
M = randn(n)*2;

b = log(M.^2);
b = b(:);

A = zeros(n^2,n); 
for i=1:n
    for j=1:n
        A((j-1)*n+i, i) = A((j-1)*n+i, i) + 1;
        A((j-1)*n+i, j) = A((j-1)*n+i, j) - 1;
    end
end
A = 2*A;

% formulating the problem as a GP 
cvx_begin
    variables y(n)
    minimize ( logsumexp_sdp(A*y + b) )
cvx_end
d = exp(y);
D = diag(d); 

% displaying results 
disp('The matrix that minimizes ||DMD^{-1}||_F is: '); 
disp(D);
disp('The minimium Frobenius norm achieved is: ');
disp(norm(D*M*inv(D),'fro'));
disp('while the Frobunius norm of the original matrix M is: ');
disp(norm(M,'fro'));

