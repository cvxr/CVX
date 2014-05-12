% Exercise 4.47: Maximum determinant PSD matrix completion
% Boyd & Vandenberghe "Convex Optimization"
% Almir Mutapcic - Jan 2006
%
% Given a symmetric matrix A in R^(n-by-n) with some entries unspecified
% we find its completion such that A is positive semidefinite and
% it has a maximum determinant out of all possible completions.
% This problem can be formulated as a log det (and det_rootn) problem.
%
% This is a numerical instance of the specified book exercise.

% problem size
n = 4;

% create and solve the problem
cvx_begin
  % A is a PSD symmetric matrix (n-by-n)
  variable A(n,n) semidefinite;

  % constrained matrix entries.
  A(1,1) == 3; %#ok
  A(2,2) == 2; %#ok
  A(3,3) == 1; %#ok
  A(4,4) == 5; %#ok
  % Note that because A is symmetric, these off-diagonal
  % constraints affect the corresponding element on the
  % opposite side of the diagonal.
  A(1,2) == .5; %#ok
  A(1,4) == .25; %#ok
  A(2,3) == .75; %#ok

  % find the solution to the problem
  maximize( log_det( A ) )
  % maximize( det_rootn( A ) )
cvx_end

% display solution
fprintf('Matrix A with maximum determinant (%g) is:\n', det(A));
disp(A)
disp('Its eigenvalues are:')
disp(eig(A))

