% CVX: Additional nonlinear functions added by CVX
%   berhu             - Reverse Huber penalty function.
%   det_root2n        - 2nth-root of the determinant of an SPD matrix.
%   det_rootn         - nth-root of the determinant of an SPD matrix.
%   geomean           - Geometric mean. prod(x).^(1/length(x))
%   huber             - Huber penalty function.
%   inv_pos           - Reciprocal of a positive quantity. pos(x).^-1
%   lambda_max        - Maximum eigenvalue of a symmetric matrix.
%   lambda_min        - Minimum eigenvalue of a symmetric matrix.
%   logsumexp_sdp     - SDP-based approximation of log(sum(exp(x))).
%   matrix_frac       - Matrix fractional function. x'*inv(Y)*x, Y s.d.p.
%   norm_largest      - Sum of the k largest magnitudes of a vector.
%   norms             - Computation of multiple vector norms.
%   norms_largest     - Computation of multiple norm_largest() norms.
%   polyenv           - Evaluate the convex or concave envelope of a polynomial.
%   polyval_trig      - Evaluate a trigonometric polynomial.
%   pos               - Positive part. max(x,0)
%   quad_form         - Quadratic form. 0.5*x'*(Q+Q')*x=real(x'*Q*x)
%   quad_over_lin     - Sum of squares over linear. sum(x.^2)./y
%   quad_pos_over_lin - Sum of squares of positives over linear. sum(pos(x).^2)./y
%   sigma_max         - Maximum singular value. max(eig(X'*X))
%   square            - Square. x.^2
%   square_abs        - Square of absolute value. abs(x).^2
%   square_pos        - Square of positive part. pos(x).^2
%   sum_largest       - Sum of the largest k values of a vector.
%   sum_smallest      - Sum of the smallest k elements of a vector.
%   sum_square        - Sum of squares. sum(x.^2)
%   sum_square_abs    - sum of squares of absolute values. sum(abs(x).^2)
%   sum_square_pos    - Sum of squares of positive parts. sum(pos(x).^2)
%   vec               - Vectorize. reshape(x,numel(x),1)

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
