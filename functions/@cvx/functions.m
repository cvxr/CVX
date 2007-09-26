% CVX: Additional functions added by CVX.
%
%   These functions have been provided to expand the variety of constraints
%   and objectives that can be specified in CVX models. But in fact, they 
%   can be used with numeric arguments *outside* of CVX as well. The help 
%   text for each of these functions contains general information about the
%   computations it performs, as well as specific information about its
%   proper use in CVX models, as dictated by its convexity/concavity and
%   monotonicity properties.
%
%   For a list of Matlab's built-in functions that have been extended to
%   provide CVX support, type "help cvx/builtins".
%
%   berhu             - Reverse Huber penalty function.
%   det_inv           - Determinant of the inverse of an SPD matrix.
%   det_root2n        - 2nth-root of the determinant of an SPD matrix.
%   det_rootn         - nth-root of the determinant of an SPD matrix.
%   geomean           - Geometric mean.
%   huber             - Huber penalty function.
%   inv_pos           - Reciprocal of a positive quantity.
%   lambda_max        - Maximum eigenvalue of a symmetric matrix.
%   lambda_min        - Minimum eigenvalue of a symmetric matrix.
%   logsumexp_sdp     - SDP-based approximation of log(sum(exp(x))).
%   matrix_frac       - Matrix fractional function.
%   norm_largest      - Sum of the k largest magnitudes of a vector.
%   norms             - Computation of multiple vector norms.
%   norms_largest     - Computation of multiple norm_largest() norms.
%   polyenv           - Convex or concave envelope of a polynomial.
%   polyval_trig      - Evaluate a trigonometric polynomial.
%   pos               - Positive part.
%   pow_pos           - Convex/concave branches of the power function.
%   pow_abs           - Absolute value raised to a fixed power.
%   quad_form         - Quadratic form.
%   quad_over_lin     - Sum of squares over linear.
%   quad_pos_over_lin - Sum of squares of positives over linear.
%   sigma_max         - Maximum singular value.
%   square            - Square.
%   square_abs        - Square of absolute value.
%   square_pos        - Square of positive part.
%   sum_largest       - Sum of the largest k values of a vector.
%   sum_smallest      - Sum of the smallest k elements of a vector.
%   sum_square        - Sum of squares.
%   sum_square_abs    - sum of squares of absolute values.
%   sum_square_pos    - Sum of squares of positive parts.
%   vec               - Vectorize.

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
