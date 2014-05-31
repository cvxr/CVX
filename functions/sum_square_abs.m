function y = sum_square_abs( varargin )

%SUM_SQUARE   Sum of squares.
%   For vectors, SUM_SQUARE_ABS(X) is the sum of the squares of the
%   absolute values of the elements of the vector; i.e., SUM(ABS(X).^2).
%
%   For matrices, SUM_SQUARE_ABS(X) is a row vector containing the
%   application of SUM_SQUARE_ABS to each column. For N-D arrays, the 
%   SUM_SQUARE_ABS operation is applied to the first non-singleton 
%   dimension of X.
%
%   SUM_SQUARE(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       If X is real, then SUM_SQUARE(X,...) is convex and nonmonotonic in
%       X. If X is complex, then SUM_SQUARE(X,...) is neither convex nor
%       concave. Thus, when used in CVX expressions, X must be affine. DIM
%       must be constant.

try
    varargin{end+1:2} = [];
    y = sum_square( varargin{:}, true );
catch exc
	cvx_throw( exc );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
