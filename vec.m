function v = vec( x )

% VEC   Vectorize.
%
% VEC(X), where X is a matrix or N-D array, returns a column
% vector containing the elements of X 'unrolled'. The elements
% are arranged in columnwise order. For example, if
% X = [ 1, 2 ; 3, 4 ], then vec(X) = [1;3;2;4].

v = reshape( x, numel( x ), 1 );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
