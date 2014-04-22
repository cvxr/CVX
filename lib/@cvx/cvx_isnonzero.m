function y = cvx_isnonzero( x, full )
if nargin < 2 || ~full,
    y = nnz( x.basis_ ) ~= 0;
else
    y = reshape( full( any( x.basis_, 1 ) ), x.size_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
