function x = cvx_fastref( x, ndxs )
x.basis_ = x.basis_(:,ndxs);
x.size_ = [size(x.basis_,2),1];

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
