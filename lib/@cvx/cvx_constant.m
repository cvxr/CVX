function y = cvx_constant( x )
y = reshape( full( x.basis_( 1, : ) ), x.size_ );
if cvx_use_sparse( y ), y = sparse( y ); end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
