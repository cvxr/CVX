function v = value( x, data )
b = x.basis_;
nb = size( b, 1 );
nx = size( data, 1 );
if nx < nb, 
    data( end + 1 : nb, : ) = NaN;
elseif nx > nb,
    data = data( 1 : nb, : );
end
v = reshape( data.' * b, x.size_ );
if any( x.size_ == 1 ) || ~cvx_use_sparse( v ),
    v = full( v );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
