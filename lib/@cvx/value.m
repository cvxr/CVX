function v = value( x, data )
global cvx___
b = x.basis_;
nb = size( b, 1 );
if nargin == 1, 
    data = cvx___.x; 
end
nx = size( data, 1 );
if nx < nb, 
    data( end + 1 : nb, : ) = NaN;
    if nargin == 1,
        cvx___.x = data;
    end
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
