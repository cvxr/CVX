function [ s, pattern ] = nnc( x )
error( cvx_verify( x ) );
ndxs = any( cvx_basis( x ), 1 );
ndxs = find( ndxs( 2 : end ) );
s = length( ndxs );
if nargin > 1,
    sz = size( x );
    s = reshape( sparse( ndxs, 1, 1, prod( sz ), 1 ), sz );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
