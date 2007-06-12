function [ ans, ndxs ] = cvx_fieldnames( x )
ans = fieldnames( x );
ndxs = true( 1, length( ans ) );
for k = 1 : length( ans ),
    str = ans{ k };
    ndxs( k ) = str( end ) ~= '_';
end
ans = ans( ndxs );
if nargout > 2,
    ndxs = find( ndxs );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
