function x = subsref( x, S )
error( nargchk( 2, 2, nargin ) );

try
    ndxs = builtin( 'subsref', reshape( 1 : prod( x.size_ ), x.size_ ), S );
catch
    error( lasterr );
end

x = cvx( size( ndxs ), x.basis_( :, ndxs ) );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
