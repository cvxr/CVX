function z = cvxprob( varargin )
try
    st = dbstack;
    depth = length( st );
    if depth <= 2, name = '';
    else name = st(3).name; end
    [ index, id ] = cvx_push( name, depth, varargin{:} );
    z = class( struct( 'index_', index, 'id_', id ), 'cvxprob' );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
