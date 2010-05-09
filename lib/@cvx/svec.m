function z = svec( x, nrm )

if nargin < 2 || isempty( nrm ) || isequal( nrm, 'fro' ),
    nrm = 2;
elseif ~isnumeric( nrm ) || length( nrm ) ~= 1 || nrm < 1,
    error( 'Second argument must be a number between 1 and Inf, or ''fro''.' );
end

if ~isreal( x ) && nrm ~= 2,
    z = vec( x );
    return
else
    [ xR, y ] = bcompress( x );
    if isempty( y ),
        z = cvx( 0 );
    else
        z = y .* norms( xR, nrm, 2 );
    end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
