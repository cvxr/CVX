function y = svec( x, nrm )
error( cvx_verify( x ) );

if nargin < 2,
    nrm = 2;
elseif ~isnumeric( nrm ) | length( nrm ) ~= 1 | all( nrm ~= [ 1, 2, Inf ] ),
    error( 'Second argument must be 1, 2, or Inf.' );
end

[ xL, xR ] = cvx_bcompress( cvx_basis( x ) );
if nrm < Inf,
    xL = diag( xL' * xL );
    if nrm == 2,
        xL = sqrt( xL );
    end
    xR = diag( xL ) * xR;
end

y = cvx( problem( x ), size( xR, 1 ), xR );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
