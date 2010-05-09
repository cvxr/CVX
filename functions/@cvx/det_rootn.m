function cvx_optval = det_rootn( X )

%DET_ROOTN   Internal cvx version.

error( nargchk( 1, 1, nargin ) );
n = size( X, 1 );
if ndims( X ) > 2,

    error( 'N-D arrays are not supported.' );

elseif size( X, 2 ) ~= n,

    error( 'Matrix must be square.' );

elseif nnz( X ) <= n && nnz( diag( X ) ) == nnz( X ),

    cvx_optval = geo_mean( diag( X ) );

elseif cvx_isconstant( X ),

    cvx_optval = cvx( det_rootn( cvx_constant( X ) ) );

elseif isreal( X ),

    cvx_begin
        variable Z(n,n) lower_triangular
        D = diag( Z );
        maximize( geo_mean( D ) );
        subject to
            [ diag( D ), Z' ; Z, X ] == semidefinite(2*n);
    cvx_end

else

    cvx_begin
        variable Z(n,n) lower_triangular complex
        D = diag( Z );
        maximize( geo_mean( real( D ) ) );
        subject to
            [ diag( D ), Z' ; Z, X ] == hermitian_semidefinite(2*n);
    cvx_end

end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
