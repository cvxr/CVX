function cvx_optval = det_inv( X, p )

%DET_INV   Internal cvx version.

error( nargchk( 1, 2, nargin ) );
n = size( X, 1 );
if ndims( X ) > 2,
    error( 'N-D arrays are not supported.' );
elseif size( X, 2 ) ~= n,
    error( 'Matrix must be square.' );
elseif nargin < 2,
    p = 1;
elseif ~isnumeric( p ) || ~isreal( p ) || numel( p ) ~=  1 || p <= 0,
    error( 'Second argument must be a positive scalar.' );
end

w = [ ones(n,1) ; p ];
if cvx_isconstant( X ),
    
    cvx_optval = cvx( det_inv( cvx_constant( X ), p ) );

elseif nnz( X ) <= n && nnz( diag( X ) ) == nnz( X ),
    
    cvx_begin
        epigraph variable y
        geo_mean( [ diag(X) ; y ], w ) >= 1; %#ok
    cvx_end

elseif isreal( X ),

    cvx_begin
        epigraph variable y
        variable Z(n,n) lower_triangular
        D = diag( Z );
        [ diag( D ), Z' ; Z, X ] == semidefinite(2*n);
        geo_mean( [ D ; y ], [], w ) >= 1; %#ok
    cvx_end

else

    cvx_begin
        epigraph variable y
        variable Z(n,n) lower_triangular complex
        D = diag( Z );
        [ diag( D ), Z' ; Z, X ] == hermitian_semidefinite(2*n);
        geo_mean( [ real( D ) ; y ], [], w ) >= 1; %#ok
    cvx_end

end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
