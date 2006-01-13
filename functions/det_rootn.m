function cvx_optval = det_rootn( X )
error( nargchk( 1, 1, nargin ) );

% DET_ROOTN    n-th root of the determinant of a symmetric matrix.
%     For a square matrix X, DET_ROOTN(X) returns
%         POW(DET(X),1/(size(X,1))
%     if X is symmetric (real) or Hermitian (complex) and positive
%     semidefinite, and -Inf otherwise.
%
%     This function can be used in many convex optimization problems that
%     call for LOG(DET(X)) instead. For example, if the objective function
%     contains nothing but LOG(DET(X)), it can be replaced with
%     DET_ROOTN(X), and the same optimal point will be produced.
%
%     Disciplined convex programming information:
%         DET_ROOTN is concave and nonmonotonic; therefore, when used in
%         CVX specifications, its argument must be affine.

if ndims( X ) > 2,
    
    error( 'N-D arrays are not supported.' );
    
elseif diff( size( X ) ) ~= 0,
    
    error( 'Matrix must be square.' );
    
elseif cvx_isconstant( X ),
    
    cvx_optval = det_rootn( cvx_constant( X ) );
    
elseif isreal( X ),

    n = size(X,1);
    cvx_begin
        variable Z(n,n) lower_triangular
        D = diag( Z );
        maximize geomean( D )
        subject to
            [ diag( D ), Z' ; Z, X ] == semidefinite(2*n);
    cvx_end
    
else,
    
    n = size(X,1);
    cvx_begin
        variable Z(n,n) lower_triangular complex
        D = diag( Z );
        maximize geomean( D )
        subject to
            [ diag( D ), Z' ; Z, X ] == hermitian_semidefinite(2*n);
    cvx_end
    
end
