function z = sparse( i, j, x, m, n )
error( nargchk( 3, 5, nargin ) );
if nargin == 4, error( 'Not enough arguments.' ); end

%
% Check sizes and indices
%

ni = prod( size( i ) );
nj = prod( size( j ) );
nk = prod( x.size_ );
nz = [ ni, nj, nk ];
if any( diff( nz( nz ~= 1 ) ) ~= 0 ),
    error( 'Vectors must be the same length.' );
elseif ( ~isnumeric( i ) & ~ischar( i ) ) | ( ~isnumeric( j ) & ~ischar( j ) ),
    error( 'Indices into a matrix must be numeric.' );
elseif any( i <= 0 ) | any( j <= 0 ) | any( i ~= floor( i ) ) | any( j ~= floor( j ) ),
    error( 'Indices into a matrix must be positive integers.' );
end
i = i( : );
j = j( : );
if nargin == 3,
    m = max( i );
    n = max( j );
elseif ( ~isnumeric( m ) & ~ischar( m ) ) | ( ~isnumeric( n ) & ~ischar( n ) ) | n < 0 | m < 0 | n ~= floor( n ) | m ~= floor( m ),
    error( 'Sparse matrix sizes must be positive integers.' );
elseif any( i > m ) | any( j > n ),
    error( 'Index exceeds matrix dimensions.' );
end

%
% Recalculate indices
%

nz = max( nz );
ij = i + m * ( j - 1 );
if length( ij ) < nz,
    ij = ij( ones( nz, 1 ), : );
end

%
% Reconstruct basis matrices
%

[ ix, jx, vx ] = find( x.basis_ );
xb = sparse( ix, ij(jx), vx, max(ix), m * n );
z = cvx( [ m, n ], xb );

%
% Verify vexity is preserved
%

v = cvx_vexity( z );
if any( isnan( v( : ) ) ),
    error( 'Disciplined convex programming error:\n    Sparse matrix construction produced invalid sums of convex and concave terms.' );
end
