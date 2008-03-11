function cvx_optval = berhu( x, M, t )

%BERHU   Internal cvx version.

%
% Check arguments
%

error( nargchk( 1, 3, nargin ) );
if ~cvx_isaffine( x ),
    error( sprintf( 'Disciplined convex programming error:\n    HUBER is nonmonotonic in X, so X must be affine.' ) );
end
if nargin < 2,
    M = 1;
elseif isa( M, 'cvx' ),
    error( 'Second argument must be numeric.' );
elseif ~isreal( M ) | any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end
if nargin < 3,
    t = 1;
elseif ~isreal( t ),
    error( 'Third argument must be real.' );
elseif cvx_isconstant( t ) & nnz( cvx_constant( t ) <= 0 ),
    error( 'Third argument must be real and positive.' );
elseif ~cvx_isconcave( t ),
    error( sprintf( 'Disciplined convex programming error:\n    HUBER is convex and nonincreasing in T, so T must be concave.' ) );
end
sz = cvx_size_check( x, M, t );
if isempty( sz ),
    error( 'Sizes are incompatible.' );
end

%
% Compute result
%

cvx_begin separable
    variables v( sz ) w( sz )
    minimize( quad_over_lin( w, t, 0 ) ./ (2*M) + v + w )
    abs( x ) <= w + v;
    v <= M * t;
    w >= 0;
cvx_end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
