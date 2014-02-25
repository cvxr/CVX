function cvx_optval = berhu_pos( x, M, t ) %#ok

%BERHU   Internal cvx version.

%
% Check arguments
%

error( nargchk( 1, 3, nargin ) ); %#ok
if ~cvx_isconvex( x ),
    error( 'Disciplined convex programming error:\n    X must be convex.' );
end
if nargin < 2,
    M = 1;
elseif isa( M, 'cvx' ),
    error( 'Second argument must be numeric.' );
elseif ~isreal( M ) || any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end
if nargin < 3,
    t = 1;
elseif ~isreal( t ),
    error( 'Third argument must be real.' );
elseif cvx_isconstant( t ) && nnz( cvx_constant( t ) <= 0 ),
    error( 'Third argument must be real and positive.' );
elseif ~cvx_isconcave( t ),
    error( 'Disciplined convex programming error:\n   T must be affine or concave.', 1 ); %#ok
end
sz = cvx_size_check( x, M, t );
if isempty( sz ),
    error( 'Sizes are incompatible.' );
end

%
% Compute result
%

v = []; w = [];
cvx_begin separable
    variable v( sz ) nonnegative
    variable w( sz ) nonnegative
    minimize( quad_over_lin( w, t, 0 ) ./ ( 2 * M ) + w + v )
    x <= w + v; %#ok
    v <= M .* t; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
