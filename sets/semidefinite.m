function x = semidefinite( sz, iscplx ) %#ok

%SEMIDEFINITE   Real symmetric positive semidefinite matrices.
%    SEMIDEFINITE(N), where N is an integer, creates a symmetric matrix
%    variable of size [N,N] and constrains it to be positive semidefinite.
%    Therefore, given the declaration
%       variable x(n,n) symmetric
%    the constraint
%       x == semidefinite(n)
%    is equivalent to
%       lambda_min(x) >= 0;
%    In fact, lambda_min is implemented in CVX using SEMIDEFINITE for
%    real matrices.
%
%    SEMIDEFINITE(SX), where SX is a valid size vector, creates an array
%    variable of size SX and constrains each subarray along the leading two
%    dimensions to be positive semidefinite. SX(1) and SX(2) must be equal.
%    Therefore, given the declaration
%       variable x(sx) symmetric
%    the constraint
%       x == semidefinite(sx)
%    is equivalent to
%       for k = 1:prod(sx(3:end)),
%          lambda_min(x(:,:,k)) >= 0;
%       end
%
%    SEMIDEFINITE(N,CPLX) and SEMIDEFINITE(SX,CPLX) create real semidefinite
%    sets if CPLX is FALSE, and complex Hermitian semidefinite sets if CPLX
%    is TRUE. The latter case is equivalent to calling the function
%    HERMITIAN_SEMIDEFINITE.
%
%   Disciplined convex programming information:
%       SEMIDEFINITE is a cvx set specification. See the user guide for
%       details on how to use sets.

%
% Check size vector
%

if ~isnumeric( sz ) || isempty( sz ) || any( sz < 0 ) || any( sz ~= floor( sz ) ),
    cvx_throw( 'First argument must be a nonnegative integer or a valid size vector.' );
elseif length( sz ) == 1,
    sz = [ sz, sz ];
elseif sz( 1 ) ~= sz( 2 ),
    cvx_throw( 'If a size vector is supplied, the first two dimensions must be equal.' );
end
if sz( 1 ) == 1,
    str = { 'nonnegative' };
else
    str = { 'semidefinite' };
end

%
% Check iscplx flag
%

if nargin > 1,
    if ( ~isnumeric( iscplx ) && ~islogical( iscplx ) ) || length( iscplx ) ~= 1,
        cvx_throw( 'Second argument must be a numeric or logical scalar.' );
    elseif iscplx && sz(1) > 1,
        str{end+1} = 'hermitian';
    end
end

%
% Construct the cone
%

cvx_begin set
    variable( 'x(sz)', str{:} );
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
