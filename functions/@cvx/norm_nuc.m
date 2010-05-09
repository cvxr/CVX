function cvx_optval = norm_nuc( X ) %#ok

%NORM_NUC   Internal cvx version.

error( nargchk( 1, 1, nargin ) );
if ndims( X ) > 2,
    error( 'norm_nuc is not defined for N-D arrays.' );
elseif ~cvx_isaffine( X ),
    error( 'Input must be affine.' );
end

%
% Construct problem
% 

[ m, n ] = size( X ); %#ok
cvx_begin sdp
    variable W1(m,m) symmetric
    variable W2(n,n) symmetric
    minimize(0.5*(trace(W1)+trace(W2)));
    [W1,X;X',W2] >= 0; %#ok
cvx_end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
