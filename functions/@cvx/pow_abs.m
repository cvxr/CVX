function y = pow_abs( x, p )

%POW_ABS   Internal cvx version.

error(nargchk(2,2,nargin)); %#ok
if ~cvx_isconstant( p ),
    error( 'Second argument must be constant.' );
elseif ~isreal( p ),
    error( 'Second argument must be real.' );
elseif any( cvx_constant( p(:) ) < 1 ),
    error( 'Second argument must be greater than or equal to one.' );
end
y = pow_cvx( x, p, 'pow_abs' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
