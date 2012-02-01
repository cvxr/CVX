function y = pow_pos( x, p )

%POW_POS   Internal cvx version.

error(nargchk(2,2,nargin));
if ~cvx_isconstant( p ),
    error( 'Second argument must be constant.' );
elseif ~isreal( p ),
    error( 'Second argument must be real.' );
end
p = cvx_constant( p );
if nnz( p < 1 ),
    error( 'Second argument must be greater than or equal to 1.\n(Use POW_P for exponents less than 1.)', 1 ); %#ok
end
y = pow_cvx( x, p, 'pow_pos' );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

