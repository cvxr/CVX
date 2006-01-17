function cvx_optval = inv_pos( x )
error( nargchk( 1, 1, nargin ) );

% INV_POS   Inverse of a positive quantity.
%     For a real scalar X, INV_POS(X) returns 1/X if X is positive, and
%     +Inf otherwise.
%
%     For matrices and N-D arrays, the function is applied to each element.
%
%     Disciplined convex programming information:
%         INV_POS is convex and nonincreasing; therefore, when used in CVX
%         specifications, its argument must be concave (or affine).

if ~isreal( x ),
    
    error( 'Input must be real.' );
    
elseif cvx_isconstant( x ),
    
    cvx_optval = inv_pos( cvx_constant( x ) );
    
elseif cvx_isconcave( x ),
    
    cvx_begin
        variable y( size( x ) )
        minimize y
        quad_over_lin( 1, x ) <= y;
    cvx_end
    
else,
    
    error( 'Disciplined convex programming error:\n    INV_POS is convex and nonincreasing, so its argument must be concave.' );
        
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

