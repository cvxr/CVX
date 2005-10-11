function z = mrdivide( x, y )
error( cvx_verify( x, y ) );

%
% Check conformance with DCP rules
%

if ~cvx_isconstant( y ),
    error( 'Division by a non-cvx_constant expression is forbidden.' );
end

%
% Check size of divisor
%

sy = size( y );
if length( sy ) > 2,
    error( 'Input arguments must be 2-D.' );
elseif sy( 1 ) ~= sy( 2 ),
    error( 'Right-hand matrix must be square.' );
end

%
% Invert divisor
%

try
    y = inv( cvx_constant( y ) );
catch
    if all( sy == 1 ),
        error( 'Division by zero.' );
    else,
        error( 'Right-hand matrix is singular.' );
    end
end

%
% Multiply
%

z = mtimes( x, y );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
