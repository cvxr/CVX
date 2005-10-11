function z = mldivide( x, y )
error( cvx_verify( x, y ) );

%
% Check conformance with DCP rules
%

if ~cvx_isconstant( x ),
    error( 'Division by a non-cvx_constant expression is forbidden.' );
end

%
% Check size of divisor
%

sx = size( x );
if length( sx ),
    error( 'Input arguments must be 2-D.' );
elseif sx( 1 ) ~= sx( 2 ),
    error( 'Left-hand matrix must be square.' );
end

%
% Invert divisor
%

try
    x = inv( cvx_constant( x ) );
catch
    if all( sx == 1 ),
        error( 'Division by zero.' );
    else,
        error( 'Left-hand matrix is singular.' );
    end
end

%
% Multiply
%

z = mtimes( x, y );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
