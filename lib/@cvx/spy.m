function spy( x, mode )
error( cvx_verify( x ) );

switch nargin,
    case 0,
        error( 'Not enough arguments.' );
    case 1,
        mode = '';
    case 2,
        if ~ischar( mode ) | size( mode, 1 ) > 1,
            error( 'Second argument must be a string.' );
        end
end

b = cvx_basis( x );
s = size( x );

switch mode,
    case { '2-d', '2-D', '2d', '2D' },
        b( find( b ) ) = 1;
        b = sum( b, 2 );
        b = reshape( b, s );
    case { '', '3-d', '3-D', '3d', '3D' },
    	p = dimension( problem( x ) ) + 1;
        if size( b, 2 ) < p,
            b( end, n ) = 0;
        end
    otherwise,
        error( [ 'Unknown spy mode: ', mode ] );
end
        
spy( b );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
