function y = sum_log( x, dim )

%SUM_LOG Internal CVX version.

global cvx___
if ~cvx___.expert,
    error( sprintf( 'Disciplined convex programming error:\n    Logarithms are not yet supported.' ) );
end

error( nargchk( 1, 2, nargin ) );
if nargin == 2,
    y = size( x, dim ) * log( geo_mean( x, dim ) );
else
    y = length( x ) * log( geo_mean( x ) );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
