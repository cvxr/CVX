function y = sum_log( x, dim )

%SUM_LOG Internal CVX version.

cvx_expert_check( 'sum_log', x );
if nargin < 2 || isempty( dim ),
	dim = cvx_default_dimension( size( x ) );
elseif ~cvx_check_dimension( dim ),
    error( 'Second argument must be a dimension.' );
end
error( nargchk( 1, 2, nargin ) ); %#ok
try
    y = size( x, dim ) * log( geo_mean( x, dim ) );
catch exc
    exc.message = strrep( exc.message, 'geo_mean', 'sum_log' );
    throw( exc );
end
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
