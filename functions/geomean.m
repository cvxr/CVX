function y = geomean( varargin )

%GEOMEAN   Geometric mean.
%   GEOMEAN(X) = GEO_MEAN(X) = PROD(X).^(1/LENGTH(X)). We have replaced this
%   function with GEO_MEAN to better match our function naming conventions.
%   Please start using it instead.

warning( 'CVX:Renamed', [ ...
    'The function "geomean" has been renamed "geo_mean". Please start\n', ...
    'using the new name. The old name will be removed in a future release.' ], 1 );

y = geo_mean( varargin{:} );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
