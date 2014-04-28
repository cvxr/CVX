function x = newvar( prob, name, siz, str, geo )

global cvx___
p = prob.index_;
pstr = cvx___.problems( p );

% Check for conflicts
if ~isempty( name ),
    if ~isvarname( name ),
        error( 'CVX:Variable', 'Invalid variable name: %s', name );
    elseif isfield( pstr.variables, name ),
        error( 'CVX:Variable', 'Duplicate variable name: %s', name );
    elseif isfield( pstr.duals, name ),
        error( 'CVX:Variable', 'Primal/dual variable name conflict: %s', name );
    end
end

% Determine structure
len = prod( siz );
siz(end+1:2) = 1;
if nargin < 4 || isempty( str ),
    dof = len;
    str = [];
else
    temp = any( str, 2 );
    dof = full( sum( temp ) );
    if dof ~= length( temp ),
        str = str( temp, : );
    end
end

% Create the variable
x = cvx_pushvar( dof );
if nargin == 5 && ~isempty( geo ) && geo,
    x = cvx_pushexp( x );
end
x = sparse( x, 1 : dof, 1 );
if ~isempty( str ), x = x * str; end
x = cvx( siz, x );

% Save
if ~isempty( name ),
    cvx___.problems( p ).variables.(name) = x;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
