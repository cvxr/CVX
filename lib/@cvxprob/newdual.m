function z = newdual( prob, name, reps )

%
% Check problem
%

global cvx___
[ p, pstr ] = verify( prob );
vars = pstr.duals;

%
% Check name
%

if isempty( name ),
    error( 'CVX:Dual', 'Anonymous dual variables are not allowed.' );
elseif ischar( name ),
    if ~isempty( name ),
        if size( name, 1 ) ~= 1,
            error( 'CVX:Dual', 'Second argument must be a string or a subscript structure array.' );
        elseif ~isvarname( name ),
            error( 'CVX:Dual', 'Invalid dual variable name: %s', name );
        elseif isfield( vars, name ),
            error( 'CVX:Dual', 'Dual variable name conflict: %s', name );
        elseif isfield( pstr.variables,name ),
            error( 'CVX:Dual', 'Primal/dual variable name conflict: %s', name );
        end
    end
end

%
% Add the variable to the problem
%

nstr = struct( 'type', '.', 'subs', name );
duals = pstr.duals;
if ~isempty( reps ),
    duals.(name) = cell( reps );
else
    duals.(name) = [];
end
pstr.duals = duals;
dvars = pstr.dvars;
z = cvxdual( p, nstr );
dvars.(name) = z;
pstr.dvars = dvars;
cvx___.problems(p) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
