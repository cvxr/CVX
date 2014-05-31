function z = cvx_pushdual( name, reps )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end

if ~isvarname( name ),
    if isempty( name ),
        cvx_throw( 'Anonymous dual variables are not allowed.' );
    elseif ~ischar( name ) || size( name, 1 ) ~= 1,
        cvx_throw( 'Second argument must be a string or a subscript structure array.' );
    else
        cvx_throw( 'Invalid dual variable name: %s', name );
    end
elseif isfield( pstr.duals, name ),
    cvx_throw( 'Dual variable name conflict: %s', name );
elseif isfield( pstr.variables,name ),
    cvx_throw( 'Primal/dual variable name conflict: %s', name );
end

duals = pstr.duals;
if ~isempty( reps ),
    duals.(name) = cell( reps );
else
    duals.(name) = [];
end
pstr.duals = duals;
dvars = pstr.dvars;
z = cvxdual( length(cvx___.problems), name );
dvars.(name) = z;
pstr.dvars = dvars;
cvx___.problems(end) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
