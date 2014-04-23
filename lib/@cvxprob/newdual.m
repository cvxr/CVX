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
if ~isempty( reps ),
    y = cell( reps );
    [ y{:} ] = deal( cvx );
    z = cell( reps );
    ndxs = cell( 1, length( reps ) - ( reps(end) == 1 ) );
    [ ndxs{:} ] = ind2sub( reps, 1 : prod( reps ) );
    ndxs = vertcat( ndxs{:} );
    nstr(2).type = '{}';
    for k = 1 : prod( reps ),
        nstr(2).subs = sprintf( '%d,', ndxs(:,k) );
        nstr(2).subs = eval( [ '{', nstr(2).subs(1:end-1), '}' ] );
        z{k} = cvxdual( p, nstr );
    end
else
    y = [];
    z = cvxdual( p, nstr );
end
vars = pstr.dvars;
vars = builtin( 'subsasgn', vars, nstr(1), z );
pstr.dvars = vars;
vars = pstr.duals;
vars = builtin( 'subsasgn', vars, nstr(1), y );
pstr.duals = vars;
cvx___.problems(p) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
