function v = newvar( prob, args )

global cvx___
p = prob.index_;
pstr = cvx___.problems( p );

name = args(1).name;
if ~isvarname( name ),
    error( 'CVX:Variable', 'Invalid variable name: %s', name );
elseif isfield( pstr.variables, name ),
    error( 'CVX:Variable', 'Variable name conflict: %s', name );
elseif isfield( pstr.duals,name ),
    error( 'CVX:Variable', 'Primal/dual variable name conflict: %s', name );
end

switch args(end).name,
    case 'epigraph_', iseh = +1; args(end) = [];
    case 'hypograph_', iseh = -1; args(end) = [];
    otherwise, iseh = 0;
end
if iseh,
    if isequal( pstr.direction, 'find' ),
        error( 'CVX:Variable', 'Epigraph/hypograph variables cannot be added to sets.' );
    elseif ~isempty( pstr.objective ),
        error( 'CVX:Variable', 'An objective has already been supplied for this problem.' );
    end
end
v = cvx_createvar( args, pstr.gp );

%
% Add variable to the problem structure
%

if iseh,
    if pstr.gp,
        vv = log( v );
    else
        vv = v;
    end
    if iseh > 0,
        dir = 'epigraph';
    else
        dir = 'hypograph';
    end
    pstr.objective = vv;
    pstr.direction = dir;
    pstr.geometric = pstr.gp;
end
pstr.variables.(name) = v;
cvx___.problems( p ) = pstr;

