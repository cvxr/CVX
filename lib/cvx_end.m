function cvx_end

global cvx___
prob = evalin( 'caller', 'cvx_problem', '[]' );
if isempty( prob ) | ~isa( prob, 'cvxprob' ),
    error( 'No cvx problem exists in this scope.' );
elseif cvx___.stack{ end } ~= prob,
    error( 'Internal cvx data corruption.' );
end
p = index( prob );

if prob.complete,

    %
    % Check the integrity of the variables
    %

    bad_vars = {};
    vars = prob.variables;
    if ~isempty( vars ),
        fn = cvx_fieldnames( vars );
        temp = struct( 'type', '.' );
        for k = 1 : length( fn ),
            temp.subs = fn{k};
            nvar = evalin( 'caller', fn{k}, 'NaN' );
            if ~isa( nvar, 'cvx' ) | id( nvar ) ~= id( subsref( vars, temp ) ),
                bad_vars{end+1} = fn{k};
            end
        end
    end
    vars = prob.duals;
    if ~isempty( vars ),
        fn = fieldnames( vars );
        for k = 1 : length( fn ),
            nvar = evalin( 'caller', fn{k}, 'NaN' );
            if ~isa( nvar, 'cvxdual' ) | ~isequal( name( nvar ), fn{k} ) | problem( nvar ) ~= prob,
                bad_vars{end+1} = fn{k};
            end
        end
    end
    if ~isempty( bad_vars ),
        temp = sprintf( ' %s', bad_vars{:} );
        temp = sprintf( 'The following cvx variable(s) have been overwritten:\n  %s\nThis is often an indication that an equality constraint was\nwritten with one equals ''='' instead of two ''==''. The model\nmust be rewritten before cvx can proceed.', temp );
        eval( 'caller', 'cvx_clear' );
        error( temp );
    end

    if cvx___.pause,
        disp( ' ' );
        input( 'Press Enter/Return to call the solver:' );
        disp( ' ' );
    end

    %
    % Compress and solve
    %

    eliminate( prob );
    solve( prob, cvx___.quiet );

    if cvx___.pause & ~cvx___.quiet,
        disp( ' ' );
        input( 'Press Enter/Return to continue:' );
        disp( ' ' );
    end

    %
    % Fill the named primal and dual variables with their numeric values
    %

    vars = prob.variables;
    if ~isempty( vars ),
        vars = apply_value( vars, prob.x );
        fn   = cvx_fieldnames( vars );
        temp = struct( 'type', '.' );
        for k = 1 : length( fn ),
            temp.subs = fn{k};
            assignin( 'caller', fn{k}, subsref( vars, temp ) );
        end
    end
    assignin( 'caller', 'cvx_optpnt', vars );
    vars = prob.duals;
    if ~isempty( vars ),
        vars = apply_value( vars, prob.y );
        fn   = cvx_fieldnames( vars );
        for k = 1 : length( fn ),
            temp.subs = fn{k};
            assignin( 'caller', fn{k}, subsref( vars, temp ) );
        end
    end
    assignin( 'caller', 'cvx_optdpt', vars );
    assignin( 'caller', 'cvx_status', prob.status );
    assignin( 'caller', 'cvx_optval', prob.result );

else,

    %
    % Determine the parent problem; and while we're at it,
    % zero out the 'vexity' of the substitution variables
    % now that the problem is finished.
    %

    p = index( prob );
    if length( cvx___.stack ) < 2,
        error( 'Internal cvx data corruption.' );
    else,
        nprob = cvx___.stack{ end - 1 };
    end
    np = index( nprob );

    %
    % Compress the subproblem
    %

    eliminate( prob );

    %
    % Allocate the variable and cone sdata
    %

    srcm = 1 : length( prob.reserved );
    dstm = length( nprob.reserved );
    dstm = [ 1, dstm + 1 : dstm + length( srcm ) - 1 ];
    map  = sparse( srcm, dstm, 1 );
    cvx___.problems( np ).reserved( dstm, : ) = cvx___.problems( p ).reserved;
    cvx___.problems( np ).vexity( dstm, : )   = cvx___.problems( p ).vexity;
    cvx___.problems( np ).x( end + 1 : dstm( end ), : ) = NaN;
    ocones = cvx___.problems( np ).cones;
    cones = prob.cones;
    for k = 1 : length( cones ),
        cone = cones( k );
        cone.indices = reshape( dstm( cone.indices ), size( cone.indices ) );
        if isempty( ocones ),
            ocones = cone;
        else,
            match = find( strcmp( { ocones.type }, cone.type ) );
            if ~isempty( match ),
                match = match( cellfun( 'size', { ocones(match).indices }, 1 ) == size( cone.indices, 1 ) );
            end
            if ~isempty( match ),
                ocones(match(1)).indices = [ ocones(match(1)).indices, cone.indices ];
            else,
                ocones = [ ocones, cone ];
            end
        end
    end
    cvx___.problems( np ).cones = ocones;

    %
    % Map the named variables
    %

    nvars = cvx___.problems( np ).variables;
    base(1).type = '.';
    base(1).subs = [ prob.name, '_' ];
    base(2).type = '{}';
    try,
        ndxs = length( subsref( nvars, base(1) ) ) + 1;
    catch,
        ndxs = 1;
    end
    base(2).subs = { ndxs };
    vars = apply_map( prob.variables, map, nprob );
    vars = cvx_collapse( vars, true, false );
    nvars = builtin( 'subsasgn', nvars, base, vars );
    cvx___.problems( np ).variables = nvars;

    %
    % Process the substitutions
    %

    subs = cvx___.problems( p ).substitutions;
    if ~isempty( subs ),
        newcnstr( nprob, subs, vars.map_, '=', true );
    end

    %
    % Process the equality constraints
    %

    if ~isempty( prob.equalities ),
        newcnstr( nprob, apply_map( prob.equalities, map, nprob ), 0, '=' );
    end

    %
    % Process the objective
    %

    assignin( 'caller', 'cvx_optpnt', cvx_collapse( vars, false, false ) );
    if ~isempty( prob.objective ),
        assignin( 'caller', 'cvx_optval', apply_map( prob.objective, map, nprob ) );
    else,
        assignin( 'caller', 'cvx_optval', 0 );
    end

    %
    % Clear primal and dual variables from parent workspace.
    %

    if isempty( prob.variables ),
        prif = {};
    else
        prif = fieldnames( prob.variables );
    end
    if isempty( prob.duals ),
        dulf = {};
    else,
        dulf = fieldnames( prob.duals );
    end
    if ~isempty( dulf ) | ~isempty( prif ),
        evalin( 'caller', sprintf( '%s ', 'clear', dulf{:}, prif{:} ) );
    end
    assignin( 'caller', 'cvx_status', 'Incorporated into parent' );

end

% Pop the problem off the stack and clear it from internal storage

cvx___.stack( prob.stackpos : end ) = [];
cvx___.problems( p : end ) = [];
evalin( 'caller', 'clear cvx_problem' );
cvx_clearpath( 1 );

% END

function obj = apply_value( obj, x )

s = size( obj );
n = prod( s );
switch class( obj ),
    case 'cell',
        for k = 1 : n,
            obj{ k } = apply_value( obj{ k }, x );
        end
    case 'struct',
        obj = cell2struct( apply_value( struct2cell( obj ), x ), fieldnames( obj ) );
    otherwise,
        sz  = size( obj );
        obj = cvx_basis( obj );
        obj = obj * x( 1 : size( obj, 2 ), : );
        obj = reshape( obj, sz );
end

function obj = apply_map( obj, map, nprob )

s = size( obj );
n = prod( s );
switch class( obj ),
    case 'cell',
        for k = 1 : n,
            obj{ k } = apply_map( obj{ k }, map, nprob );
        end
    case 'struct',
        obj = cell2struct( apply_map( struct2cell( obj ), map, nprob ), fieldnames( obj ) );
    otherwise,
        obj = cvx_basis( obj );
        obj = obj * map( 1 : size( obj, 2 ), : );
        obj = cvx( nprob, s, obj );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
