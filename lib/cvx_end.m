function cvx_end

global cvx___
prob = evalin( 'caller', 'cvx_problem', '[]' );
if isempty( prob ) | ~isa( prob, 'cvxprob' ),
    error( 'No cvx problem exists in this scope.' );
elseif index( prob ) ~= length( cvx___.problems ),
    error( 'Internal cvx data corruption.' );
end
p = index( prob );
pstr = cvx___.problems( p );

if isempty( pstr.objective ) & isempty( pstr.variables ) & isempty( pstr.duals ) & nnz( pstr.t_variable ) == 1,

    warning( 'Empty cvx model; no action taken.' );
    evalin( 'caller', 'pop( cvx_problem, ''none'' )' );

elseif pstr.complete & nnz( pstr.t_variable ) == 1,

    %
    % Check the integrity of the variables
    %

    if isempty( pstr.variables ),
        fn1 = cell( 0, 1 );
        vv1 = fn1;
    else
        [ fn1, ndx1 ] = cvx_fieldnames( pstr.variables );
        vv1 = struct2cell( pstr.variables );
        vv1 = vv1(ndx1);
    end
    if isempty( pstr.duals ),
        fn2 = cell( 0, 1 );
        vv2 = fn2;
    else
        [ fn2, ndx2 ] = cvx_fieldnames( pstr.dvars );
        vv2 = struct2cell( pstr.dvars );
        vv2 = vv2(ndx2);
    end
    i1 = cvx_ids( vv1{:}, vv2{:} );
    i2 = sprintf( '%s,', fn1{:}, fn2{:} );
    i2 = evalin( 'caller', sprintf( 'cvx_ids( %s )', i2(1:end-1) ) );
    if any( i1 ~= i2 ),
        vv = [ fn1 ; fn2 ];
        temp = sprintf( ' %s', vv{i1~=i2} );
        temp = sprintf( 'The following cvx variable(s) have been overwritten:\n  %s\nThis is often an indication that an equality constraint was\nwritten with one equals ''='' instead of two ''==''. The model\nmust be rewritten before cvx can proceed.', temp );
        eval( 'caller', 'cvx_clear' );
        error( temp );
    end

    %
    % Pause
    %

    if cvx___.pause,
        disp( ' ' );
        input( 'Press Enter/Return to call the solver:' );
        disp( ' ' );
    end

    %
    % Compress and solve
    %

    clear pstr
    solve( prob, cvx___.quiet );
    pstr = cvx___.problems( p );

    %
    % Pause again!
    %

    if cvx___.pause & ~cvx___.quiet,
        disp( ' ' );
        input( 'Press Enter/Return to continue:' );
        disp( ' ' );
    end

    %
    % Copy the variable data to the workspace
    %

    assignin( 'caller', 'cvx_optpnt', pstr.variables );
    assignin( 'caller', 'cvx_optdpt', pstr.duals );
    assignin( 'caller', 'cvx_status', pstr.status );
    assignin( 'caller', 'cvx_optval', pstr.result );

    %
    % Compute the numerical values and clear out
    %

    evalin( 'caller', 'pop( cvx_problem, ''value'' )' );

else

    %
    % Determine the parent problem
    %

    if length( cvx___.problems ) < 2,
        error( 'Internal cvx data corruption.' );
    end
    np = p - 1;

    %
    % Move the variable structure into the parent
    %

    vars = cvx_collapse( pstr.variables, true, false );
    dvars = cvx_collapse( pstr.duals, true, false );
    if ~isempty( vars ) | ~isempty( dvars ),
        base(1).type = '.';
        base(1).subs = [ pstr.name, '_' ];
        base(2).type = '{}';
    end
    if ~isempty( vars ),
        nvars = cvx___.problems( np ).variables;
        base(2).subs = { eval( 'length(subsref(nvars,base(1)))+1', '1' ) };
        nvars = builtin( 'subsasgn', nvars, base, vars );
        cvx___.problems( np ).variables = nvars;
    end
    if ~isempty( dvars ),
        nvars = cvx___.problems( np ).duals;
        base(2).subs = { eval( 'length(subsref(nvars,base(1)))+1', '1' ) };
        nvars = builtin( 'subsasgn', nvars, base, dvars );
        cvx___.problems( np ).duals = nvars;
    end

    %
    % Merge the touch information
    %

    v = cvx___.problems( np ).t_variable;
    v = v | pstr.t_variable( 1 : size( v, 1 ), : );
    cvx___.problems( np ).t_variable = v;

    %
    % Process the objective and optimal point, converting to pure
    % epigraph/hypograph form if necessary
    %

    assignin( 'caller', 'cvx_optpnt', cvx_collapse( vars, false, false ) );
    assignin( 'caller', 'cvx_optdpt', cvx_collapse( dvars, false, false ) );
    x = prob.objective;
    if isempty( x ),

        assignin( 'caller', 'cvx_optval', 0 );
        temp = length( cvx___.problems( p ).t_variable ) + 1 : length( cvx___.readonly );
        cvx___.readonly( temp ) = cvx___.readonly( temp ) - 1;

    else

        switch pstr.direction,
            case 'minimize',
                force = false;
                os = +1;
            case 'epigraph',
                force = true;
                os = +1;
            case 'maximize',
                force = false;
                os = -1;
            case 'hypograph',
                force = true;
                os = -1;
        end
        
        if ~force,
            x = sparsify( x, 'objective' );
        end
        xB = cvx_basis( x );
        [ r, c ] = find( xB );
        t = r ~= 1; r = r( t );
        cvx___.canslack( r ) = true;

        %
        % Set the vexity flags
        %

        cvx___.vexity( r, : ) = os;
        cvx___.readonly( r, : ) = np;
        
        %
        % Convert to geometric form if necessary
        %

        tt = pstr.geometric;
        if any( tt ),
            if all( tt ),
                x = exp( x );
            else
                x( tt ) = exp( x( tt ) );
            end
        end

        assignin( 'caller', 'cvx_optval', x );
    end

    %
    % Set the status and clear the problem from internal storage
    %

    clear pstr
    assignin( 'caller', 'cvx_status', 'Incorporated' );
    evalin( 'caller', 'pop( cvx_problem, ''none'' )' );

end

if isempty( cvx___.problems ) & cvx___.profile,
    profile off;
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
