function cvx_end

%CVX_END  Completes a cvx specification.
%   CVX_BEGIN marks the end of a new cvx model, and instructs cvx to
%   complete its processing. For standard, complete models, cvx will send
%   a transformed version of the problem to a solver to obtain numeric
%   results, and replace references to cvx variables with numeric values.

global cvx___
prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    error( 'No cvx problem exists in this scope.' );
elseif index( prob ) ~= length( cvx___.problems ),
    error( 'Internal cvx data corruption.' );
end
p = index( prob );
pstr = cvx___.problems( p );

if isempty( pstr.objective ) && isempty( pstr.variables ) && isempty( pstr.duals ) && nnz( pstr.t_variable ) == 1,

    warning( 'CVX:EmptyModel', 'Empty cvx model; no action taken.' );
    evalin( 'caller', 'cvx_pop( cvx_problem, ''none'' )' );

elseif pstr.complete && nnz( pstr.t_variable ) == 1,

    %
    % Check the integrity of the variables
    %

    if isempty( pstr.variables ),
        fn1 = cell( 0, 1 );
        vv1 = fn1;
    else
        fn1  = fieldnames( pstr.variables );
        ndxs = horzcat( fn1{:} );
        ndxs = ndxs( cumsum( cellfun( 'length', fn1 ) ) ) ~= '_';
        fn1  = fn1( ndxs );
        vv1  = struct2cell( pstr.variables );
        vv1  = vv1(ndxs);
    end
    if isempty( pstr.duals ),
        fn2 = cell( 0, 1 );
        vv2 = fn2;
    else
        fn2  = fieldnames( pstr.duals );
        ndxs = horzcat( fn2{:} );
        ndxs = ndxs( cumsum( cellfun( 'length', fn2 ) ) ) ~= '_';
        fn2  = fn2( ndxs );
        vv2  = struct2cell( pstr.dvars );
        vv2  = vv2(ndxs);
    end
    i1 = cvx_ids( vv1{:}, vv2{:} );
    i2 = sprintf( '%s,', fn1{:}, fn2{:} );
    i2 = evalin( 'caller', sprintf( 'cvx_ids( %s )', i2(1:end-1) ) );
    tt = i1 ~= i2;
    if any( tt ),
        vv = [ fn1 ; fn2 ];
        evalin( 'caller', 'cvx_clear' );
        temp = sprintf( ' %s', vv{tt} );
        error( 'The following cvx variable(s) have been overwritten:\n  %s\nThis is often an indication that an equality constraint was\nwritten with one equals ''='' instead of two ''==''. The model\nmust be rewritten before cvx can proceed.', temp ); %#ok
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

    solve( prob );
    pstr = cvx___.problems( p );

    %
    % Pause again!
    %

    if cvx___.pause && ~cvx___.quiet,
        disp( ' ' );
        input( 'Press Enter/Return to continue:' );
        disp( ' ' );
    end

    %
    % Copy the variable data to the workspace
    %

    if numel( pstr.objective ) > 1,
        if strfind( pstr.status, 'Solved' ),
            pstr.result = value( pstr.objective );
            if pstr.geometric, pstr.result = exp( pstr.result ); end
        else
            pstr.result = pstr.result * ones(size(pstr.objective));
        end
    end
    assignin( 'caller', 'cvx_optpnt',  pstr.variables );
    assignin( 'caller', 'cvx_optdpt',  pstr.duals );
    assignin( 'caller', 'cvx_status',  pstr.status );
    assignin( 'caller', 'cvx_optval',  pstr.result );
    assignin( 'caller', 'cvx_slvitr',   pstr.iters );
    assignin( 'caller', 'cvx_slvtol',     pstr.tol );
    
    %
    % Compute the numerical values and clear out
    %

    evalin( 'caller', 'cvx_pop( cvx_problem, ''value'' )' );

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
    if ~isempty( vars ) || ~isempty( dvars ),
        pname = [ pstr.name, '_' ];
        if ~isempty( vars ),
            try
                ovars = cvx___.problems(np).variables.(pname);
            catch
                ovars = {};
            end
            ovars{end+1} = vars;
            cvx___.problems(np).variables.(pname) = ovars;
        end
        if ~isempty( dvars ),
            try
                ovars = cvx___.problems(np).duals.(name);
            catch
                ovars = {};
            end
            ovars{end+1} = dvars;
            cvx___.problems(np).duals.(pname) = ovars;
        end
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

    assignin( 'caller', 'cvx_optpnt', cvxtuple( cvx_collapse( vars, false, false ) ) );
    assignin( 'caller', 'cvx_optdpt', cvxtuple( cvx_collapse( dvars, false, false ) ) );
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
        [ r, c ] = find( xB ); %#ok
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

    assignin( 'caller', 'cvx_status', 'Incorporated' );
    evalin( 'caller', 'cvx_pop( cvx_problem, ''none'' )' );

end

assignin( 'caller', 'cvx_cputime', cputime - pstr.cputime );
    
if isempty( cvx___.problems ) && cvx___.profile,
    profile off;
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
