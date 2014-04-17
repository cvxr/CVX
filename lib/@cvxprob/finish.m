function finish( prob )

global cvx___
try

p = prob.index_;
np = length( cvx___.problems );
if np < p,
    if np > 1, prob = cvx___.problems(1).self; mode = 'reset'; end
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
end
pstr = cvx___.problems( p );
if pstr.self ~= prob,
    prob = cvx___.problems(1).self; mode = 'reset';
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
elseif np > p,
    % Cruft left over from an interrupted model.
    pop( cvx___.problems(p+1).self, 'reset' );
end

if isempty( pstr.objective ) && isempty( pstr.variables ) && isempty( pstr.duals ) && nnz( pstr.t_variable ) == 1,

    pop( prob, 'none' );
    warning( 'CVX:EmptyModel', 'Empty cvx model; no action taken.' );

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
        vv2  = vv2( ndxs );
    end
    fn1 = [ fn1 ; fn2 ];
    i1  = cvx_ids( vv1{:}, vv2{:} );
    i2  = sprintf( '%s,', fn1{:} );
    try
        i2 = evalin( 'caller', sprintf( 'cvx_ids( %s )', i2(1:end-1) ) );
    catch
        i2 = zeros(1,numel(fn1));
        for k = 1 : length(fn1),
            try
                i2(k) = evalin( 'caller', sprintf( 'cvx_ids( %s )', fn1{k} ) );
            catch
            end
        end
    end
    if any( i1 ~= i2 ),
        temp = sprintf( ' %s', fn1{ i1 ~= i2 } );
        error( 'CVX:CorruptModel', 'The following cvx variable(s) have been cleared or overwritten:\n  %s\nThis is often an indication that an equality constraint was\nwritten with one equals ''='' instead of two ''==''. The model\nmust be rewritten before CVX can solve it.', temp );
    end

    clear pstr;
    solve( prob );
    pstr = cvx___.problems( p );
    
    if numel( pstr.objective ) > 1 && ~isempty(pstr.result),
        if strfind( pstr.status, 'Solved' ),
            pstr.result = value( pstr.objective );
            if pstr.geometric, pstr.result = exp( pstr.result ); end
        else
            pstr.result = pstr.result * ones(size(pstr.objective));
        end
    end

    assignin( 'caller', 'cvx_status',  pstr.status );
    assignin( 'caller', 'cvx_optval',  pstr.result );
    assignin( 'caller', 'cvx_optbnd',  pstr.bound );
    assignin( 'caller', 'cvx_slvitr',  pstr.iters );
    assignin( 'caller', 'cvx_slvtol',  pstr.tol );
    assignin( 'caller', 'cvx_cputime', cputime - pstr.cputime );
    evalin( 'caller', 'pop( cvx_problem, ''value'' )' );
    
elseif np < 2,

    prob = cvx___.problems(1).self; mode = 'reset';
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );

else
    
    np = p - 1;
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

    x = pstr.objective;
    if isempty( x ),

        assignin( 'caller', 'cvx_optval', 0 );
        temp = length( pstr.t_variable ) + 1 : length( cvx___.readonly );
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

        %
        % Set the vexity flags
        %

        cvx___.canslack( r, : ) = true;
        cvx___.readonly( r, : ) = np;
        cvx___.classes( r, : ) = int8( cvx___.classes( r, : ) + 3 * os );
        
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
        pop( prob, 'none' );
    end

end

evalin( 'caller', 'clear cvx_problem' );

catch exc

    if ~isempty( mode ), pop( prob, 'reset' ); end
    evalin( 'caller', 'clear cvx_problem' );
    rethrow( exc );
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
