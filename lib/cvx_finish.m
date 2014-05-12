function cvx_finish

global cvx___
try
    p = length(cvx___.problems);
    pstr = cvx___.problems(p);
catch
    error( 'CVX:NoModel', 'No CVX model is present.' );
end

if length(cvx___.classes) == pstr.checkpoint(1) && ...
   length(cvx___.equalities) == pstr.checkpoint(2) && ...
   length(cvx___.cones) == pstr.checkpoint(3) && ...
   isempty(pstr.objective),
    
    warning( 'CVX:EmptyModel', 'Empty CVX model; no action taken.' );

elseif pstr.complete,

    %
    % Check the integrity of the variables
    %

    vars = pstr.variables;
    if isempty( vars ),
        fn1 = cell( 0, 1 );
        vv1 = fn1;
    else
        fn1  = fieldnames( vars );
        ndxs = cellfun( @(x)x(end)~='_', fn1 );
        fn1  = fn1( ndxs );
        vv1  = struct2cell( vars );
        vv1  = vv1(ndxs);
    end
    if isempty( pstr.duals ),
        fn2 = cell( 0, 1 );
        vv2 = fn2;
    else
        fn2  = fieldnames( pstr.duals );
        ndxs = cellfun( @(x)x(end)~='_', fn2 );
        fn2  = fn2( ndxs );
        vv2  = struct2cell( pstr.dvars );
        vv2  = vv2( ndxs );
    end
    fn1 = [ fn1 ; fn2 ];
    i1  = cellfun( @cvx_id, [ vv1 ; vv2 ] )';
    i2  = evalin( 'caller', [ 'cellfun( @cvx_id, {', sprintf( '%s ', fn1{:} ), '} )' ] );
    if any( i1 ~= i2 ),
        temp = sprintf( ' %s', fn1{ i1 ~= i2 } );
        error( 'CVX:CorruptModel', 'The following cvx variable(s) have been cleared or overwritten:\n  %s\nThis is often an indication that an equality constraint was\nwritten with one equals ''='' instead of two ''==''. The model\nmust be rewritten before CVX can solve it.', temp );
    end

    cvx_solve;
    pstr = cvx___.problems(p);

    if numel( pstr.objective ) > 1 && ~isempty( pstr.result ),
        if strfind( pstr.status, 'Solved' ),
            pstr.result = value( pstr.objective );
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

elseif p < 2,

    cvx_pop( 0 );
    error( 'CVX:IncompleteModel', 'Cannot construct a set object outside of a CVX model.' );

else

    np = p - 1;
    nstr = cvx___.problems(np);
    if nstr.gp ~= pstr.gp,
        error( 'CVX:MixedGP', 'Cannot embed a convex model into a geometric model, nor vice versa.' );
    end
    cp = pstr.checkpoint(4);
    if nstr.checkpoint(4) > cp
        cvx___.problems(np).checkpoint(4) = cp;
        cvx___.problems(np).complete = false;
    end
    vars = cvx_collapse( pstr.variables, true, false );
    dvars = cvx_collapse( pstr.duals, true, false );
    if ~isempty( vars ) || ~isempty( dvars ),
        pname = [ pstr.name, '_' ];
        if ~isempty( vars ),
            try
                ovars = nstr.variables.(pname);
            catch
                ovars = {};
            end
            ovars{end+1} = vars;
            nstr.variables.(pname) = ovars;
        end
        if ~isempty( dvars ),
            try
                ovars = nstr.duals.(name);
            catch
                ovars = {};
            end
            ovars{end+1} = dvars;
            nstr.(pname) = ovars;
        end
    end
    x = pstr.objective;
    if isempty( x ),
        assignin( 'caller', 'cvx_optval', 0 );
    else
        persistent mapsgn %#ok
        if isempty(mapsgn),
            mapsgn = ( cvx_remap( 'positive', 'p_nonconst' ) - cvx_remap( 'negative', 'n_nonconst' ) )';
        end    
        os = 1 - 2 * strcmp( pstr.direction, 'maximize' );
        [ r, c ] = find( cvx_basis( x ) ); %#ok
        if pstr.geometric,
            r = cvx___.logarithm( r, : );
        end
        cx = cvx___.classes( r, : );
        cx = int8( 9 + mapsgn(cx) + os * ( 3 + ( cx == 2 ) ) );
        cvx___.classes( r, : ) = cx;
        if pstr.geometric,
            cvx_pushexp( r, true );
        end
        assignin( 'caller', 'cvx_optval', x );
    end

end

pstr.finished = true;
cvx___.problems(end) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
