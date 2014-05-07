function disp( prob, prefix )

if nargin < 2, prefix = ''; end

global cvx___
[ pn, p ] = verify( prob, false ); %#ok

if isempty( p.variables ),
    nvars = 0;
else
    nvars = length( fieldnames( p.variables ) );
end
if isempty( p.duals ),
    nduls = 0;
else
    nduls = length( fieldnames( p.duals ) );
end
persistent naff ylog
if isempty( naff ),
    naff = ~cvx_remap( 'affine' );
    ylog = cvx_remap( 'l_valid_' );
end
neqns  = cellfun( @(z)size(z,2), cvx___.equalities );
nineqs = neqns .* ~cvx___.inequality;
neqns  = sum( neqns - nineqs );
nineqs = sum( nineqs );
if isempty( p.name ) || strcmp( p.name, 'cvx_' ),
    nm = '';
else
    nm = [ p.name, ': ' ];
end

nt  = length( cvx___.classes );
fv = length( p.t_variable );
qv = fv + 1 : nt;
tt = p.t_variable;
ni = nnz( tt ) - 1;
nv = nt - fv + ni;
tt( qv ) = true;
gfound = nnz(ylog(cvx___.classes(tt)));
cfound = false;
for k = 1 : length( cvx___.cones ),
    if any( any( tt( cvx___.cones( k ).indices ) ) ),
        cfound = true;
        break;
    end
end

if all( [ numel( p.objective ), nv, nvars, nduls, neqns, nineqs, cfound, gfound ] == 0 ),
    disp( [ prefix, nm, 'empty CVX model' ] );
else
    if ( p.gp ),
        ptype =' geometric ';
    elseif ( p.sdp ),
        ptype = ' semidefinite ';
    else
        ptype = ' ';
    end
    if isempty( p.objective ),
        tp = 'feasibility';
    else
        switch p.direction,
            case 'minimize',  tp = 'minimization';
            case 'epigraph',  tp = 'epigraph minimization';
            case 'hypograph', tp = 'hypograph maximization';
            case 'maximize',  tp = 'maximization';
        end
        if numel( p.objective ) > 1,
            sz = sprintf( '%dx', size( p.objective ) );
            tp = [ sz(1:end-1), '-objective ', tp ];
        end
    end
    disp( [ prefix, nm, 'CVX', ptype, tp, ' model' ] );
    mlen = 0;
    if nvars > 0
        [ namep, sizep ] = cvx_dispvar( p.variables, '', false );
        mlen = max([mlen,max(cellfun(@numel,namep))]);
    end
    if nduls > 0
        [ named, sized ] = cvx_dispvar( p.duals, '', true );
        mlen = max([mlen,max(cellfun(@numel,named))]);
    end
    mlen = mlen + 3;
    spc = ' '; spc(1,2:mlen) = ' ';
    if nvars > 0,
        disp( [ prefix, 'variables: ' ] );
        for k = 1 : numel( namep ),
            disp( [ prefix, '   ', namep{k}, spc(1:mlen-length(namep{k})), sizep{k} ] );
        end
    end
    if nduls > 0,
        disp( [ prefix, 'dual variables: ' ] );
        for k = 1 : numel( named ),
            disp( [ prefix, '   ', named{k}, spc(1:mlen-length(named{k})), sized{k} ] );
        end
    end
    if neqns > 0 || nineqs > 0,
        disp( [ prefix, 'linear constraints:' ] );
        if neqns > 0,
            if neqns > 1, plural = 'ies'; else plural = 'y'; end
            fprintf( 1, '%s   %d equalit%s\n', prefix, neqns, plural );
        end
        if nineqs > 0,
            if nineqs > 1, plural = 'ies'; else plural = 'y'; end
            fprintf( 1, '%s   %d inequalit%s\n', prefix, nineqs, plural );
        end
    end
    if cfound || gfound,
        disp( [ prefix, 'nonlinearities:' ] );
        if gfound > 0,
            if gfound > 1, plural = 's'; else plural = ''; end
            fprintf( 1, '%s   %d exponential pair%s\n', prefix, gfound, plural );
        end
        if cfound,
            for k = 1 : length( cvx___.cones ),
                ctyp = cvx___.cones( k ).type;
                ndxs = cvx___.cones( k ).indices;
                ndxs = ndxs( :, any( reshape( tt( ndxs ), size( ndxs ) ), 1 ) );
                if ~isempty( ndxs ),
                    isint = isequal( ctyp(1:2), 'i_' );
                    if isint,
                        ncones = numel( ndxs );
                    elseif isequal( cvx___.cones( k ).type, 'nonnegative' ),
                        ncones = 1;
                        csize = numel(  ndxs  );
                    else
                        [ csize, ncones ] = size( ndxs );
                    end
                    if ncones == 1, plural = ''; else plural = 's'; end
                    if isint,
                        fprintf( 1, '%s   %d %s variable%s\n', prefix, ncones, ctyp(3:end), plural );
                    else
                        fprintf( 1, '%s   %d order-%d %s cone%s\n', prefix, ncones, csize, ctyp, plural );
                    end
                end
            end
        end
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
