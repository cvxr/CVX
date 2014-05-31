function cvx_display( pn, prefix )
if nargin < 2, prefix = '    '; end

global cvx___
p = cvx___.problems( pn );

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
nt = length( cvx___.classes );
ne = length( cvx___.equalities );
if ne <= p.checkpoint(2),
    neqns = 0;
    nineqs = 0;
    touched = false(nt,1);
else
    neqns  = cellfun( @(z)size(z,2), cvx___.equalities );
    nineqs = neqns(:) .* ~cvx___.inequality(:);
    neqns  = neqns( p.checkpoint(2) + 1 : end );
    nineqs = nineqs( p.checkpoint(2) + 1 : end );
    touched = cellfun( @(z)[any(z,2);false(nt-size(z,1),1)], ...
        cvx___.equalities(p.checkpoint(2)+1:end), 'UniformOutput', false );
    touched = any( [ touched{:} ], 2 );
end
touched( p.checkpoint(1) + 1 : end ) = true;
nineqs = sum( nineqs );
neqns  = sum( neqns ) - nineqs;
cfound = false;
for k = p.checkpoint(3) + 1 : length( cvx___.cones ),
    ndxs = cvx___.cones(k).indices;
    qq = any(reshape(touched(ndxs),size(ndxs)),1);
    if any( qq ),
        touched(ndxs(:,qq)) = true;
        cfound = true;
        break;
    end
end
nv = nnz(touched) - touched(1);
np = nnz(touched(1:p.checkpoint(1))) - touched(1);
gfound = nnz(ylog(cvx___.classes(touched)));

if isempty( p.name ) || strcmp( p.name, 'cvx_' ),
    nm = '';
else
    nm = [ p.name, ': ' ];
end
if ( p.gp ),
    ptype =' geometric ';
elseif ( p.sdp ),
    ptype = ' semidefinite ';
elseif abs( p.direction ) > 1,
    if p.direction > 0, 'log-convex ';
    else ptype = 'log-concave '; end
else
    ptype = ' ';
end
if isempty( p.objective ),
    tp = 'feasibility';
else
    if p.direction < 0,
        tp = 'maximization';
    else
        tp = 'minimization';
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
    if np,
        disp( [ prefix, '    + ', sprintf( '%d parent variables', np ) ] );
    end
elseif np,
    disp( [ prefix, sprintf( 'variables: %d parent variables', np ) ] );
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
            ndxs = ndxs( :, any( reshape( touched( ndxs ), size( ndxs ) ), 1 ) );
            if ~isempty( ndxs ),
                if size( ndxs, 1 ) == 1,
                    isvar = true;
                    ncones = numel( ndxs );
                else
                    isvar = false;
                    [ csize, ncones ] = size( ndxs );
                end
                if ncones == 1, 
                    plural = ''; 
                else 
                    plural = 's'; 
                end
                if isvar,
                    fprintf( 1, '%s   %d %s variable%s\n', prefix, ncones, ctyp, plural );
                else
                    fprintf( 1, '%s   %d order-%d %s cone%s\n', prefix, ncones, csize, ctyp, plural );
                end
            end
        end
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
