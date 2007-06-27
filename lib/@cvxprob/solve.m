function solve( prob )

global cvx___
p = index( prob );
quiet = cvx___.problems(p).quiet;
[ At, cones, objsize, sgn, Q, P, esrc ] = eliminate( prob );
nobj = prod( objsize );

if nobj == 0,
    error( 'Shouldn''t have an empty objective here.' );
elseif nobj > 1 & ~cvx___.problems( p ).separable,
    error( 'Non-separable multiobjective problems are not supported.' );
elseif ~isempty( esrc ),
    error( 'Non-GP exp() and log() constraints are not supported.' );
end

c = At( :, 1 : nobj );
At( :, 1 : nobj ) = [];
d = c( 1, : );
c( 1, : ) = [];
[ n1, m ] = size( At );
if n1 < 1,
    b = zeros( m, 1 );
else
    b = - At( 1, : ).';
    At( 1, : ) = [];
end
n = n1 - 1;
lambda = ones( nobj, 1 );
tot_ineqs = 0;
for k = 1 : length( cones ),
    cones(k).indices = cones(k).indices - 1;
    tot_ineqs = tot_ineqs + length(cones(k).indices);
end

%
% Ferret out the degenerate and overdetermined problems
%

found = true;
tt = ( b' ~= 0 ) & ~any( At, 1 );
infeas = any( tt );
if m > n & n > 0,
    
    %
    % Overdetermined problem
    %
    
    x      = NaN * ones( n, 1 );
    y      = NaN * ones( m, 1 );
    value  = NaN * ones( objsize );
    status = 'Overdetermined';
    estr = sprintf( 'Overdetermined equality constraints detected.\n   CVX cannot solve this problem; but it is likely infeasible.' );
    if ~quiet,
        disp( estr );
    else
        warning( estr );
    end
    pval = NaN;
    dval = NaN;

elseif n ~= 0 & ~infeas & ( any( b ) | any( c ) ),
        
    %
    % Call solver
    %
    
    prob = cvx___.problems( p );
    solv = prob.solver;
    lsolv = lower(solv);
    prec = prob.precision;
    spacer = '-';
    spacer = spacer(:,ones(1,60));
    sfunc  = [ 'cvx_solve_', lsolv ];
    if ~quiet,
        disp( ' ' );
        disp( sprintf( 'Calling %s: %d variables, %d equality constraints', solv, n, m ) );
        disp( spacer );
    end
    opath = path;
%    try
        path( [ getfield( cvx___.path.solvers, lsolv ), opath ] );
        if cvx___.profile,
            profile off
        end
        [ x, y, status ] = feval( sfunc, At, b, c * lambda, sgn, cones, quiet, prec );
        if cvx___.profile,
            profile resume
        end
%    catch
%        path( opath );
%        rethrow( lasterror );
%    end
    switch status,
    case { 'Solved', 'Inaccurate/Solved' },
        value = sgn * reshape( c' * x + d', objsize );
        pval = 1;
        dval = 1;
    case { 'Infeasible', 'Inaccurate/Infeasible' },
        value = sgn * Inf * ones( objsize );
        pval = NaN;
        dval = 0;
    case { 'Unbounded', 'Inaccurate/Unbounded' },
        value = -sgn * Inf * ones( objsize );
        pval = 0;
        dval = NaN;
    otherwise,
        value = NaN * ones( objsize );
        pval = NaN;
        dval = NaN;
    end
    if ~quiet,
        disp( spacer );
    end
    
elseif infeas,
    
    %
    % Infeasible
    %
    
    if ~quiet,
        disp( 'Trivial infeasibilities detected; solution determined analytically.' );
    end
    status = 'Infeasible';
    x = NaN * ones( n, 1 );
    b( ~tt ) = 0;
    y = - b / ( b' * b );
    value = sgn * Inf * ones( objsize );
    pval = NaN;
    dval = 0;
    
else
    
    %
    % The origin is optional
    %
    
    if ~quiet,
        disp( 'Homogeneous problem detected; solution determined analytically.' );
    end
    status = 'Solved';
    x = zeros( n, 1 );
    y = zeros( m, 1 );
    value = sgn * reshape( d, objsize );
    pval = 1;
    dval = 1;
    
end

gvec = cvx___.problems( p ).geometric;
if nnz( gvec ),
    value( gvec ) = exp( value( gvec ) );
end
value = full( value );
if ~quiet,
    disp( sprintf( 'Status: %s', status ) );
    if length( value ) == 1,
        disp( sprintf( 'Optimal value (cvx_optval): %+g', value ) );
    else
        disp( sprintf( 'Optimal value (cvx_optval): (multiobjective)' ) );
    end
end

%
% Push the results into the master CVX workspace
%

global cvx___
cvx___.x = full( Q * [ pval ; x ] );
cvx___.y = full( P * [ dval * lambda ; y ] );
cvx___.problems( p ).result = value;
cvx___.problems( p ).status = status;
if nnz( cvx___.exponential ),
    esrc = find( cvx___.exponential );
    edst = cvx___.exponential( esrc );
    cvx___.x( edst ) = exp( cvx___.x( esrc ) );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
