function solve( prob )

global cvx___
p = index( prob );
pr = cvx___.problems(p);
nobj = numel(pr);
if nobj > 1 & ~pr.separable,
    error( 'Non-separable multiobjective problems are not supported.' );
end
quiet = cvx___.problems(p).quiet;
[ At, cones, sgn, Q, P, esrc, edst, dualized ] = eliminate( prob, true );
if ~isempty( esrc ),
    error( 'Non-GP exp() and log() constraints are not supported.' );
end
clear pr

dbca = At;
c = At( :, 1 );
At( :, 1 ) = [];
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
for k = 1 : length( cones ),
    cones(k).indices = cones(k).indices - 1;
end

%
% Ferret out the degenerate and overdetermined problems
%

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
        if dualized,
            disp( sprintf( 'Note: for improved efficiency, %s is solving the dual problem.', solv ) );
        end
        disp( spacer );
    end
    opath = path;
    try
        path( [ getfield( cvx___.path.solvers, lsolv ), opath ] );
        if cvx___.profile,
            profile off
        end
        [ x, y, status ] = feval( sfunc, At, b, c, sgn, cones, quiet, prec );
        if cvx___.profile,
            profile resume
        end
    catch
        path( opath );
        rethrow( lasterror );
    end
    switch status,
    case { 'Solved', 'Inaccurate/Solved' },
        value = sgn * ( c' * x + d' );
        pval = 1;
        dval = 1;
    case { 'Infeasible', 'Inaccurate/Infeasible' },
        value = sgn * Inf;
        pval = NaN;
        dval = 0;
    case { 'Unbounded', 'Inaccurate/Unbounded' },
        value = -sgn * Inf;
        pval = 0;
        dval = NaN;
    otherwise,
        value = NaN;
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

if dualized,
    switch status,
        case 'Infeasible', status = 'Unbounded';
        case 'Unbounded',  status = 'Infeasible';
        case 'Inaccurate/Infeasible', status = 'Inaccurate/Unbounded';
        case 'Inaccurate/Unbounded',  status = 'Inaccurate/Infeasible';
    end
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
x = full( Q * [ pval ; x ] );
y = full( P * [ dval ; y ] );
if dualized,
    cvx___.x = y;
    cvx___.y = x(2:end);
else
    cvx___.x = x;
    cvx___.y = y(2:end);
end
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
