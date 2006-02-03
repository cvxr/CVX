function solve( prob, quiet )
if nargin < 2, quiet = false; end

global cvx___
prob = index( prob );
p = cvx___.problems( prob );

objsize = size( p.objective );
nobj = prod( objsize );
if nobj > 1 & ~p.separable,
    error( 'Non-separable, multiobjective problems are not supported.' );
end

n = length( p.reserved );
A = cvx_basis( p.equalities );
if size( A, 2 ) < n, 
    if size( A, 1 ) == 0,
        A = sparse( [], [], [], 0, n );
    else,
        A( :, n ) = 0;
    end
end
if isempty( p.objective ),
    sign = +1;
    c = zeros( n - 1, 1 );
    d = 0;
else,
    c = cvx_basis( p.objective );
    if size( c, 2 ) < n,
        c( :, n ) = 0;
    end
    switch p.direction,
        case 'minimize', sign = +1;
        case 'maximize', sign = -1;
    end
    d = c( :, 1 ).'; 
    c = c( :, 2 : end ).';
end
lambda = ones( max( nobj, 1 ), 1 );

[ m, n ] = size( A );
n = n - 1;
x = zeros( n, 1 );
y = zeros( m, 1 );

b = - A( :, 1 ); 
A = + A( :, 2 : end );

if m > n & n > 0,
    
    x( : ) = NaN;
    y( : ) = NaN;
    value  = NaN * ones( objsize );
    status = 'Overdetermined';
    warning( 'Overdetermined equality constraints; problem is likely infeasible.' );
    pval = NaN;
    dval = NaN;

elseif n == 0,
    
    %
    % No variables
    %

    if any( b ~= 0 ),

        status = 'Infeasible';
        y = - P * sign( b );
        value = sign * Inf * ones( objsize );
        pval = 1;
        dval = 0;

    else,

        status = 'Solved';
        value = reshape( d, objsize );
        pval = 1;
        dval = 1;

    end

else,
    
    [ value, x, y, status ] = cvx_solve_sedumi( A, b, c * lambda, d * lambda, sign, p.cones, quiet );
    switch status,
    case { 'Solved', 'Inaccurate/Solved' },
        if nobj > 1, value = reshape( c' * x + d', objsize ); end
        pval = 1;
        dval = 1;
    case { 'Infeasible', 'Inaccurate/Infeasible' },
        if nobj > 1, value = sign * Inf * ones( objsize ); end
        pval = 1;
        dval = 0;
    case { 'Unbounded', 'Inaccurate/Unbounded' },
        if nobj > 1, value = -sign * Inf * ones( objsize ); end
        pval = 0;
        dval = 1;
    otherwise,
        if nobj > 1, value = NaN * ones( objsize ); end
        pval = NaN;
        dval = NaN;
    end

end

p.x = full( [ pval ; x ] );
p.y = full( [ dval ; y ] );
p.result = full( value );
p.status = status;
cvx___.problems( prob ) = p;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
