function solve( prob, quiet )
if nargin < 2, quiet = false; end

global cvx___
prob = index( prob );
p = cvx___.problems( prob );

if length( p.objective ) > 1,
    error( 'Multiobjective problems are not supported.' );
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
    if any( size( c ) < [ 1, n ] ), 
        c( 1, n ) = 0; 
    end
    switch p.direction,
        case 'minimize', sign = +1;
        case 'maximize', sign = -1;
    end
    d = sign * c( :, 1 ); 
    c = sign * c( :, 2 : end )';
end

[ m, n ] = size( A );
n = n - 1;
x = zeros( n, 1 );
y = zeros( m, 1 );

b = - A( :, 1 ); 
A = + A( :, 2 : end );

if m > n & n > 0,
    
    x( : ) = NaN;
    y( : ) = NaN;
    value  = NaN;
    status = 'Overdetermined equality constraints; problem is likely infeasible';
    warning( 'Overdetermined equality constraints; problem is likely infeasible.' );

elseif n == 0,
    
    %
    % No variables
    %

    if any( b ~= 0 ),

        status = 'Infeasible';
        y = - P * sign( b );
        value = Inf;

    else,

        status = 'Solved';
        value = d;

    end

else,
    
    [ value, x, y, status ] = cvx_solve_sedumi( A, b, c, d, p.cones, quiet );

end

p.x = full( [ 1 ; x ] );
p.y = full( [ 1 ; y ] );
p.result = full( sign * value );
p.status = status;
cvx___.problems( prob ) = p;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
