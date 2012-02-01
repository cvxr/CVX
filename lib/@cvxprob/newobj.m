function newobj( prob, dir, x )
error( nargchk( 3, 3, nargin ) );

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end
global cvx___
p = index( prob );
if ~isempty( cvx___.problems( p ).objective ),
    error( 'An objective has already been supplied for this problem.' );
end

%
% Check direction
%

if ~ischar( dir ) || size( dir, 1 ) ~= 1 || ~any( strcmpi( dir, { 'minimize', 'maximize' } ) ),
    error( 'The second argument must be either "minimize" or "maximize".' );
end

%
% Check objective expression
%

if ~isreal( x ),
    error( 'Expressions in objective functions must be real.' );
elseif isempty( x ),
    warning( 'CVX:EmptyObjective', 'Empty objective.' );
end

%
% Store the objective
%

persistent remap
if isempty( remap ),
    remap = cvx_remap( 'log-valid' ) & ~cvx_remap( 'constant' );
end
vx = remap( cvx_classify( x ) );
if any( vx ),
    if all( vx ),
        x = log( x );
    else
        x( vx ) = log( x( vx ) );
    end
end
if isa( x, 'cvx' ),
    zndx = any( cvx_basis( x ), 2 );
    v = cvx___.problems( p ).t_variable;
    cvx___.problems( p ).t_variable = v | zndx( 1 : size( v, 1 ), : );
end
cvx___.problems( p ).objective = x;
cvx___.problems( p ).direction = dir;
cvx___.problems( p ).geometric = vx;

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
