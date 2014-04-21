function newobj( prob, dir, x )

persistent remap_min remap_max remap_log
if isempty( remap_max ),
    remap_min = cvx_remap( 'convex', 'l_convex' );
    remap_max = cvx_remap( 'concave', 'l_concave' );
    remap_log = cvx_remap( 'l_valid' ) & ~cvx_remap( 'constant' );
end

%
% Check problem
%

global cvx___
[ p, pstr ] = verify( prob );
if ~isempty( pstr.direction ),
	if isequal( dir, 'find' ),
        error( 'CVX:ArgError', 'Objective functions cannot be added to sets.' );
    else
	    error( 'CVX:ArgError', 'An objective has already been supplied for this problem.' );
	end
end

%
% Check direction
%

if ~ischar( dir ) || size( dir, 1 ) ~= 1,
    error( 'CVX:ArgError', 'The second argument must be a string.' );
end

%
% Check objective expression
%

if ~isa( x, 'cvx' ) && ~isa( x, 'double' ) && ~isa( x, 'sparse' ),
    error( 'CVX:ArgError', 'Cannot accept an objective of type ''%s''.', class( arg ) );
elseif ~isreal( x ),
    error( 'CVX:ArgError', 'Expressions in objective functions must be real.' );
elseif isempty( x ),
    warning( 'CVX:EmptyObjective', 'Empty objective.' );
end
cx = cvx_classify( x(:) );
switch dir,
    case { 'minimize', 'minimise' }, vx = remap_min( cx ); dir = 'minimize';
    case { 'maximize', 'maximise' }, vx = remap_max( cx ); dir = 'maximize';
    otherwise, error( 'CVX:ArgError', 'Invalid objective type: %s', dir );
end
if ~all( vx ),
    cvx_dcp_error( dir, 'unary', cvx_subsref( x, vx == 0 ) );
end

%
% Store the objective
%

vx = remap_log( cx );
if any( vx ),
    if all( vx ),
        x = log( x );
    else
        x( vx ) = log( x( vx ) );
    end
end
if isa( x, 'cvx' ),
    zndx = any( cvx_basis( x ), 2 );
    v = pstr.t_variable;
    pstr.t_variable = v | zndx( 1 : size( v, 1 ), : );
end
pstr.objective = x;
pstr.direction = dir;
pstr.geometric = vx;
cvx___.problems( p ) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
