function cvx_pushobj( dir, x )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    error( 'CVX:NoModel', 'No CVX model is present.' );
end

%
% Check problem
%

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

persistent remap_min remap_max remap_log remap_gp
if isempty( remap_max ),
    remap_min = cvx_remap( 'convex', 'l_convex' );
    remap_max = cvx_remap( 'concave', 'l_concave' );
    remap_log = cvx_remap( 'l_valid' ) & ~cvx_remap( 'constant' );
    remap_gp  = cvx_remap( 'g_valid' );
end
if ~isa( x, 'cvx' ) && ~isa( x, 'double' ) && ~isa( x, 'sparse' ),
    error( 'CVX:ArgError', 'Cannot accept an objective of type ''%s''.', class( arg ) );
elseif ~isreal( x ),
    error( 'CVX:ArgError', 'Expressions in objective functions must be real.' );
elseif isempty( x ),
    warning( 'CVX:EmptyObjective', 'Empty objective.' );
end
cx = cvx_classify( x );
switch dir,
    case { 'minimize', 'minimise' }, vx = remap_min( cx ); dir = 'minimize';
    case { 'maximize', 'maximise' }, vx = remap_max( cx ); dir = 'maximize';
    otherwise, error( 'CVX:ArgError', 'Invalid objective type: %s', dir );
end
if ~all( vx ),
    if ~epi, type = 'dir'; %#ok
    elseif dir > 0, type = 'epigraph';
    else type = 'hypograph'; end
    cvx_dcp_error( type, 'unary', cvx_fastref( x, vx == 0 ) );
end

%
% Store the objective
%

vx = remap_log( cx ) + remap_gp( cx );
if any( vx ),
    if any( diff( vx ) ),
        error( 'CVX:LogMix', 'Invalid mix of logarithm and non-logarithmic objectives.' );
    end
    x = log( x );
    vx = vx(1);
end
pstr.objective = x;
pstr.direction = dir;
pstr.geometric = vx;
if pstr.n_variable > 1 && pstr.complete,
    x = cvx_basis( x );
    nv = min(pstr.n_variable,size(x,1));
    if nnz(x(2:nv,:)), pstr.complete = false; end
end
cvx___.problems( end ) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
