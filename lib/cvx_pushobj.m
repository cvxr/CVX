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

persistent remap_min remap_max
if isempty( remap_max ),
    remap_min = cvx_remap( { { 'constant' } ; { 'l_convex' } ; { 'convex' } } );
    remap_max = cvx_remap( { { 'constant' } ; { 'l_concave' } ; { 'concave' } } );
end
if ~isa( x, 'cvx' ) && ~isnumeric( x ),
    error( 'CVX:ArgError', 'Cannot accept an objective of type ''%s''.', class( arg ) );
elseif ~isreal( x ),
    error( 'CVX:ArgError', 'Expressions in objective functions must be real.' );
elseif isempty( x ),
    warning( 'CVX:EmptyObjective', 'Empty objective.' );
end
cx = cvx_classify( x );
switch dir,
    case { 'minimize', 'minimise' }, 
        vx = remap_min( :, cx ); 
        dir = 'minimize';
    case { 'maximize', 'maximise' }, 
        vx = remap_max( :, cx ); 
        dir = 'maximize';
    otherwise, 
        error( 'CVX:ArgError', 'Invalid objective type: %s', dir );
end
avx = all( vx, 2 );
if ~any( avx ),
    if ~epi, type = 'dir'; %#ok
    elseif dir > 0, type = 'epigraph';
    else type = 'hypograph'; end
    vx = any( vx, 1 );
    if all( vx ),
        error( 'CVX:LogMix', 'Invalid mix of logarithmic and linear objectives.' );
    else
        cvx_dcp_error( type, 'unary', cvx_fastref( x, ~vx ) );
    end
end

%
% Store the objective
%

pstr.objective = cvx( x );
pstr.direction = dir;
pstr.geometric = avx(2) & ( ~avx(1) | pstr.gp );
if pstr.checkpoint(1) > 1 && pstr.complete,
    x = cvx_basis( x );
    nv = min(pstr.checkpoint(1),size(x,1));
    if nnz(x(2:nv,:)), pstr.complete = false; end
end
cvx___.problems( end ) = pstr;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
