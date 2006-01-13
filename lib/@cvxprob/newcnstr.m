function newcnstr( prob, x, y, op, cheat )
error( nargchk( 4, 5, nargin ) );
if nargin < 5, cheat = false; end

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end
global cvx___
p = index( prob );

%
% Check arguments
%

cx = isnumeric( x ) | isa( x, 'cvx' );
if ~cx,
    x  = cvx_collapse( x, false, true );
    cx = isnumeric( x ) | isa( x, 'cvx' );
end
cy = isnumeric( y ) | isa( y, 'cvx' );
if ~cy,
    y  = cvx_collapse( y, false, true );
    cy = isnumeric( y ) | isa( y, 'cvx' );
end
sx = size( x );
sy = size( y );
if ~cx | ~cy,
    if ~iscell( x ) | ~iscell( y ),
        error( 'Cannot form a valid CVX constraint from this expression.' );
    elseif ~isequal( sx, sy ),
        error( 'The left- and right-hand sides have incompatible sizes.' );
    elseif op(1) ~= '=',
        error( 'Only equality constraints may use composite forms.' );
    else,
        nx = prod( sx );
        for k = 1 : nx,
            newcnstr( prob, x{k}, y{k}, op, cheat );
        end
        return
    end
else,
    [ prob, x, y ] = cvx_operate( prob, x, y );
end

%
% Check for a dual reference
%

dx = getdual( x );
if ~isempty( dx ),
    dy = getdual( y );
    if ~isempty( dy ),
        error( [ 'Two dual variables found: "', dx, '","', dy, '"' ] );
    end
else,
    dx = getdual( y );
end
if ~isempty( dx ),
    duals = cvx___.problems( p ).duals;
    try,
        dual = builtin( 'subsref', duals, dx );
    catch,
        nm = cvx_subs2str( dx );
        error( [ 'Dual variable "', nm(2:end), '" has not been declared.' ] );
    end
    if ~isempty( dual ),
        nm = cvx_subs2str( dx );
        error( [ 'Dual variable "', nm(2:end), '" already in use.' ] );
    end
end

%
% Handle the SDP case
%

mx = sx( 1 ) > 1 & sx( 2 ) > 1;
my = sy( 1 ) > 1 & sy( 2 ) > 1;
if op(1) ~= '=' & cvx___.problems( p ).sdp & ( mx | my ),
    if sx( 1 ) ~= sx( 2 ) | sy( 1 ) ~= sy( 2 ),
        error( 'SDP constraints must be square.' );
    elseif ~isequal( sx, sy ) & ( mx | cvx_isnonzero( x ) ) & ( my | cvx_isnonzero( y ) ),
        error( 'Left- and right-hand sides have incompatible sizes.' );
    elseif op( 1 ) == '>',
        x = pluslikeoper( x, y, 'minus', cheat );
        y = semidefinite( size( x ), ~isreal( x ) );
    else,
        y = pluslikeoper( y, x, 'minus', cheat );
        x = semidefinite( size( y ), ~isreal( y ) );
    end    
    op = '==';
end

%
% Detect redundancies and conflicts and act accordingly
%

z = pluslikeoper( x, y, 'minus', cheat );
zS = size( z );
zB = cvx_basis( z );
dim = length( cvx___.problems( p ).reserved );
[ zL, zR, zI ] = cvx_bcompress( zB );
if op( 1 ) == '=',
    temp = zR ~= 0;
    temp = temp( :, 1 ) & ( temp( :, 1 ) == sum( temp, 2 ) );
    if any( temp ),
        error( sprintf( 'Trivially infeasible equality constraints detected (e.g., 1 = 0).\n   Trivially infeasible constraints are not permitted.\n    See the user guide for more details.' ) );
    end
else,
    if any( zL( : ) < 0 ),
        if length( op ) < 2,
            error( sprintf( 'Conflicting strict inequalities detected (e.g., a < b & b < a )\n   Trivially infeasible constraints are not permitted.\n   See the user guide for more details.' ) );
        else,
            error( sprintf( 'Conflicting non-strict inequalities detected (e.g, a <= b & b <= a )\n   These must be converted to equalities in order to proceed.\n   See the user guide for more details.' ) );
        end
        temp = zR ~= 0;
        temp = temp( :, 1 ) == sum( temp, 2 );
        if any( temp ),
            if op( 1 ) == '<',
                temp = temp & ( zR( :, 1 ) > 0 );
            else,
                temp = temp & ( zR( :, 1 ) < 0 );
            end
            if any( temp ),
                error( sprintf( 'Trivally infeasible inequality constraints detected (e.g., 1 < 0).\n   Trivially infeasible constraints are not permitted.\n   See the user guide for more details.' ) );
            end
        end
    end
end
%
% Add slacks
%

nZ = size( zR, 1 );
if op( 1 ) ~= '=',
    ndxs = cvx_vexity( cvx( prob, nZ, zR ) ) == 0;
    if any( ndxs ),
        if all( ndxs ),
            v = [];
        else,
            ndxs = find( ndxs );
            temp = length( ndxs );
            v = sparse( ndxs, 1 : temp, 1, nZ, temp );
        end
        v = cvx_basis( newslack( prob, nZ, v ) );
        if op( 1 ) == '>', v = -v; end
        zR = [ zR, v( :, size( zR, 2 ) + 1 : end ) ];
    end
end

%
% Add the equalities
%

oeqs = length( cvx___.problems( p ).equalities ) + 1;
neqs = oeqs + nZ - 1;
if oeqs == 1,
    cvx___.problems( p ).equalities = cvx( prob, nZ, zR );
else,
    t1 = cvx_basis( cvx___.problems( p ).equalities ); 
    c1 = size( t1, 2 );
    c2 = size( zR, 2 );
    if c1 < c2, t1( end, c2 ) = 0; end
    if c2 < c1, zR( end, c1 ) = 0; end
    cvx___.problems( p ).equalities = cvx( prob, neqs, [ t1 ; zR ] );
    clear t1 t2
end
cvx___.problems( p ).x = [];
cvx___.problems( p ).y = [];

%
% Create the dual
%

if ~isempty( dx ),
    zI = zI * sparse( 1 : size( zI, 2 ), oeqs + 1 : neqs + 1, 1 );
    cvx___.problems( p ).duals = builtin( 'subsasgn', duals, dx, cvx( prob, zS, zI ));
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
