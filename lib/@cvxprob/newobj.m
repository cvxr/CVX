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

if ~ischar( dir ) | size( dir, 1 ) ~= 1 | ~any( strcmpi( dir, { 'minimize', 'maximize' } ) ),
    error( 'The second argument must be either "minimize" or "maximize".' );
end

%
% Check objective expression
%

if ~cvx_isvalid( x ),
    error( 'Objective expression must be a valid cvx object.' );
elseif ~isreal( x ),
    error( 'Expressions in objective functions must be real.' );
elseif isempty( x ),
    warning( 'Empty objective.' );
end
    
%
% Convert to pure epigraph/hypograph form
%

sz = size( x );
switch dir,
   case 'minimize', nm = 'epi_'; dirn = +1;
   case 'maximize', nm = 'hyp_'; dirn = -1;
end
temp = newvar( prob, nm, sz, cvx_bcompress( cvx_basis( x ), true ) );
newcnstr( prob, temp, x, '=' );
cvx___.problems( p ).objective = temp;
[ r, c ] = find( cvx_basis( temp ) );
cvx___.problems( p ).vexity( c ) = dirn;
cvx___.problems( p ).vexity( 1 ) = 0;
cvx___.problems( p ).direction   = dir;
cvx___.problems( p ).x = [];
cvx___.problems( p ).y = [];

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
