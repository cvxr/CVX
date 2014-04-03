function x = cvx_accept_cvxccv( x )
global cvx___
if ~isa( x, 'cvx' ), return; end
t = cvx_vexity( x ) > 0;
if nnz( t ) == 0, return; end
prob = cvx___.problems(end).self;
tt = t(:) ~= 0;
src = x( tt );
dst = newtemp( prob, size( src ) );
x( tt ) = dst;
t = nonzeros( t );
newcnstr( prob, t.*src(:), t.*dst(:), '<=' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
