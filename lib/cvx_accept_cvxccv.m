function x = cvx_accept_cvxccv( x )
global cvx___
persistent cmap
if isempty( cmap ),
	cmap = [0,1,-1]*cvx_remap( { 'affine' ; 'convex' ; 'concave' } );
end
t = cmap(cvx_classify(x));
if ~any(t(:)), return; end
prob = cvx___.problems(end).self;
tt = t ~= 0;
src = x( tt );
dst = newtemp( prob, size( src ) );
x( tt ) = dst;
t = t( tt );
newcnstr( prob, t .* src(:), t .* dst(:), '<=' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
