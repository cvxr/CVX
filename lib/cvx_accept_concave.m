function x = cvx_accept_concave( x )
global cvx___
persistent cmap
if isempty( cmap ),
	cmap = cvx_remap('concave') & ~cvx_remap('affine');
end
tt = cmap(cvx_classify(x));
if ~any(tt(:)), return; end
prob = cvx___.problems(end).self;
src = x( tt );
dst = newtemp( prob, size( src ) );
x( tt ) = dst;
newcnstr( prob, src(:), dst(:), '>=' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
