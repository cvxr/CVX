function y = cvx_isconcave( x, full )
persistent remap
if isempty( remap ),
	remap = cvx_remap('concave');
end
y = remap( cvx_classify( x ) );
if nargin < 2 || ~full,
	y = all( y );
else
	y = reshape( y, x.size_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
