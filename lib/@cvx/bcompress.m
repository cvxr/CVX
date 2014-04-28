function [ xR, x, sx ] = bcompress( x, mode, nsrt )
if nargin < 2 || isempty( mode ),
    mode = 'full';
end
if nargin < 3 || isempty( nsrt ),
    nsrt = 0;
end
sx = x.size_;
xb = x.basis_;
if nargout <= 1,
    xR = cvx_bcompress( xb, mode, nsrt );
else
    [ xR, xb ] = cvx_bcompress( xb, mode, nsrt );
    x = cvx( size( xb, 2 ), xb );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
