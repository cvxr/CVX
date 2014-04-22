function y = cvx_readlevel( x )

global cvx___
b = x.basis_;
s = size( b );
[ r, c ] = find( b );
y = max( sparse( r, c, cvx___.readonly( r ), s(1), s(2) ), [], 1 );
y = reshape( full( y ), x.size_ );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
