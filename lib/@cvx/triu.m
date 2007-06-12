function x = triu( x, k )
if nargin < 2, k = 0; end

%
% Check inputs
%

s = x.size_;
if length( s ) > 2,
    error( 'The first argument must be 2-D.' );
elseif ~isnumeric( k ) | length( k ) ~= 1,
    error( 'The second argument must be an integer scalar.' );
end

%
% Determine the indices of the zeroed-out elements
%

rows = [ 1 : s( 1 ) ]';
rows = rows( :, ones( 1, s( 2 ) ) );
cols = [ 1 : s( 2 ) ];
cols = cols( ones( 1, s( 1 ) ), : );
ndxs = rows >= cols + k;
if ~any( ndxs ), return; end

%
% Zero them out
%

b = x.basis_;
b( :, ndxs ) = 0;
x = cvx( s, b );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
