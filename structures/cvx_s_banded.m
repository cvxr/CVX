function y = cvx_s_banded( m, n, lower, upper )
%CVX_S_BANDED (U,L)-banded matrices.

if nargin < 3,
    error( 'Bandwidth arguments missing.' );
elseif ~isnumeric( lower ) || length( lower ) ~= 1 || lower < 0 || lower ~= floor( lower ),
    error( 'Bandwidth arguments must be nonnegative integers.' );
elseif nargin < 4, 
    upper = lower;
elseif ~isnumeric( upper ) || length( upper ) ~= 1 || upper < 0 || upper ~= floor( upper ),
    error( 'Bandwidth arguments must be nonnegative integers.' );
end

c  = 0 : n - 1;
c  = c( ones( 1, m ), : );
r  = ( 0 : m - 1 )';
r  = r( :, ones( 1, n ) );
temp = r - c;
temp = temp <= lower & temp >= -upper;
r  = r( temp );
c  = c( temp );
y = sparse( 1 : length( r ), r + m * c + 1, 1, length( r ), m * n );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.



