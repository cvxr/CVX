function y = cvx_sign( x, dir )
if isreal( x ),
    y = sign( x );
else
    y = sign( real( x ) ) .* ( imag( x ) == 0 );
end
if nargin == 2,
    y( x == 0 ) = dir;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
