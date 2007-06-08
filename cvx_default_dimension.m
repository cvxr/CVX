function y = cvx_default_dimension( sx )

% CVX_DEFAULT_DIMENSION   Default dimension for SUM, MAX, etc. 
%
% CVX_DEFAULT_DIMENSION( SX ) selects the dimension index DIM that would be used
% in a call to such functions like SUM, MAX, and so forth if none is supplied.
% SX is the size of the matrix being considered.
%
% For example suppose size(X) = [1,3,4]; then SUM(X) would sum over dimension
% 2; and CVX_DEFAULT_DIMENSION([1,3,4]) would return 2.

y = find( sx ~= 1 );
if isempty( y ), 
    y = 1; 
else
    y = y( 1 ); 
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
