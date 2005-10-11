function H = hankel(c,r)
error( nargchk( 1, 2, nargin ) );

%
% Check arguments
%

error( cvx_verify( c ) );
if nargin < 2,
    r = zeros(size(c));
else,
    error( cvx_verify( r ) );
    [ prob, c, r ] = cvx_operate( [], c, r );
    temp = cvx_subref( r, 1 ) - cvx_subref( c, prod(size(c)) );
    if ~cvx_isconstant( temp ) | cvx_constant( temp ) ~= 0,
        warning('MATLAB:hankel:AntiDiagonalConflict',['Last element of ' ...
               'input column does not match first element of input row. ' ...
               '\n         Column wins anti-diagonal conflict.'])
    end
end

%
% Compute indices and construct data vector
% 

r = vec( r );
c = vec( c );
nc = length( c );
nr = length( r );
x = [ c ; cvx_subref( r, 2 : nr, 1 ) ];

%
% Construct matrix
%

cidx = [ 1 : nc ]';
ridx = 0 : nr - 1;
H    = cidx(:,ones(nr,1)) + ridx(ones(nc,1),:);
H    = reshape( cvx_subref( x, H( : ) ), size( H ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
