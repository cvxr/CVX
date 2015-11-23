function xi = cvx_invert_structure( x, compact )

%CVX_INVERT_STRUCTURE Compute a right-inverse of a structure mapping.
%Assumption: each column has no more than one nonzero element.

if nargin == 2 && compact,
    
    [ m, n ] = size( x );
    [ ii, jj, vv ] = find( x );
    [ is, nn ] = sort( ii );
    dd = [ true ; diff(is) ~= 0 ];
    nd = nn(dd);
    xi = sparse( 1 : size(x,1), jj(nd), 1.0 ./ conj( vv(nd) ), m, n );
    
else
    
    [ ii, jj, vv ] = find( x );
    tmp = 1.0 ./ full(sparse(ii,1,vv.*conj(vv),size(x,1),1));
    xi = sparse( ii, jj, vv .* tmp(ii), size(x,1), size(x,2) );
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
