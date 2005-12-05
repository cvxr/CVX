function [ xL, xR, xI ] = cvx_bcompress( x, magonly, nsrt )
error( nargchk( 1, 3, nargin ) );
if nargin < 2, magonly = false; end
if nargin < 3, nsrt = 0; end

%
% Quick exit for all-zero matrices
%
    
[ m, n ] = size( x );
if nnz( x ) == 0,
    xL = sparse( [], [], [], m, 0 );
    if nargout > 1,
        xR = sparse( [], [], [], 0, n );
        xI = xL;
    end
    return
end

%
% Separate real and imaginary parts
%

iscplx = ~isreal( x );
if iscplx,
    x = cvx_c2r( x, 1 );
    m = m * 2;
end

[ ndxs, scls ] = cvx_bcompress_mex( sparse( x' ), magonly, nsrt );
temp = scls ~= 0; 
scls( temp ) = 1.0 ./ scls( temp );
xL = sparse( 1 : m, ndxs, scls, m, m );
t2 = any( xL, 1 );
xL = xL( :, t2 );
    
if nargout > 1,
   if all( t2 ),
       xR = x;
       xI = xL;
   else,
       xR = x( t2, : );
       if nargout > 2,
           dof = size( xL, 2 );
           xI = xL * sparse( 1 : dof, 1 : dof, 1.0 ./ diag( xL' * xL ), dof, dof );
       end
   end
end

if iscplx,
    xL = cvx_r2c( xL, 1 );
    if nargout > 2,
        xI = cvx_r2c( xI, 1 );
    end
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
