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

global cvx___
if m == 1,

    xL = 1;
    xR = x;
    xI = 1;
    return
    
elseif cvx___.has_mex,
    
    [ ndxs, scls ] = cvx_bcompress_mex( sparse( x' ), magonly, nsrt );
    
else,

    %
    % Normalize rows by the sign of the first non-zero column
    %

    xsgni = sum( abs( x ), 2 );
    approx = 1;
    if ~magonly,
        [ c, r, v ] = find( x' );
        temp = r( v( : ) < 0 & [ 1 ; diff( r( : ) ) ] ~= 0 );
        xsgni( temp ) = - xsgni( temp );
    end
    xsgn = 1.0 ./ ( xsgni + ( xsgni == 0 ) );
    xR = spdiags( xsgn, 0, m, m ) * x;

    %
    % Sort the rows and determine the unique ones
    %

    nzrows = find( xsgni ~= 0 );
    switch length( nzrows ),
        case 0,
            d = [];
            dc = [];
        case 1,
            d = true;
            dc = 1;
        otherwise,
            nzcols = find( any( xR, 1 ) );
            for k = nzcols( end : -1 : 1 ),
               [ v, ind ] = sort( xR( nzrows, k ) );
               nzrows = nzrows( ind );
            end
            xRL = xR( nzrows, : );
            xRR = xRL( 2 : end, : );
            xRL( end, : ) = [];
            d = abs( xsgni( nzrows ) );
            d = sum( abs( xRL - xRR ), 2 ) > eps * ( d( 1 : end - 1 ) + d( 2 : end ) );
            d = full( [ true; d ] );
            dc = cumsum( d );
            if 1,
                % This produces the exact same ordering as the mex file,
                % but it is mathematically unnecessary
                [ nzrows, ind ] = sort( nzrows );
                [ v, ind ] = sort( dc( ind ) );
                nzrows = nzrows( ind );
            end
    end
    
    ndxs_src       = 1 : m;
    ndxs           = ndxs_src;
    temp           = nzrows( d );
    ndxs( nzrows ) = ndxs( temp( dc ) );
    scls = xsgni( ndxs_src ) .* xsgn( ndxs );
    
end

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
