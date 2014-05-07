function [ xR, x ] = cvx_bcompress( x, mode, num_sorted )

if nargin < 3 || isempty( num_sorted ),
    num_sorted = 0;
end

if nargin < 2 || isempty( mode ),
    mode = 0;
else
    switch mode,
        case 'full',      mode = 0;
        case 'magnitude', mode = 1;
        case 'none',      mode = 2;
        otherwise,        error( [ 'Invalid normalization mode: ', mode ] );
    end
end

[ m, n ] = size( x ); %#ok
iscplx = ~isreal( x );
if iscplx,
    xR = real( x );
    xI = imag( x );
    xB = xR & xI;
    xB(xB) = abs(xI(xB)) < 8 * eps * abs(xR(xB));
    xI(xB) = 0;
    x = [ xR, xI ];
    n = n * 2;
end

[ ndxs, scls ] = cvx_bcompress_mex( sparse( x ), mode, num_sorted );
xR = sparse( ndxs, 1 : n, scls, n, n );
t2 = any( xR, 2 );
xR = xR( t2, : );

if nargout > 1 && ~all( t2 ),
    x = x( :, t2 );
end

if iscplx,
    xR = xR(:,1:end/2) + 1j * xR(:,end/2+1:end);
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
