function [ xR, x, elims ] = cvx_bcompress( x, mode, num_sorted, num_ended )

[ m, n ] = size( x );
if nargin < 4 || isempty( num_ended ),
    num_ended = n;
end

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
        otherwise,        cvx_throw( [ 'Invalid normalization mode: ', mode ] );
    end
end

iscplx = ~isreal( x );
if iscplx,
    n = n * 2;
    num_sorted = num_sorted * 2;
    num_ended = num_ended * 2;
    if issparse( x ),
        [rr,cc,vi] = find( x );
        vr = real( vi );
        vi = imag( vi );
        vi(abs(vi)<8*eps*abs(vr)) = 0;
        x = sparse( [ rr ; rr ], [ 2 * cc - 1 ; 2 * cc ], [ vr ; vi ], m, n );
    else
        xR = real( x );
        xI = imag( x );
        xB = xR & xI;
        xB(xB) = abs(xI(xB)) < 8 * eps * abs(xR(xB));
        xI(xB) = 0;
        x = reshape( [ xR ; xI ], m, n );
    end
end

[ ndxs, scls ] = cvx_bcompress_mex( sparse( x ), mode, num_sorted, num_ended );
xR = sparse( ndxs, 1 : n, scls, n, n );
t2 = any( xR, 2 );
xR = xR( t2, : );

if nargout > 1 && ~all( t2 ),
    x = x( :, t2 );
end

if iscplx,
    xR = xR(:,1:2:end) + 1j * xR(:,2:2:end);
end

if nargout == 3,
    elims = find(ndxs~=(1:length(ndxs)));
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
