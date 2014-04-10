function w = cvx_sign( x, dz )

global cvx___
if 1,
    if nargin < 2, dz = 0; end
    w = cvx_sign_mex( x.basis_, cvx___.sign, dz );
    w = reshape( w, x.size_ );
else
    s = x.size_; %#ok
    if ~all( s ), w = zeros( s ); return; end
    b = x.basis_;
    nb = size( b, 1 );
    s = cvx___.sign;
    n = numel( s );
    if nb < n,
        s = s( 1 : nb, 1 );
    elseif n < nb,
        s( nb, 1 ) = 0;
        cvx___.sign = s;
    end
    w = ~any( b( ~s, : ), 1 );
    if ~isreal( b ), 
        w = w & ~any( imag( b ), 1 ); 
    end
    b = b( s ~= 0, w );
    if ~isempty( b ),
        q = full( nonzeros(s).' * b );
        q( abs( q ) ~= sum( abs( b ), 1 ) ) = 0;
        w(w) = sign(q);
    else
        w = +w;
    end
    if nargin == 2,
        w( ~any( x.basis_, 1 ) ) = dz;
    end
    w = reshape( w, x.size_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
