function [ psd, X, Z ] = cvx_check_psd( X, mode )
err = X - X';
X   = 0.5 * ( X + X' );
if norm( err, 'fro' ) > 8 * eps * norm( X, 'fro' ),
    psd = false;
elseif nargin == 3 && isequal( mode, 'sym' ),
    psd = true;
elseif nargin < 3 || ~isequal( mode, 'eig' ),
    [ Z, p ] = chol( X );
    psd = p == 0;
else
    Z = eig( full( X ) );
    psd = all( Z >= 0 );
end
