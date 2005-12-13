% MIN_PHASE_SPECTRAL_FACT   Minimal phase spectral factorization
% An SDP example from the SeDuMi documentation.

echo on

b = [2; 0.2; -0.3];
n = length( b );
cvx_begin
    variable X( n, n ) symmetric
    variable sums( n )
    dual variable y
    minimize( ( n - 1 : -1 : 0 ) * diag( X ) );
    for k = 0 : n - 1,
        sums( k + 1 ) == sum( diag( X, k ) );
    end
    y : sums == b;
    X == semidefinite( n );
cvx_end

echo off
