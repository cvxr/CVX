% Closest Toeplitz SDP search.
%
% An SDP example from the SeDuMi documentation, this script finds
% a Toeplitz Hermitian SDP matrix that is closest, in the Frobenius
% sense, to a given matrix.

echo on

P = [ 4, 1+2*j,3-j;1-2*j,3.5,0.8+2.3*j;3+j,0.8-2.3*j,4 ];
n = size( P, 1 );
cvx_begin
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z == hermitian_semidefinite( n ) : Q;
cvx_end

echo off

