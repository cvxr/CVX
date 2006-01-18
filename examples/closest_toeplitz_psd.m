% Closest Toeplitz SDP search.
%
% This script finds a Toeplitz Hermitian PSD matrix that is closest to a
% given Hermitian matrix, as measured by the Frobenius norm. That is, for
% a given matrix P, it solves:
%    minimize   || Z - P ||_F
%    subject to Z >= 0
%
% Adapted from an example provided in the SeDuMi documentation. Notice
% the use of SDP mode to simplify the semidefinite constraint.

% The data. P is Hermitian, but is neither Toeplitz nor PSD.
P = [ 4,     1+2*j,     3-j       ; ...
      1-2*j, 3.5,       0.8+2.3*j ; ...
      3+j,   0.8-2.3*j, 4         ];
  
% Construct and solve the model
n = size( P, 1 );
cvx_begin sdp
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z >= 0 : Q;
cvx_end

% Display resuls
disp( 'The original matrix, P: ' );
disp( P )
disp( 'The optimal point, Z:' );
disp( Z )
disp( 'min( eig( Z ) ) (should be nonnegative):' );
disp( min( eig( Z ) ) )
disp( 'The optimal value, || Z - P ||_F:' );
disp( norm( Z - P, 'fro' ) );
disp( 'Complementary slackness: Z * Q, should be near zero:' );
disp( Z * Q )
