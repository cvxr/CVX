function y = entr( x )

%ENTR   Scalar entropy.
%   ENTR(X) returns an array of the same size as X with the unnormalized
%   entropy function applied to each element:
%                { -X.*LOG(X) if X > 0,
%      ENTR(X) = { 0          if X == 0,
%                { -Inf       otherwise.
%   If X is a vector representing a discrete probability distribution, then
%   SUM(ENTR(X)) returns its entropy.
%
%   Disciplined convex programming information:
%       ENTR(X) is concave and nonmonotonic in X. Thus when used in CVX
%       expressions, X must be real and affine. Its use will effectively 
%       constrain X to be nonnegative: there is no need to add an
%       additional X >= 0 to your model in order to enforce this.

global cvx___
if isa( x, 'cvx' ) & ~cvx___.expert,
    error( sprintf( 'Disciplined convex programming error:\n    Entropy is not yet supported.' ) );
end

error(nargchk(1,1,nargin));
y = -rel_entr( x, 1 );
