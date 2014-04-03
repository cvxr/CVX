function x = pdom( x )

% PDOM    Positive domain.
%    For constant X, S, PDOM(X,S)=X if X>=0, and S*Inf if X<0. If S is not
%    supplied, then S=+1 is assumed.
%
%    Disciplined convex programming information:
%        PDOM(X) is intended to be used primarily as a helper function for
%        other CVX functions, and it has some strange behavior.
%        --- For arguments known to be nonnegative, PDOM(X)=X. That is,
%            PDOM(X) has no effect.
%        --- For affine and concave arguments, PDOM(X)=X, with the added
%            constraints that X is nonnegative. If X is known to be 
%            nonpositive, a warning is issued, since this means that the
%            only feasible value of X is X=0.
%        --- Complex affine and unsigned convex X are rejected.
%        This last rule makes sense if you consider that PDOM(X,-1) is 
%        concave and increasing and PDOM(X,1) is convex and sign-monotonic.
%        Thus PDOM(X) is, in a sense, a "permissive" fusion of PDOM(X,-1)
%        and PDOM(X,1), accepting any argument that would be accepted by
%        either specific form of the function.

if ~isreal( x ), error( 'Argument must be real.' ); end
if nargin < 2, s = 1; end
x( x < 0 ) = s * Inf;

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
