function y = conj( x )

%Disciplined convex/geometric programming information for CONJ:
%   CONJ(X) imposes no convexity restrictions on its arguments. However,
%   since CONJ(X)=X when X is real, it is only useful for complex
%   affine expressions.

y = cvx( x.size_, conj( x.basis_ ) );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
