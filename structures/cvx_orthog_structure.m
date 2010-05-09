function xi = cvx_orthog_structure( x )

%CVX_ORTHOG_STRUCTURE   Computes the orthogonal basis to a matrix structure.

[m,n] = size(x);

% Compute the orthogonal projection matrix
xi = x';
xi = speye(n)-xi*((x*xi)\x);

% Reduced row echelon form
xi = cvx_cleanup_structure(xi);

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
