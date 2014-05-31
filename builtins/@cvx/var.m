function y = var( x, w, dim )

%VAR    Internal cvx version.

if nargin < 3, dim = []; end
if nargin < 2, w = []; end

y = std( x, w, dim, true );

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
