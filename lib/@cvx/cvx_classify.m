function v = cvx_classify( x )

global cvx___
v = cvx_classify_mex( x.basis_, cvx___.classes );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
