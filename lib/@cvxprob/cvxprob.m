function w = cvxprob( w )
global cvx___
if nargin,
    w = class( w, 'cvxprob' );
else
    w = cvx___.pobj;
    np = length(cvx___.problems);
    w.id_ = cvx___.problems(np).id;
    w.index_ = np;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
