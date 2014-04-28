function y = exp_nc( x )

y = cvx( x.size_, cvx_getexp( x.basis_ ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
