function y = cvx_value( x )
global cvx___
y = cellfun( @cvx_value, x, 'UniformOutput', false );

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
