function y = vertcat( varargin )

%   Disciplined convex/geometric programming information for VERTCAT:
%      VERTCAT imposes no convexity restrictions on its arguments.

y = cat( 1, varargin{:} );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
