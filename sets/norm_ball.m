function x = norm_ball( x, varargin ) %#ok

%NORM_BALL   Norm ball.
%   NORM_BALL( sz, ... ) returns a variable of size sz, say 'x', that is
%   constrained to satisfy NORM( x, ... ) <= 1. Any syntactically valid
%   and _convex_ use of the NORM() function has a direct analog in
%   NORM_BALL. The convex requirement specifically excludes, then, all
%   instances of NORM( x, p ) where p < 1.
%
%   See NORM for more detaills.
%
%   Disciplined convex programming information:
%       NORM_BALL is a cvx set specification. See the user guide for
%       details on how to use sets.

cvx_begin set
    variable x( size(x) )
    norm( x, varargin{:} ) <= 1; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
