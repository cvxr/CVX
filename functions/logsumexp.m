function y = logsumexp( varargin )

%LOGSUMEXP    log(sum(exp(x))).
%   LOGSUMEXP(X) = LOG_SUM_EXP(X) = LOG(SUM(EXP(X)). We have replaced this
%   function with LOG_SUM_EXP to better match our function naming
%   conventions. Please start using it instead.

warning( 'CVX:Renamed', [ ...
    'The function "logsumexp" has been renamed "log_sum_exp". Please start\n', ...
    'using the new name. The old name will be removed in a future release.' ], 1 );

y = log_sum_exp( varargin{:} );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
