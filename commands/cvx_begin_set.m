function cvx_begin_set( varargin )

% CVX_BEGIN_SET    Starts a new cvx CVX specification.
%
% The command
%    CVX_BEGIN_SET
% can be used to mark the beginning of a set definition---a cvx feasibility
% problem intended to describe a set for use in other models. See the files
% in the cvx subdirectory sets/ for examples. It has actually been
% deprecated in favor of
%    CVX_BEGIN SET
% (two separate words). However, it will continue to be available for back-
% compatability reasons, so you are free to use either.

if ~iscellstr( varargin ), error( 'Arguments must be strings.' ); end
assignin( 'caller', 'cvx_problem', cvxprob( 'set', varargin{:} ) );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
