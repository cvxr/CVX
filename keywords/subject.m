function subject( varargin )

%SUBJECT Implements the "subject to" keyword.
%   The keyword
%      SUBJECT TO
%   is a "no-op"---that is, it has no functional value. It is provided
%   solely to allow CVX models to read more closely to their mathematical
%   counterparts; e.g.
%      MINIMIZE( <objective> )
%      SUBJECT TO
%           <constraint1>
%           ...
%   It may be omitted without altering the model in any way.

if ~iscellstr( varargin ),
    error( 'SUBJECT TO must be used in command mode.' );
elseif nargin ~= 1 || ~strcmpi( varargin{1}, 'to' ),
    error( 'Syntax: subject to' );
elseif ~isa( evalin( 'caller', 'cvx_problem', '[]' ), 'cvxprob' ),
    error( 'SUBJECT TO can only be used within a CVX model.' );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
