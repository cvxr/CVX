function expressions( varargin )

%EXPRESSIONS Declares one or more CVX expression holders.
%   EXPRESSIONS x1 x2 x3 ..., where x1, x2, x3, etc. are valid
%   variable names, declares multiple cvx expression holders. It is
%   exactly equivalent to issuing a separate EXPRESSION command
%   for each x1, x2, x3, ...
%        
%   EXPRESSIONS allows the declaration of vector, matrix, and
%   array variables. 
%
%   For more information about expression holders, see the help for 
%   EXPRESSION or the CVX user guide.
%
%   Examples:
%      expressions x y z
%
%   See also EXPRESSION.

if ~iscellstr( varargin ),
    error( 'EXPRESSIONS must be used in command mode.' );
end
for k = 1 : nargin,
    evalin( 'caller', [ 'expression ', varargin{k} ] );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
