function varargout = expressions( varargin )

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

if nargin < 1,
    cvx_throw( 'Incorrect syntax for EXPRESSIONS. Type HELP EXPRESSIONS for details.' );
elseif nargout && nargout ~= nargin,
    cvx_throw( 'Incorrect number of output arguments.' );
elseif ~iscellstr( varargin ),
    cvx_throw( 'All arguments must be strings.' );
end

global cvx___

try

    evalin( 'caller', 'cvx_verify' );
    cvx___.args = { varargin, nargin, [] };
    args = evalin( 'caller', 'cvx_parse' );
    for k = 1 : nargin,
        v = cvx_pushexpr( args(k) );
        if nargout,
            varargout{k} = v; %#ok
        else
            assignin( 'caller', args(k).name, v );
        end
    end
    
catch exc
    
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
    
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
