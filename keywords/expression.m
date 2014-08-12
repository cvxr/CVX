function vout = expression( varargin )

%EXPRESSION Declares a single CVX object for storing subexpressions.
%   EXPRESSION x
%   where x is a valid MATLAB variable nm, declares a scalar expression 
%   holder for the current cvx model. Like a variable, an expression holder
%   can be used in constraints and objectives, according to the DCP ruleset.
%   However, unlike a variable, an expression holder is initialized to
%   zero, because the intent is for it to hold intermediate computations.
%
%   EXPRESSION x(n1,n2,...,nk)
%   declares a vector, matrix, or array expression holder with dimensions
%   n1, n2, ..., nk, each of which must be nonnegative integers. The value
%   of the expression holder is initialized to zero.
%
%   Examples:
%      variable x y
%      expression z
%      z = 2 * x - y;
%
%   See also EXPRESSIONS.

if nargin > 1,
    cvx_throw( 'Too many input arguments.\nTrying to declare multiple expression holders? Use the EXPRESSIONS keyword instead.', 1 );
end

global cvx___

try

    evalin( 'caller', 'cvx_verify' );
    cvx___.args = { varargin, 1, [] }; %#ok
    args = evalin( 'caller', 'cvx_parse' );
    v = cvx_pushexpr( args );
    if nargout,
        vout = v;
    else
        assignin( 'caller', args(1).name, v );
    end
        
catch exc
    
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
