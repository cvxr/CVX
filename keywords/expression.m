function varargout = expression( nm, varargin )

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

prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    error( 'A cvx problem does not exist in this scope.' );
elseif nargin > 1,
    error( 'Too many input arguments.\nTrying to declare multiple expression holders? Use the EXPRESSIONS keyword instead.', 1 ); %#ok
end

%
% Step 1: separate the name from the parenthetical
%

xt = find( nm == '(' );
if isempty( xt ),
    x.name = nm;
    x.size = [1,1];
elseif nm( end ) ~= ')',
    error( 'Invalid expression specification: %s', nm );
else
    x.name = nm( 1 : xt( 1 ) - 1 );
    x.size = nm( xt( 1 ) + 1 : end - 1 );
end
if ~isvarname( x.name ),
    error( 'Invalid expression specification: %s', nm );
elseif x.name( end ) == '_',
    error( 'Invalid expression specification: %s\n   Variables ending in underscores are reserved for internal use.', nm );
elseif exist( [ 'cvx_s_', x.name ], 'file' ) == 2,
    error( [ 'Invalid expression specification: %s\n', ...
        '   The name "%s" is reserved as a matrix structure modifier,\n', ...
        '   which can be used only with the VARIABLE keyword.' ], nm, x.name );
end
tt = evalin( 'caller', x.name, '[]' );
if isa( tt, 'cvxobj' ) && cvx_id( tt ) >= cvx_id( prob ),
    error( 'Invalid expression specification: %s\n   Name already used for another CVX object.', nm );
end

%
% Step 2: Parse the size. In effect, all we do here is surround what is
% replace the surrounding parentheses with square braces and evaluate. All
% that matters is the result is a valid size vector. In particular, it
% need no be a simple comma-delimited list.
%

if ischar( x.size ),
    x.size = evalin( 'caller', [ '[', x.size, '];' ], 'NaN' );
    [ temp, x.size ] = cvx_check_dimlist( x.size, true );
    if ~temp,
        error( 'Invalid expression specification: %s\n   Dimension list must be a vector of finite nonnegative integers.', nm );
    end
end

%
% Step 3. Initialize
%

v = cvx( x.size, [] );
if nargout > 0,
    varargout{1} = v;
else
    assignin( 'caller', x.name, v );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
