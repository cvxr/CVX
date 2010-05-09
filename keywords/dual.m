function varargout = dual( varargin )

%DUAL Declares one or more dual variables.
%
%   DUAL VARIABLE x
%   where x is a valid MATLAB variable name, declares a dual variable
%   for the current cvx problem. A dual variable with that name is added
%   to the problem and a cvxdual object with that name is created in
%   the current worksapce. An error is generated if a cvx problem is not
%   in the current workspace.
%
%   Note that a dual variable initially has no size; i.e., SIZE(x) returns
%   [0,0]. This is because the size of a dual variable is determined by the
%   constraint to which it is attached. To attach a dual variable to a
%   constraint, use the colon notation as follows:
%
%      variable x(n)
%      dual variable z
%      A * x == b : z
%
%   In this example, the dual variable z is attached to the equality
%   constraint A * x == b. Assuming SIZE(A) = [m,n] and SIZE(B) = [m,1],
%   the size of z will be set at [m,1] by this operation.
%
%   DUAL VARIABLES x1 x2 ... xk
%   where x1, x2, ..., xk are valid MATLAB variable names, declares
%   multiple dual variables.
%
%   DUAL VARIABLE x{s1,s2,...,sk}
%   DUAL VARIABLES x1{s1,s2,...,sk} x2{...
%   creates a cell array of dual variables, of dimension [s1,s2,..,sk].
%   This is useful if you want to create a sequence of related constraints,
%   using a FOR loop for example, and assign a separate dual variable to
%   each one; for example:
%
%       variable x(n)
%       dual variables z{10}
%       for k = 1 : 10,
%           A(:,:,k) * x <= b(:,k) : z{k}
%       end
%
%   Assuming SIZE(A) = [m,n,10] and SIZE(B) = [m,10], this loop assigns
%   each of the dual variables z{k} to a constraint and sets their sizes
%   to [m,1]. Note that the elements of the dual variable cell array need
%   not end up the same size.
%
%   See also VARIABLE, VARIABLES.

if ~evalin( 'caller', 'exist(''cvx_problem'',''var'')', '0' ),
    error( 'A cvx problem does not exist in this scope.' );
elseif ~iscellstr( varargin ),
    error( 'All arguments must be strings.' );
end

prob = evalin( 'caller', 'cvx_problem' );
if strcmp( varargin{1}, 'variable' ),
    error( nargchk( 2, 2, nargin ) );
    varargin(1) = [];
elseif strcmp( varargin{1}, 'variables' ),
    error( nargchk( 2, Inf, nargin ) );
    varargin(1) = [];
end

nargs = length( varargin );
if nargout > 0,
    error( nargoutchk( nargs, nargs, nargout ) );
end

for k = 1 : nargs,
    nm = varargin{k};
    xt = find( nm == '{' );
    if isempty( xt ),
        x.name = nm;
        x.size = [];
    elseif nm(end) ~= '}',
        error( 'Invalid dual variable specification: %s', nm );
    else
        x.name = nm( 1 : xt( 1 ) - 1 );
        x.size = nm( xt( 1 ) + 1 : end - 1 );
    end
    if ~isvarname( x.name ),
        error( 'Invalid dual variable specification: %s', nm );
    elseif x.name( end ) == '_',
        error( 'Invalid dual variable specification: %s\n   Variables ending in underscores are reserved for internal use.', nm );
    end
    if ischar( x.size ),
        x.size = evalin( 'caller', [ '[', x.size, '];' ], 'NaN' );
        [ temp, x.size ] = cvx_check_dimlist( x.size, true );
        if ~temp,
            error( 'Invalid variable specification: %s\n   Dimension list must be a vector of finite nonnegative integers.', nm );
        end
    end
    temp = newdual( prob, x.name, x.size );
    if nargout > 0,
        varargout{k} = temp;
    end
    assignin( 'caller', x.name, temp );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
