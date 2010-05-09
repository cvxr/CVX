function varargout = variable( nm, varargin )

%VARIABLE Declares a single CVX variable with optional matrix structure.
%   VARIABLE x
%   where x is a valid MATLAB variable nm, declares a scalar
%   variable for the current cvx problem. A variable with that
%   name is added to the problem, and a cvx object with that
%   name is created in the current workspace. An error is
%   generated if a cvx problem isn't in the current workspace.
%
%   VARIABLE x(n1,n2,...,nk)
%   declares a vector, matrix, or array variable with dimensions
%   n1, n2, ..., nk, each of which must be positive integers.
%
%   VARIABLE x(n1,n2,...,nk) mod1 mod2 mod3 ... modp
%   declares a vector, matrix, or array with structure. The
%   modifiers mod1, mod2, ... can each be one of the following:
%       complex   symmetric   skew-symmetric   hermitian
%       skew-hermitian   toeplitz   hankel   upper-hankel
%       lower-triangular   upper-triangular   tridiagonal
%       diagonal   lower-bidiagonal   upper-bidiagonal
%   Appropriate combinations of these modifiers can be chosen
%   as well. All except "complex" require that the matrix be
%   square. If an N-D (N>2) array is specified, then the matrix
%   structure is applied to each 2-D "slice" of the array.
%
%   Examples:
%      variable x(100,100) symmetric tridiagonal
%      variable z(10,10,10)
%      variable y complex
%
%   See also VARIABLES, DUAL, DUALS.

prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    error( 'A cvx problem does not exist in this scope.' );
end
global cvx___
p = index( prob );

%
% Step 1: separate the name from the parenthetical
%

xt = find( nm == '(' );
if isempty( xt ),
    x.name = nm;
    x.size = [1,1];
elseif nm( end ) ~= ')',
    error( 'Invalid variable specification: %s', nm );
else
    x.name = nm( 1 : xt( 1 ) - 1 );
    x.size = nm( xt( 1 ) + 1 : end - 1 );
end
if ~isvarname( x.name ),
    error( 'Invalid variable specification: %s', nm );
elseif x.name( end ) == '_',
    error( 'Invalid variable specification: %s\n   Variables ending in underscores are reserved for internal use.', nm );
elseif exist( [ 'cvx_s_', x.name ], 'file' ) == 2,
    error( [ 'Invalid variable specification: %s\n', ...
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
        error( 'Invalid variable specification: %s\n   Dimension list must be a vector of finite nonnegative integers.', nm );
    end
end

%
% Step 3. Parse the structure.
%

nmods = length( varargin );
filt = true( 1, nmods );
islin = false;
isgeo = false;
isepi = false;
ishypo = false;
modifiers = '';
for k = 1 : length( varargin ),
    strs = varargin{k};
    if isempty( strs ),
        continue;
    elseif ~ischar( strs ) || size( strs, 1 ) ~= 1,
        error( 'Matrix structure modifiers must be strings.' );
    end
    switch strs,
        case 'geometric_',
            isgeo = true;
            filt( k ) = false;
            continue;
        case 'linear_',
            islin = true;
            filt( k ) = false;
            continue;
        case 'epigraph_',
            isepi = true;
            filt( k ) = false;
            continue;
        case 'hypograph_',
            ishypo = true;
            filt( k ) = false;
            continue;
    end
    valid = true;
    xt = find( strs == '(' );
    if isempty( xt ),
        nm = strs;
    elseif strs( end ) ~= ')',
        valid = false;
    else
        nm = strs( 1 : xt(1) - 1 );
        strx.name = nm;
        strx.args = evalin( 'caller', [ '{', strs(xt(1)+1:end-1), '}' ], 'NaN' );
        if ~iscell( strx.args ),
            valid = false;
        end
        varargin{k} = strx;
    end
    if valid,
        valid = isvarname( nm ) & exist( [ 'cvx_s_', nm ], 'file' ) == 2;
    end
    if valid,
        modifiers = [ modifiers, ' ', strs ];
    elseif isvarname( nm ),
        error( [ 'Invalid matrix structure modifier: %s\n', ...
               'Trying to declare multiple variables? Use the VARIABLES keyword instead.' ], strs );
    else
        error( 'Invalid matrix structure modifier: %s', strs );
    end
    modifiers = [ modifiers, ' ', strs ];
end
if ~isempty( varargin ),
    str = cvx_create_structure( x.size, varargin{filt} );
    if isempty( str ),
        error( 'Incompatible structure modifiers:%s', modifiers );
    end
else
    str = [];
end
if isgeo && islin,
    error( 'GEOMETRIC and LINEAR keywords cannot be used simultaneously.' );
end
if isepi && ishypo,
    error( 'EPIGRAPH and HYPOGRAPH keywords cannot be used simultaneously.' );
end
geo = isgeo | ( ~islin & cvx___.problems( p ).gp );
v = newvar( prob, x.name, x.size, str, geo );
if isepi || ishypo,
    if geo, vv = log( v ); else vv = v; end
    if isepi, dir = 'epigraph'; else dir = 'hypograph'; end
    cvx___.problems( p ).objective = vv;
    cvx___.problems( p ).direction = dir;
    cvx___.problems( p ).geometric = geo;
end
if nargout > 0,
    varargout{1} = v;
else
    assignin( 'caller', x.name, v );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
