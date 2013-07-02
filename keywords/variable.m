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

global cvx___
prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    error( 'No CVX model exists in this scope.' );
elseif isempty( cvx___.problems ) || cvx___.problems( end ).self ~= prob,
    error( 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
end

%
% Step 1: separate the name from the parenthetical, verify the name
%

toks = regexp( nm, '^\s*([a-zA-Z]\w*)\s*(\(.*\))?\s*$', 'tokens' );
if isempty( toks ),
    error( 'Invalid variable specification: %s', nm );
end
toks = toks{1};
x.name = toks{1};
x.size = toks{2};
if x.name(end) == '_',
    error( 'Invalid variable specification: %s\n   Variables ending in underscores are reserved for internal use.', nm );
end

%
% Step 2: Parse the size. In effect, all we do here is surround what is
% replace the surrounding parentheses with square braces and evaluate. All
% that matters is the result is a valid size vector. In particular, it
% need to be a simple comma-delimited list.
%

if isempty( x.size ),
	x.size = [];
else
    try
        x.size = evalin( 'caller', [ '[', x.size(2:end-1), '];' ] );
    catch exc
        throw( MException( exc.identifier, exc.message ) );
    end
    [ temp, x.size ] = cvx_check_dimlist( x.size, true );
    if ~temp,
        error( 'Invalid variable specification: %s\n   Dimension list must be a vector of finite nonnegative integers.', nm );
    end
end

%
% Step 3. Parse the structure.
%

str = struct( 'name', {}, 'args', {}, 'full', {} );
for k = 1 : length( varargin ),
	strs = varargin{k};
	if ~isempty( strs ),
		if ~ischar( strs ) || size( strs, 1 ) ~= 1,
			error( 'Matrix structure modifiers must be strings.' );
		end
		toks = regexp( strs, '^\s*([a-zA-Z]\w*)\s*(\(.*\))?\s*$', 'tokens' );
		if isempty( toks ),
			error( 'Invalid structure specification: %s\n', strs );
		end
		toks = toks{1};
		strx.name = toks{1};
		if isempty( toks{2} ),
			strx.args = {};
		else
			strx.args = evalin( 'caller', [ '{', toks{2}(2:end-1), '}' ] );	
		end
		strx.full = strs;
		str(end+1) = strx; %#ok
	end
end

islin = false;
isgeo = false;
isepi = false;
ishypo = false;
isnneg = false;
n_itypes = 0;
itype = '';
str = [];
for k = 1 : length( varargin ),
    strs = varargin{k};
    if isempty( strs ),
        continue;
    elseif ~ischar( strs ) || size( strs, 1 ) ~= 1,
        error( 'Matrix structure modifiers must be strings.' );
    end
    switch strs,
        case 'geometric_',  isgeo  = true;
        case 'linear_',     islin  = true;
        case 'epigraph_',   isepi  = true;
        case 'hypograph_',  ishypo = true;
        case 'integer',     n_itypes = n_itypes + 1; itype = 'i_integer';
        case 'binary',      n_itypes = n_itypes + 1; itype = 'i_binary';
        case 'nonnegative', isnneg = true;
        otherwise,
            if any( strs == '(' ) && strs(end) == ')',
                xt = find( strs == '(' );
                strx.name = strs(1:xt(1)-1);
                strx.args = evalin( 'caller', [ '{', strs(xt(1)+1:end-1), '}' ] );
                strx.full = strs;
                str{end+1} = strx; %#ok
            else
                str{end+1} = strs; %#ok
            end
    end
end
if ~isempty(str),
    try
        str = cvx_create_structure( x.size, str{:} );
    catch exc
        throw( MException( exc.identifier, sprintf( '%s\nTrying to declare multiple variables? Use the VARIABLES keyword instead.', exc.message ) ) );
    end
end
if isgeo && islin,
    error( 'GEOMETRIC and LINEAR keywords cannot be used simultaneously.' );
end
if isepi && ishypo,
    error( 'EPIGRAPH and HYPOGRAPH keywords cannot be used simultaneously.' );
end
if n_itypes,
    if isgeo,
        error( 'Integer variables cannot be use be GEOMETRIC.' );
    elseif ~islin && cvx___.problems( end ).gp,
        error( 'Integer variables cannot be used in geometric programs.' );
    elseif n_itypes > 1,
        error( 'At most one integer keyword may be specified.' );
    end
end

geo = isgeo || ( ~islin && cvx___.problems( end ).gp );
v = newvar( prob, x.name, x.size, str, geo );
if isepi || ishypo,
    if geo, vv = log( v ); else vv = v; end
    if isepi, dir = 'epigraph'; else dir = 'hypograph'; end
    cvx___.problems( end ).objective = vv;
    cvx___.problems( end ).direction = dir;
    cvx___.problems( end ).geometric = geo;
end
if itype,
    [ tx, dummy ] = find( cvx_basis( v ) ); %#ok
    newnonl( prob, itype, tx(:)' );
    cvx___.canslack( tx ) = false;
end
if isnneg && ~geo,
    if ~n_itypes,
        [ tx, dummy ] = find( cvx_basis( v ) );  %#ok
        newnonl( prob, 'nonnegative', tx(:)' );
    elseif ~isequal( itype, 'i_binary' ),
        newcnstr( prob, v, 0, '>=' );
    end
end
if nargout > 0,
    varargout{1} = v;
else
    assignin( 'caller', x.name, v );
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
