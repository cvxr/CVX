function varargout = variable( varargin )

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
try
    
cvx___.args = { varargin, 1, [] };
[ prob, name, args ] = evalin( 'caller', 'cvx_parse' );
[ p, pstr ] = verify( prob );

%
% Parse the structure
%

isepi  = false;
ishypo = false;
isnneg = false;
issemi = false;
asnneg = pstr.gp;
itype  = '';
if nargin > 1,
    try
        [ str, itypes ] = cvx_create_structure( varargin, name, args );
    catch exc
        error( exc.identifier, exc.message );
    end
    n_itypes = 0;
    for k = 1 : length( itypes ),
        strs = itypes{k};
        switch strs,
            case 'epigraph_',    isepi  = true;
            case 'hypograph_',   ishypo = true;
            case 'integer',      n_itypes = n_itypes + 1; itype = 'i_integer';
            case 'binary',       n_itypes = n_itypes + 1; itype = 'i_binary';
            case 'nonnegative',  isnneg = ~pstr.gp; asnneg = true;
            case 'semidefinite', issemi = true;
            case 'nonnegative_', asnneg = true;
        end
    end
    if isepi && ishypo,
        error( 'CVX:Structure', 'EPIGRAPH and HYPOGRAPH keywords cannot be used simultaneously.' );
    end
    if issemi && pstr.gp,
        error( 'CVX:Structure', 'SEMIDEFINITE keywords cannot be used in geometric programs.' );
    end
    if asnneg && ~isreal( str ),
        error( 'CVX:Structure', 'Complex variables cannot also be nonnegative.' );
    end
    if n_itypes,
        if pstr.gp,
            error( 'CVX:Structure', 'Integer variables cannot be used in geometric programs.' );
        elseif isepi || ishypo,
            error( 'CVX:Structure', 'Integer variables cannot be used as epigraphs or hypograph variables.' );
        elseif n_itypes > 1,
            error( 'CVX:Structure', 'At most one integer keyword may be specified.' );
        end
    end
else
    str = [];
end

%
% Create the variables
%

tx = [];
v = newvar( prob, name{1}, args{1}, str, pstr.gp );
if isepi || ishypo,
    if pstr.gp, vv = log( v ); else vv = v; end
    if isepi, dir = 'epigraph'; else dir = 'hypograph'; end
    cvx___.problems( p ).objective = vv;
    cvx___.problems( p ).direction = dir;
    cvx___.problems( p ).geometric = pstr.gp;
end
if itype,
    [ tx, dummy ] = find( cvx_basis( v ) ); %#ok
    newnonl( prob, itype, tx(:)' );
    cvx___.canslack( tx ) = false;
    if ~asnneg && ~strcmp( itype, 'i_binary' ), tx = []; end
end
if issemi,
    if xsize(1) > 1,
        isnneg = false;
        newcnstr( prob, v, 0, '>=', true );
        vv = reshape( v, xsize(1) * xsize(2), prod(xsize(3:end)) );
        vv = vv(1:xsize(1)+1:end,:);
        [ tx, dummy ] = find( cvx_basis( vv ) ); %#ok
        asnneg = false;
    else
        isnneg = true;
    end
end
if isnneg,
    newcnstr( prob, v, 0, '>=', false );
    asnneg = true;
end
if asnneg,
    [ tx, dummy ] = find( cvx_basis( v ) ); %#ok
end
if ~isempty( tx ),
    cvx_setnneg( cvx( numel(tx), sparse(tx,1:numel(tx),1) ) );
end
if nargout > 0,
    varargout{1} = v;
else
    assignin( 'caller', name{1}, v );
end

catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
