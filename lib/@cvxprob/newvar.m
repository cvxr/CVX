function y = newvar( prob, name, siz, str, geo )

global cvx___
[ p, pstr ] = verify( prob );

%
% Check name
%

if isempty( name ),
    nstr = [];
elseif ischar( name ),
    nstr = struct( 'type', '.', 'subs', name );
else
    nstr = name;
    name = cvx_subs2str( name, [ 1, 0, 1 ], 1 );
    name(1) = [];
end

%
% Retrieve an existing variable, and check for conflicts
%

vars = pstr.variables;
if ~isempty( nstr ),
    try
        y = builtin( 'subsref', vars, nstr );
    catch
        y = [];
    end
    if nargin == 2,
        if ~isempty( y ), return; end
        error( 'CVX:Variable', [ 'Unknown variable: ', name ] );
    elseif ~isempty( y ),
        error( 'CVX:Variable', [ 'Duplicate variable name: ', name ] );
    end
elseif nargin == 2,
    error( 'CVX:Variable', 'Second argument must be a non-empty string or subscript structure array.' );
else
    y = []; %#ok
end

%
% Check for conflict with dual variable
%

if ~isempty( nstr ) && isfield( pstr.duals, nstr(1).subs ),
    error( 'CVX:Variable', 'Primal/dual variable name conflict: %s', nstr(1).subs );
end

%
% Quick exit for retrieval mode
%

if nargin == 2,
    y = vars;
    try
        y = subsref( y, name );
        return
    catch
        error( 'CVX:Variable', [ 'Unknown variable: ', nstr ] );
    end
end

%
% Create the variable
%

if isa( siz, 'cvx' ),

    %
    % Out of an existing object
    %

    y = newsubst( prob, siz );

else

    %
    % Check structure
    %

    siz = cvx_get_dimlist( { siz }, 'default', [ 1, 1 ] );
    len = prod( siz );
    if nargin < 4 || isempty( str ),
        dof = len;
        str = [];
    elseif ~isnumeric( str ) || ndims( str ) > 2 || size( str, 2 ) ~= len, %#ok
        error( 'CVX:Variable', 'Fourth argument must be a valid structure matrix.' );
    elseif nnz( str ) == 0,
        error( 'CVX:Variable', 'Structure matrix cannot be identically zero.' );
    else
        temp = any( str, 2 );
        dof = full( sum( temp ) );
        if dof ~= length( temp ),
            str = str( temp, : );
        end
    end

    %
    % Geometric flag
    %

    if nargin < 5 || isempty( geo ),
        geo = false;
    elseif ~isnumeric( geo ) && ~islogical( geo ) || length( geo ) ~= 1,
        error( 'CVX:Variable', 'Fifth argument must be true or false.' );
    end

    %
    % Allocate the raw variable data
    %

    geo  = any( geo( : ) );
    omax = length( cvx___.classes );
    nmax = omax + dof;
    ndim = omax + 1 : nmax;
    cvx___.classes ( ndim, 1 ) = int8(9);
    cvx___.canslack( ndim, 1 ) = true;
    cvx___.readonly( ndim, 1 ) = p;
    if geo,
        gdim = nmax + 1 : nmax + dof;
        cvx___.classes( gdim, 1 ) = int8(16);
        cvx___.canslack( gdim, 1 ) = true;
        cvx___.readonly( gdim, 1 ) = p;
        cvx___.logarithm( gdim, 1 ) = ndim';
        cvx___.exponential( ndim, 1 ) = gdim';
        cvx___.exponential( nmax + dof, 1 ) = 0;
        ndim = gdim;
    end
    cvx___.x = [];
    cvx___.y = [];

    %
    % Create the variable object
    %

    if dof == 0,
        str2 = sparse( 1, 0 );
    else
        str2 = sparse( ndim, 1 : dof, 1 );
        if ~isempty( str ),
            str2 = str2 * str;
        end
    end
    y = cvx( siz, str2 );

end

%
% If the variable is named, save it in the problem structure
%

if ~isempty( nstr ),
    try
        cvx___.problems( p ).variables = builtin( 'subsasgn', vars, nstr, y );
    catch
        error( 'CVX:Variable', [ 'Invalid variable name: ', name ] );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
