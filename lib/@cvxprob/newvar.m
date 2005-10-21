function y = newvar( prob, name, siz, str )
error( nargchk( 2, 4, nargin ) );

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end
p = index( prob );
global cvx___

%
% Check name
%

if ischar( name ),
    if ~isempty( name ),
        if size( name, 1 ) ~= 1,
            error( 'Second argument must be a string or a subscript structure array.' );
        elseif ~isvarname( name ),
            error( sprintf( 'Invalid variable name: %s', name ) );
        elseif isfield( cvx___.problems( p ).variables, name ),
            error( sprintf( 'Variable already exists: %s', name ) );
        end
        nstr = struct( 'type', '.', 'subs', name );
    else,
        nstr = [];
    end
elseif ~isstruct( name ),
    error( 'Second argument must be a string or a subscript structure array.' );
else
    nstr = name;
    name = cvx_subs2str( name, [ 1, 0, 1 ], 1 );
    name(1) = [];
    if ~isequal( nstr(1).type, '.' ),
        error( 'Invalid subscript structure: first element must be a field reference.' );
    end
end

%
% Quick exit for retrieval mode
%

if nargin == 2,
    y = cvx___.problems( p ).variables;
    try
        y = subsref( y, name );
        return
    catch,
        error( [ 'Unknown variable: ', nstr ] );
    end
end

%
% Quick exit for creating variables out of existing data
%

if isa( siz, 'cvx' ) & problem( siz ) == prob,
    cvx___.problems( p ).variables = builtin( 'subsasgn', cvx___.problems( p ).variables, nstr, siz );
    y = siz;
    return
end

%
% Creating variables from new data
%

vars = cvx___.problems( p ).variables;
if ~isempty( nstr ),
    if isfield( cvx___.problems( p ).duals, nstr(1).subs ),
        error( sprintf( 'Primal/dual variable name conflict: %s', nstr(1).subs ) );
    elseif isfield( vars, nstr(1).subs ) & isa( subsref( vars, nstr(1) ), 'cvx' ),
        error( [ 'Variable name conflict: ', nstr(1).subs ] );
    end
end

%
% Check size
%

[ temp, siz ] = cvx_check_dimlist( siz, false );
if ~temp,
    if cvx_check_dimlist( siz, true ),
        error( 'Invalid size vector (cannot be empty).' );
    else,
        error( 'Invalid size vector.' );
    end
end

%
% Check structure
%

len = prod( siz );
if nargin < 4 | isempty( str ),
    dof = len;
    str = [];
elseif ~isnumeric( str ) | ndims( str ) > 2 | size( str, 1 ) ~= len,
    error( 'Fourth argument must be a valid structure matrix.' );
elseif nnz( str ) == 0,
    error( 'Structure matrix cannot be identically zero.' );
else,
    temp = any( str, 1 );
    dof = full( sum( temp ) );
    if dof ~= length( temp ),
        str = str( :, temp );
    end
end

%
% Add the variable to the problem
%

dims = length( cvx___.problems( p ).reserved );
dims = dims + 1 : dims + dof;
str2 = sparse( 1 : dof, dims, 1, dof, dims(end) );
if ~isempty( str ), str2 = str * str2; end
y = cvx( prob, siz, str2, dof );
if ~isempty( nstr ),
    cvx___.problems( p ).variables = builtin( 'subsasgn', vars, nstr, y );
end
cvx___.problems( p ).reserved( dims, 1 ) = 0;
cvx___.problems( p ).vexity(   dims, 1 ) = 0;
cvx___.problems( p ).x = [];
cvx___.problems( p ).y = [];

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
