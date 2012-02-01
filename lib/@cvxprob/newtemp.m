function z = newtemp( prob, model, base )
global cvx___
error( nargchk( 1, 3, nargin ) );

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end

%
% Check model
%

if nargin < 2 || isempty( model ),
    ismod = 0;
    siz   = [ 1, 1 ];
    str   = sparse( 1, 1, 1 );
elseif isa( model, 'cvx' ),
    ismod = 1;
    siz = size( model );
    str = bcompress( model );
elseif cvx_check_dimlist( model ),
    ismod = 0;
    siz   = model( : ).';
    str   = 1 : prod( siz );
    str   = sparse( str, str, 1 );
else
    error( 'Second argument must be a cvx object or a non-empty size vector.' );
end

%
% Check base
%

if nargin < 3 || isempty( base ),
    base( 1 ).type = '.';
    base( 1 ).subs = 'temp_';
    base( 2 ).type = '{}';
    try
        ndx = length(cvx___.problems( p ).variables.temp_);
    catch
        ndx = 0;
    end
    base( 2 ).subs = { ndx + 1 };
end

%
% Create temporary variable
%

len = length( cvx___.vexity );
z = newvar( prob, base, siz, str );

%
% Set vexity
%

if ismod,
    model = cvx_vexity( model );
    cvx___.vexity( len + 1 : end ) = sign( str * model( : ) );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
