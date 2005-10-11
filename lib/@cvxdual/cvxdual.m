function y = cvxdual( prob, name )
error( nargchk( 2, 2, nargin ) );

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end
global cvx___
p = index( prob );

%
% Check name
%

if ~isvarname( name ),
    error( 'Second arugment must be a valid MATLAB variable name.' );
elseif isfield( cvx___.problems( p ).duals, name ),
    y = subsref( cvx___.problems( p ).duals, struct( 'type', '.', 'subs', name ) );
    return
elseif isfield( cvx___.problems( p ).variables, name ),
    error( [ 'Primal/dual variable name conflict: "', name, '"' ] );
end

%
% Create the object
%

y = class( struct( 'name_', name ), 'cvxdual', cvxobj( prob ) );
cvx___.problems( p ).duals = builtin( 'subsasgn', cvx___.problems( p ).duals, struct( 'type', '.', 'subs', name ), cvx( prob, [0,0], [] ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
