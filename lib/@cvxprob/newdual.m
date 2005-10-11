function z = newdual( prob, name )
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
    error( [ 'Second argument must be a valid MATLAB variable name.' ] );
elseif isfield( cvx___.problems( p ).variables, name ),
    error( [ 'Dual variable duplicates a primal variable name: ', name ] );
elseif isfield( cvx___.problems( p ).duals, name ),
    error( [ 'Duplicate dual variable name: ', name ] );
end

%
% Creating the dual variable structure
%

z = cvxdual( prob, name );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
