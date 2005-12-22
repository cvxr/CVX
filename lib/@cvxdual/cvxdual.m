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
% Create the object
%

y = class( struct( 'name_', name ), 'cvxdual', cvxobj( prob ) );
        
% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
