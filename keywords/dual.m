function dual( varargin )

if ~evalin( 'caller', 'exist(''cvx_problem'',''var'')', '0' ),
    error( 'A cvx problem does not exist in this scope.' );
else,
    prob = evalin( 'caller', 'cvx_problem' );
    p = index( prob );
end

if any( strcmpi( varargin{1}, { 'variable', 'variables' } ) ),
    varargin(1) = [];
end

for k = 1 : length(varargin),
    name = varargin{k};
    if ~isvarname( name ),
        error( [ 'Invalid dual variable specification: ', name ] );
    end
    assignin( 'caller', name, newdual( prob, name ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
