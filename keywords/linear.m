function linear( varargin )
error( nargchk( 1, Inf, nargin ) );

%LINEAR Declare multiple cvx variables.
%   LINEAR VARIABLE x
%   LINEAR VARIABLES x1 x2 x3 ...
%   LINEAR VARIABLE x(n1,n2,...,nk) ...
%   The LINEAR keyword is used to declare linear variables in
%   a CVX model that has been declared as a geometric program
%   using the CVX_BEGIN GP command. In SDPs and normal CVX
%   models, the LINEAR keyword has no practical effect.
%
%   See also GEOMETRIC, VARIABLE, VARIABLES, DUAL, DUALS.

if ~iscellstr( varargin ),
    error( 'LINEAR must be used in command mode.' );
end

if nargin < 1,
    k1 = 1;
elseif isequal( varargin{1}, 'variable' ),
    evalin( 'caller', sprintf( '%s ', 'variable', varargin{2:end}, ' linear_' ) );
elseif isequal( varargin{1}, 'variables' ),
    for k = 2 : nargin,
        evalin( 'caller', [ 'variable ', varargin{k}, ' linear_' ] );
    end
else
    for k = 1 : nargin,
        evalin( 'caller', [ 'variable ', varargin{k}, ' linear_' ] );
    end
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
