function maximize( varargin )

%MAXIMIZE Specifiies a concave (or affine) objective to be maximized.

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isempty( cvx_problem ) || ~isa( cvx_problem, 'cvxprob' ),
    error( 'A cvx problem does not exist in this scope.' );
end
if nargin < 1,
    error( 'Objective expression missing.' );
elseif iscellstr( varargin ),
	x = sprintf( '%s,', varargin{:} );
    x = evalin( 'caller', x(1:end-1) );
elseif nargin > 1,
    error( 'Too many input arguments.' );
else
    x = varargin{1};
end
try
    newobj( cvx_problem, 'maximize', x );
catch exc
    rethrow( exc )
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
