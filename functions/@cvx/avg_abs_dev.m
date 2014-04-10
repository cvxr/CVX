function y = avg_abs_dev( varargin )

%AVG_ABS_DEV    Internal cvx version.

persistent remap funcs
if isempty( remap ),
	params.map = cvx_remap( 'affine' );
	params.funcs = { @avg_abs_dev_1 };
	params.zero = NaN;
	params.reduce = true;
	params.reverse = false;
	params.dimarg = 2;
	params.name = 'avg_abs_dev';
end

try
    if nargin < 2, dim = []; end
    y = reduce_op( params, varargin{:} );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function y = avg_abs_dev_1( x )
[nx,nv] = size(x); %#ok
% In theory we could just say y = mean(abs(x-mean(x))). However, by
% adding an extra variable we preserve sparsity.
cvx_begin
    variable y( 1, nv );
    y == sum( x ) / nx; %#ok
    minimize( sum( abs( x - repmat(y,[nx,1]) ) ) / nx );
cvx_end
y = cvx_optval;

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
