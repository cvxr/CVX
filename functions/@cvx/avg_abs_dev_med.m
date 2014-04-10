function y = avg_abs_dev_med( varargin )

%AVG_ABS_DEV_MED    Internal cvx version.

persistent params
if isempty( params ),
	params.map = cvx_remap( 'affine' );
	params.funcs = { @avg_abs_dev_med_1 };
	params.zero = NaN;
	params.reduce = true;
	params.reverse = false;
	params.dimarg = 2;
	params.name = 'avg_abs_dev_med';
end

try
    y = reduce_op( params, varargin{:} );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function y = avg_abs_dev_med_1( x )
y = [];
[nx,nv] = size(x); %#ok
cvx_begin
    variable y( 1, nv );
    minimize( sum( abs( x - repmat(y,[nx,1]) ) ) / nx );
cvx_end
y = cvx_optval;

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
