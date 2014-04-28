function y = avg_abs_dev( varargin )

%AVG_ABS_DEV   Average absolute deviation (about the mean).
%   For vectors, AVG_ABS_DEV(X) is the average absolute deviation of X
%   about its mean; that is, AVG_ABS_DEV(X)=MEAN(ABS(X-MEAN(X))). For
%   matrices, AVG_ABS_DEV(X) is a row vector containing the average
%   absolute deviation of each column. For N-D arrays, AVG_ABS_DEV(X)
%   is the average absolute deviation of the elements along the first
%   non-singleton dimension of X.
%
%   AVG_ABS_DEV(X,DIM) performs the computation along the dimension DIM. 
%
%   See also AVG_ABS_DEV_MED.
%
%   Disciplined convex programming information:
%       AVG_ABS_DEV is convex and nonmontonic in X. 
%       Therefore, X must be affine.
%      

persistent params
if isempty( params ),
	params.map = cvx_remap( { 'real' ; 'r_affine' } );
    params.funcs = { @avg_abs_dev_cnst, @avg_abs_dev_aff };
    params.constant = 1;
	params.zero = NaN;
	params.reduce = true;
	params.reverse = false;
	params.name = 'avg_abs_dev';
	params.dimarg = 2;
end

try
    y = cvx_reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = avg_abs_dev_cnst( x )
y = mean( abs( bsxfun( @minus, x, mean( x, 1 ) ) ), 1 );

function y = avg_abs_dev_aff( x )
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