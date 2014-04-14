function y = prod( x, dim )

%   Disciplined geometric programming information for PROD:
%      PROD(X) and PROD(X,DIM) are vectorized versions of multiplication
%      so in most cases it would be incompatible with DCPs. Therefore it
%      has not been implemented to support DCPs. DGPs however support
%      products more liberally. When PROD is used in a DCP, elements in
%      each subvector must satisfy the corresponding combination rule
%      for multiplication (see TIMES). For example, suppose that X looks
%      like this:
%         X = [ log-convex log-concave log-affine  ;
%               log-affine log-concave log-concave ]
%      Then PROD(X,1) would be permittted, but PROD(X,2) would not, 
%      because the top row contains the product of log-convex and 
%      log-concave terms, in violation of the DGP ruleset.

persistent params
if isempty( params ),
	params.map      = cvx_remap( { 'constant' ; 'l_convex' ; 'l_concave' } );
	params.funcs    = { @prod_1, @prod_2, @prod_2 };
	params.zero     = 1;
	params.reduce   = true;
	params.reverse  = true;
	params.constant = 1;
	params.dimarg   = 2;
	params.name     = 'prod';
end

try
    y = cvx_reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function x = prod_1( x )
x = prod( x, 2 );

function x = prod_2( x )
s = x.size_;
if s(2) ~= 1,
	x = log( x );
    b = reshape( x.basis_, [], s(2) );
    b = sum( b, 2 );
    x = exp( cvx( s, reshape( b, [], s(1) ) ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
