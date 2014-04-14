function y = cumprod( varargin )

%   Disciplined geometric programming information for PROD:
%      CUMPROD(X,DIM) is a vectorized version of multiplication,
%      so in most cases it would be incompatible with DCPs. Therefore it
%      has not been implemented to support DCPs. DGPs however support
%      products more liberally. When PROD is used in a DCP, elements in
%      each subvector must satisfy the corresponding combination rule
%      for multiplication (see TIMES). For example, suppose that X looks
%      like this:
%         X = [ log-convex log-concave log-affine  ;
%               log-affine log-concave log-concave ]
%      Then CUMPROD(X,1) would be legal, but CUMPROD(X,2) would not, 
%      because the top row contains the product of log-convex and 
%      log-concave terms, in violation of the DGP ruleset.

persistent params
if isempty( params ),
	params.map      = cvx_remap( { 'constant' ; 'l_convex' ; 'l_concave' } );
	params.funcs    = { @cumprod_1, @cumprod_2, @cumprod_2 };
	params.constant = 1;
	params.zero     = 1;
	params.reduce   = false;
	params.reverse  = true;
	params.dimarg   = 2;
	params.name     = 'cumprod';
end

try
    y = reduce_op( params, varargin );
catch exc
	if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
	else rethrow( exc ); end
end

function x = cumprod_1( x )
x = builtin( 'prod', x, 1 );

function x = cumprod_2( x )
s = x.size_;
if s(2) ~= 1,
	x = log( x );
    b = reshape( x.basis_, [], s(2) );
    b = cumsum( b, 2 );
    x = exp( cvx( s, reshape( b, [], prod(s) ) ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
