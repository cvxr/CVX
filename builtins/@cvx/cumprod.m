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

persistent P
if isempty( P ),
	P.map      = cvx_remap( { 'constant' ; 'l_convex' ; 'l_concave' } );
	P.funcs    = { @cumprod_1, @cumprod_2, @cumprod_2 };
	P.constant = 1;
	P.zero     = 1;
	P.reduce   = false;
	P.reverse  = true;
	P.dimarg   = 2;
	P.name     = 'cumprod';
end
y = reduce_op( P, varargin );

function x = cumprod_1( x )
x = builtin( 'cumprod', x, 2 );

function x = cumprod_2( x )
x = exp( cumsum( log( x ), 2 ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
