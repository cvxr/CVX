function y = prod( varargin )

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

persistent P
if isempty( P ),
	P.map      = cvx_remap( { 'constant' ; 'l_convex' ; 'l_concave' } );
	P.funcs    = { @prod_1, @prod_2, @prod_2 };
	P.zero     = 1;
	P.reduce   = true;
	P.reverse  = true;
	P.constant = 1;
	P.dimarg   = 2;
	P.name     = 'prod';
end
y = cvx_reduce_op( P, varargin{:} );

function x = prod_1( x )
x = builtin( 'prod', x, 2 );

function x = prod_2( x )
x = exp( sum( log( x ), 2 ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
