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

persistent remap,
if isempty( remap ),
    remap = remap( { 'constant' ; 'l_convex' ; 'l_concave' } );
    funcs = { @prod_1, @prod_2 };
end

try
    if nargin < 2, dim = []; end
    y = reduce_op( 'prod', funcs, remap, 1, true, false, x, dim );
catch exc
	if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
	else rethrow( exc ); end
end

function x = prod_1( x )
if x.size_(1) ~= 1,
    x = cvx( prod( cvx_constant( x ), 1 ) );
end

function x = prod_2( x )
if x.size_(1) ~= 1,
    x = exp( sum( log( x ), 1 ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
