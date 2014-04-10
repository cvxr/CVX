function y = sum( varargin )

%   Disciplined convex/geometric programming information for SUM:
%      SUM(X) and SUM(X,DIM) is a vectorized version of addition. So 
%      when SUM is used in a DCP or DGP, elements in each subvector must
%      satisfy the corresponding combination rules for addition (see 
%      PLUS). For example, suppose that X looks like this:
%         X = [ convex concave affine  ;
%               affine concave concave ]
%      Then SUM(X,1) would be permittted, but SUM(X,2) would not, 
%      because the top row contains the sum of convex and concave terms,
%      in violation of the DCP ruleset. For DGPs, addition rules dictate
%      that the elements of X must be log-convex or log-affine.

persistent params
if isempty( params ),
    params.map     = cvx_remap( { 'affine' ; 'convex' ; 'concave' } );
    params.funcs   = { @sum_1, @sum_1, @sum_1 };
    params.zero    = 0;
    params.reduce  = true;
    params.reverse = true;
    params.dimarg  = 2;
    params.name    = 'sum';
end

try
    y = reduce_op( params, varargin{:} );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function x = sum_1( x )
s = x.size_;
if s(1) > 1,
    b = x.basis_;
    [m,n] = size(b);
    b = reshape( b, m*s(1), s(2) );
    b = sum( b, 2 );
    b = reshape( b, m, s(1) );
    x = cvx( s, b );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
