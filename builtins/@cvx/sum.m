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

persistent P
if isempty( P ),
    P.map      = {};
    P.funcs    = { @sum_1 };
    P.reduce   = true;
    P.reverse  = true;
    P.name     = 'sum';
    P.constant = [];
    P.dimarg   = 2;
end
y = cvx_reduce_op( P, varargin{:} );

function x = sum_1( x )
s = x.size_;
b = reshape( x.basis_, [], s(2) );
b = sum( b, 2 );
x = cvx( s(1), reshape( b, [], s(1) ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
