function y = cumsum( varargin )

%Disciplined convex/geometric programming information for SUM:
%   CUMSUM(X) and CUMSUM(X,DIM) are vectorized forms of addition. So 
%   when CUMSUM is used in a DCP or DGP, elements in each subvector 
%   must satisfy the corresponding combination rules for addition (see
%   PLUS). For example, suppose that X looks like this:
%      X = [ convex concave affine  ;
%            affine concave concave ]
%   Then CUMSUM(X,1) would be permittted, but CUMSUM(X,2) would not, 
%   because the top row contains the sum of convex and concave terms, in
%   violation of the DCP ruleset. For DGPs, addition rules dictate that
%   the elements of X must be log-convex or log-affine.

persistent P
if isempty( P ),
    P.map      = cvx_remap( { 'constant' ; 'affine' ; 'convex' ; 'concave' } );
    P.funcs    = { @cumsum_1, @cumsum_2, @cumsum_2 };
    P.zero     = 0;
    P.reduce   = false;
    P.reverse  = true;
    P.dimarg   = 2;
    P.constant = 1;
    P.name     = 'cumsum';
end
y = cvx_reduce_op( P, varargin{:} );

function x = cumsum_1( x )
x = builtin( 'cumsum', x, 2 );

function x = cumsum_2( x )
s = x.size_;
if s(2) ~= 1,
    b = reshape( x.basis_, [], s(2) );
    b = cumsum( b, 2 );
    x = cvx( s, reshape( b, [], prod(s) ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
