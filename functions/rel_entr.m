function z = rel_entr( x, y )

%REL_ENTR   Scalar relative entropy.
%   REL_ENTR(X,Y) returns an array of the same size as X+Y with the 
%   relative entropy function applied to each element:
%                      { X.*LOG(X./Y) if X >  0 & Y >  0,
%      REL_ENTR(X,Y) = { 0            if X == 0 & Y >= 0,
%                      { +Inf         otherwise.
%   X and Y must either be the same size, or one must be a scalar. If X and
%   Y are vectors, then SUM(REL_ENTR(X,Y)) returns their relative entropy.
%   If they are PDFs (that is, if X>=0, Y>=0, SUM(X)==1, SUM(Y)==1) then
%   this is equal to their Kullback-Liebler divergence SUM(KL_DIV(X,Y)).
%   -SUM(REL_ENTR(X,1)) returns the entropy of X.
%
%   Disciplined convex programming information:
%       REL_ENTR(X,Y) is convex in both X and Y, nonmonotonic in X, and
%       nonincreasing in Y. Thus when used in CVX expressions, X must be
%       real and affine and Y must be concave. The use of REL_ENTR(X,Y) in
%       an objective or constraint will effectively constrain both X and Y 
%       to be nonnegative, hence there is no need to add additional
%       constraints X >= 0 or Y >= 0 to enforce this.

cvx_expert_check( 'rel_entr', x, y );

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { { 'negative' }, { 'any' } }, ...
        { { 'any' }, { 'nonpositive' } }, ...
        { { 'real' } }, ...
        { { 'r_affine' }, { 'concave' } }, ...
        [0,0,1,2] );
    P.funcs = { @rel_entr_1, @rel_entr_2 };
    P.constant = 1;
    P.name = 'rel_entr';
    P.test = @rel_entr_test;
end
z = cvx_binary_op( P, x, y );

function z = rel_entr_1( x, y )   
z  = x .* builtin( 'log', x ./ y );

function z = rel_entr_2( x, y ) %#ok
sz = max( numel(y), numel(x) );
cvx_begin
    epigraph variable z( sz )
    { -z, x, y } == exponential( sz ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
