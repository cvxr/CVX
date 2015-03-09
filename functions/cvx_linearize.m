function y = cvx_linearize( x, do_abs ) %#ok

% CVX_LINEARIZE    Linearize.
% CVX_LINEARIZE is used primarily within CVX functions to efficiently
% implement certain monotonic functions.
% For Y = CVX_LINEARIZE(X):
%    For real affine X, Y = X. 
%    For convex X, Y is linear, and satisfies Y >= X.
%    For concave X, Y is linear, and satisfies Y <= X.
% For Y = CVX_LINEARIZE(X,1):
%    For affine X, Y = X. 
%    For complex affine X, Y is linear, and satisfies Y >= abs(X).
%    For nonnegative convex X, Y is linear, and satisfies Y >= X.
%    For nonpositive concave X, Y is linear, and satisfies Y >= -X.
% Under normal usage, CVX users should never have to use this function.

persistent Pstd Pabs Pcplx
if nargin < 2, do_abs = 'std'; end
switch do_abs,
    case 'std',
        if isempty(Pstd),
            Pstd.map = cvx_remap( { 'r_affine' }, { 'convex' }, ...
                { 'concave' }, [2,3,4] );
            Pstd.funcs = { [], @lin_aff, @lin_cvx, @lin_ccv };
        end
        P = Pstd;
    case 'abs',
        if isempty(Pabs),
            Pabs.map = cvx_remap( { 'complex' }, { 'r_affine' }, ...
                { 'c_affine', 'p_convex', 'n_concave' } );
            Pabs.funcs = { @lin_abs, @lin_aff, @lin_absn };
        end
        P = Pabs;
    case 'cplx',
        if isempty(Pcplx),
            Pcplx.map = cvx_remap( { 'affine' }, { 'convex' }, ...
                { 'concave' }, [2,3,4] );
            Pcplx.funcs = { [], @lin_aff, @lin_cvx, @lin_ccv };
        end
        P = Pcplx;
    otherwise,
        error( 'Invalid mode: %s', do_abs );
end
y = cvx_unary_op( P, x );

function y = lin_aff( x )
y = x;

function y = lin_abs( x )
y = abs(x);	

function y = lin_cvx( x ) %#ok
cvx_begin set
    variable y(size(x))
    x <= y; %#ok
cvx_end

function y = lin_ccv( x ) %#ok
cvx_begin set
    variable y(size(x))
    x >= y; %#ok
cvx_end

function y = lin_absn( x ) %#ok
cvx_begin set
    variable y(size(x))
    abs(x) <= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
