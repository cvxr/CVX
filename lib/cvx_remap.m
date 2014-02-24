function r = cvx_remap( varargin )

%CVX_REMAP   CVX expression type map generator.
%   This is an internal function used to help filter CVX expressions by type.

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - nonpositive concave
% 6  - concave
% 7  - real affine
% 8  - convex
% 9  - nonnegative convex
% 10 - complex affine
% 11 - log concave
% 12 - log affine
% 13 - log convex monomial
% 14 - log convex posynomial
% 15 - invalid

persistent remap_big remap_str
if isempty( remap_str ),
    remap_str = { ...
        'negative', 'zero', 'positive', 'complex', 'nonnegative', 'nonzero', 'nonpositive', 'real', 'constant', ...
        'non-constant', 'np-concave' 'concave', 'affine', 'convex', 'nn-convex', 'real-affine', 'complex-affine', 'non-affine', ...
        'log-concave', 'log-affine', 'log-convex', 'log-valid', 'monomial', 'posynomial', ...
        'valid', 'invalid' ...
    };
    remap_big = [ ...
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % negative
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % zero
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % positive
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % complex
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % nonnegative
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % nonzero (real)
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % nonpositive
        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % real
        1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % constant
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0; ... % non-constant
        1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ... % nonpositive concave
        1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; ... % concave
        1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0; ... % affine
        1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0; ... % convex
        0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0; ... % nonnegative convex
        1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; ... % real-affine
        0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; ... % complex-affine
        0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0; ... % non-affine
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0; ... % log-concave
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; ... % log-affine
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0; ... % log-convex
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0; ... % log-valid
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0; ... % monomial
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0; ... % posynomial
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0; ... % valid
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; ... % invalid
     ];
     remap_str = remap_str(:);
     [ remap_str, ndx ] = sort( remap_str );
     remap_big = remap_big( ndx, : );
end

[ c, ndx ] = sort( [ remap_str ; varargin(:) ] );
d = strcmp(c(1:end-1),c(2:end));
r = +any( remap_big( ndx(d), : ), 1 );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
