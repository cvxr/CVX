function v = cvx_class( x, needsign, needreal, needzero )
if nargin < 2, needsign = false; end
if nargin < 3, needreal = false; end
if nargin < 4, needzero = needsign; end

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - negative concave
% 6  - concave
% 7  - positive concave
% 8  - negative affine
% 9  - real affine
% 10 - positive affine
% 11 - negative convex
% 12 - convex
% 13 - positive convex
% 14 - complex affine
% 15 - log concave
% 16 - log affine
% 17 - log convex (monomial)
% 18 - log convex (posynomial)
% 19 - gp monomial (log affine)
% 20 - ggp monomial (log convex)
% 21 - gp posynomial (log convex)
% 22 - invalid
% 23 - real constant (negative or positive)
% 24 - constant (real or complex)
% 25 - affine (real or complex)

if isempty( x ),
    v = 'empty';
    return
end
persistent remap_s remap_r remap_z strs
if isempty( strs ),
    strs = { 
        'negative constant', 'zero', 'positive constant', 'complex constant', ...
        'negative concave', 'concave', 'positive concave', ...
        'negative affine', 'real affine', 'positive affine', ...
        'negative convex', 'convex', 'positive convex', 'complex affine', ...
        'log concave', 'log affine', 'log convex', 'log convex', ...
        'monomial', 'posynomial', 'posynomial', ...
        'invalid', ...
        'constant', 'affine' };
    remap_s = [23, 2,23, 4, 6, 6, 6, 9, 9, 9,12,12,12,14,15,16,17,18,19,20,21,22,23,24,25];
    remap_r = [ 1, 2, 3,24, 6, 6, 6, 9,25, 9,12,12,12,25,15,16,17,18,19,20,21,12,24,24,25];
    remap_z = [ 1,24, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]; 
end
x = cvx_classify( x );
if ~needsign,
    x = remap_s( x );
end
if ~needreal,
    x = remap_r( x );
end
if ~needzero,
    x = remap_z( x );
end
v = sparse( double(x), 1, 1, numel(strs), 1 ) ~= 0;
v = strs( v );
if length( v ) == 1,
    v = v{1};
else
    v = sprintf( '%s/', v{:} );
    v = [ 'mixed ', v(1:end-1) ];
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
