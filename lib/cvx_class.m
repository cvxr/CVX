function v = cvx_class( x, needsign, needreal, needzero )
if nargin < 2, needsign = false; end
if nargin < 3, needreal = false; end
if nargin < 4, needzero = needsign; end

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
% ---
% 16 - constant
% 17 - affine
% 18 - real constant

if isempty( x ),
    v = 'empty';
    return
end
persistent remap_s remap_r remap_z strs
if isempty( strs ),
    remap_s = [18,2,18,4,6,6,7,8,8,10,11,12,13,14,15];
    remap_r = [1,2,3,16,5,6,17,8,9,17,11,12,13,14,15,16,17,16];
    remap_z = [1,16,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
    strs = { 'negative constant', 'zero', 'positive constant', 'complex constant', ...
             'nonpositive concave', 'concave', 'real affine', 'convex', 'nonnegative convex', 'complex affine', ...
             'log-concave', 'log-affine', 'log-convex', 'log-convex', ...
             'invalid', 'constant', 'affine', 'real constant' };
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
v = sparse( x, 1, 1, 18, 1 ) ~= 0;
if nnz( v ) ~= v( 2 ),
    v( 2 ) = false;
end
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
