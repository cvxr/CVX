function v = cvx_class( x, needsign, needreal )
if nargin < 3,
    if nargin < 2, needsign = false; end
    needreal = false;
end

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - concave
% 6  - real affine
% 7  - convex
% 8  - complex affine
% 9  - log concave
% 10 - log affine
% 11 - log convex monomial
% 12 - log convex posynomial
% 13 - invalid

if isempty( x ),
    v = 'empty';
    return
end
persistent remap_1 remap_2 strs
if isempty( strs ),
    remap_1 = [ 14,2,14,14,5,6,7,15,9,10,11,12,13 ];
    remap_2 = [ 16,2,16,4,5,6,7,8,9,10,11,12,13 ];
    strs = { 'negative constant', 'zero', 'positive constant', 'complex constant', ...
             'concave', 'real affine', 'convex', 'complex affine', ...
             'log-concave', 'log-affine', 'log-convex', 'log-convex', ...
             'invalid', 'constant', 'affine', 'real constant' };
end
x = cvx_classify( x );
if ~needsign,
    x = remap_1( x );
elseif ~needreal,
    x = remap_2( x );
end
v = sparse( x, 1, 1, 16, 1 ) ~= 0;
if nnz( v ) ~= v( 2 ),
    v( 2 ) = false;
end
v = strs(find( v ));
if length( v ) == 1,
    v = v{1};
else
    v = sprintf( '%s/', v{:} );
    v = [ 'mixed ', v(1:end-1) ];
end
