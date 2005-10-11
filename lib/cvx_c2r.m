function x = cvx_c2r( x, dim )

%
% Quick exit for real quantities
%

if isreal( x ), return; end

%
% Determine expansion dimension
%

sx = size( x );
if nargin < 2,
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim <= 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension.' );
end
sx = [ sx, ones( 1, dim - length( sx ) ) ];
nd = length( sx );

%
% Permute if necessary
%

perm = [];
if any( sx( 1 : dim - 1 ) ~= 1 ),
    perm = [ dim, 1 : dim - 1, dim + 1 : nd ];
    x = permute( x, perm );
    sx = sx( perm );
    dim = 1;
end

%
% Perform expansion
%

x = x( : ).';
sx( dim ) = 2 * sx( dim );
x = reshape( [ real( x ) ; imag( x ) ], sx );

%
% Reverse permute if necessary
%

if ~isempty( perm ),
    x = ipermute( x, perm );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

