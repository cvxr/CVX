function cones = cvx_pushcone( cones, ctype, indices, ndims )
sx = size( indices );
[ indices, cx ] = find( cvx_basis( indices ) ); %#ok
if nargin < 4, ndims = 1; end
indices = reshape( indices, [], prod(sx(ndims+1:end)) );
cones = cvx_pushcone( cones, ctype, indices );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
