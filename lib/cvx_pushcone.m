function cvx_pushcone( ctype, indices, ndims )
global cvx___
isnneg = isequal( ctype, 'nonnegative' );
if isa( indices, 'cvx' ),
    sx = size( indices );
    [ indices, cx ] = find( cvx_basis( indices ) ); %#ok
    if ~isnneg
        if nargin < 3, ndims = 1; end
        indices = reshape( indices, [], prod(sx(ndims+1:end)) );
    end
end
cones = cvx___.cones;
if isnneg,
    indices = indices(:)';
    if ~isempty(cones) && isequal( cones(1).type, 'nonnegative' ),
        cones(1).indices = [ cones(1).indices, indices ];
    else
        ncone = struct( 'type', ctype, 'indices', indices );
        if isempty( cones )
            cones = ncone;
        else
            cones = [ ncone, cones ];
        end
    end
elseif isempty( cones ),
    cones = struct( 'type', ctype, 'indices', indices );
else
    match = find( strcmp( { cones.type }, ctype ) );
    nlsiz = size( indices, 1 );
    match = match( cellfun( 'size', { cones(match).indices }, 1 ) == nlsiz );
    if isempty( match ),
        cones = [ cones, struct( 'type', ctype, 'indices', indices ) ];
    else
        match = match(1);
        cones(match).indices = [ cones(match).indices, indices ];
    end
end
cvx___.cones = cones;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
