function newnonl( prob, ncones, arg )
global cvx___
verify( prob );
cones = cvx___.cones;
if nargin == 3,
    ncones = struct( 'type', ncones, 'indices', arg );
end
for k = 1 : length( ncones ),
    ncone = ncones(k);
    if isequal( ncone.type, 'nonnegative' ),
        ncone.indices = ncone.indices(:)';
    end
    if isempty( cones ),
        cones = ncone;
    else
        match = find( strcmp( { cones.type }, ncone.type ) );
        if ~isempty( match ),
            nlsiz = size( ncone.indices, 1 );
            match = match( cellfun( 'size', { cones(match).indices }, 1 ) == nlsiz );
            if isempty( match ),
                cones = [ cones, ncone ]; %#ok
            else
                match = match(1);
                cones(match).indices = [ cones(match).indices, ncone.indices ];
            end
        else
            cones = [ cones, ncone ]; %#ok
        end
    end
end
cvx___.cones = cones;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
