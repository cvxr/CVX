function cones = cvx_pushcone( cones, ctype, indices )
if isempty( indices ),
    return
end
gcone = islogical( cones );
if gcone
    global cvx___ %#ok
    cones = cvx___.cones;
end
nn = size(indices,1);
if any( strcmp( ctype, { 'nonnegative', 'integer', 'binary' } ) ),
    indices = indices(:)';
    nn = 1;
elseif nn == 1,
    ctype = 'nonnegative';
end
match = [];
if ~gcone && ~isempty( cones ),
    match = find( strcmp( { cones.type }, ctype ) );
    if nn > 1
        match = match(cellfun('size',{cones(match).indices},1)==nn);
    end
end
if ~isempty( match ),
    cones(match(1)).indices = [ cones(match).indices, indices ];
    if length(match) > 1, cones(match(2:end)) = []; end
    cvx___.cones = cones;
else
    if nn == 1,
        slacks = 1;
    else
        slacks = zeros(nn,1);
        switch ctype,
            case 'nonnegative',
                slacks = 1;
            case 'exponential', 
                slacks = [-1;0;1];
            case 'lorentz',         
                slacks = [zeros(nn-1,1);1];
            case 'rotated_lorentz',
                if nn <= 2, 
                    cones = cvx_pushcone( cones, 'nonnegative', indices(:)' );
                    if gcone, cvx___.cones = cones; end
                    return; 
                end
                slacks = [zeros(nn-2,1);1;1];
            case 'semidefinite';
                q = round(0.5*(sqrt(8*nn+1)-1));
                slacks(cumsum([1,q:-1:2]),:) = 1;
            case 'hermitian_semidefinite',
                q = round(sqrt(nn));
                slacks(cumsum([1,2*q-1:-2:2]),:) = 1;
        end
    end
    ncone = struct( 'type', ctype, 'indices', indices, 'slacks', slacks );
    if isempty( cones ), 
        cones = ncone;
    elseif nn == 1 && ~gcone,
        cones = [ ncone, cones ];
    else
        cones = [ cones, ncone ];
    end
end
if gcone,
    cvx___.cones = cones;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
