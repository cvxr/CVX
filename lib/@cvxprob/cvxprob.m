function z = cvxprob( x )

global cvx___

switch nargin,

case 0,

    name = '';
    st = dbstack;
    isver7 = isfield( st, 'file' );
    for k = 2 : length( st ),
        name = st(k).name;
        if ~isver7,
            if ispc, ds = '\'; else, ds = '/'; end
            name( 1 : max( find( name == ds ) ) ) = [];
            tt = find( name == '(' );
            if ~isempty( tt ),
                name = name( tt + 1 : end - 1 );
            else,
                tt = find( name == '.' );
                name = name( 1 : tt( 1 ) - 1 );
            end
        end
        if any( strcmp( name, { 'cvx_begin', 'cvx_begin_set', 'cvx_create_problem' } ) ),
            name = '';
        else,
            break;
        end
    end

    if ~isvarname( name ), name = 'cvx_'; end
    z = class( struct( 'dummy', [] ), 'cvxprob', cvxobj( length( cvx___.problems ) + 1 ) );
    temp = struct( ...
        'name',          name,   ...
        'complete',      true,   ...
        'sdp',           false,  ...
        'locked',        false,  ...
        'variables',     [],     ...
        'duals',         [],     ...
        'direction',     '',     ...
        'objective',     cvx( z, [ 0, 1 ], [] ), ...
        'equalities',    cvx( z, [ 0, 1 ], [] ), ...
        'cones',         struct( 'type', {}, 'indices', {} ), ...
        'reserved',      true, ...
        'vexity',        0,    ...
        'substitutions', [],   ...
        'x',             [],   ...
        'y',             [],   ...
        'result',        [],   ...
        'status',        'unsolved', ...
        'stackpos',      length( cvx___.stack ) + 1, ...
        'depth',         length( dbstack ) - 1, ...
        'self',          z );
    if isempty( cvx___.problems ),
        cvx___.problems = temp;
    else,
        cvx___.problems( end + 1 ) = temp;
    end
    cvx___.stack{ temp.stackpos } = z;

case 1,

    switch class( x ),
        case 'cvxobj',
            z = cvx___.problems( index( x ) ).self;
        case 'cvxprob',
            z = x;
        otherwise,
            if isempty( cvx___.stack ),
                error( 'There is no problem in progress.' );
            else,
                z = cvx___.stack{ end };
            end
    end
    
end        

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
