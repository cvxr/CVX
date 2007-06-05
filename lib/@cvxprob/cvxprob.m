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
            if ispc, ds = '\'; else ds = '/'; end
            name( 1 : max( find( name == ds ) ) ) = [];
            tt = find( name == '(' );
            if ~isempty( tt ),
                name = name( tt + 1 : end - 1 );
            else
                tt = find( name == '.' );
                name = name( 1 : tt( 1 ) - 1 );
            end
        end
        switch k,
            case 2,
                if ~strcmp( name, 'cvx_create_problem' ),
                    error( 'Creating a cvxprob object manually is not permitted.' );
                else
                    name = '';
                end
            case 3,
                if ~any( strcmp( name, { 'cvx_begin', 'cvx_begin_set' } ) ),
                    error( 'Calling cvx_create_problem manually is not permitted.' );
                else
                    name = '';
                end
            otherwise,
                break;
        end
    end

    if ~isempty( cvx___.problems ),
        ndx = [ cvx___.problems.depth ];
        ndx = min( find( ndx >= length( st ) - 3 ) );
        if ~isempty( ndx ),
            pop( cvx___.problems( ndx ).self, 'reset' );
        end
    end

    if ~isempty( cvx___.problems ),
        nprec  = cvx___.problems( end ).precision;
        ngprec = cvx___.problems( end ).gptol;
        nrprec = cvx___.problems( end ).rat_growth;
        nsolv  = cvx___.problems( end ).solver;
    else
        nprec  = cvx___.precision;
        ngprec = cvx___.gptol;
        nrprec = cvx___.rat_growth;
        nsolv  = cvx___.solver;
    end

    if ~isvarname( name ), name = 'cvx_'; end
    z = class( struct( 'index_', length( cvx___.problems ) + 1 ), 'cvxprob', cvxobj );
    nres = length( cvx___.reserved );
    neqs = length( cvx___.equalities );
    temp = struct( ...
        'name',          name,   ...
        'complete',      true,   ...
        'sdp',           false,  ...
        'gp',            false,  ...
        'separable',     false,  ...
        'locked',        false,  ...
        'precision',     nprec,  ...
        'solver',        nsolv,  ...
        'gptol',         ngprec, ...
        'rat_growth',    nrprec, ...
        't_variable',    logical( sparse( nres, 1 ) ), ...
        'n_equality',    0,          ...
        'n_linform',     0,          ...
        'n_uniform',     0,          ...
        'variables',     [],         ...
        'duals',         [],         ...
        'dvars',         [],         ...
        'direction',     '',         ...
        'geometric',     [],         ...
        'objective',     [],         ...
        'status',        'unsolved', ...
        'result',        [],         ...
        'depth',         length( st ) - 3, ...
        'self',          z );
    temp.t_variable( 1 ) = true;
    if isempty( cvx___.problems ),
        cvx___.problems = temp;
    else
        cvx___.problems( end + 1 ) = temp;
    end

case 1,

    if ~isequal( x, 'current' ),
        error( 'Argument must be the string ''current''.' );
    elseif isempty( cvx___.problems ),
        error( 'There is no problem in progress.' );
    else
        z = cvx___.problems( end ).self;
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
