function z = cvxprob( x )

global cvx___

switch nargin,

case 0,

    name = '';
    st = dbstack;
    depth = length( st ) - 3;

    if ~isempty( cvx___.problems ),
        ndx = find( [ cvx___.problems.depth ] >= depth );
        if ~isempty( ndx ),
            cvx_pop( cvx___.problems( ndx(1) ).self, 'reset' );
            if ndx == 1, cvx_setpath( 1 ); end
        end
    end
    
    if ~isempty( cvx___.problems ),
        nprec  = cvx___.problems( end ).precision;
        ngprec = cvx___.problems( end ).gptol;
        nrprec = cvx___.problems( end ).rat_growth;
        nsolv  = cvx___.problems( end ).solver;
        nquiet = cvx___.problems( end ).quiet;
    else
        nprec  = cvx___.precision;
        ngprec = cvx___.gptol;
        nrprec = cvx___.rat_growth;
        nsolv  = cvx___.solver;
        nquiet = cvx___.quiet;
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
        'quiet',         nquiet,  ... 
        'rat_growth',    nrprec, ...
        't_variable',    logical( sparse( nres, 1 ) ), ...
        'n_equality',    length(cvx___.equalities), ...
        'n_linform',     length(cvx___.linforms), ...
        'n_uniform',     length(cvx___.uniforms), ...
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

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
