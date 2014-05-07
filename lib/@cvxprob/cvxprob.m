function z = cvxprob( varargin )

global cvx___
cvx_global

%
% Clear out any old problems left at this depth. These will be here due to
% an error while constructing a CVX model, or due to the user deciding to
% start over when constructing a model
%

st = dbstack;
depth = length( st ) - 2;
if length(st) <= 2,
    name = '';
else
    name = st(3).name;
end
if ~isempty( cvx___.problems ) && cvx___.problems(end).depth >= depth,
    ndx = find( [ cvx___.problems.depth ] >= depth, 1, 'first' );
    cvx_pop( ndx, true, true );
end

%
% Grab the latest defaults to place in the new problem
%

if ~isempty( cvx___.problems ),
    temp = cvx___.problems( end );
    nprec  = temp.precision;
    npflag = temp.precflag;
    nrprec = temp.rat_growth;
    nsolv  = temp.solver;
    nquiet = temp.quiet;
else
    nprec  = cvx___.precision;
    npflag = cvx___.precflag;
    nrprec = cvx___.rat_growth;
    selected = cvx___.solvers.selected;
    nsolv  = struct( 'index', selected, 'settings', { cvx___.solvers.list(selected).settings } );
    nquiet = cvx___.quiet;
end

%
% Construct the object
%

z = class( struct( 'index_', length( cvx___.problems ) + 1, 'id_', cvx_id( cvx ) ), 'cvxprob' );
temp = struct( ...
    'name',          name,   ...
    'complete',      true,   ...
    'sdp',           false,  ...
    'gp',            false,  ...
    'separable',     false,  ...
    'locked',        false,  ...
    'precision',     nprec,  ...
    'precflag',      npflag, ...
    'solver',        nsolv,  ...
    'quiet',         nquiet, ... 
    'cputime',       cputime, ...
    'rat_growth',    nrprec, ...
    't_variable',    sparse( 1, 1, 1, length( cvx___.classes ), 1 ), ...
    'n_equality',    length( cvx___.equalities ), ...
    'variables',     [],         ...
    'duals',         [],         ...
    'dvars',         [],         ...
    'direction',     '',         ...
    'geometric',     [],         ...
    'objective',     [],         ...
    'status',        'unsolved', ...
    'result',        [],         ...
    'bound',         [],         ...
    'iters',         Inf,        ...
    'tol',           Inf,        ...
    'cleared',       false,      ...
    'clearmode',     'clear',    ...
    'depth',         depth, ...
    'self',          z );

%
% Process the argument strings
%

for k = 1 : nargin,
    mode = varargin{k};
    if ~ischar( mode ),
        error( 'CVX:ArgError', 'Arguments must be strings.' );
    end
    switch lower( mode ),
        case 'quiet',
            temp.quiet = true;
        case 'set',
            temp.complete  = false;
            temp.direction = 'find';
        case 'sdp',
            if temp.gp,
                error( 'CVX:ArgError', 'The GP and SDP modifiers cannot be used together.' );
            end
            temp.sdp = true;
        case 'gp',
            if temp.sdp,
                error( 'CVX:ArgError', 'The GP and SDP modifiers cannot be used together.' );
            end
            temp.gp = true;
            cvx___.geometric = true;
            if cvx___.expert == 0,
                cvx___.expert = -1;
            end
        case 'separable',
            temp.separable = true;
        otherwise,
            error( 'CVX:ArgError', 'Invalid CVX problem modifier: %s', mode );
    end
end

%
% Add the problem to the stack
%

if isempty( cvx___.problems ),
    cvx___.problems = temp;
    cvx_setpath(1);
else
    cvx___.problems( end + 1 ) = temp;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
