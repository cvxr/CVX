function [ z, id ] = cvx_push( name, depth, varargin )

global cvx___
try
    np = length( cvx___.problems );
catch
    cvx_global
    np = 0;
end

%
% Clear out any old problems left at this depth. These will be here due to
% an error while constructing a CVX model, or due to the user deciding to
% start over when constructing a model
%

if np && cvx___.problems(np).depth >= depth,
    ndx = find( [ cvx___.problems.depth ] >= depth, 1, 'first' );
    cvx_pop( ndx, cvx___.problems(ndx).id );
    np = ndx - 1;
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

id = cvx___.id + 1;
temp = struct( ...
    'name',          name,   ...
    'complete',      true,   ...
    'finished',      false,  ...
    'sdp',           false,  ...
    'gp',            false,  ...
    'separable',     false,  ...
    'precision',     nprec,  ...
    'precflag',      npflag, ...
    'solver',        nsolv,  ...
    'quiet',         nquiet, ... 
    'cputime',       cputime, ...
    'rat_growth',    nrprec, ...
    'checkpoint',    [ length( cvx___.classes ), length( cvx___.equalities ), length( cvx___.cones ), length( cvx___.classes ) ], ...
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
    'depth',         depth,      ...
    'id',            id );

%
% Process the argument strings
%

for k = 1 : nargin - 2,
    mode = varargin{k};
    if ~ischar( mode ),
        error( 'CVX:ArgError', 'Arguments must be strings.' );
    end
    switch lower( mode ),
        case 'quiet',
            temp.quiet = true;
        case 'set',
            if isempty( cvx___.problems ),
                error( 'CVX:ArgError', 'Cannot construct a set object outside of a CVX model.' );
            end
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

cvx___.id = id;
if np,
    z = np + 1;
    cvx___.problems( z ) = temp;
else
    cvx___.problems = temp;
    cvx_setpath(1);
    z = 1;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
