function cvx_push( name, depth, args )

tstart = tic;
cstart = cputime;
global cvx___
persistent tstruct

%
% Grab the latest defaults to place in the new problem
%

np = length( cvx___.problems );
if np
    pstr = cvx___.problems( np );
    if pstr.depth >= depth,
        p = find( [ cvx___.problems.depth ] >= depth, 1, 'first' );
        id = cvx___.problems( p ).id;
        cvx_pop( p, id );
        np = p - 1;
        if np, pstr = cvx___.problems( np ); end
    end
end

%
% Construct the object
%

id = cvx___.id + 1;
if isempty( tstruct ),
    q = [];
    tstruct = struct( ...
        'name',          '',    ...
        'finished',      false, ...
        'sdp',           false, ...
        'gp',            false, ...
        'dualize',       0,     ...
        'precision',     q,     ...
        'precflag',      q,     ...
        'solver',        q,     ...
        'settings',      q,     ...
        'quiet',         q,     ... 
        'cputime',       q,     ...
        'tictime',       q,     ...
        'checkpoint',    q,     ...
        'variables',     q,     ...
        'duals',         q,     ...
        'dvars',         q,     ...
        'direction',     Inf,   ...
        'objective',     q,     ...
        'depth',         0,     ...
        'id',            0 );
end
temp = tstruct;
temp.name    = name;
temp.id      = id;
temp.depth   = depth;
temp.cputime = cstart;
temp.tictime = tstart;
if np
    temp.precision  = pstr.precision;
    temp.precflag   = pstr.precflag;
    temp.solver     = pstr.solver;
    temp.settings   = pstr.settings;
    temp.quiet      = pstr.quiet;
    ncheck = length( cvx___.classes );
    temp.checkpoint = [ ncheck, length( cvx___.equalities ), ...
        length( cvx___.cones ), ncheck ];
else
    temp.precision  = cvx___.precision;
    temp.precflag   = cvx___.precflag;
    temp.solver     = cvx___.solvers.selected;
    temp.settings   = cvx___.solvers.list(temp.solver).settings;
    temp.quiet      = cvx___.quiet;
    temp.checkpoint = [ 1, 0, 0, 1 ];
end

%
% Process the argument strings
%

for k = 1 : length(args),
    mode = args{k};
    if ~ischar( mode ),
        cvx_throw( 'Arguments must be strings.' );
    end
    switch lower( mode ),
        case 'dualize',
            temp.dualize = 1;
        case 'dualize_off',
            temp.dualize = -1;
        case 'quiet',
            temp.quiet = true;
        case 'set',
            if isempty( cvx___.problems ),
                cvx_throw( 'Cannot construct a set object outside of a CVX model.' );
            end
            temp.direction = NaN;
        case 'sdp',
            if temp.gp,
                cvx_throw( 'The GP and SDP modifiers cannot be used together.' );
            end
            temp.sdp = true;
        case 'gp',
            if temp.sdp,
                cvx_throw( 'The GP and SDP modifiers cannot be used together.' );
            end
            temp.gp = true;
            if cvx___.expert == 0,
                cvx___.expert = -1;
            end
        otherwise,
            cvx_throw( 'Invalid CVX problem modifier: %s', mode );
    end
end

%
% Add the problem to the stack
%

cvx___.id = id;
if np,
    cvx___.problems(np+1) = temp;
else
    cvx___.problems = temp;
    cvx_setpath(1);
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
