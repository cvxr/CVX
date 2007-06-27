function cvx_setpath( arg )

%CVX_SETPATH   Sets the cvx path.
%   CVX_SETPATH adds the internal cvx directories to Matlab's path so that the
%   CVX system can find the functions that they contain. There is no reason to 
%   call this function during normal use of CVX; it is done automatically as
%   needed. However, if you are debugging CVX, calling this function can help to
%   insure that breakpoints stay valid.

% Create the global cvx data structure
global cvx___
if isempty( cvx___ ),
    ver = version;
    temp = find( ver == '.' );
    if length( temp ) > 1,
        ver( temp( 2 ) : end ) = [];
    end
    ver = eval( ver, 'NaN' );
    pstr = struct( ...
        'string', '',    ...
        'solvers','',    ...
        'formed', false, ...
        'active', false, ...
        'hold',   false );
    cvx___ = struct( ...
        'path',         pstr,  ...
        'problems',     [],    ...
        'mversion',     ver,   ...
        'id',           0,     ...
        'pause',        false, ...
        'quiet',        false, ...
        'profile',      false, ...
        'gptol',        0.001, ...
        'solver',       'SDPT3', ...
        'precision',    [ eps^0.5, eps^0.25 ], ...
        'rat_growth',   10, ...
        'reserved',     1, ...
        'geometric',    sparse( 1, 1 ), ...
        'logarithm',    sparse( 1, 1 ), ...
        'exponential',  sparse( 1, 1 ), ...
        'vexity',       0, ... % sparse( 1, 1 ), ...
        'canslack',     false, ...
        'readonly',     0,  ...
        'equalities',   [], ...
        'needslack',    logical( zeros( 0, 1 ) ), ...
        'linforms',     [], ...
        'linrepls',     [], ...
        'uniforms',     [], ...
        'unirepls',     [], ...
        'cones',        struct( 'type', {}, 'indices', {} ), ...
        'x',            zeros( 0, 1 ), ...
        'y',            zeros( 0, 1 ) );
end

% Set the hold flag
if nargin == 0,
    cvx___.path.hold = true;
end

% Determine the path string to add
if ~cvx___.path.formed,
    opath = matlabpath;
    s = which( 'cvx_begin' );
    if isempty( s ),
        error( 'Cannot set the cvx path unless cvx_begin is in the current path.' );
    end
    if ispc,
        fs = '\';
        ps = ';';
    else
        fs = '/';
        ps = ':';
    end
    temp = strfind( s, fs );
    s( temp(end-1) + 1 : end ) = [];
    solvers = { 'sdpt3/Solver', 'sdpt3/SolverHSD', 'sdpt3/Solver/Mexfun', 'sdpt3/Linsysolver/spchol', 'sdpt3/Linsysolver/MA47', 'sedumi' };
    subs = solvers;
    miss_solv = 0;
    if cvx___.mversion >= 6.5,
        subs{end+1} = 'keywords';
        subs{end+1} = 'sets';
    end
    npath = '';
    spaths = [];
    for k = 1 : length( subs ),
        tsub = subs{k};
        tt = tsub == '/';
        if any( tt ),
            base = tsub(1:min(find(tt))-1);
            tsub(tt) = fs;
        else
            base = tsub;
        end
        temp = [ s, tsub ];
        if exist( temp, 'dir' ),
            temp2 = [ temp, ps ];
            ndxs = strfind( opath, temp2 );
            if ~isempty( ndxs ),
                opath( ndxs(1) : ndxs(1) + length(temp2) - 1 ) = [];
                matlabpath( opath );
            end
            if k > length(solvers),
                npath = [ npath, temp2 ];
            elseif isempty( spaths ),
                spaths = struct( base, temp2 ); 
            elseif isfield( spaths, base ),
                spaths = setfield( spaths, base, [ getfield( spaths, base ), temp2 ] );
            else
                spaths = setfield( spaths, base, temp2 );
            end
        elseif k > length(solvers),
            error( sprintf( [ ...
                'Cannot find the required cvx subdirectory: %s\n', ...
                'The cvx distribution is corrupt; please reinstall.' ], temp ) );
        elseif isfield( spaths, base ),
            error( sprintf( [ ...
                'The cvx solver directory %s is incomplete.\n', ...
                'The cvx distribution is corrupt; please reinstall.' ], temp ) );
        end
    end
    cvx___.path.solvers = spaths;
    cvx___.path.string = npath;
end

% Add the string to the path, if it's not already there
if ~cvx___.path.active,
    if ~isempty( cvx___.path.string ),
        matlabpath( [ cvx___.path.string, matlabpath ] );
    end
end

if ~isa( cvx___.equalities, 'cvx' ),
    cvx___.equalities = cvx( [1,0], [] );
    cvx___.linforms   = cvx( [1,0], [] );
    cvx___.linrepls   = cvx( [1,0], [] );
    cvx___.uniforms   = cvx( [1,0], [] );
    cvx___.unirepls   = cvx( [1,0], [] );
end

cvx___.path.formed = true;
cvx___.path.active = true;

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
