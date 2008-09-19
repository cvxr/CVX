%CVX_SETPATH   Create CVX's global internal data structure, if needed.
%   CVX_GLOBAL creates a hidden structure CVX needs to do its work. It is
%   harmless for the user to call it, but it is also useless to do so.

global cvx___
if isempty( cvx___ ),

    %
    % Determine the version number, used to select the proper version-specific
    % behavior in those cases where back-compatability fails.
    %

    isoct = exist('OCTAVE_VERSION');
    ver = version;
    temp = find( ver == '.' );
    if length( temp ) > 1,
        ver( temp( 2 ) : end ) = [];
    end
    ver = eval( ver, 'NaN' );
    if ~isoct & ver >= 7.1,
        warning( 'off', 'MATLAB:dispatcher:ShadowedMEXExtension' );
    end

    %
    % Construct the path control structure
    %

    pstr = struct(        ...
        'string', '',     ...
        'solvers','',     ...
        'sactive', '',    ...
        'active',  false, ...
        'hold',    false );

    %
    % Construct the primary data structure
    %

    cvx___ = struct( ...
        'path',         pstr,  ...
        'problems',     [],    ...
        'octave',       isoct, ...
        'mversion',     ver,   ...
        'hcellfun',     isoct | ver > 7.1, ... 
        'id',           0,     ...
        'pause',        false, ...
        'quiet',        false, ...
        'profile',      false, ...
        'expert',       false, ...
        'gptol',        0.001, ...
        'solver',       'SDPT3', ...
        'precision',    [ eps^0.5, eps^0.5, eps^0.25 ], ...
        'rat_growth',   10, ...
        'reserved',     1, ...
        'geometric',    sparse( 1, 1 ), ...
        'logarithm',    sparse( 1, 1 ), ...
        'exponential',  sparse( 1, 1 ), ...
        'vexity',       0, ... % sparse( 1, 1 ), ...
        'nan_used',     false, ...
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


    %
    % Construct the path string
    %

    opath = matlabpath;
    s = which( 'cvx_begin' );
    if isempty( s ),
        clear cvx___
        error( 'CVX has not been set up properly. Please re-run cvx_setup.' );
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
    subs = { 'sdpt3', 'sdpt3/Solver', 'sdpt3/HSDSolver', 'sdpt3/Solver/Mexfun', 'sdpt3/Linsysolver/spchol' };
    if ~cvx___.octave & strcmp( mexext, 'mexw32' ) & cvx___.mversion >= 7.5,
        subs{end+1} = 'sedumi/mexw32';
    end
    subs{end+1} = 'sedumi';
    nsolver = length( subs );
    miss_solv = 0;
    if cvx___.octave | cvx___.mversion >= 7.0,
        subs{end+1} = 'keywords';
        subs{end+1} = 'sets';
    end
    npath = '';
    spaths = [];
    needupd = false;
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
                if k > nsolver,
                    cvx___.path.active = true;
                elseif isempty( cvx___.path.sactive ) & strcmpi( cvx___.solver, base ),
                    cvx___.path.sactive = base;
                else
                    opath( ndxs(1) : ndxs(1) + length(temp2) - 1 ) = [];
                end
                needupd = true;
            end
            if k > nsolver,
                npath = [ npath, temp2 ];
            elseif isempty( spaths ),
                spaths = struct( base, temp2 );
            elseif isfield( spaths, base ),
                spaths = setfield( spaths, base, [ getfield( spaths, base ), temp2 ] );
            else
                spaths = setfield( spaths, base, temp2 );
            end
        elseif k > nsolver,
            error( sprintf( [ ...
                'Cannot find the required cvx subdirectory: %s\n', ...
                'The cvx distribution is corrupt; please reinstall.' ], temp ) );
        elseif isfield( spaths, base ),
            error( sprintf( [ ...
                'The cvx solver directory %s is incomplete.\n', ...
                'The cvx distribution is corrupt; please reinstall.' ], temp ) );
        end
    end
    if exist( 'mosekopt', 'file' ) == 3,
        spaths.mosek = '';
    end
    cvx___.path.solvers = spaths;
    cvx___.path.string  = npath;
    if needupd,
        if cvx___.path.active,
            opath = [ cvx___.path.string, opath ];
        end
        if ~isempty( cvx___.path.sactive ),
            opath = [ getfield( spaths, cvx___.path.sactive ), opath ];
        end
        s = warning('off');
        matlabpath(opath);
        warning(s);
    end
    
    %
    % Create the global cvx objects
    %

    cvx___.equalities = cvx( [1,0], [] );
    cvx___.linforms   = cvx( [1,0], [] );
    cvx___.linrepls   = cvx( [1,0], [] );
    cvx___.uniforms   = cvx( [1,0], [] );
    cvx___.unirepls   = cvx( [1,0], [] );
    
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
