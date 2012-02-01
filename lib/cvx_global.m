%CVX_SETPATH   Create CVX's global internal data structure, if needed.
%   CVX_GLOBAL creates a hidden structure CVX needs to do its work. It is
%   harmless for the user to call it, but it is also useless to do so.

global cvx___
if isempty( cvx___ ),

    %
    % Determine the version number, used to select the proper version-specific
    % behavior in those cases where back-compatability fails.
    %

    isoct = exist( 'OCTAVE_VERSION', 'var' );
    ver = version;
    ver(ver=='.') = ' ';
    ver = sscanf(ver,'%d');
    ver = ver(1) + 0.01 * ( ver(2) + 0.01 * ver(3) );
    if ~isoct && ver >= 7.01,
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
        'hcellfun',     isoct | ver > 7.01, ... 
        'id',           0,     ...
        'pause',        false, ...
        'quiet',        false, ...
        'profile',      false, ...
        'expert',       false, ...
        'gptol',        0.001, ...
        'solver',       'sedumi', ...
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
        'needslack',    false( 0, 1 ) , ...
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
    subs = { 'sedumi/pre7.5', 'sedumi', 'sdpt3', 'sdpt3/Solver', 'sdpt3/HSDSolver', 'sdpt3/Solver/Mexfun/pre7.5', 'sdpt3/Solver/Mexfun', 'keywords', 'sets' };
    if cvx___.octave,
        smap = [ 0, 1, 1, 1, 1, 0, 1, +1, +1 ];
    elseif cvx___.mversion < 7.00,
        smap = [ 1, 1, 0, 0, 0, 0, 0, -1, -1 ];
    elseif cvx___.mversion < 7.03,
        smap = [ 1, 1, 0, 0, 0, 0, 0, +1, +1 ];
    elseif cvx___.mversion < 7.05,
        smap = [ 1, 1, 1, 1, 1, 1, 1, +1, +1 ];
    else
        smap = [ 0, 1, 1, 1, 1, 0, 1, +1, +1 ];
    end
    npath = '';
    spaths = [];
    needupd = false;
    miss_solv = 0;
    nsolv = length(subs) - 2;
    for k = 1 : nsolv + 2,
        tsub = subs{k};
        tt = tsub == '/';
        if k > nsolv,
            base = '';
        elseif any(tt),
            base = tsub(1:min(find(tt))-1);
        else
            base = tsub;
        end
        tsub(tt) = fs;
        temp = [ s, tsub ];
        if exist( temp, 'dir' ),
            temp2 = [ temp, ps ];
            ndxs = strfind( opath, temp2 );
            if ~isempty( ndxs ),
                if smap(k) > 0,
                    if k > nsolv,
                        cvx___.path.active = true;
                    elseif isempty( cvx___.path.sactive ) && strcmpi( cvx___.solver, base ),
                        cvx___.path.sactive = base;
                    end
                end
                if smap(k) >= 0,
                    opath( ndxs(1) : ndxs(1) + length(temp2) - 1 ) = [];
                    needupd = true;
                end
            end
            if smap(k) > 0,
                if isempty( base ),
                    npath = [ npath, temp2 ];
                elseif isempty( spaths ),
                    spaths = struct( base, temp2 );
                elseif isfield( spaths, base ),
                    spaths.(base) = [ spaths.(base), temp2 ];
                else
                    spaths.(base) = temp2;
                end
            end
        elseif isempty(base) || smap(k) && ~isempty( spaths ) && isfield( spaths, base ),
            error( [ ...
                'Cannot find the required cvx subdirectory: %s\n', ...
                'The cvx distribution is corrupt; please reinstall.' ], temp );
        end
    end
    cvx___.path.solvers = spaths;
    cvx___.path.string  = npath;
    if needupd,
        if cvx___.path.active,
            opath = [ cvx___.path.string, opath ];
        end
        if ~isempty( cvx___.path.sactive ),
            opath = [ spaths.(cvx___.path.sactive), opath ];
        end
        s = warning('off');
        matlabpath(opath);
        warning(s);
    end
    
    %
    % Create the global cvx objects
    %

    cvx___.equalities = cvx( [0,1], [] );
    cvx___.linforms   = cvx( [0,1], [] );
    cvx___.linrepls   = cvx( [0,1], [] );
    cvx___.uniforms   = cvx( [0,1], [] );
    cvx___.unirepls   = cvx( [0,1], [] );
    
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
