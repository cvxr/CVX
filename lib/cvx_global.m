function cvx_global( prefs )

%CVX_GLOBAL   Create CVX's global internal data structure, if needed.
%   CVX_GLOBAL creates a hidden structure CVX needs to do its work. It is
%   harmless for the user to call it, but it is also useless to do so.

global cvx___ 
if isfield( cvx___, 'problems' ),
    return
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the global data structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1,
    try 
        if strncmp( computer, 'PC', 2 ), fs = '\'; else fs = '/'; end
        prefs = load( [ prefdir, fs, 'cvx_prefs.mat' ] );
    catch %#ok
        error( 'CVX:Preferences', 'Could not load saved CVX preferences; cannot continue. Please run CVX_SETUP.' );
    end
end
if ~isempty( prefs.license ),
    if ~usejava('jvm'),
        warning( 'CVX:Licensing', cvx_error( ...
            'The CVX licensing system requires the Java VM. The professional features of CVX are disabled.', ...
            [66,75], false, '', true ) );
        prefs.license = [];
    elseif ~exist( 'cvx_license', 'file' ), 
        warning( 'CVX:Licensing', cvx_error( ...
            'The CVX licensing system cannot be found. The professional features of CVX are disabled.', ...
            [66,75], false, '', true ) );
        prefs.license = [];
    else
        try
            prefs.license = cvx_license( prefs.license );
            if prefs.license.days_left < 0,
                [ dummy, errmsg ] = cvx_license(''); %#ok
                errmsg = sprintf( '%s\n', errmsg{2:end} );
                warning( 'CVX:Licensing', 'A problem was found with the current CVX license:\n%sThe professional features of CVX are disabled.\nPlease correct the issues and re-run CVX_SETUP.', errmsg );
                prefs.license = [];
            end
        catch errmsg
            prefs.license = [];
            errmsg = cvx_error( errmsg, 67, false, '    ' );
            warning( 'CVX:Licensing', 'The CVX licensing system encountered an unexpected error:\n%sThe professional features of CVX are disabled.\nPlease correct the issue or contact CVX Research support.', errmsg );
        end
    end
end

cvx___ = struct( ...
    'expert',       prefs.expert,     ...
    'precision',    prefs.precision,  ...
    'precflag',     prefs.precflag,   ...
    'rat_growth',   prefs.rat_growth, ...
    'path',         prefs.path,       ...
    'license',      prefs.license,    ...
    'solvers',      prefs.solvers,    ...
    'problems',     [],    ...
    'id',           0,     ...
    'pause',        false, ...
    'quiet',        false, ...
    'profile',      false, ...
    'reserved',     1, ...
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
temp = cvx( [0,1], [] );
cvx___.equalities = temp;
cvx___.linforms   = temp;
cvx___.linrepls   = temp;
cvx___.uniforms   = temp;
cvx___.unirepls   = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run each shim to connect/reconnect the solvers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cur_d = pwd;
solvers = cvx___.solvers.list;
nsolv = length(solvers);
nrej = 0;
for k = 1 : length(solvers),
    tsolv = solvers(k);
    try
        cd(tsolv.spath);
        tsolv.warning = '';
        tsolv = feval(tsolv.sname,tsolv);
    catch errmsg
        errmsg = cvx_error( errmsg, 63, false, '    ' );
        if isempty( tsolv.name ),
            tsolv.name = [ tsolv.spath, tsolv.sname ];
        end
        tsolv.error = sprintf( 'unexpected error:\n%s', errmsg );
    end
    if ~isempty(tsolv.error),
        nrej = nrej + 1;
    end
    solvers(k) = tsolv;
end
cvx___.solvers.list = solvers;
cd( cur_d );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If any solvers have errors, force the user to re-run cvx_setup. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nrej,
    reject = {};
    reject_lic = {};
    for k = 1 : nsolv,
        if ~isempty( solvers(k).error ),
            if isequal( solvers(k).error, 'A CVX Professional license is required.' ),
                reject_lic{end+1} = solvers(k).name; %#ok
            else
                errmsg = [ solvers(k).name, ': ', solvers(k).error ];
                reject{end+1} = cvx_error( errmsg, 67, false, '    ' ); %#ok
            end
        end
    end
    if ~isempty( reject_lic ),
        reject_lic = sprintf( '%s ', reject_lic{:} );
        warning( 'CVX:SolverErrors', 'The following solvers are are disabled due to licensing issues: %s', reject_lic );
    end
    if ~isempty( reject ),
        reject = sprintf( '%s', reject{:} );
        warning( 'CVX:SolverErrors', 'The following errors were issued when initializing the solvers:\n%sPlease check your installation and re-run CVX_SETUP.\nThese solvers are unavailable for this session.%s', reject );
    end
    if nrej == length(solvers ),
        clear global cvx___
        error( 'CVX:SolverErrors', 'All solvers were disabled due to errors or license issues.\nPlease re-run CVX_SETUP and, if necessary, contact CVX Research for support.' );
    elseif ~isempty(solvers(cvx___.solvers.map.default).error),
        for k = 1 : nsolv,
            if isempty(solvers(k).error),
                cvx___.solvers.map.default = k;
                cvx___.solvers.selected = k;
                break
            end
        end
        warning( 'CVX:SolverErrors', 'The default solver has temporarily been changed to %s.', solvers(k).name );
    end
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
