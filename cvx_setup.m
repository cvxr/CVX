function cvx_setup( license_file )

% CVX_SETUP   Sets up and tests the cvx distribution.
%    This function is to be called any time CVX is installed on a new machine,
%    to insure that the paths are set properly and the MEX files are compiled.

% Clear out the global CVX structure
global cvx___
cvx___ = [];
nret = false;
oldpath = '';
line = '---------------------------------------------------------------------------'; 

try 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get version and portability information %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    expected = MException( 'CVX:Expected', '' );
    unexpected = MException( 'CVX:Unexpected', '' );
    if nargin < 1, license_file = []; end
    [ nver, isoctave, fs, ps, mpath, problem ] = cvx_version( license_file ); %#ok
    if problem, throw(expected); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the CVX paths %
    %%%%%%%%%%%%%%%%%%%%%%%%

    oldpath = cvx_startup( false );
    cpath = struct( 'string', '', 'active', false, 'solver', 0, 'hold', false );
    subs = strcat( [ mpath, fs ], { 'keywords', 'sets' } );
    cpath.string = sprintf( [ '%s', ps ], subs{:} );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Search for saved preferences %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf( 'Saved preferences...' ); nret = true;
    pfile = [ prefdir, fs, 'cvx_prefs.mat' ];
    oprefs = [];
    if exist( pfile, 'file' ),
        oprefs = load( pfile );
        fprintf( 'found.\n' ); nret = false;
    else
        fprintf( 'not found; defaults created.\n' ); nret = false;
    end
    prefs = struct( 'expert', false, 'precision', [eps^0.5,eps^0.5,eps^0.25], ...
                    'precflag', 'default', 'rat_growth', 10, ...
                    'path', cpath, 'license', [], 'solvers', [] );
    selected = 'sdpt3';
    if ~isempty( oprefs ),
        try prefs.expert = oprefs.expert; catch end %#ok
        try prefs.precision = oprefs.precision; catch end %#ok
        try prefs.precflag = oprefs.precflag; catch end %#ok
        try prefs.rat_growth = oprefs.rat_groth; catch end %#ok
        try selected = prefs.solvers.names{prefs.solvers.map.default}; catch end %#ok
    end
    prefs.path = cpath;
    prefs.license = cvx___.license;
    prefs.solvers = struct( 'selected', 0, 'active', 0, 'list', [], 'names', {{}}, 'map', struct( 'default', 0 ) );
    prefs.solvers.list = struct( 'name', {}, 'error', {}, 'warning', {}, 'dualize', {}, 'path', {}, 'check', {}, 'solve', {}, 'settings', {}, 'spath', {}, 'sname', {} );
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Search for solvers %
    %%%%%%%%%%%%%%%%%%%%%%

    fprintf( 'Searching for solvers...' ); nret = true;
    shimpath = [ mpath, fs, 'shims', fs ];
    solvers = dir( shimpath );
    solvers = { solvers(~[solvers.isdir]).name };
    solvers = solvers( ~cellfun( @isempty, regexp( solvers, '\.(m|p)$' ) ) );
    solvers = unique( cellfun( @(x)x(1:end-2), solvers, 'UniformOutput', false ) );
    solvers = struct( 'name', '', 'error', '', 'warning', '', 'dualize', '', 'path', '', 'check', [], 'solve', [], 'settings', [], 'sname', solvers, 'spath', shimpath );
    solver2 = which( 'cvx_solver_shim', '-all' );
    for k = 1 : length(solver2),
        tsolv = solver2{k};
        ndxs = find(tsolv==fs,1,'last');
        solvers(end+1).spath = tsolv(1:ndxs-1); %#ok
        solvers(end).sname = tsolv(ndxs+1:end);
    end
    cur_d = pwd;
    nsolv = length(solvers);
    nrej = 0; nwarn = 0;
    for k = 1 : length(solvers),
        tsolv = solvers(k);
        try
            cd(tsolv.spath);
            tsolv = feval( tsolv.sname, tsolv );
        catch errmsg
            errmsg = cvx_error( errmsg, 63, false, '  ' );
            if isempty( tsolv.name ),
                tsolv.name = [ tsolv.spath, tsolv.sname ];
            end
            tsolv.error = sprintf( 'unexpected error:\n%s', errmsg );
        end
        nrej = nrej + ~isempty(tsolv.error);
        nwarn = nwarn + ~isempty(tsolv.warning);
        solvers(k) = tsolv;
    end
    cd( cur_d );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Process solver errors %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    plurals = { 's', '', 's' };
    nwarn = sum(~cellfun(@isempty,{solvers.warning}));
    fprintf( '%d shim%s found.\n', nsolv, plurals{min(nsolv+1,3)} ); nret = false;
    if nsolv,
        fprintf( '%d solver%s initialized:\n', nsolv-nrej, plurals{min(nsolv-nrej+1,3)} );
        for k = 1 : nsolv,
            if isempty( solvers(k).error ),
                if strcmpi( solvers(k).name, selected ),
                    fprintf( '    %s (default)\n', solvers(k).name );
                else
                    fprintf( '    %s\n', solvers(k).name );
                end
            end
        end
    end
    if nrej,
        fprintf( '%d solver%s skipped:\n', nrej, plurals{min(nrej+1,3)} );
        for k = 1 : nsolv,
            if ~isempty( solvers(k).error ),
                cvx_error( [ solvers(k).name, ': ', solvers(k).error ], 67, false, '    ' );
            end
        end
    end
    if nwarn,
        fprintf( '%d solver%s issued warnings:\n', nwarn, plurals{min(nwarn+1,3)} );
        for k = 1 : nsolv,
            if ~isempty( solvers(k).warning ) && isempty( solvers(k).error ),
                cvx_error( [ solvers(k).name, ': ', solvers(k).warning ], 67, false, '    ' );
            end
        end
    end
    if nrej == nsolv,
        fprintf( [ ... 
            'No valid solvers were found. This suggests a corrupt installation. Please\n', ...
            'try re-installing the files and re-running cvx_setup. If the same error\n', ...
            'occurs, please contact CVX support.\n' ] );
        throw(unexpected);
    end
    solvers = solvers(cellfun(@isempty,{solvers.error}));
    prefs.solvers.list  = solvers;
    prefs.solvers.names = { solvers.name };
    prefs.solvers.map   = struct( 'default', 0 );
    for k = 1 : length(solvers),
        prefs.solvers.map.(lower(solvers(k).name)) = k;
        if strcmpi( solvers(k).name, selected ),
            prefs.solvers.map.default = k;
            prefs.solvers.selected = k;
        end
    end
    if prefs.solvers.selected == 0,
        prefs.solvers.selected = 1;
        prefs.solvers.map.default = 1;
        fprintf( [ ...
            'WARNING: The default solver %s is missing; %s has been selected as a\n', ...
            '    new default. If this was unexpected, try re-running cvx_setup.\n' ], ...
            selected, prefs.solvers.list(prefs.solvers.selected).name, lower(prefs.solvers.list(prefs.solvers.selected).name) );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the global data structure %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cvx_global( prefs );
    if isempty( cvx___.solvers.list ),
        throw(unexpected);
    end
    
    %%%%%%%%%%%%%%%%%%%%
    % Save preferences %
    %%%%%%%%%%%%%%%%%%%%
    
    fprintf( 'Saving updated preferences...' ); nret = true;
    try
        cvx_save_prefs( true );
        fprintf( 'done.\n' ); nret = false;
    catch errmsg
        fprintf( 'unexpected error:\n' ); nret = false;
        cvx_error( errmsg, 67, true, '    ' );
        fprintf( 'Please attempt to correct this error and re-run cvx_setup. If you cannot,\n' );
        fprintf( 'you will be forced to re-run cvx_setup every time you start MATLAB.\n' );
    end
    if ~isempty( oldpath ),
        fprintf( 'Saving updated path...' ); nret = true;
        if savepath,
            fprintf('failed. (see below)\n');
        else
            fprintf('done.\n');
            oldpath = [];
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Test the distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf( 'Testing with a simple model...' ); nret = true;
    m = 16; n = 8;
    A = randn(m,n);
    b = randn(m,1);
    cvx_begin('quiet')
        variable('x(n)');
        minimize( norm(A*x-b,1) );
    cvx_end
    fprintf( 'done!\n' ); nret = false;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Quick instructions on changing the solver %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if length( cvx___.solvers.list ) > 1,
        fprintf( '%s\n', line );
        fprintf( 'To change the default solver, type "cvx_solver <solver_name>".\n')
        fprintf( 'To save this change for future sessions, type "cvx_save_prefs".\n' );
        fprintf( 'Please consult the users'' guide for more information.\n' );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instruct the user to save the path %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty( oldpath )
        need_upr = false;
        need_disclaim = true;
        user_path = which( 'startup.m' );
        if isempty( user_path ),
            user_path = userpath;
            if length(user_path) <= 1,
                need_upr = true;
                user_path = system_dependent('getuserworkfolder', 'default');
                if ~isempty( user_path ),
                    if isempty( strfind( user_path, [ fs, 'MATLAB' ] ) ),
                        user_path = [ user_path, fs, 'MATLAB' ];
                    end
                end
            elseif user_path(end) == ps,
                user_path(end) = '';
            end
            if ~isempty( user_path ),
                user_file = [ user_path, fs, 'startup.m' ];
            else
                user_file = '';
            end
        else
            user_file = user_path;
            user_path = user_path(1:end-10);
        end
        fprintf( '%s\n', line );
        fprintf('NOTE: the MATLAB path has been changed to point to the CVX distribution. To\n' );
        fprintf('use CVX without having to re-run CVX_SETUP every time MATLAB starts, you\n' );
        fprintf('will need to save this path permanently. This script attempted to do this\n' ); 
        fprintf('for you, but failed---likely due to UNIX permissions restrictions.\n' );
        if exist( user_file, 'file' ),
            fprintf( 'To solve the problem, edit the file\n    %s\nand add the following line to the end of the file:\n', user_file ); 
            fprintf( '    run %s%scvx_startup.m\n', mpath, fs );
        elseif exist( user_path, 'dir' ),
            fprintf( 'To solve the problem, create a new file\n    %s\ncontinaing the following line:\n', user_file );
            fprintf( '    run %s%scvx_startup.m\n', mpath, fs );
        elseif ~isempty( user_path ),
            need_upr = false;
            fprintf( 'To solve the problem, perform the following steps:\n' );
            fprintf( '    1) Create a new directory %s:\n           !mkdir -p %s\n', user_path, user_path );
            fprintf( '    2) Create a new startup.m file in that directory:\n           edit %s\n', user_file );
            fprintf( '    3) Add the following line to that file, and save:\n           run %s%scvx_startup.m\n', mpath, fs );
            fprintf( '    4) Execute the commands:\n           userpath reset; startup\n' );
        else
            fprintf( 'To solve the problem, create a startup.m file containing the line:\n' );
            fprintf( '    run %s%scvx_startup.m\n', mpath, fs );
            fprintf( 'Consult the MATLAB documentation for the proper location for that file.\n' );
            need_disclaim = false;
        end
        if need_upr,
            fprintf( 'Finally, execute the following MATLAB commands:\n    userpath reset; startup\n' );
        end
        if need_disclaim,
            fprintf( 'Please consult the MATLAB documentation for more information about the\n' );
            fprintf( 'startup.m file and its proper placement and usage.\n' );
        end
    end

catch errmsg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Restore the environment in the event of an error %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    unexpected = false;
    if nret, fprintf( '\n' ); end
    switch errmsg.identifier,
        case { 'CVX:Expected', 'CVX:Licensing' },
            if ~isempty( errmsg.message ),
                cvx_error( errmsg, 67, 'ERROR: ', '    ' );
            end
        case 'CVX:Unexpected',
            unexpected = true;
        otherwise,
            cvx_error( errmsg, 67, 'UNEXPECTED ERROR: ', '    ' );
            unexpected = true;
    end
    if ~isempty( oldpath ),
        path( oldpath );
    end
    clear global cvx___
    if unexpected,
        fprintf( 'Please report this error to support, and include entire output of\n' );
        fprintf( 'CVX_SETUP in your support request.\n' );
    end

end

fprintf( '%s\n\n', line );

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
