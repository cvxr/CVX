function cvx_setup( varargin )

% CVX_SETUP   Sets up and tests the cvx distribution.
%    This function is to be called any time CVX is installed on a new machine,
%    to insure that the paths are set properly and the MEX files are compiled.

disp( ' ' );
dd = pwd;

% Clear out the global CVX structure
global cvx___
cvx___ = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get version and portability information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numeric version
ver = version;
ver(ver=='.') = ' ';
ver = sscanf(ver,'%d');
ver = ver(1) + 0.01 * ( ver(2) + 0.01 * ver(3) );

% Octave?
isoctave = exist( 'OCTAVE_VERSION', 'var' );
if isoctave,
    error( 'Sorry, CVX does not yet run under octave.' );
    newext = 'mex';
    needLarge = false;
    if isnan( ver ) || ver < 3.0202,
        error( 'CVX requires octave 3.2.2 or later.' );
    end
else
    newext = mexext;
    needLarge = strcmp( newext(end-1:end), '64' );
    if isnan( ver ) || ver < 6.05 || ( ( ver < 7.03 ) && needLarge ),
        error( [ 'CVX requires 32-bit MATLAB 6.5 or later\n', ...
            'or 64-bit MATLAB 7.3 or later.' ], 1 );
    end
    if ver >= 7.01,
        warning( 'off', 'MATLAB:dispatcher:ShadowedMEXExtension' );
    end
end

% File path separators
if ispc,
    fs = '\';
    ps = ';';
else
    fs = '/';
    ps = ':';
end

% Mex locations
fullRecompile = any( strcmp( varargin, '-force' ) );
usePre75 = ~isoctave && ver < 7.05;

%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the CVX paths %
%%%%%%%%%%%%%%%%%%%%%%%%

try
    mpath = mfilename('fullpath');
catch
    dbs = dbstack;
    mpath = dbs.name;
end
temp = strfind( mpath, fs );
mpath( temp(end) : end ) = [];
if ispc, mpath = lower(mpath); end
solvers = { 'sdpt3', 'sedumi' };
rmpaths = { 'sets', 'keywords', 'builtins', 'commands', 'functions', 'lib', [ 'lib', fs, 'pre7.5' ], 'structures' };
for k = 1 : length(rmpaths),
    rmpaths{k} = [ mpath, fs, rmpaths{k} ];
end
rmpaths{end+1} = mpath;
addpaths = rmpaths;
if ~usePre75,
    addpaths(7) = [];
end
delepaths = {};
if isoctave || ver >= 7.0,
    addpaths(1:2) = [];
end
missing = {};
for k = 1 : length(addpaths),
    if ~exist( addpaths{k}, 'dir' ),
        missing{end+1} = addpaths{k};
    end
end
msolv = 0;
for k = 1 : length(solvers),
    if ~exist( [ mpath, fs, solvers{k} ], 'dir' ),
        missing{end+1} = solvers{k};
        msolv = msolv + 1;
    end
end
if length(missing) > msolv || msolv == length(solvers),
    error( [ 'The following directories in the CVX distribution are missing:\n  ', ...
             sprintf( ' %s', missing{:} ), ...
             '\nPlease reinstall the distribution and try again.' ], 1 );
end
oldpath = path;
if ispc, oldpath = lower(oldpath); end
newpath = [ ps, oldpath, ps ];
for k = 1 : length(rmpaths),
    ndxs = strfind( newpath, [ ps, rmpaths{k}, ps ] );
    if ~isempty( ndxs ),
        delepaths{end+1} = rmpaths{k};
        len = length( rmpaths{k} ) + 1;
        for j = 1 : length( ndxs ),
            newpath( ndxs(j) + 1 : ndxs(j) + len ) = [];
            ndxs = ndxs - len;
        end
    end
end
for k = 1 : length(addpaths),
    newpath = [ ps, addpaths{k}, newpath ];
end
newpath = newpath(2:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the CVX MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd( [ mpath, fs, 'lib' ] );
mexfiles = dir( '*.c' );
mexfiles = { mexfiles.name };
if ~fullRecompile,
    has_mex = 1;
    for k = 1 : length( mexfiles ),
        str = mexfiles{k};
        if exist( str(1:end-2), 'file' ) ~= 3,
            has_mex = 0;
            break;
        end
    end 
end
if fullRecompile || ~has_mex,
    disp( 'Attempting to generate the CVX MEX files...' );
    disp( '-------------------------------------------------------------' );
    if isoctave,
        mexcmd = {};
    else
        mexcmd = { '-O' };
    end
    if needLarge,
        mexcmd{end+1} = '-largeArrayDims';
    end
    has_mex = 1;
    for k = 1 : length( mexfiles ),
        str = mexfiles{k};
        try
            cmd = sprintf( '%s ', 'mex', mexcmd{:}, str );
            fprintf( 1, '%s\n', cmd );
            mex( mexcmd{:}, str );
        catch
            has_mex = 0;
        end
    end
    if has_mex,
        disp( '-------------------------------------------------------------' );
    else
        disp( 'ERROR: One or more of cvx''s required MEX files could not be generated.' );
        disp( 'This is likely because your MEX system has not been properly configured.' );
        disp( 'Please consult the MATLAB documentation for details on how to do this.' );
    end
end
cd( dd );
if ~has_mex,
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the SeDuMi MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sedpath = [ mpath, fs, 'sedumi' ];
if exist( sedpath, 'dir' ),
    has_mex = 0;
    if ~fullRecompile,
        has_mex = ~isempty( dir( [ sedpath, fs, '*.', newext ] ) );
        if ~has_mex && ~isoctave && ver < 7.05,
            has_mex = ~isempty( dir( [ sedpath, fs, 'pre7.5', fs, '*.', newext ] ) );
        end
    end
    if ~has_mex,
        cd( sedpath );
        try
            disp( 'Generating the SeDuMi MEX files.' );
            disp( '-------------------------------------------------------------' );
            install_sedumi(1);
        catch
            disp( '-------------------------------------------------------------' );
            disp( 'SeDuMi was NOT built successfully. Please try CVX on a supported' );
            disp( 'platform, or manually run the ''install_sedumi'' command in the' );
            disp( 'sedumi/ subdirectory to try and find and correct the error.' );
            disp( 'Once the error has been fixed, please re-run cvx_setup.' );
            cd( dd );
            return
        end
        disp( '-------------------------------------------------------------' );
        cd( dd );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the SDPT3 MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sedpath = [ mpath, fs, 'sdpt3' ];
mexpath = [ sedpath, fs, 'Solver', fs, 'Mexfun' ];
if ~isoctave && ver < 7.03
    disp( 'Warning: SDPT3 is no longer supported for this version of MATLAB.' );
    disp( 'SeDuMi will be the only solver available. To use SDPT3, please use' );
    disp( 'Matlab version 7.3 or later.' );
elseif exist( sedpath, 'dir' ),
    has_mex = 0;
    if ~fullRecompile,
        has_mex = ~isempty( dir( [ mexpath, fs, '*.', newext ] ) );
        if ~has_mex && ~isoctave && ver < 7.05,
            has_mex = ~isempty( dir( [ mexpath, fs, 'pre7.5', fs, '*.', newext ] ) );
        end
    end
    if ~has_mex,
        cd( sedpath );
        try
            disp( 'Generating the SDPT3 MEX files.' );
            Installmex;
        catch
            disp( '-------------------------------------------------------------' );
            disp( 'SDPT3 was NOT built successfully. Please try CVX on a supported' );
            disp( 'platform, or manually run the ''Installmex'' command in the' );
            disp( 'sdpt3/ subdirectory to try and find and correct the error.' );
            disp( 'Once the error has been fixed, please re-run cvx_setup.' );
            cd( dd );
            return
        end
        disp( '-------------------------------------------------------------' );
        cd( dd );
    end
end

path(newpath);
disp( 'Testing the cvx distribution. If this script aborts with' );
disp( 'an error, please report the error to the authors.' );
disp( '-------------------------------------------------------------' );

m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);
xls = A \ b;
cvx_begin
    cvx_quiet(true);
    variable('x(n)');
    minimize( norm(A*x-b) );
cvx_end
try
    cvx_clearpath;
end

if norm( x - xls ) > 0.01 * norm( x ),
    err = norm( x - xls ) / norm( x );
    disp( '-------------------------------------------------------------' );
    fprintf( 1, 'cvx differs from native Matlab by %g%%\n', 100 * err );
    disp( 'Unexpected numerical errors were found when solving the test problem.' );
    disp( 'Please report this to the authors.' );
    disp( ' ' );
    return
else
    disp( 'No errors! cvx has been successfully installed.' );
    disp( ' ' );
end

if ~strcmp(newpath,oldpath),
    disp( 'NOTE: The MATLAB path has been updated to point to the cvx distribution.' );
    disp( 'In order to use cvx regularly, you must save this new path definition.' );
    switch computer,
        case 'PCWIN',
        case 'MAC',
        case 'MACI',
            disp( 'To accomplish this, type the command' );
            if ver >= 7,
                disp( '    savepath' );
                disp( 'at the MATLAB prompt. Alternatively, type the command' );
            end
            disp( '    pathtool' );
            disp( 'at the MATLAB prompt, which brings up the ''Set Path'' dialog. Press' );
            disp( 'the ''Save'' button, and then the ''Close'' button.' );
        otherwise,
            disp( 'To accomplish this, add these lines to your startup.m file:' );
            for k = 1 : length( delepaths ),
                disp( [ '    rmpath ', delepaths{k} ] );
            end
            for k = length( addpaths ) : -1 : 1,
                disp( [ '    addpath ', addpaths{k} ] );
            end
            disp( 'Consult the MATLAB documentation if necessary.' );
    end
    disp( ' ' );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
