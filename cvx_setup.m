function cvx_setup( varargin )

% CVX_SETUP   Sets up and tests the cvx distribution.
%    This function is to be called any time CVX is installed on a new machine,
%    to insure that the paths are set properly and the MEX files are compiled.

disp( ' ' );
dd = pwd;

global cvx___
cvx___ = [];

ver = version;
verm = 0;
temp = find( ver == '.' );
if length( temp ) > 1,
	verm = ver( temp( 2 ) + 1 : end );
    ver( temp( 2 ) : end ) = [];
end
needLarge = 0;
ver = eval( ver, 'NaN' );
try
    mpath = mfilename('fullpath');
catch
    dbs = dbstack;
    mpath = dbs.name;
end
isoctave = exist('OCTAVE_VERSION');
if isoctave,
	if isnan( ver ) | ver < 3.0 | ( ver == 3.0 & verm < 1 ),
		error( 'CVX requires octave 3.0.1 or later.' );
	end
else
	newext = mexext;
        needLarge = strcmp( newext(end-1:end), '64' );
	if isnan( ver ) | ver < 6.1 | ( ( ver < 7.3 ) & needLarge ),
	    error( sprintf( ...
		[ 'CVX requires 32-bit MATLAB 6.1 or later,\n', ...
          '          or 64-bit MATLAB 7.3 or later.' ] ) );
    end
    if ver >= 7.1,
        warning( 'off', 'MATLAB:dispatcher:ShadowedMEXExtension' );
    end
end
if ispc,
    fs = '\';
    ps = ';';
else
    fs = '/';
    ps = ':';
end
temp = strfind( mpath, fs );
mpath( temp(end) : end ) = [];
rmpaths = { 'sets', 'keywords', 'builtins', 'commands', 'functions', 'lib', 'structures', 'matlab6' };
needpaths = {};
delepaths = {};
if isoctave | ver >= 7.0,
    checkpaths = rmpaths(1:end-1);
    addpaths = rmpaths(3:end-1);
elseif ver >= 6.5,
    checkpaths = rmpaths(1:end-1);
    addpaths = rmpaths(1:end-1);
else
    checkpaths = rmpaths;
    addpaths = rmpaths;
end
if ~isoctave & strcmp( mexext, 'mexw64' ) & ver < 7.5,
    skip_sedumi = 1;
    solvers = { 'sdpt3' };
else
    skip_sedumi = 0;
    solvers = { 'sdpt3', 'sedumi' };
end
for k = 1 : length(addpaths),
	addpaths{k} = [ mpath, fs, addpaths{k} ];
end
addpaths{end+1} = mpath;
missing = {};
for k = 1 : length(checkpaths),
    temp = [ mpath, fs, checkpaths{k} ];
    if ~exist( temp, 'dir' ),
        missing{end+1} = checkpaths{k};
    end
end
msolv = 0;
for k = 1 : length(solvers),
    if ~exist( [ mpath, fs, solvers{k} ], 'dir' ),
        missing{end+1} = solvers{k};
        msolv = msolv + 1;
    end
end
if length(missing) > msolv | msolv == length(solvers),
    error( sprintf( ...
        [ 'The following directories in the CVX distribution are missing:\n', ...
          '  %s\n', ...
          'Please reinstall the distribution and try again.' ], ...
          sprintf( ' %s', missing{:} ) ) );
end
needupd = 0;
newpath = [ ps, path, ps ];
if ispc,
    newpath2 = lower(newpath);
    mpath2 = lower(mpath);
else
    newpath2 = newpath;
    mpath2 = mpath;
end
for k = 1 : length(rmpaths),
    temp = [ ps, mpath2, fs, rmpaths{k}, ps ];
    ndxs = strfind( newpath2, temp );
    if isempty( ndxs ),
        temp = [ ps, rmpaths{k}, ps ];
        ndxs = strfind( newpath2, temp );
    end
    if ~isempty( ndxs ),
        needupd = 1;
        delepaths{end+1} = temp(2:end-1);
        len = length( temp ) - 1;
        for j = 1 : length( ndxs ),
            newpath( ndxs(j) + 1 : ndxs(j) + len ) = [];
            newpath2( ndxs(j) + 1 : ndxs(j) + len ) = [];
            ndxs = ndxs - len;
        end
    end
end
for k = 1 : length(addpaths),
    temp = [ ps, addpaths{k}, ps ];
    if ispc,
        temp = lower( temp );
    end
    ndxs = strfind( newpath2, temp );
    if isempty( ndxs ),
        needupd = 1;
        needpaths{end+1} = temp(2:end-1);
        newpath2 = [ temp(1:end-1), newpath2 ];
        newpath = [ ps, addpaths{k}, newpath ];
    end
end
newpath = newpath(2:end-1);
path(newpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the CVX MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newext = mexext;
if strcmp( newext, 'mexglx' ) & ver >= 7.2,
    fullRecompile = true;
else
    fullRecompile = any( strcmp( varargin, '-force' ) );
end
isw32 = strcmp( newext, 'mexw32' );
isNewLinux = strcmp( newext, 'mexglx' ) & ~isoctave & ver >= 7.5;
if isoctave,
	mexcmd = {};
else
	mexcmd = { '-O' };
end
if needLarge,
    mexcmd{end+1} = '-largeArrayDims';
end
libpath = [ mpath, fs, 'lib' ];
try
    cd( libpath );
    mexfiles = dir( '*.c' );
    mexfiles = { mexfiles.name };
    has_mex = 1;
    for k = 1 : length( mexfiles ),
        str = mexfiles{k};
        if ~exist( [ str(1:end-1), newext ] ),
            if ~isw32 | ~exist( [ str(1:end-1), 'dll' ] ),
                has_mex = 0;
                break;
            end
        end
    end
catch
    disp( 'The cvx distribution seems to be incomplete: the lib/ directory is' );
    disp( 'missing. Please re-unpack the distribution and try again.' );
    disp( ' ' );
    cd( dd );
    return
end
if ~has_mex | fullRecompile,
    disp( 'Attempting to generate the CVX MEX files...' );
    disp( '-------------------------------------------------------------' );
    has_mex = 1;
    for k = 1 : length( mexfiles ),
        str = mexfiles{k};
        try
            disp( sprintf( '%s ', 'mex', mexcmd{:}, str ) );
            mex( mexcmd{:}, str );
        catch
            has_mex = 0;
        end
    end
    if ~has_mex,
        disp( 'ERROR: One or more of cvx''s required MEX files could not be generated.' );
        disp( 'This is likely because your MEX system has not been properly configured.' );
        disp( 'Please consult the MATLAB documentation for details on how to do this.' );
        disp( ' ' );
        cd( dd );
        return
    else
        disp( '-------------------------------------------------------------' );
    end
end
cd( dd );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the SeDuMi MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newext = mexext;
isw32 = strcmp( newext, 'mexw32' );
sedpath = [ mpath, fs, 'sedumi' ];
if ~skip_sedumi & exist( sedpath, 'dir' ),
    if ~fullRecompile,
        mexfiles = dir( [ sedpath, fs, '*.', newext ] );
        if isempty( mexfiles ) & ispc & isw32,
            mexfiles = dir( [ sedpath, fs, '*.dll' ] );
        end
    end
    if fullRecompile | isempty( mexfiles ),
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
if exist( sedpath, 'dir' ),
    if ~fullRecompile,
        mexfiles = dir( [ sedpath, fs, 'Solver', fs, 'Mexfun', fs, '*.', newext ] );
        if isempty( mexfiles ) & ispc & isw32,
            newext = 'dll';
            mexfiles = dir( [ sedpath, fs, 'Solver', fs, 'Mexfun', fs, '*.', newext ] );
        end
    end
    if fullRecompile | isempty( mexfiles ),
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
    disp( sprintf( 'cvx differs from native Matlab by %g%%', 100 * err ) );
    disp( 'Unexpected numerical errors were found when solving the test problem.' );
    disp( 'Please report this to the authors.' );
    disp( ' ' );
    return
else
    disp( 'No errors! cvx has been successfully installed.' );
    disp( ' ' );
end

if skip_sedumi,
    disp( 'NOTE: SeDuMi 1.2 does not work with this version of Matlab.' );
    disp( 'The only solver available for this platform is SDPT3. If you wish' );
    disp( 'to use SeDuMi, upgrade to Matlab 7.5 or later.' );
    disp( ' ' );
end

if needupd,
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
            for k = 1 : length( needpaths ),
                disp( [ '    addpath ', needpaths{k} ] );
            end
            disp( 'Consult the MATLAB documentation if necessary.' );
    end
    disp( ' ' );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
