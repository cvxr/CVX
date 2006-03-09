function cvx_setup

% CVX_SETUP   Sets up and tests the cvx distribution.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

disp( ' ' );
dd = cd;

global cvx___
cvx___ = [];

ver = version;
temp = find( ver == '.' );
if length( temp ) > 1,
    ver( temp( 2 ) : end ) = [];
end
ver = eval( ver, 'NaN' );
if isnan( ver ) | ver < 6.1,
    error( sprintf( 'cvx has not yet been configured to handle your MATLAB version (%s).\nPlease contact the authors to inquire about future support plans.', version ) ) 
elseif ver < 7,
    mpath = dbstack;
    mpath = mpath(1);
    mpath = mpath.name;
else,
    mpath = dbstack( '-completenames' );
    mpath = mpath(1);
    mpath = mpath.file;
end
if ispc, 
    fs = '\'; 
    ps = ';';
else, 
    fs = '/'; 
    ps = ':';
end
temp = strfind( mpath, fs );
mpath( temp(end) : end ) = [];

rmpaths = { 'lib', 'functions', 'sets', 'structures', 'sedumi', 'doc' };
addpaths = { mpath };
needpaths = {};
if ver < 7,
    addpaths{end+1} = [ mpath, fs, 'keywords' ];
else,
    rmpaths{end+1} = 'keywords';
end
if ver < 6.5,
    addpaths{end+1} = [ mpath, fs, 'matlab6' ];
else,
    rmpaths{end+1} = 'matlab6';
end
needupd = 0;
missing = 0;
newpath = [ ps, matlabpath, ps ]; 
if ispc,
    newpath2 = lower(newpath);
    mpath2 = lower(mpath);
else,
    newpath2 = newpath;
    mpath2 = mpath;
end
for k = 1 : length(rmpaths),
    temp = [ ps, mpath2, fs, rmpaths{k}, ps ];
    ndxs = strfind( newpath2, temp );
    if ~isempty( ndxs ),
        needupd = 1;
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
        newpath2 = [ temp(1:end-1), newpath2 ];
        newpath = [ ps, addpaths{k}, newpath ];
        if ~exist( temp(2:end-1), 'dir' ),
            missing = 1;
        end
    end
end
if missing,
    disp( 'The cvx distribution seems to be incomplete; one or more of the required' );
    disp( 'directories is missing. Please re-unpack the distribution and try again.' );
    disp( ' ' )
    return
end
matlabpath(newpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the SeDuMi MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newext = mexext;
isw32 = strcmp( newext, 'mexw32' );
sedpath = [ mpath, fs, 'sedumi' ];
mexfiles = dir( [ sedpath, fs, '*.', newext ] );
if isempty( mexfiles ) & ispc & isw32,
    mexfiles = dir( [ sedpath, fs, '*.dll' ] );
end
if isempty( mexfiles ),
    try,
        cd( sedpath );
    catch,
        disp( 'The cvx distribution seems to be incomplete: the lib/ directory is' );
        disp( 'missing. Please re-unpack the distribution and try again.' );
        disp( ' ' );
        cd( dd );
        return
    end
    try,
        disp( 'Running the SeDuMi MEX compilation script in 5 seconds.' );
        pause(5);
        install_sedumi;
    catch,
        disp( '-------------------------------------------------------------' );
        disp( 'SeDuMi was NOT built successfully. Please try CVX on a supported' );
        disp( 'platform, or manually run the ''install_sedumi'' command in the' );
        disp( 'sedumi/ subdirectory to try and find and correct the error.' );
    end
    disp( '-------------------------------------------------------------' );
    cd( dd );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile the CVX MEX files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

libpath = [ mpath, fs, 'lib' ];
try,
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
catch,
    disp( 'The cvx distribution seems to be incomplete: the lib/ directory is' );
    disp( 'missing. Please re-unpack the distribution and try again.' );
    disp( ' ' );
    cd( dd );
    return
end
if ~has_mex,
    disp( 'Attempting to generate the CVX MEX files...' );
    disp( '-------------------------------------------------------------' );
    has_mex = 1;
    for k = 1 : length( mexfiles ),
        str = mexfiles{k};
        try,
            disp( [ 'mex -O ', str(1:end-1), newext ] );
            mex( '-O', str );
        catch,
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
    else,
        disp( '-------------------------------------------------------------' );
    end    
end
cd( dd );

disp( 'Testing the cvx distribution. If this script aborts with' );
disp( 'an error, please report the error to the authors.' );
disp( '-------------------------------------------------------------' );

m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);
xls = A \ b;
s = cvx_quiet( 1 );
cvx_begin
    variable('x(n)');
    minimize( norm(A*x-b) );
cvx_end
cvx_quiet( s );
try,
    cvx_clearpath;
end

if norm( x - xls ) > 0.01 * norm( x ),
    disp( '-------------------------------------------------------------' );
    disp( sprintf( 'cvx differs from native Matlab by %g%%', 100 * err ) );
    disp( 'Unexpected numerical errors were found when solving the test problem.' );
    disp( 'Please report this to the authors.' );
    disp( ' ' );
    return
else,
    disp( 'No errors! cvx has been successfully installed.' );
    disp( ' ' );
end

if needupd,
    disp( 'NOTE: The MATLAB path has been updated to point to the cvx distribution.' );
    disp( 'In order to use cvx regularly, you must save this new path definition.' );
    if ispc,
        disp( 'To accomplish this, type the command' );
        if ver < 7,
            disp( '    pathtool' );
            disp( 'at the MATLAB prompt, which brings up the ''Set Path'' dialog. Press' );
            disp( 'the ''Save'' button, and then the ''Close'' button.' );
        else,
            disp( '    savepath' );
            disp( 'at the MATLAB prompt.' );
        end
    else,
        disp( 'To accomplish this, add these lines to your startup.m file:' );
        for k = 1 : length( addpaths ),
            disp( [ '    addpath ', addpaths{k} ] );
        end
        disp( 'Consult the MATLAB documentation if necessary.' );
    end
    disp( ' ' );
end
