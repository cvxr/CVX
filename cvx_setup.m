function cvx_setup

% CVX_SETUP   Sets up and tests the cvx distribution.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

disp( ' ' );

ver = str2num( version( '-release' ) );
if ver < 14,
    mpath = dbstack;
    mpath = mpath(1);
    mpath = mpath.name;
else,
    mpath = dbstack( '-completenames' );
    mpath = mpath(1);
    mpath = mpath.file;
end
if ispc, fs = '\'; else, fs = '/'; end
temp = strfind( mpath, fs );
mpath( temp(end) : end ) = [];

addpaths = { mpath };
needpaths = {};
if ver < 14,
    addpaths{end+1} = [ mpath, fs, 'keywords' ];
end
if ver < 13,
    addpaths{end+1} = [ mpath, fs, 'matlab6' ];
end
s = lastwarn;
lastwarn('');
oldpath = path;
addpath( addpaths{:}, '-begin' );
if ~isempty( lastwarn ),
    disp( 'The cvx distribution seems to be incomplete; one or more of the required' );
    disp( 'directories is missing. Please re-unpack the distribution and try again.' );
    disp( ' ' )
    return
end
lastwarn(s);
needupd = strcmp(oldpath,path) == 0;

disp( 'Testing the cvx distribution. Please consult the documentation if an' );
disp( 'error occurs at this point.' );
disp( '----------------------------------------------------------------------' );

m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);

xls = A \ b;

s = cvx_quiet( true );
cvx_begin
    variable('x(n)');
    minimize( norm(A*x-b) );
cvx_end
cvx_quiet( s );
try,
    cvx_clearpath;
end

if norm( x - xls ) > 0.01 * norm( x ),
    disp( sprintf( 'cvx differs from native Matlab by %g%%', 100 * err ) );
    disp( 'Unexpected numerical errors were found when solving the test problem.' );
    disp( 'Please report this to the authors.' );
    disp( ' ' );
    return
end

disp( 'No errors! cvx has been successfully installed.' );
disp( ' ' );

global cvx___
if ~cvx___.has_mex,
    disp( 'WARNING: The MEX file ''lib/cvx_bcompress_mex'' was not found. It' );
    disp( 'is not strictly necessary, but its presence greatly improves the' );
    disp( 'performance of cvx. Consult the cvx documentation for detials on' );
    disp( 'how to create this MEX file for your platform.' );
    disp( ' ' );
end

if strcmp( oldpath, path ) == 0,
    disp( 'NOTE: The MATLAB path has been updated to point to the cvx distribution.' );
    disp( 'In order to use cvx regularly, you must save this new path definition.' );
    if ispc,
        disp( 'To accomplish this, type the command' );
        if ver < 14,
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
