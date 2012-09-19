function prevpath = cvx_startup( quiet )

%CVX_STARTUP   Quietly add CVX to your MATLAB path (for startup).
%    Running CVX_STARTUP upon startup ensures that CVX is properly included
%    in your MATLAB path.
%
%    On Mac and PC systems, this function is not necessary, because you can
%    simply run "savepath" to save the CVX path changes permanently, or use
%    the "pathtool" command to do the same. In fact, CVX_SETUP attempts to
%    run "savepath" as part of its setup process.
%
%    On UNIX systems, the MATLAB installation is usually installed in such
%    a manner that individual users cannot change the global MATLAB path.
%    This can be circumvented by making the PATHDEF.M file world-writable;
%    but if the system serves multiple users, this may not be desirable.
%
%    CVX_STARTUP provides a simple solution to this problem. Simply add the
%    following line to the end of your local STARTUP.M file:
%        run <full_path_to_CVX>/cvx_startup
%    Note that using the "run" command *and* providing the full path to the
%    CVX installation are critical. With this change, CVX_STARTUP.M will be
%    called every time the user runs MATLAB.
%
%    Note that CVX_STARTUP is *not* a substitute for running CVX_SETUP.
%    Please run CVX_SETUP first when installing CVX, and *then* add the
%    CVX_STARTUP line if instructed to do so.

if nargin < 1 || isempty( quiet ), 
    quiet = false; 
end
if ~quiet,
    fprintf( 'Setting CVX paths...' );
end
if strncmp(computer,'PC',2), fs = '\'; ps = ';'; else fs = '/'; ps = ':'; end
mpath = mfilename('fullpath');
temp = strfind( mpath, fs );
mpath = mpath( 1 : temp(end) - 1 );
addpaths = { 'builtins', 'commands', 'functions', 'lib', 'structures' };
addpaths = strcat( [ mpath, fs ], addpaths );
addpaths{end+1} = mpath;
prevpath = matlabpath;
oldpath = textscan( prevpath, '%s', 'Delimiter', ps );
oldpath = oldpath{1}(:)';
other_homes = regexprep( which( 'cvx_setup.m', '-all' ), '.cvx_setup\.m$', '' );
other_homes( strcmp( other_homes, mpath ) ) = [];
ndxs = false(size(oldpath));
for k = 0 : length(other_homes),
    if k, tpath = other_homes{k}; else tpath = mpath; end
    plen = length(tpath);
    tndxs = strncmp( tpath, oldpath, plen );
    tndxs(tndxs) = cellfun(@(x)length(x)<=plen||x(plen+1)==fs,oldpath(tndxs));
    ndxs = ndxs | tndxs;
end
newpath = horzcat( addpaths, oldpath(~ndxs) );
npath = sprintf( [ '%s', pathsep ], newpath{:} );
npath = npath(1:end-1);
matlabpath(npath);
if ~quiet,
    if isequal( newpath, oldpath ),
        fprintf( 'already set!\n' );
    else
        fprintf( 'done.\n' );
    end
    if ~isempty( other_homes ),
        fprintf( 'WARNING: other CVX installations were found in your MATLAB path:\n' );
        fprintf( '    %s\n', other_homes{:} );
        fprintf( 'They have been removed to prevent conflicts.\n' );
    end
end
if nargout,
    if isequal( npath, prevpath ),
        prevpath = '';
    end
else
    clear prevpath
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
