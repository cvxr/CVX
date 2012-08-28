function [ nver, isoctave, fs, ps, mpath, problem ] = cvx_version( license_file )

% CVX_VERSION   Returns version and environment information for CVX.
%
%    When called with no arguments, CVX_VERSION prints out version and
%    platform information that is needed when submitting CVX bug reports.
%
%    This function is also used internally to return useful variables that
%    allows CVX to adjust its settings to the current environment.

% Numeric version
global cvx___
nver = version;
nver(nver=='.') = ' ';
nver = sscanf(nver,'%d');
nver = nver(1) + 0.01 * ( nver(2) + 0.01 * nver(3) );

% Matlab / Octave flag
isoctave = exist( 'OCTAVE_VERSION', 'var' );

% File and path separators
if ispc, 
    fs = '\'; 
    ps = ';'; 
else
    fs = '/'; 
    ps = ':';
end
if nargout <= 4 && nargout, 
    return; 
end
    
% Install location
mpath = mfilename('fullpath');
temp = strfind( mpath, fs );
mpath = mpath( 1 : temp(end) - 1 );

if nargout <= 5 && nargout, return; end

% Version and build
cvx_ver   = '2.0 (beta)';
cvx_bld   = ':::BUILD:::';
cvx_bdate = ':::BDATE:::';
cvx_ddate = ':::DDATE:::';
cvx_dbld  = ':::DBUILD:::';
line = '---------------------------------------------------------------------------';
fprintf( '\n' );
disp( line );
fprintf( 'CVX, version %-13s                     (c) 2012, CVX Research, Inc.\n', cvx_ver );
fprintf( 'Software for Disciplined Convex Programming\n' );
disp( line );
fprintf( 'Version info:\n' );
fprintf( '    Code: build %s, %s\n', cvx_bld, cvx_bdate );
fprintf( '    Documentation: build %s, %s\n', cvx_dbld, cvx_ddate );
fprintf( 'Installation info:\n    Path: %s\n', mpath );
if usejava('jvm'),
    os_name = char(java.lang.System.getProperty('os.name'));
    os_arch = char(java.lang.System.getProperty('os.arch'));
    os_version = char(java.lang.System.getProperty('os.version'));
    fprintf('    OS: %s %s ver. %s\n', os_name, os_arch, os_version );
end
if isoctave,
    fprintf( '    GNU Octave %s on %s\n', version, computer );
else
    verd = ver('MATLAB');
    fprintf( '    MATLAB version %s %s\n', verd.Version, verd.Release );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for valid version %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem = true;
if isoctave,
    fprintf( 'Sorry, CVX does not yet run under octave.\n' );
elseif nver < 7.05,
    fprintf( 'CVX requires MATLAB 7.5 or later.\n' );
else
    problem = false;
end
if nargout == 0,
    clear nver
end
if problem,
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify file and directory lists from the manifest %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen( [ mpath, fs, 'MANIFEST' ], 'r' );
if fid > 0,
    fprintf( 'Verfying CVX directory contents:' );
    manifest = textscan( fid, '%s' );
    manifest = manifest{1};
    fclose( fid );
    newman = get_manifest( mpath, fs );
    if ~isequal( manifest(:), newman(:) ),
        missing = setdiff( manifest, newman );
        additional = setdiff( newman, manifest );
        if ~isempty( missing ) || ~isempty( additional ),
            if ~isempty( missing ),
                fprintf( '\n    WARNING: The following files/directories are missing:\n' );
                isdir = cellfun(@(x)x(end)==fs,missing);
                missing_d = missing(isdir);
                missing_f = missing(~isdir);
                while ~isempty( missing_d ),
                    mdir = missing_d{1};
                    ss = strncmp( missing_d, mdir, length(mdir) );
                    tt = strncmp( missing_f, mdir, length(mdir) );
                    fprintf( '        %s%s%s + %d files, %d directories\n', mpath, fs, mdir, nnz(tt), nnz(ss) );
                    missing_d(ss) = [];
                    missing_f(tt) = [];
                end
                for k = 1 : min(length(missing_f),10),
                    fprintf( '        %s%s%s\n', mpath, fs, missing_f{k} );
                end
                if length(missing_f) > 10,
                    fprintf( '        (and %d more files)\n', length(missing_f) - 10 );
                end
                disp( '    These missing may prevent CVX from operating properly.'  );
            end
            if ~isempty( additional ),
                if isempty( missing ), fprintf( '\n' ); end
                disp( '    WARNING: The following extra files/directories were found:' );
                isdir = cellfun(@(x)x(end)==fs,additional);
                additional_d = additional(isdir);
                additional_f = additional(~isdir);
                while ~isempty( additional_d ),
                    mdir = additional_d{1};
                    ss = strncmp( additional_d, mdir, length(mdir) );
                    tt = strncmp( additional_f, mdir, length(mdir) );
                    fprintf( '        %s%s%s + %d files, %d directories\n', mpath, fs, mdir, nnz(tt), nnz(ss) );
                    additional_d(ss) = [];
                    additional_f(tt) = [];
                end
                for k = 1 : min(length(additional_f),10),
                    disp( [ '        ', mpath, fs, additional_f{k} ]  );
                end
                if length(additional_f) > 10,
                    fprintf( '        (and %d more files)\n', length(additional_f) - 10 );
                end
                disp( '    These files may alter the behavior of CVX in unsupported ways.' );
            end
        else
            fprintf( '\n    All files found.\n' );
        end
    else
        fprintf( '\n    All files found.\n' );
    end
elseif 1,
    fid = fopen( [ mpath, fs, 'MANIFEST' ], 'w' ); %#ok
    if fid,
        newman = get_manifest( mpath, fs );
        fprintf( fid, '%s\n', newman{:} );
        fclose( fid );
    end
else    
    fprintf( 'Manifest missing; cannot verify file structure.\n' ) ;
end

%%%%%%%%%%%%%%%%
% License file %
%%%%%%%%%%%%%%%%

try
    if nargin < 1 || isempty( license_file ),
        cvx___.license = cvx_license;
    elseif isequal( license_file, '*clear*' ),
        fprintf( 'Clearing license information.\n' );
        cvx___.license = [];
    else
        cvx___.license = cvx_license( license_file );
    end
catch %#ok
    cvx___.license = [];
end
disp( line );
if nargout == 0,
    fprintf( '\n' );
end

function newman = get_manifest( mpath, fs )
dirs   = {};
files  = {};
nfiles = dir( mpath );
ndir   = '';
dndx   = 0;
pat    = '^\.|~$';
while true,
    isdir  = [ nfiles.isdir ];
    nfiles = { nfiles.name };
    tt     = cellfun( @isempty, regexp( nfiles, pat ) );
    isdir  = isdir(tt);
    nfiles = nfiles(tt);
    ndirs  = nfiles(isdir);
    if ~isempty(ndirs),
        dirs = horzcat( dirs, strcat(strcat(ndir,sort(ndirs)),fs) ); %#ok
    end
    nfiles = nfiles(~isdir);
    if ~isempty(nfiles),
        files = horzcat( files, strcat(ndir,sort(nfiles)) ); %#ok
    end
    if length( dirs ) == dndx, break; end
    dndx = dndx + 1;
    ndir = dirs{dndx};
    nfiles = dir( [ mpath, fs, ndir ] );
end
newman = horzcat( dirs, files );
if ~isequal( fs, '/' ),
    newman = strrep( newman, fs, '/' );
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

