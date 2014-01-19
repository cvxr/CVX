function [ fs, ps, mpath, mext, nver, isoctave, problem ] = cvx_version( license_file )

% CVX_VERSION   Returns version and environment information for CVX.
%
%    When called with no arguments, CVX_VERSION prints out version and
%    platform information that is needed when submitting CVX bug reports.
%
%    This function is also used internally to return useful variables that
%    allows CVX to adjust its settings to the current environment.

global cvx___

% File and path separators
comp = computer;
if strncmp( comp, 'PC', 2 ), 
    fs = '\'; 
    ps = ';'; 
else
    fs = '/'; 
    ps = ':';
end
    
% Install location
mpath = mfilename('fullpath');
temp = strfind( mpath, fs );
mpath = mpath( 1 : temp(end) - 1 );

% MEX extension

mext = mexext;

if nargout <= 4 && nargout, 
    return; 
end

% Numeric version
nver = version;
nver(nver=='.') = ' ';
nver = sscanf(nver,'%d');
nver = nver(1) + 0.01 * ( nver(2) + 0.01 * nver(3) );

% Matlab / Octave flag
isoctave = exist( 'OCTAVE_VERSION', 'builtin' );

if nargout <= 6 && nargout, 
    return; 
end

% Version and build
java_version = 0;
cvx_ver = '2.0 (beta)';
cvx_bld = '999';
cvx_bdate = '9999-99-99 99:99:99';
cvx_ddate = '9999-99-99 99:99:99';
cvx_dbld = '999';
line = '---------------------------------------------------------------------------';
fprintf( '\n%s\nCVX, version %-13s                     (c) 2012, CVX Research, Inc.\n', line, cvx_ver );
fprintf( 'Software for Disciplined Convex Programming\n%s\n', line );
fprintf( 'Version info:\n' );
fprintf( '    Code: build %s, %s\n', cvx_bld, cvx_bdate );
fprintf( '    Documentation: build %s, %s\n', cvx_dbld, cvx_ddate );
fprintf( 'Installation info:\n    Path: %s\n', mpath );
if isoctave,
    fprintf( '    GNU Octave %s on %s\n', version, computer );
else
    verd = ver('MATLAB');
    fprintf( '    MATLAB version: %s %s\n', verd.Version, verd.Release );
    if usejava( 'jvm'),
        os_name = char(java.lang.System.getProperty('os.name'));
        os_arch = char(java.lang.System.getProperty('os.arch'));
        os_version = char(java.lang.System.getProperty('os.version'));
        java_version = char(java.lang.System.getProperty('java.version'));
        fprintf('    OS: %s %s version %s\n', os_name, os_arch, os_version );
        fprintf('    Java version: %s\n', java_version );
        try
            ndxs = strfind( java_version, '.' );
            java_version = str2double( java_version(1:ndxs(2)-1) );
            if java_version < 1.6,
                fprintf('       WARNING: full support for CVX Professional licenses\n' );
                fprintf('       requres Java version 1.6.0 or later. Please upgrade.\n' );
            end
        catch
        end
    else
        fprintf( '    Architecture: %s\n', computer );
        fprintf( '    Java version: disabled\n' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for valid version %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem = true;
if isoctave && nver < 3.08,
	fprintf( '%s\nCVX requires Octave 3.8 or later.\n%s\n', line, line );
elseif ~isoctave && nver < 7.08 && strcmp( comp(end-1:end), '64' ),
    fprintf( '%s\nCVX requires MATLAB 7.8 or later (7.5 or later on 32-bit platforms).\n' , line, line );
elseif ~isoctave && nver < 7.05,
    fprintf( '%s\nCVX requires MATLAB 7.5 or later (7.8 or later on 64-bit platforms).\n' , line, line );
else
	problem = false;
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
    if ~isequal( manifest, newman ),
        missing = setdiff( manifest, newman );
        additional = setdiff( newman, manifest );
        if ~isempty( missing ) || ~isempty( additional ),
            if ~isequal( fs, '/' ),
                missing = strrep( missing, '/', fs );
                additional = strrep( additional, '/', fs );
            end
            if ~isempty( missing ),
                fprintf( '\n    WARNING: The following files/directories are missing:\n' );
                isdir = cellfun(@(x)x(end)==fs,missing);
                missing_d = missing(isdir);
                missing_f = missing(~isdir);
                while ~isempty( missing_d ),
                    mdir = missing_d{1};
                    ss = strncmp( missing_d, mdir, length(mdir) );
                    tt = strncmp( missing_f, mdir, length(mdir) );
                    fprintf( '        %s%s%s + %d files, %d subdirectories\n', mpath, fs, mdir, nnz(tt), nnz(ss) - 1 );
                    missing_d(ss) = [];
                    missing_f(tt) = [];
                end
                for k = 1 : min(length(missing_f),10),
                    fprintf( '        %s%s%s\n', mpath, fs, missing_f{k} );
                end
                if length(missing_f) > 10,
                    fprintf( '        (and %d more files)\n', length(missing_f) - 10 );
                end
                fprintf( '    These omissions may prevent CVX from operating properly.\n'  );
            end
            if ~isempty( additional ),
                if isempty( missing ), fprintf( '\n' ); end
                fprintf( '    WARNING: The following extra files/directories were found:\n' );
                isdir = cellfun(@(x)x(end)==fs,additional);
                issedumi = cellfun(@any,regexp( additional, [ '^sedumi.*[.]', mexext, '$' ] ));
                additional_d = additional(isdir&~issedumi);
                additional_f = additional(~isdir&~issedumi);
                additional_s = additional(issedumi);
                while ~isempty( additional_d ),
                    mdir = additional_d{1};
                    ss = strncmp( additional_d, mdir, length(mdir) );
                    tt = strncmp( additional_f, mdir, length(mdir) );
                    fprintf( '        %s%s%s + %d files, %d subdirectories\n', mpath, fs, mdir, nnz(tt), nnz(ss) - 1 );
                    additional_d(ss) = [];
                    additional_f(tt) = [];
                end
                for k = 1 : min(length(additional_f),10),
                    fprintf( '        %s%s%s\n', mpath, fs, additional_f{k} );
                end
                if length(additional_f) > 10,
                    fprintf( '        (and %d more files)\n', length(additional_f) - 10 );
                end
                fprintf( '    These files may alter the behavior of CVX in unsupported ways.\n' );
                if ~isempty( additional_s ),
                	fprintf( '    ERROR: obsolete versions of SeDuMi MEX files were found:\n' );
	                for k = 1 : length(additional_s),
	                    fprintf( '        %s%s%s\n', mpath, fs, additional_f{k} );
	            	end
	            	fprintf( '    These files are now obsolete, and must be removed to ensure\n' );
	            	fprintf( '    that SeDuMi operates properly and produces sound results.\n' );
	            	if ~problem,
		            	fprintf( '    Please remove these files and re-run CVX_SETUP.\n' );
		            	problem = true;
		            end
	            end
            end
        else
            fprintf( '\n    No missing files.\n' );
        end
    else
        fprintf( '\n    No missing files.\n' );
    end
else    
    fprintf( 'Manifest missing; cannot verify file structure.\n' ) ;
end
if ( ~exist( [ mpath, fs, 'lib', fs, 'cvx_eliminate_mex.', mext ], 'file' ) || ...
     ~exist( [ mpath, fs, 'lib', fs, 'cvx_bcompress_mex.', mext ], 'file' ) ) && ~problem,
    fprintf( '    ERROR: one or more MEX files for this platform are missing.\n' );
    fprintf( '    These files end in the suffix ".%s". CVX will not operate\n', mext );
    fprintf( '    without these files. Please visit\n' );
    fprintf( '        http://cvxr.com/cvx/download\n' );
    fprintf( '    And download a distribution targeted for your platform.\n' );
    problem = true;
end

%%%%%%%%%%%%%%%%
% License file %
%%%%%%%%%%%%%%%%

cvx___.license = [];
exception = '';
if java_version >= 1.6 && exist( 'cvx_license', 'file' ),
    try
        if nargin < 1, license_file = ''; end
        cvx___.license = cvx_license( license_file );
    catch exception
    end
elseif ~isoctave && exist( 'cvx_license', 'file' ),
    fprintf( 'CVX Professional disabled; requires the Java virtual machine.\n' );
end
fprintf( '%s\n', line );
if nargout == 0,
    clear fs
    fprintf( '\n' );
end
if ~isempty( exception ),
    rethrow( exception )
end

function newman = get_manifest( mpath, fs )
dirs   = {};
files  = {};
nfiles = dir( mpath );
ndir   = '';
dndx   = 0;
pat2   = '^\.|~$|';
pat    = '^\.|~$|^cvx_license.[md]at$|^doc$|^examples$';
while true,
    isdir  = [ nfiles.isdir ];
    nfiles = { nfiles.name };
    tt     = cellfun( @isempty, regexp( nfiles, pat ) ); pat = pat2;
    isdir  = isdir(tt);
    nfiles = nfiles(tt);
    ndirs  = nfiles(isdir);
    if ~isempty(ndirs),
        dirs = horzcat( dirs, strcat(strcat(ndir,ndirs), fs ) ); %#ok
    end
    nfiles = nfiles(~isdir);
    if ~isempty(nfiles),
        files = horzcat( files, strcat(ndir,nfiles) ); %#ok
    end
    if length( dirs ) == dndx, break; end
    dndx = dndx + 1;
    ndir = dirs{dndx};
    nfiles = dir( [ mpath, fs, ndir ] );
end
[tmp,ndxs1] = sort(upper(dirs)); %#ok
[tmp,ndxs2] = sort(upper(files)); %#ok
newman = horzcat( dirs(ndxs1), files(ndxs2) );
if ~isequal( fs, '/' ),
    newman = strrep( newman, fs, '/' );
end
newman = newman(:);

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
