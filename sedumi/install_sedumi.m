function install_sedumi( varargin )

%SeDuMi installation script
%
% Heavy rewrite by Michael C. Grant for SeDuMi 1.34
% Now detects if binaries are already 
%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA

need_rebuild = any( strcmp( varargin, '-rebuild' ) );

targets64={...
    'bwblkslv.c sdmauxFill.c sdmauxRdot.c',...
    'choltmpsiz.c',...
    'cholsplit.c',...
    'dpr1fact.c auxfwdpr1.c sdmauxCone.c  sdmauxCmp.c sdmauxFill.c sdmauxScalarmul.c sdmauxRdot.c blkaux.c',...
    'symfctmex.c symfct.c',...
    'ordmmdmex.c ordmmd.c',...
    'quadadd.c',...
    'sqrtinv.c sdmauxCone.c',...
    'givensrot.c auxgivens.c sdmauxCone.c',...
    'urotorder.c auxgivens.c sdmauxCone.c sdmauxTriu.c sdmauxRdot.c',...
    'psdframeit.c reflect.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c sdmauxScalarmul.c',...
    'psdinvjmul.c reflect.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c sdmauxScalarmul.c blkaux.c',...
    'bwdpr1.c sdmauxCone.c sdmauxRdot.c',...
    'fwdpr1.c auxfwdpr1.c sdmauxCone.c sdmauxScalarmul.c',...
    'fwblkslv.c sdmauxScalarmul.c',...
    'qblkmul.c sdmauxScalarmul.c',...
    'blkchol.c blkchol2.c sdmauxFill.c sdmauxScalarmul.c',...
    'vecsym.c sdmauxCone.c',...
    'qrK.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c',...
    'finsymbden.c sdmauxCmp.c',...
    'symbfwblk.c',...
    'whichcpx.c sdmauxCone.c',...
    'ddot.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c',...
    'makereal.c sdmauxCone.c sdmauxCmp.c',...
    'partitA.c sdmauxCmp.c',...
    'getada1.c sdmauxFill.c',...
    'getada2.c sdmauxCone.c sdmauxRdot.c sdmauxFill.c',...
    'getada3.c spscale.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c sdmauxCmp.c',...
    'adendotd.c sdmauxCone.c',...
    'adenscale.c',...
    'extractA.c',...
    'vectril.c sdmauxCone.c sdmauxCmp.c',...
    'qreshape.c sdmauxCone.c sdmauxCmp.c',...
    'sortnnz.c sdmauxCmp.c',...
    'iswnbr.c',...
    'incorder.c',...
    'findblks.c sdmauxCone.c sdmauxCmp.c',...
    'invcholfac.c triuaux.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c sdmauxScalarmul.c blkaux.c',...
    };

fs = filesep;
mpath = mfilename('fullpath');
mpath = mpath( 1 : max([1,strfind(mpath,fs)]) - 1 );
ISOCTAVE = exist('OCTAVE_VERSION','builtin');
if ISOCTAVE, prog = 'Octave'; else prog = 'Matlab'; end
COMPUTER = computer;
mext = mexext;

%
% We don't want to rebuild the binaries if they're already present, unless
% the user has specifically asked for a rebuild. Frankly, we can't
% guarantee that rebuilding will work.
%

if ISOCTAVE,
    page_output_immediately( true, 'local' );
end

line = '---------------------------------------------------------------------------';
fprintf( '\n%s\nSeDuMi installation script\n   Directory: %s\n   %s %s on %s\n%s\n', ...
    line, mpath, prog, version, COMPUTER, line );

if ~need_rebuild,
    fprintf( 'Looking for existing binaries...' );
    mdir = '';
    if ISOCTAVE,
        if ispc,
            mdir = 'o_win32';
        elseif ismac,
            mdir = 'o_mac32';
        elseif isunix && any( strfind( computer, 'linux' ) ),
            mdir = 'o_lin32';
        end
        if ~isempty(mdir) && strncmpi( COMPUTER, 'x86_64', 6 ),
            mdir(end-1:end) = '64';
        end
        if ~exist( [ mpath, fs, mdir ], 'dir' ),
            mdir = '';
        end
    end
    nfound = [ 0, 0 ];
    for k = 1 : length(targets64),
        targ = targets64{k};
        targ = [ targ(1:min(strfind(targ,'.'))), mext ];
        if exist( [ mpath, fs, targ ], 'file' ),
            nfound(1) = nfound(1) + 1;
        elseif ~isempty(mdir) && exist( [ mpath, fs, mdir, targ ], 'file' ),
            nfound(2) = nfound(2) + 1;
        end
    end
    if sum(nfound) == 0,
        fprintf( 'none found; building...\n' );
        need_rebuild = true;
    elseif sum(nfound) < length(targets64),
        fprintf( 'incomplete set found.\n' );
        disp( line );
        disp( 'Some of the binaries for this platform were found, but some' );
        disp( 'were missing as well. This may mean your download was corrupt;' );
        disp( 'consider downloading and unpacking SeDuMi again. Otherwise, to' );
        disp( 'try rebuilding the MEX files yourself, run this command:' );
        disp( '    install_sedumi -rebuild' );
        fprintf( '%s\n\n', line );
        return;
    else
        fprintf( 'found!\n' );
        fprintf( '   If for some reason you need to rebuild the binaries, use this command:\n' );
        fprintf( '      install_sedumi -rebuild\n' );
    end
else
    nfound = [1,0];
end

if need_rebuild,
    disp( 'Attempting to recompile the SeDuMi binaries:' );
    % Note the use of 0.01 here. That's because version 7 had more than 10
    % minor releases, so 7.10-7.14 need to be ordered after 7.01-7.09.
    VERSION  = [1,0.01]*sscanf(version,'%d.%d');
    IS64BIT  = ~ISOCTAVE & strcmp(COMPUTER(end-1:end),'64');
    if ispc,
        flags = {'-DPC'};
    elseif isunix,
        flags = {'-DUNIX'};
    end
    libs = {};
    if ISOCTAVE,
        % Octave has mwSize and mwIndex hardcoded in mex.h as ints.
        % There is no definition for mwSignedIndex so include it here.  
        % This means that Octave ignores the -largeArrayDims flag.
        if VERSION < 3.08,
            flags{end+1} = '-DmwSignedIndex=int';
        end
        libs{end+1} = '-lblas';
    else
        flags{end+1} = '-O';
        if IS64BIT && ( VERSION >= 7.03 ),
            flags{end+1} = '-largeArrayDims';
        elseif VERSION < 7.03,
            flags{end+1} = '-DmwIndex=int';
            flags{end+1} = '-DmwSize=int';
            flags{end+1} = '-DmwSignedIndex=int';
        end
        if VERSION >= 7 && ispc,
            if IS64BIT, dirval = 'win64'; else dirval = 'win32'; end
            libdir = [ matlabroot, fs, 'external', fs, 'lib', fs, dirval, fs ];
            if exist( [ libdir, 'microsoft' ], 'dir' ),
                libdir = [ libdir, 'microsoft' ];
                found = true;
            elseif exist( [ libdir, 'msvc60' ], 'dir' ),
                libdir = [ libdir, 'msvc60' ];
                found = true;
            elseif exist( [ libdir, 'lcc' ], 'dir' ),
                libdir = [ libdir, 'lcc' ];
                found = true;
            end
            if found,
                libs{end+1} = [ '-L"', libdir, '"' ];
            end
        end
        if VERSION >= 7.05,
            libs{end+1} = '-lmwblas';
        else
            libs{end+1} = '-lmwlapack';
        end
    end
    libs = sprintf( ' %s', libs{:} );
    flags = sprintf( ' %s', flags{:} );
    olddir = pwd;
    cd( mpath );
    failed = false;
    fprintf( 'Template: mex%s <sources>%s\n', flags, libs );
    for i=1:length(targets64),
        targ = targets64{i};
        mfile = [ targ(1:min(strfind(targ,'.'))), mext ];
        temp = [ 'mex ', flags, ' ', targets64{i}, libs ];
        fprintf( '   %s: %s\n', mfile, targ );
        eval( temp, 'failed=true;' ); %#ok
    end
    cd( olddir );
    if failed,
        fprintf( 'At least one compilation failure occurred.\n' );
        nfound = [0,0];
    else
        fprintf( 'Compilation successful.\n' );
        nfound = [1,0];
    end
end

if any(nfound),
    disp( line );
    fprintf( 'Adding SeDuMi to the %s path:\n', prog );
    ps = pathsep;
    pp = [ ps, path, ps ];
    already = true;
    fprintf( '   Base directory...' );
    if ~any(strfind(pp,[ps,mpath,ps])),
        already = false;
        pp = [ pp, mpath, ps ];
        fprintf( 'added.\n' );
    else
        fprintf( 'already there.\n' );
    end
    if nfound(2),
        fprintf( '   Binaries directory...' );
        if nfound(2) && ~any(strfind(pp,[ps,mpath,fs,mdir,ps])),
            already = false;
            pp = [ pp, mpath, fs, mdir, ps ];
            fprintf( 'added.\n' );
        else
            fprintf( 'already there.\n' );
        end
    end
    fprintf( '   Conversion directory...' );
    if ~any(strfind(pp,[ps,mpath,fs,'conversion',ps])),
        already = false;
        pp = [ pp, mpath, fs, 'conversion', ps ];
        fprintf( 'added.\n' );
    else
        fprintf( 'already there.\n' );
    end
    fprintf( '   Examples directory...'  );
    if ~any(strfind(pp,[ps,mpath,fs,'examples',ps])),
        already = false;
        pp = [ pp, mpath, fs, 'examples', ps ];
        fprintf( 'added.\n' );
    else
        fprintf( 'already there.\n' );
    end
    if ~already,
        path(pp);
        fprintf( 'Please save the %s path if you want to use SeDuMi from any directory.\n', prog );
    end
    disp( line );
    disp('SeDuMi has been succesfully installed.' );
    disp( 'For more information, type "help sedumi" or see the user guide.')
else
    disp( line );
    disp( 'SeDuMi was not successfully installed.' );
    disp( 'Please attempt to correct the errors and try again.' );
end

fprintf('%s\n\n',line);

