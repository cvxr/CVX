function install_sedumi

%SeDuMi installation script
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
    'ddot.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c',...
    'partitA.c sdmauxCmp.c',...
    'getada1.c sdmauxFill.c',...
    'getada2.c sdmauxCone.c sdmauxRdot.c sdmauxFill.c',...
    'getada3.c spscale.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c sdmauxCmp.c',...
    'adendotd.c sdmauxCone.c',...
    'adenscale.c',...
    'extractA.c',...
    'sortnnz.c sdmauxCmp.c',...
    'iswnbr.c',...
    'incorder.c',...
    'findblks.c sdmauxCone.c sdmauxCmp.c',...
    'invcholfac.c triuaux.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c sdmauxScalarmul.c blkaux.c',...
    };

disp( 'Building SeDuMi binaries...' )
ISOCTAVE = exist('OCTAVE_VERSION','builtin');
COMPUTER = computer;
% Note the use of 0.01 here. That's because version 7 had more than 10
% minor releases, so 7.10-7.14 need to be ordered after 7.01-7.09.
VERSION  = [1,0.01]*sscanf(version,'%d.%d');
IS64BIT  = ~ISOCTAVE & strcmp(COMPUTER(end-1:end),'64');
mexprog  = 'mex';
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
    flags{end+1} = '-DmwSignedIndex=int';
    libs{end+1} = '-lblas';
else
    if nargin > 1 && ~isempty(endpath),
        flags{end+1} = '-outdir';
        flags{end+1} = endpath;
    end
    flags{end+1} = '-O';
    if IS64BIT && ( VERSION >= 7.03 ),
        flags{end+1} = '-largeArrayDims';
    elseif VERSION < 7.03,
        flags{end+1} = '-DmwIndex=int';
        flags{end+1} = '-DmwSize=int';
        flags{end+1} = '-DmwSignedIndex=int';
    end
    if VERSION >= 7,
        if VERSION >= 7.05, libval = 'blas'; else libval = 'lapack'; end
        if IS64BIT, dirval = 'win64'; else dirval = 'win32'; end
        libdir = [ matlabroot, '\extern\lib\', dirval, '\microsoft' ];
        if exist( [ libdir, '\msvc60' ], 'file' ),
            libdir = [ libdir, '\msvc60' ];
        end
        libs{end+1} = [ '-L"', libdir, '"' ];
        libs{end+1} = [ '-lmw', libval ];
    end
end
libs = sprintf( ' %s', libs{:} );
flags = sprintf( ' %s', flags{:} );
for i=1:length(targets64)
    temp =  [ mexprog, flags, ' ', targets64{i}, libs ];
    disp( temp );
    eval( temp );
end
disp( 'Done!' )
if nargin < 1,
    if ISOCTAVE,
	disp('Adding SeDuMi to the Octave path')
    else
    disp('Adding SeDuMi to the Matlab path')
    end
    path(path,pwd);
    cd conversion
    path(path,pwd);
    cd ..
    cd examples
    path(path,pwd);
    cd ..
    if ISOCTAVE
        disp('Please save the Octave path if you want to use SeDuMi from any directory.'); 
        disp('To do this type savepath at the Octave prompt.');
    else
    disp('Please save the Matlab path if you want to use SeDuMi from any directory.');
    disp('Go to File/Set Path and click on Save.');
    end
    disp('SeDuMi has been succesfully installed. For more information type help sedumi or see the User guide.')
end
