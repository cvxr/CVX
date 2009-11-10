% PRELP  Loads and preprocesses LP from an MPS file.
%
% > [A,b,c,lenx,lbounds] = PRELP('problemname')
%    The above command results in an LP in standard form,
%    - Instead of specifying the problemname, you can also use PRELP([]), to
%    get the problem from the file /tmp/default.mat.
%    - Also, you may type PRELP without any input arguments, and get prompted
%    for a name.
%
% MINIMIZE  c'*x SUCH THAT  A*x = b AND  x>= 0.
%
%   So, you can solve it with SeDuMi:
% > [x,y,info] = SEDUMI(A,b,c);
%
% After solving, post-process it with
% > [x,objp] = POSTPROCESS(x(1:lenx),lbounds).
%
% REMARK  x(lenx+1:length(x)) will contain upper-bound slacks.
%
% IMPORTANT works only if LIPSOL is installed on your system.
%
% See also sedumi, getproblem, postprocess (LIPSOL), frompack, lipsol.

function [A,b,c,lenx,lbounds,times] = prelp(pname)
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
% 02110-1301, USA%


global OUTFID
global Ubounds_exist

if ~exist('loadata','file') || ~exist('preprocess','file')
    error('To use PRELP, you need to have LIPSOL installed.')
end

%--------------------------------------------------
% LOAD LP PROBLEM INTO MEMORY
%--------------------------------------------------
if (nargin == 0)
    pname = input('Enter problem name: ','s');
end

t0 = cputime;
[A,b,c,lbounds,ubounds,BIG] = loadata(pname);
times(1) = cputime - t0;

%--------------------------------------------------
% PREPROCESS LP PROBLEM
% NB: Y.Zhang's preprocess returns lbounds for post-
% processing; the pre-processed problem has x>=0.
%--------------------------------------------------
t0 = cputime;
[A,b,c,lbounds,ubounds,FEASIBLE] = ...
    preprocess(A,b,c,lbounds,ubounds,BIG);
if ~FEASIBLE
    fprintf('\n');
    if isempty(OUTFID)
        return;
    end;
    msginf =  'Infeasibility detected in preprocessing';
    fprintf(OUTFID, [pname '  0   ' msginf '\n']);
    return;
end;
%[A,b,c,ubounds] = scaling(A,b,c,ubounds);
%--------------------------------------------------
% INSERT UBOUND- CONSTRAINTS IN THE A-MATRIX
%--------------------------------------------------
b = full(b); c = full(c);
[m,lenx] = size(A);
if Ubounds_exist
    nub = nnz(ubounds);
    A= [ A sparse(m,nub); sparse(1:nub,find(ubounds),1,nub,lenx) speye(nub) ];
    b = [b; nonzeros(ubounds)];
    c = [c; zeros(nub,1)];
end
%--------------------------------------------------
% LOCATE DENSE COLUMNS
%--------------------------------------------------
%checkdense(A);
times(2) = cputime - t0;
