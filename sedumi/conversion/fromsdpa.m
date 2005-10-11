% FROMSDPA   Reads SDP problem from sparse SDPA-formatted input file.
%    [At,b,c,K] = fromsdpa(fname) produces conic optimization problem
%    data, as can be used by sedumi. "fname" contains a full pathname
%    to the SDPA 4.10-formatted file. (SDPA is a stand-alone solver for
%    semidefinite programming, by Fujisawa, Kojima and Nakata.  It is
%    used as a standard format in the collection SDPLIB by Borchers.)
%
%    To read and solve the problem "arch0", you may type
%
%    [At,b,c,K, perm] = fromsdpa('arch0.dat-s');
%    [x,y,info] = sedumi(At,b,c,K);
%
%    The above 2 lines assume that arch0.dat-s is somewhere in your MATLAB
%    search path, that it is not compressed, and that you know the extension
%    'dat-s'.  To alleviate these conditions, you may like to use the script
%    GETPROBLEM.
%
% SEE ALSO SeDuMi, getproblem, frompack, prelp.

function [At,b,c,K,perm] = fromsdpa(fname)

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

fid = fopen(fname,'r');
if fid == -1
    error('File not found.')
end
m = fscanf(fid,'%d',1);
while(isempty(m))
    fgetl(fid); m = fscanf(fid,'%d',1);
end
sdpN = fscanf(fid,'%d',1);
Ks = fscanf(fid,'%d',sdpN);
Ks = abs(Ks);
perm = find(Ks == 1);
K.l = length(perm);
perm2 = find(Ks>1);
K.s = Ks(perm2);
perm = [perm; perm2];
if length(perm) ~= sdpN
    error('PSD orders should be positive integers')
end
invperm = zeros(sdpN,1);
invperm(perm) = 1:sdpN;
b = fscanf(fid,'%f',m);
if(isempty(b))
    b1 = fscanf(fid,'{%f,',1);
    if(isempty(b1))
        error('Invalid sdpa file format\n')
    end
    b = [b1; fscanf(fid,'%f,',m-2); fscanf(fid,'%f}',1)];
end
if length(b) ~= m
    error('Invalid sdpa file format\n')
end
% ------------------------------------------------------------
% Get data in "ikrcf" format: A(constraint i, block k) has
%    entry "f" at position (r,c).
% ------------------------------------------------------------
E = fscanf(fid,'%d %d %d %d %f',[5. inf]);
% ------------------------------------------------------------
% Split E into C, A1, A2, ..., Am
% ------------------------------------------------------------
[xir,xjc] = sdpasplit(E(1,:),m+1);
% ------------------------------------------------------------
% Get c-vector. We've to change sign, because SDPA assumes "max"
% instead of "min".
% ------------------------------------------------------------
c = sparse([],[],[],sum(Ks.^2),1);
knz = 1;
if ~isempty(xir)
    if xir(1) == 0
        c = -sdpa2vec( E(2:5,xjc(1):xjc(2)-1),  K, invperm);
        knz = 2;
    end
end
% ------------------------------------------------------------
% Construct each constraint i=1:m
% ------------------------------------------------------------
N = length(c);
At = sparse([],[],[],N,m,2*size(E,2));
for i = knz:length(xir)
    At(:,xir(i)) = sdpa2vec( E(2:5,xjc(i):xjc(i+1)-1),  K, invperm);
end
fclose(fid);
