function Apart = extractA(At,Ajc,blk0,blk1,blkstart,blkstart2) %#ok
% Apart = extractA(At,Ajc,blk0,blk1,blkstart[,blkstart2])
%
% EXTRACTA  Fast alternative to
%  Apart = At(blkstart(1):blkstart(2)-1,:).
%  Instead of blkstart(2), it takes "blkstart2" (if supplied) or
%  size(At,1)+1 (if neither blkstart(2) nor blkstart2) are available.
%
%  Extract submatrix of
%  A with subscripts "blkstart(1):blkstart(2)-1" and has its nonzeros per column
%  bounded bij Ajc1:Ajc2. If Ajc1 = [] (Ajc2=[]) then start at column start
%  (end at column end) of A.
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
%
%  See also partitA

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

%Indicate to the user Matlab cannot find the SeDuMi binaries
sedumi_binary_error();