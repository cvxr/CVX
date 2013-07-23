function Lden = symbcholden(L,dense,DAt)
% Lden = symbcholden(L,dense,DAt)
%
% SYMBCHOLDEN  Creates Lden.{LAD, perm,dz, sign, first}
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
%
% See also sedumi, dpr1fact

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

% ------------------------------------------------------------
% Symbolic forward Cholesky of dense columns, in order
% [LP, Q-blk, Q-norm, Q-tr]
% ------------------------------------------------------------
i1 = dense.l + 1;
i2 = i1 + length(dense.q);
LAD = [symbfwblk(L,dense.A(:,1:i1-1)), symbfwblk(L,DAt.denq),...
    symbfwblk(L,dense.A(:,i2:end)),symbfwblk(L,dense.A(:,i1:i2-1))];
% ------------------------------------------------------------
% Incremental ordering heuristic, excluding the Lorentz-trace cols
% ------------------------------------------------------------
[perm, dz] = incorder(LAD(:,1:length(dense.cols)));
% ------------------------------------------------------------
% Insert the trace cols with a "-1"-factor just after corresponding
% Lorentz-block columns
% ------------------------------------------------------------
Lden = finsymbden(LAD,perm,dz,i1);