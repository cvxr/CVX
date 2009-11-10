%                                     [L.L, L.d, L.skip, L.add] = blkchol(L,X,pars,absd)
% BLKCHOL Fast block sparse Cholesky factorization.
%    The sparse Cholesky factor will be placed in the fields L.L, L.d;
%    the symbolic factorization fields remain unchanged.
%    On input, L should be the symbolic factorization structure of X,
%    as created by SYMBCHOL.
%    Performs LDL' factorization of X(L.perm,L.perm).
%
%   There are important differences with CHOL(X(L.perm,L.perm))':
%
%   -  BLKCHOL uses the supernodal partition L.XSUPER,
%    to use dense linear algebra on dense subblocks.
%    This explains the performance benefit of BLKCHOL over CHOL.
%
%   -  BLKCHOL never fails.
%
%   - To solve "X*y = b", use
%    >> [L.L,L.d,L.skip,L.add] = blkchol(symbchol(X),X);
%    >> L.d(find(L.skip)) = inf;
%    >> y = sparbwslv(L, sparfwslv(L,b) ./ L.d);
%
% For numerical reasons, "blkchol" may skip unstable pivots,
% or add on the diagonal. Such pivots are then listed in L.{skip,add}.
%
% See also symbchol, sparfwslv, sparbwslv, [symbfact, symmmd, chol].

function [LL, Ld, Lskip, Ladd] = blkchol(L,X,pars,absd) %#ok

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
 %

error('Build the SeDuMi binaries by typing make at the command prompt.');