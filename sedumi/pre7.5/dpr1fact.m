%                               [Lden,L.d] = dpr1fact(x, d, Lsym, smult, maxu)
% DPR1FACT  Factor d[iag] p[lus] r[ank] 1:
%    [Lden,L.d] = dpr1fact(x, d, Lsym, smult, maxu)
%    Computes fi and d such that
%       diag(d_IN) + x*diag(smult)*x' =
%(PI_{i=1}^n L(p_OUT^i,beta_i)) * diag(d_OUT) * (PI_{i=1}^n L(p_OUT^i,beta_i))'
%    where L(p,beta) = eye(n) + tril(p*beta',-1).
%
% Lden.dopiv(k) = 1 if p(:,k) has been reordered, with permutation in
% Lden.pivperm.
% We reorder if otherwise |p(i,k)*beta(j,k)| > maxu.
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
%
% See also fwdpr1,bwdpr1,sedumi


function [Lden,Ld] = dpr1fact(x, d, Lsym, smult, maxu)
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
error('At OS prompt, type "make" to create cholTool mex-files.')