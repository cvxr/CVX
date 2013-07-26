function mineig = minpsdeig(x,K)
% mineig = minpsdeig(x,K)
%
% MINPSDEIG  Computes the smallest spectral coefficients of x w.r.t. K
% Uses an iterative method if the matrix is large, takes the minimum of all
% the eigenvalues if the matrix is small.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

% This file is part of SeDuMi 1.3 by Imre Polik
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

% disp('The SeDuMi binaries are not installed.')
% disp('In Matlab, launch "install_sedumi" in the folder you put the SeDuMi files.')
% disp('For more information see the file Install.txt.')
% error(' ')

Ks = K.s;
if isempty(Ks),
    mineig = [];
    return
end
Kq = Ks .* Ks;
nr = K.rsdpN;
nc = length(Ks);
xi = length(x) - sum(Kq) - sum(Kq(nr+1:end));
eigv = zeros(nc,1);
for i = 1 : nc,
    ki = Ks(i);
    qi = Kq(i);
    XX = x(xi+1:xi+qi);
    xi = xi + qi;
    if i > nr,
        XX = XX + 1j*x(xi+1:xi+qi);
        xi = xi + qi;
    end
    XX = reshape(XX,ki,ki);
    XX = XX + XX';
    if ki > 500,
        if nnz(XX) < 0.1 * numel(XX), XX = sparse(XX); end
        [v,eigv(i),flag] = eigs(XX,1,'SA',struct('issym',true)); %#ok
        if flag, eigv(i) = min(eig(XX)); end
    else
        eigv(i) = min(eig(XX));
    end
end
mineig = min(eigv) / 2;
