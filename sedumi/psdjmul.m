function z = psdjmul(x,y, K)
% z = psdmul(x,y, K)
%
% PSDMUL  for full x,y. Computes (XY+YX)/2
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

Ks = K.s;
if isempty(Ks),
    z = [];
    return
end
Kq = Ks .* Ks;
nr = K.rsdpN;
nc = length(Ks);
N  = K.N - K.sblkstart(1) + 1;
z  = zeros(N,1);
xi = length(x) - N;
zi = 0;
for i = 1 : nc,
    ki = Ks(i);
    qi = Kq(i);
    XX = x(xi+1:xi+qi);
    YY = y(xi+1:xi+qi);
    xi = xi + qi;
    if i > nr,
        XX = XX + 1j * x(xi+1:xi+qi);
        YY = YY + 1j * y(xi+1:xi+qi);
        xi = xi + qi;
    end
    XX = reshape(XX,ki,ki);
    YY = reshape(YY,ki,ki);
    ZZ = XX * YY;
    ZZ = 0.5 * ( ZZ + ZZ' );
    z(zi+1:zi+qi) = real(ZZ);
    zi = zi + qi;
    if i > nr,
        z(zi+1:zi+qi) = imag(ZZ);
        zi = zi + qi;
    end
end
