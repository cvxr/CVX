%                          [y, ddotx, Dx, xTy] = PopK(d,x,K,lpq)
% POPK  Implements the quadratic operator for symmetric cones K.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [y, ddotx, Dx, xTy] = PopK(d,x,K,lpq)
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
if nargin < 4
    lpq = 0;
end
% LP / Lorentz
i1 = K.mainblks(1); 
i2 = K.mainblks(2);
y = [d.l .* x(1:i1-1); -d.det .* x(i1:i2-1); qblkmul(d.det,x,K.qblkstart)];
ddotx = (d.q1).*x(i1:i2-1) + ddot(d.q2,x,K.qblkstart);
% PSD:
Dx = psdscale(d,x,K);
if lpq == 0
    y = [y; psdscale(d,Dx,K,1)];               % Include PSD-part
end
% xTy:
if nargout >= 4
    xTy = x(1:K.lq)'*y(1:K.lq) + sum(ddotx.^2) + sum(Dx.^2);
end