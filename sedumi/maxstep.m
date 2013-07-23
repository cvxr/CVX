function tp = maxstep(dx,x,auxx,K)
% tp = maxstep(dx,x,auxx,K)
%
% MAXSTEP  Computes maximal step length to the boundary of the cone K.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

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
% ------------------------------------------------------------
% LP:
% ------------------------------------------------------------
mindx = min(dx(1:K.l) ./ x(1:K.l));
% ------------------------------------------------------------
% Lorentz: compute min(eig(y)) with y:=D(x)\dx. We have
% lab1+lab2 = tr y = x'Jdx / detx,  and
% (lab2-lab1)^2 = (tr y)^2 - 4 det y = [(x'Jdx)^2 - tdetx*tdetdx] / detx^2
% ------------------------------------------------------------
if ~isempty(K.q)
    ix = K.mainblks;
    reltr = x(ix(1):ix(2)-1).*dx(ix(1):ix(2)-1)...
        - ddot(x(ix(2):ix(3)-1),dx,K.qblkstart);
    norm2 = reltr.^2 - tdet(dx,K).*auxx.tdet;
    if norm2 > 0
        norm2 = sqrt(norm2);
    end
    mindxq = min( (reltr - norm2)./auxx.tdet);
    mindx = min(mindx, mindxq);
end
% ------------------------------------------------------------
% PSD:
% ------------------------------------------------------------
if ~isempty(K.s)
    reldx = psdinvscale(auxx.u, dx,K);
    mindxs = minpsdeig(reldx, K);
    mindx = min(mindx, mindxs);
end
tp = 1 / max(-mindx, 1E-16);