function Lsd = sdfactor(L,Lden, dense,DAt, d,v,y, At,c,K,R,y0,pars) %#ok
% [Lsd,Rscl] = sdfactor(L,Lden, dense,DAt, d,v,y, At,K,R,y0,pars)
%
% SDFACTOR  Factor self-dual embedding
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
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

% ------------------------------------------------------------
% Compute D*Rc
% ------------------------------------------------------------
Lsd.DRc = [sqrt(d.l).*R.c(1:K.l); asmDxq(d,R.c,K); psdscale(d,R.c,K)];
% ------------------------------------------------------------
% Solve AP(d)A' * ysd = y0*Rb + AD(y0*DRc-2v)
% and let xsd = (y0*DRc-2v) - DA'*ysd.
% ------------------------------------------------------------
[Lsd.y,Lsd.x,Lsd.kcg,Lsd.b] = wrapPcg(L,Lden,At,dense,d, DAt,K, y0*R.b,...
    y0*Lsd.DRc - 2*v, pars.cg, min(1,y0)*R.maxRb);
% ------------------------------------------------------------
% Now let ysd -= y and xsd += v.
% It follows then that AD*xsd = x0 * b + Lsd.b, where Lsd.b = numerical error.
% ------------------------------------------------------------
Lsd.y = Lsd.y - y;
Lsd.x = Lsd.x + v;
% ------------------------------------------------------------
% Compute denom = norm(xsd)^2
% ------------------------------------------------------------
Lsd.denom = norm(Lsd.x)^2 + Lsd.b'*Lsd.y;