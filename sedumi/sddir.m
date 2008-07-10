%                                 [dx,dy,dz,dy0, err] = sddir(L,Lden,Lsd,p,...
%                                     d,v,vfrm,At,DAt,dense, R,K,y,y0,b, pars)
% SDDIR  Direction decomposition for Ye-Todd-Mizuno self-dual embedding.
%   Here, p is the direction p = dx+dz. If p=[], then assume p=-v,
%   the "affine scaling" direction.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [dx,dy,dz,dy0, err] = sddir(L,Lden,Lsd,pv,...
    d,v,vfrm,At,DAt,dense, R,K,y,y0,b, pars,pMode)
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

% ------------------------------------------------
% dy0 = v'*pv,    ADp = A*D(d)*pv
% ------------------------------------------------
switch pMode,
    case 1,
        % Spectral values w.r.t. vfrm
        dy0 = (vfrm.lab'*pv) / R.b0;
        pv = frameit(pv,vfrm.q,vfrm.s,K);       % Let p = VFRM * p.
    case 2,
        % Affine scaling
        dy0 = -y0;
        pv = -v;
    case 3,
        dy0 = (v'*pv) / R.b0;
end
% ------------------------------------------------------------
% Solve  AP(d)A' * ypv = dy0 * Rb + AD(dy0*DRc - pv)
% and let xpv = (dy0*DRc - pv) - DA'*ypv.
% This yields AD*xpv = err.b-dy0*Rb.
% ------------------------------------------------------------
[dy,dx,err.kcg,err.b] = wrapPcg(L,Lden,At,dense,d, DAt,K,...
    dy0*R.b, dy0*Lsd.DRc - pv, pars.cg,min(1,y0) * R.maxRb);
% ------------------------------------------------------------
% Solve denom * rdx0 = y0*(DRc'*xpv+Rb'*ypv)-(err.b'*y)
% Here, rdx0 = deltax0/x0.
% ------------------------------------------------------------
rdx0 = (y0*(Lsd.DRc'*dx+R.b'*dy) - err.b'*y) / Lsd.denom;
% ------------------------------------------------------------
% Set dy -= rdx0 * ysd
%     dx = rdx0 * xsd - dx
% error dRb := rdx0 * Lsd.b - err.b
% ------------------------------------------------------------
dy = dy - rdx0 * Lsd.y;
dx = rdx0 * Lsd.x - dx;             % Now, AD*dx = deltax0*b + dy0*Rb...
err.b = rdx0 * Lsd.b - err.b;       % ... + err.b
err.maxb = norm(err.b,inf);
dx(1) = rdx0 * v(1);
dz = pv - dx;
