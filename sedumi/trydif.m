function [t,wr,w] = trydif(t,wrIN,wIN, x,z, pars,K)
% [t,wr,w] = trydif(t,wrIN,wIN, x,z, pars,K)
%
% TRYDIF  Tries feasibility of differentiated step length w.r.t.
%  wide region and its neighborhood.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi, stepdif

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

thetaSQR = pars.theta^2;
ix = K.mainblks;
% --------------------------------------------------
% Let w = D(xM)*zM, compute lambda(w).
% --------------------------------------------------
% LORENTZ
w.tdetx = tdet(x,K);
w.tdetz = tdet(z,K);
detxz = w.tdetx .* w.tdetz / 4;
if isempty(K.q)
    lab2q = zeros(0,1);
else
    halfxz = (x(ix(1):ix(2)-1).*z(ix(1):ix(2)-1)...
        + ddot(x(ix(2):ix(3)-1),z,K.qblkstart)) / 2;
    tmp = halfxz.^2 - detxz;
    if tmp > 0
        lab2q = halfxz + sqrt(tmp);
    else
        lab2q = halfxz;
    end
end
% PSD
w.ux = psdfactor(x,K);
w.s = psdscale(w.ux,z,K);
% ALL:
w.lab = [x(1:K.l).*z(1:K.l); detxz ./ lab2q; lab2q; psdeig(w.s,K)];
[wr.delta,wr.h,wr.alpha] = iswnbr(w.lab, thetaSQR);
wr.desc = wrIN.desc;       % always descent direction.
if wr.delta > pars.beta
    t = 0; 
    w = wIN; 
    wr = wrIN;
end