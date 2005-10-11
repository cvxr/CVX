% z = veccomplex(x,cpx,K)
% **********  INTERNAL FUNCTION OF SEDUMI **********

function z = veccomplex(x,cpx,K)

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

z = zeros(length(x) - cpx.dim,1);
dimflqr = K.f+K.l+sum(K.q)+sum(K.r);
rsel = 1:dimflqr;
nfx = length(cpx.f)+length(cpx.x);
imsel = 1:nfx;
imsel = imsel' + [zeros(length(cpx.f),1); vec(cpx.x)];
rsel(imsel) = 0;
rsel = find(rsel);
z(1:dimflqr-nfx) = x(rsel);
z(cpx.f) = z(cpx.f) + sqrt(-1) * x(imsel(1:length(cpx.f)));
z(cpx.x) = z(cpx.x) + sqrt(-1) * x(imsel(length(cpx.f)+1:end));
% ----------------------------------------
% PSD:
% ----------------------------------------
hfirstk = 1+dimflqr + sum(K.s(1:K.rsdpN).^2);     % start of Hermitian blocks
zfirstk = 1+dimflqr - nfx;           % start of psd blocks in z
firstk = 1+dimflqr;                  % start if psd blocks in x
k = 1; sk = 1;
for knz = 1:length(cpx.s)
    newk = cpx.s(knz);
    if newk > k
        len = sum(K.s(sk:sk+newk-k-1).^2);
        z(zfirstk:zfirstk+len-1) = x(firstk:firstk+len-1);  % copy real PSD blocks
        zfirstk = zfirstk + len;
        firstk = firstk + len;
        sk = sk + newk - k;                           % handled newk-k real blocks
    end
    nksqr = K.s(K.rsdpN + knz)^2;                   % insert a Hermitian block
    z(zfirstk:zfirstk+nksqr-1) = x(hfirstk:hfirstk+nksqr-1) ...
        + i * x(hfirstk+nksqr: hfirstk+2*nksqr-1);
    zfirstk = zfirstk + nksqr;
    hfirstk = hfirstk + 2*nksqr;
    k = newk + 1;                                   % handled up to block newk.
end
if sk <= K.rsdpN                                  % copy remaining real blocks
    len = sum(K.s(sk:K.rsdpN).^2);
    z(zfirstk:zfirstk+len-1) = x(firstk:firstk+len-1);
end