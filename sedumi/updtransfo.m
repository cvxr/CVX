%                                      [d,vfrm] = updtransfo(x,z,w, dIN,K)
% UPDTRANSFO  Updated the Nesterov-Todd transformation using a
%  numerically stable method.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function  [d,vfrm] = updtransfo(x,z,w, dIN,K)
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

% ------------------------------------------------------------
% PSD:
% Given w = D(Xscl)Zscl, compute spec-factor Q*WLAB*Q' = W
% ------------------------------------------------------------
[wlab,q] = psdeig(w.s,K);
w.lab(K.l+2*length(K.q) + 1:end) = wlab;
% ------------------------------------------------------------
% lambda(v) = sqrt(lambda(w))
% ------------------------------------------------------------
vfrm.lab = sqrt(w.lab);
% ------------------------------------------------------------
% LP : d.l = dIN.l .* (x ./ z)
% ------------------------------------------------------------
d.l = dIN.l .* (x(1:K.l) ./ z(1:K.l));
% ------------------------------------------------------------
% Lorentz:
% Auxiliary: s := sqrt(det(x)./det(z))
%   chi = (x[k] + s(k)*Jz[k]) / (lab1(v)+lab2(v))
%   psi = (x[k] - s(k)*Jz[k]) / (lab2(v)-lab1(v))
% Scale: detd = detdIN .* s and d = D(dIN)*chi
% ------------------------------------------------------------
if isempty(K.q)
    d.det = zeros(0,1); d.q1 = zeros(0,1); d.q2 = zeros(0,1);
    d.auxdet = zeros(0,1); d.auxtr = zeros(0,1);
    vfrm.q = zeros(0,1);
else
    i1 = K.mainblks(1); i2 = K.mainblks(2); nq = i2 - i1; j3 = i2+nq-1;
    s = sqrt(w.tdetx ./ w.tdetz);
    d.det = dIN.det .* s;
    psi1 = s.*z(i1:i2-1); psi2 = qblkmul(s,z,K.qblkstart);       %s * z
    tmp = vfrm.lab(i1:i2-1) + vfrm.lab(i2:j3);
    chi1 = (x(i1:i2-1)+psi1)./tmp;
    chi2 = qblkmul(1./tmp, x(i2:K.lq)-psi2,K.qblkstart);
    psi1 = x(i1:i2-1)-psi1;
    psi2 = x(i2:K.lq)+psi2;
    dq = asmDxq(dIN,[chi1;chi2],K);                     % d = D(dIN)*chi
    d.q1 = dq(1:nq);
    d.q2 = dq(nq+1:end);
    d.auxdet = sqrt(2*d.det);
    d.auxtr = sqrt(2)*(d.q1 + d.auxdet);
    alpha = (dIN.q1 .* psi1 + ddot(dIN.q2,psi2,K.qblkstart)) ./ d.auxtr;
    tmp = 2*sqrt(s);
    psi1 = (psi1 - alpha .* chi1)./tmp;
    psi2 = psi2 - qblkmul(alpha,chi2,K.qblkstart);
    psi2 = qblkmul(1./tmp, psi2,K.qblkstart);
    gamma = (sqrt(2)*psi1+alpha) ./ dIN.auxtr;
    tmp = vfrm.lab(i2:j3) - vfrm.lab(i1:i2-1);
    tmp(tmp == 0) = 1;                                %avoid division by zero
    psi2 = psi2 + qblkmul(gamma,dIN.q2,K.qblkstart);
    vfrm.q = qblkmul(1./tmp, psi2, K.qblkstart);
end
% ------------------------------------------------------------
% D= QUD'*QUD, where QUD(:,udIN.perm) = diag(1./sqrt(vlab))*Q'*ux*ud.
% Let vinv = diag(1./sqrt(vlab))*Q'.
% ------------------------------------------------------------'
d.u = triumtriu(w.ux, dIN.u, K);             % ux * ud
[d.u,d.perm,gjc,g] = urotorder(d.u,K, 1.1, dIN.perm);  % stable reordering
q = givensrot(gjc,g,q,K);    % ROTATE Q accordingly: (G*Q)'*(G*ux*ud).'
vinv = sqrtinv(q,vfrm.lab,K);
% ------------------------------------------------------------
% QR-FACTORIZE: Qv * VINV = R
% Then the new ud is simply R*ux*ud, which is upper triangular.
% ------------------------------------------------------------
[vfrm.s, r] = qrK(vinv,K);
d.u = triumtriu(r, d.u, K);