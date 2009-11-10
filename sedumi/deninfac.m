%            [Lden, Ld] = deninfac(symLden, L,dense,DAt, d, absd, qblkstart,pars)
% DENINFAC
% Uses pars.maxuden as max. allowable |L(i,k)|. Otherwise num. reordering.
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
%
% See also sedumi

function [Lden, Ld] = deninfac(symLden, L,dense,DAt, d, absd, qblkstart,pars)
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
% INITIALIZE: check the numerical tolerance parameters in "pars",
%   and supply with defaults if necessary.
% ------------------------------------------------------------
if nargin < 7
    pars.maxuden = 10;                  % Stability of product-form factors
else
    if ~isfield(pars,'maxuden')
        pars.maxuden = 10;                % Stability of product-form factors
    end
end
% --------------------------------------------------
% Scaling: let Ad = Aden * Dden, as [LP, QBLK, QNORM, QTR]
% with QBLK=DAt.denq, and the others dense.A(:,j)*dj,
% with dj = sqrt(d.l(k)) or dj=sqrt(d.det(k)).
% NOTE that A*P(d)*A' = ADA + Ad * Ad'.
% --------------------------------------------------
if ~isempty(dense.cols)
    i1 = dense.l+1; i2 = i1 + length(dense.q);
    Ad = [dense.A(:,1:i1-1), DAt.denq, dense.A(:,i2:end), dense.A(:,i1:i2-1)];
    smult = [d.l(dense.cols(1:dense.l)); ones(length(dense.q),1); ...
        adenscale(dense,d,qblkstart); -d.det(dense.q)];
    % --------------------------------------------------
    % Correct for sparse part, i.e. let LAD = L\Ad
    % (this also handles L.perm)
    % --------------------------------------------------
    LAD = sparfwslv(L,Ad, symLden.LAD);
    clear Ad;
    % --------------------------------------------------
    % Factor the dense columns as "diag+rank-1"=diag(L.d)+LAD(:,k)*LAD(:,k)'
    % --------------------------------------------------'
    [Lden, Ld] = dpr1fact(LAD, L.d, symLden, smult, pars.maxuden);
    Lden.dz = symLden.dz;
    Lden.first = symLden.first;
    Lden.perm = symLden.perm;
    clear LAD;
else
    % --------------------------------------------------
    % If no dense columns
    % --------------------------------------------------
    Lden.betajc = 0;
    Ld = L.d;
    Lden.rowperm = 1:length(L.d);
end
% ------------------------------------------------------------
% If still Ld(i) = 0, then set to something positive. Otherwise
% we cannot do PCG.
% ------------------------------------------------------------
skip = find(L.skip);
if ~isempty(skip)
    dtol = pars.canceltol * absd(L.perm(skip));
    dtol = max(dtol, pars.abstol);
    Ld(skip(Ld(skip) <= dtol)) = 1;
end