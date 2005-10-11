%                                           y = asmDxq(d, x, K [, ddotx])
% ASMDXQ  Assemble y = D(d)x for x in Lorentz part of K.
% [y,t] = AasmDxq(d, x, K [, ddotx]) then y[k]+t(k)*d[k] = D(dk)xk.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [y,t] = asmDxq(d, x, K, ddotx)
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

if isempty(K.q)
    y = zeros(0,1);
    t = zeros(0,1);
else
    % ------------------------------------------------------------
    % Let i1, i2 such that x(i1:i2-1) = "x1", i.e. Lorentz trace part.
    % ------------------------------------------------------------
    if length(x) >= K.lq
        i1 = K.mainblks(1); i2 = K.mainblks(2);
    else
        i1 = 1; i2 = length(K.q)+1;
    end
    t = x(i1:i2-1);
    if nargin < 4
        ddotx = d.q1.*t + ddot(d.q2,x,K.qblkstart);
    end
    % --------------------------------------------------
    % Since d^{1/2} = (d+sqrt(det d)*iota) / trace(d^{1/2}),
    % and d.auxtr = trace(d^{1/2})^2, d.auxdet = sqrt(2*det d),
    % We have P(d)^{1/2}x = t*d + t*auxdet*e_1 - sqrt(det d)*Jx,
    % where t = (d'*x + x(1)*auxdet)/auxtr.
    % --------------------------------------------------
    t = (ddotx + t.* d.auxdet) ./ d.auxtr;      % old t = x1
    sdet = sqrt(d.det);
    y = [t.*(d.auxdet) - sdet .* x(i1:i2-1); qblkmul(sdet,x,K.qblkstart)];
    if nargout < 2
        y = y + [t.*d.q1; qblkmul(t,d.q2,K.qblkstart)];
    end
end