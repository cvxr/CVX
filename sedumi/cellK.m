%                                                          xcell = cellK(x,K)
% CELLK  Stores SeDuMi cone K-vector in cell-array format.
%
% On output xcell.f and xcell.l are the free and >=0 components,
% xcell.q{k}, xcell.r{k} and xcell.s{k} contain the Lorentz,
% Rotated Lorentz, and PSD-components, resp.
% xcell.s{k} is a K.s(k) x K.s(k) matrix.
%
% See also eigK, eyeK, sedumi.

function xcell = cellK(x,K)
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

% Remark: the cell-array scheme is inspired by SDPT3.
% Unlike in SDPT3 though, the cells are grouped in
% fields like xcell.q and xcell.s.
if nargin < 2
    error('cellK needs 2 input arguments')
end
if ~isstruct(K)
    error('K should be a structure')
end
i = 1;
if isfield(K,'f')
    xcell.f = x(i:i+K.f-1);
    i = i + K.f;
end
if isfield(K,'l')
    xcell.l = x(i:i+K.l-1);
    i = i + K.l;
end
% ------------------------------------------------------------
% Lorentz: this only works OUTSIDE sedumi, not for internal storage
% ------------------------------------------------------------
if isfield(K,'q')
    for k = 1: length(K.q)
        xcell.q{k} = x(i:i+K.q(k)-1);
        i = i + K.q(k);
    end
end
if isfield(K,'r')
    for k = 1: length(K.r)
        xcell.r{k} = x(i:i+K.r(k)-1);
        i = i + K.r(k);
    end
end
% ------------------------------------------------------------
% PSD-blocks:
% ------------------------------------------------------------
if isfield(K,'s')
    if ~isfield(K,'rsdpN')
        for k = 1: length(K.s)
            xcell.s{k} = mat(x(i:i+K.s(k)^2-1));
            i = i + K.s(k)^2;
        end
    else
        % ------------------------------------------------------------
        % Only applicable for internal storage: complex-as-real-storage
        % ------------------------------------------------------------
        for k = 1: K.rsdpN
            xcell.s{k} = mat(x(i:i+K.s(k)^2-1));
            i = i + K.s(k)^2;
        end
        j = sqrt(-1);
        for k = K.rsdpN+1:length(K.s)
            xcell.s{k} = mat(x(i:i+K.s(k)^2-1));
            i = i + K.s(k)^2;
            xcell.s{k} = xcell.s{k} + j * mat(x(i:i+K.s(k)^2-1));
            i = i + K.s(k)^2;
        end
    end
end