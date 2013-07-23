function [dense,Adotdden] = getdense(At,Ablkjc,K,pars)
% [dense,Adotdden] = getdense(At,Ablkjc,K,pars)
%
% GETDENSE  Creates dense.{l,cols,q}.
%   Try to find small proportion of the cone primitives that appear
%   in a large proportion of the primal constraints.
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
% Initialize: colnz = #nz-constraints per LP/LOR variable,
%  h = MAX(5, max{#nz-constraints | PSD-blocks})
% "5" is taken as a number that should not be beyond the denq-quantile.
% ------------------------------------------------------------
NORMDEN = 5;
colnz = sum(spones(extractA(At,Ablkjc,0,3,1,K.lq+1)),2);
h = max([NORMDEN; sum(findblks(At,Ablkjc,3,[],K.sblkstart),2)]);
[N,m] = size(At);
% ------------------------------------------------------------
% Replace colnz-entries for the Lorentz-trace ("x1") variables by
% nnz-constraints for that Lorentz block. Namely, these variables
% must be removed if the Lorentz block (d[k]*d[k]'*a[k]) is removed.
% ------------------------------------------------------------
i1 = K.mainblks(1); i2 = K.mainblks(2);       % Bounds Lorentz trace.
Ablkq = extractA(At,Ablkjc,1,2,i1,i2);
if i1 < i2          % MATLAB 5.0, 5.1 cannot handle spones(empty)
    Ablkq2 = findblks(At,Ablkjc,2,3,K.qblkstart); % Tag if any nonzero in block
    Ablkq = spones(spones(Ablkq) + Ablkq2);
    colnz(i1:i2-1) = sum(Ablkq,2);             % And store at Lorentz trace pos.
end
% ------------------------------------------------------------
% FIND THE denq-quantile for DENSE COLUMNS:
% If e.g. pars.denq = 0.75 and pars.denf = 20:
%    Find 75% quantile spquant. anything denser than 20*spquant is
%    tagged as dense.
% NB: spquant is chosen such that all columns with nnz <= h are to left of it.
%   This is because we don't allow PSD-block removal.
% ------------------------------------------------------------
bigcolnz = colnz(colnz > h);              % h must be left anyhow
denqN = ceil(pars.denq * length(colnz)) - (N-length(bigcolnz));
if denqN < 1                % denqN is quantile on 1:length(bigcolnz) basis
    spquant = h;              % Take h as threshold
else
    ordnzs = sort(bigcolnz);
    spquant = ordnzs(denqN);  % Take spquant > h as threshold
end
% ------------------------------------------------------------
% Find LP cols, LORENTZ cols, and LORENTZ blocks with >denf*spquant nzs
% ------------------------------------------------------------
dense.cols = find(colnz > pars.denf * spquant);     %dense LP, Q-blk, Q-cols
dense.q = find(colnz(i1:i2-1) > pars.denf * spquant);        %dense Q-blks
dense.l = length(find(dense.cols < i1));            %dense LP cols.
% ----------------------------------------
% The number of dense columns should be small relative to
% the number of constraints - otherwise we better do 1 full Chol.
% ----------------------------------------
if length(dense.cols) > m / 2;
    dense.l = 0; dense.cols = []; dense.q = [];
end
if isempty(dense.q)
    Adotdden = sparse(m,0);         % = Ablkq(dense.q,:)' block struct.
else
    % --------------------------------------------------
    % HANDLE DENSE LORENTZ BLOCKS
    % Let Adotdden lists the nz-constraints for each dense q-block.
    % --------------------------------------------------
    Adotdden = Ablkq(dense.q,:)';
end