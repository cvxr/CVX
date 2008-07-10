%                                                          L = symbchol(X)
% SYMBCHOL Symbolic block sparse Cholesky factorization.
%   L = symbchol(X) returns a structure L that can be used
%   by the efficient block sparse Cholesky solver SPARCHOL.
%   The fields in L have the following meaning:
%
%   L.perm   - Multiple minimum degree ordering.
%
%   L.L      -  Sparse lower triangular matrix, has sparsity structure
%     of Cholesky factor of X(L.perm,L.perm).
%
%   L.xsuper - Supernode partition. Supernode jsup consists of
%     the nodes   L.xsuper(jsup) : L.xsuper(jsup)-1.
%
%   L.split  - Splitting of supernodes. Recommends to split supernode
%     in blocks of sizes   L.split(xsuper(jsup):L.xsuper(jsup)-1).
%
%   L.tmpsiz - Quantity used by SPARCHOL, to allocated enough working
%     storage.
%
%   L = symbchol(X,cachsz) optimizes L.split for a computer cache
%     of size CACHSZ * 1024 byte. Default cachsz = 512.
%
% See also sparchol, sparfwslv, sparbwslv, symbfact, symmmd, chol.

function L = symbchol(X,cachsz)

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

% ----------------------------------------
% Enter here the cache-size in KB, for shaping
% optimal dense blocks of floats.
% ----------------------------------------
if ~issparse(X)
    error('X should be a sparse symmetric matrix')
end
if nargin < 2
    cachsz = 512;
end
% ----------------------------------------
% Compute multiple minimum degree ordering. 
% If the matrix is actually dense we don't bother.
% ----------------------------------------
if spars(X)<1
    perm = ordmmdmex(X);
    L = symfctmex(X,perm);
else
    L.perm=(1:size(X,1))';
    L.L=sparse(tril(ones(size(X))));
    L.xsuper=[1;size(X,1)+1];
end
% ----------------------------------------
% Symbolic Cholesky factorization structures, stored in L.
% ----------------------------------------
L.tmpsiz = choltmpsiz(L);
L.split = cholsplit(L,cachsz);