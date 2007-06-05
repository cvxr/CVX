%   L = symbchol(X)
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
%     of size CACHSZ * 1024 byte. Default cachsz = 16.
%
% SEE ALSO sparchol, sparfwslv, sparbwslv, symbfact, symmmd, chol.

function [L,flag] = symbchol(X,cachsz)

 %  
 %   This file is part of SeDuMi 1.05
 %   Copyright (C) 2001 Jos F. Sturm
 %     Dept. Econometrics & O.R., Tilburg University, the Netherlands.
 %     Supported by the Netherlands Organization for Scientific Research (NWO).
 %   Affiliation SeDuMi 1.03 and 1.04Beta (2000):
 %     Dept. Quantitative Economics, Maastricht University, the Netherlands.
 %   Affiliations up to SeDuMi 1.02 (AUG1998):
 %     CRL, McMaster University, Canada.
 %     Supported by the Netherlands Organization for Scientific Research (NWO).
 % 
 %   This program is free software; you can redistribute it and/or modify
 %   it under the terms of the GNU General Public License as published by
 %   the Free Software Foundation; either version 2 of the License, or
 %   (at your option) any later version.
 % 
 %   This program is distributed in the hope that it will be useful,
 %   but WITHOUT ANY WARRANTY; without even the implied warranty of
 %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 %   GNU General Public License for more details.
 % 
 %   You should have received a copy of the GNU General Public License
 %   along with this program; if not, write to the Free Software
 %   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 %

% ----------------------------------------
% Enter here the cache-size in KB, for shaping
% optimal dense blocks of floats.
% ----------------------------------------
 if ~issparse(X)
   error('X should be a sparse symmetric matrix')
 end
 if nargin < 2
   cachsz = 16;
 end
% ----------------------------------------
% Compute multiple minimum degree ordering
% (we may also use MATLAB's symmmd()).
% ----------------------------------------

 perm = mexordmmd(X);

% ----------------------------------------
% Symbolic Cholesky factorization structures, stored in L.
% ----------------------------------------

 [L,flag] = mexsymbfct(X,perm);
 if ~flag
    L.tmpsiz = choltmpsiz(L);
    L.split = cholsplit(L,cachsz);
 end

