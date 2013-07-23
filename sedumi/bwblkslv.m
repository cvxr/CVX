function y = bwblkslv(L,b) %#ok
% BWBLKSLV Solves block sparse upper-triangular system.
%    y = bwblkslv(L,b) yields the same result as
%              y(L.perm,:) = L.L'\b
%    However, BWBLKSLV is faster than the built-in operator "\",
%    because it uses dense linear algebra and loop-unrolling on
%    supernodes.
%
%    Typical use, with X sparse m x m positive definite and b is m x n:
%            L = sparchol(symbchol(X),X);
%            y = bwblkslv(L,fwblkslv(L,b));
%    Then y solves X*y=b.
%
% See also symbchol, fwblkslv, mldivide, mrdivide

% This file is part of CholTool 1.00
% Copyright (C) 1998 Jos F. Sturm
% CRL, McMaster University, Canada.
% Supported by the Netherlands Organization for Scientific Research (NWO).
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
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%Indicate to the user Matlab cannot find the SeDuMi binaries
sedumi_binary_error();