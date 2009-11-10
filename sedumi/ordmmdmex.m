% perm = ordmmdmex(adjncy)
%   Computes multiple-minimum-degree permutation, for sparse
%   Cholesky. Adjncy is  a sparse symmetric matrix; its diagonal
%   is irrelevant.
%
%   Invokes SPARSPAK-A Release III.
%
% **********  INTERNAL FUNCTION OF CHOLTOOL  **********

function perm = ordmmdmex(adjncy) %#ok

%
%   This file is part of CholTool 1.00
%   Copyright (C) 1998 Jos F. Sturm
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
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

error('At OS prompt, type "make" to create cholTool mex-files.')