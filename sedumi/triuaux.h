/* ************************************************************
   HEADER triuaux.h
   For use with mex-files in self-dual-minimization package.

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


   ************************************************************ */

#include "blksdp.h"
#ifndef TRIUAUX
#define TRIUAUX

void matperm(double *y, const double *x, const mwIndex *perm, const mwIndex n);
void invmatperm(double *y, const double *x, const mwIndex *perm, const mwIndex n);
void realltxl(double *y, const double *l, const double *x,
              const mwIndex n, double *xl);
void prpiltxl(double *y,double *ypi, const double *ld,const double *ldpi,
              const double *x,const double *xpi,
              const mwIndex n, double *xld);
void realutxu(double *y, const double *u, const double *x,
              const mwIndex n, double *xu);
void prpiutxu(double *y, double *ypi, const double *u, const double *upi,
              const double *x, const double *xpi,
              const mwIndex n, double *xu);
void invutxu(double *y, const double *u, const double *x,
              const mwIndex m, double *xu);
void prpiinvutxu(double *y,double *ypi, const double *u,const double *upi,
             const double *x,const double *xpi, const mwIndex m, double *xu);
void invltxl(double *y, const double *l, const double *x,
              const mwIndex m, double *xl);
void prpiinvltxl(double *y,double *ypi, const double *l,const double *lpi,
                 const double *x,const double *xpi, const mwIndex m, double *xl);
void utmulx(double *y, const double *u, const double *xu, const mwIndex n);
void prpiutmulx(double *y, double *ypi, const double *u, const double *upi,
                const double *xu, const double *xupi, const mwIndex n);
void psdscaleK(double *y, const double *ud, const mwIndex *perm, const double *x,
               const coneK cK, bool transp, double *fwork);
void scaleK(double *y, double *dmult, const double *d, const double *ud,
            const double *qdetd, const mwIndex *perm, const double *x,
            const coneK cK, bool transp, const double *invdx,
            double *fwork);
#endif
