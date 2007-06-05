/* ma47mod.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

#include "f2c.h"

integer i_sign(a,b) integer *a, *b;
{
integer x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}

integer i_dnnt(x) doublereal *x;
{
double floor();
return( (*x)>=0 ?
	floor(*x + .5) : -floor(.5 - *x) );
}

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b371 = 1.;
static doublereal c_b372 = 0.;

/* ******************************************************************* */
/* COPYRIGHT (c) 1993 Council for the Central Laboratory */
/*                    of the Research Councils */
/* All rights reserved. */

/* None of the comments in this Copyright notice between the lines */
/* of asterisks shall be removed or altered in any way. */

/* This Package is intended for compilation without modification, */
/* so most of the embedded comments have been removed. */

/* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE */
/* Licence, see http://hsl.rl.ac.uk/archive/cou.html */

/* Please note that for an HSL ARCHIVE Licence: */

/* 1. The Package must not be copied for use by any other person. */
/*    Supply of any part of the library by the Licensee to a third party */
/*    shall be subject to prior written agreement between AEA */
/*    Hyprotech UK Limited and the Licensee on suitable terms and */
/*    conditions, which will include financial conditions. */
/* 2. All information on the Package is provided to the Licensee on the */
/*    understanding that the details thereof are confidential. */
/*3. All publications issued by the Licensee that include results obtained*/
/*    with the help of one or more of the Packages shall acknowledge the */
/*    use of the Packages. The Licensee will notify the Numerical Analysis */
/*    Group at Rutherford Appleton Laboratory of any such publication. */
/* 4. The Packages may be modified by or on behalf of the Licensee */
/*    for such use in research applications but at no time shall such */
/*    Packages or modifications thereof become the property of the */
/*    Licensee. The Licensee shall make available free of charge to the */
/*    copyright holder for any purpose all information relating to */
/*    any modification. */
/* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any */
/*    direct or consequential loss or damage whatsoever arising out of */
/*    the use of Packages by the Licensee. */
/* ******************************************************************* */

/* ######DATE 24 May 1993 */
/* 4 Feb 1994. Modified to limit searches for structured pivots under */
/* control of ICNTL(4). */
/* July 1994. Code added to avoid thrashing when there are no */
/* nondefective rows of minimal row count. */
/* 7/12/94. Some defaults changed following further testing. Minor */
/*     defects (mostly in comments) remedied. */
/* 27/3/01. In MA47OD, test for stability of a tile pivot as a pair */
/*   of 1x1 pivots modified to allow for the off-diagonal entry of */
/*   the pivot block. */
/* Subroutine */ int ma47ad_(n, ne, irn, jcn, iw, liw, keep, icntl, rinfo, 
	info)
integer *n, *ne, *irn, *jcn, *iw, *liw, *keep, *icntl;
doublereal *rinfo;
integer *info;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer leaf, node, mark, perm, iwfr, lrow, i, k, ldiag;
    extern /* Subroutine */ int ma47gd_(), ma47hd_(), ma47jd_(), ma47kd_(), 
	    ma47ld_(), ma47md_(), ma47nd_();
    static integer nodes, count, sizes, lp, mp, lw, nv, father, ipe, map;

    /* Parameter adjustments */
    --info;
    --rinfo;
    --icntl;
    --keep;
    --iw;
    --jcn;
    --irn;

    /* Function Body */
    lp = icntl[1];
    mp = icntl[2];
    ldiag = icntl[3];
    for (i = 1; i <= 14; ++i) {
	info[i] = 0;
/* L10: */
    }
    if (ldiag >= 3 && mp > 0) {
/*        WRITE (MP,'(//A,I6,A,I7,A,I7)') ' Entering MA47AD with N =',
N, */
/*     +    '  NE =',NE,'  LIW =',LIW */
/*        WRITE (MP,'(A,6I6)') ' ICNTL(1:6) =', (ICNTL(I),I=1,6) */
	k = min(9,*ne);
	if (ldiag == 4) {
	    k = *ne;
	}
	if (k > 0) {
/*          WRITE (MP,'(A/(3(I6,A,2I6,A)))') ' Matrix entries:', 
*/
/*     +      (I,': (',IRN(I),JCN(I),')',I=1,K) */
/*          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .' */
	}
	k = min(10,*n);
	if (ldiag == 4) {
	    k = *n;
	}
	if (icntl[4] == 1 && k > 0) {
/*          WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I
=1,K) */
/*          IF (K.LT.N) WRITE (MP,'(7X,A)') '     . . .' */
	}
    }
    if (*n < 1 || (*n + *n + *n) / 3 != *n) {
	goto L20;
    }
    if (*ne < 1) {
	goto L30;
    }
    lw = *liw - (*n << 2) - 4;
    ipe = lw + 1;
    count = ipe + *n + 1;
    nv = count + *n + 1;
    leaf = nv + *n + 1;
    perm = 1;
    lrow = perm + *n;
    node = lrow + *n;
    mark = node + *n;
    father = mark + *n;
    sizes = father + *n;
    map = sizes + 2;
    if (icntl[4] != 1) {
	if (lw < (*ne << 1) + *n) {
	    goto L60;
	}
	ma47gd_(n, ne, &irn[1], &jcn[1], &iw[1], &iw[ipe], &iw[count], &iw[nv]
		, &keep[mark], &iwfr, &icntl[1], &info[1]);
	keep[sizes] = iw[count];
	ma47hd_(n, &iw[ipe + 1], &iw[1], &lw, &iwfr, &iw[count], &iw[nv], &
		keep[father], &keep[lrow], &keep[node], &keep[mark], &iw[leaf]
		, &keep[perm], &icntl[1], &info[1]);
	if (info[1] == -3) {
	    goto L50;
	}
    } else {
	if (lw < *ne + *n) {
	    goto L40;
	}
	ma47jd_(n, ne, &irn[1], &jcn[1], &keep[1], &iw[1], &iw[ipe], &iw[
		count], &iw[nv], &keep[mark], &iwfr, &icntl[1], &info[1]);
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    keep[father + i - 1] = iw[nv + i];
/* L12: */
	}
	if (info[1] < 0) {
	    goto L80;
	}
	ma47kd_(n, &iw[ipe + 1], &iw[1], &lw, &iwfr, &keep[1], &iw[count], &
		keep[father], &keep[lrow], &iw[nv], &keep[node], &keep[mark], 
		&iw[leaf], &icntl[1], &info[1]);
	if (info[1] == -3) {
	    goto L50;
	}
	if (info[1] < 0) {
	    goto L90;
	}
    }
    ma47ld_(n, &keep[father], &keep[lrow], &keep[node], &iw[count], &iw[nv], &
	    iw[leaf], &keep[mark], &keep[perm], &nodes, &icntl[1]);
    keep[sizes] = info[3];
    keep[sizes + 1] = nodes;
    ma47nd_(n, ne, &irn[1], &jcn[1], &keep[map], &keep[lrow], &keep[perm], &
	    iw[count], &iw[ipe]);
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	iw[keep[map + k - 1] + *liw - *ne] = jcn[k];
/* L15: */
    }
    i__1 = *ne - keep[sizes];
    i__2 = *liw - nodes - *n - 1;
    ma47md_(n, &i__1, &iw[nodes + *n + 2], &i__2, &keep[lrow], &keep[perm], &
	    nodes, &iw[1], &keep[node], &iw[nodes + 1], &keep[father], &icntl[
	    5], &info[1], &rinfo[1]);
    goto L100;
L20:
    info[1] = -1;
/*      IF (LDIAG.GT.0 .AND. LP.GT.0) WRITE (LP,'(A,I3/A,I8)') */
/*     +    ' **** Error return from MA47AD ****  INFO(1) =',INFO(1), */
/*     +    ' N has value',N */
    goto L100;
L30:
    info[1] = -2;
/*      IF (LDIAG.GT.0 .AND. LP.GT.0) WRITE (LP,'(A,I3/A,I10)') */
/*     +    ' **** Error return from MA47AD ****  INFO(1) =',INFO(1), */
/*     +    ' NE has value',NE */
    goto L100;
L40:
    info[2] = *ne + *n * 5 + 4;
    goto L70;
L50:
    info[2] = *liw + info[2];
    goto L70;
L60:
    info[2] = (*ne << 1) + *n * 5 + 4;
L70:
    info[1] = -3;
/*      IF (LDIAG.GT.0 .AND. LP.GT.0) WRITE (LP,'(A,I3/A,I10,A,I10)') */
/*     +    ' **** Error return from MA47AD ****  INFO(1) =',INFO(1), */
/*     +    ' LIW is too small. It must be increased from',LIW, */
/*     +    ' to at least',INFO(2) */
    return 0;
/*   80 IF (LDIAG.GT.0 .AND. LP.GT.0) WRITE (LP,'(A,I3/A/I10,A)') */
/*     +    ' **** Error return from MA47AD ****  INFO(1) =',INFO(1), */
/*     +    ' Invalid permutation supplied in KEEP. Component',INFO(2), */
/*     +    ' is faulty.' */
L80:
    return 0;
/*   90 IF (LDIAG.GT.0 .AND. LP.GT.0) WRITE (LP,'(A,I3/A/I10,A)') */
/*     +    ' **** Error return from MA47AD ****  INFO(1) =',INFO(1), */
/*     +    ' Invalid pivot sequence supplied in KEEP. Component',INFO(2),
 */
/*     +    ' is negative but the next component is not negative.' */
L90:
    return 0;
L100:
    if (ldiag >= 3 && mp > 0) {
/*        WRITE (MP,'(/A,I7,A,2F7.0/A,5I8/(14X,5(I8)))') */
/*     +    ' Leaving MA47AD with  INFO(1) =',INFO(1),'  RINFO(1:2) ='
, */
/*     +    RINFO(1),RINFO(2),' INFO(2:14) = ', (INFO(I),I=2,14) */
	if (info[1] >= 0) {
	    k = min(10,*ne);
	    if (ldiag == 4) {
		k = *ne;
	    }
/*          WRITE (MP,'(A/(10I6))') ' Column indices:', (JCN(K),K=
1,K) */
/*          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .' */
	    k = min(9,*n);
	    if (ldiag == 4) {
		k = *n;
	    }
/*          WRITE (MP,9000) ' KEEP(1:',N,') =',' (Permutation)', 
*/
/*     +      (KEEP(I),I=1,K) */
/*          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .' */
/* 9000     FORMAT (A,I7,A/A/ (12X,10I6)) */
/*          WRITE (MP,9000) ' KEEP(N+1:N+',N,') =', */
/*     +      ' (No. of entries in permuted rows)', (KEEP(I),I=N+1
,N+K) */
/*          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .' */
/*          WRITE (MP,9000) ' KEEP(2N+1:2N+',N,') =', */
/*     +      ' (Tree nodes at which variables eliminated)', */
/*     +      (KEEP(I),I=2*N+1,2*N+K) */
/*          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .' */
/* Computing MIN */
	    i__1 = 10, i__2 = keep[*n * 5 + 2];
	    k = min(i__1,i__2);
/*          WRITE (MP,9000) ' KEEP(3N+1:3N+',KEEP(5*N+2),') =', */
/*     +      ' (Markowitz costs at tree nodes)', (KEEP(I),I=3*N+1
,3*N+K) */
/*          IF (K.LT.KEEP(5*N+2)) WRITE (MP,'(16X,A)') ' . . .' */
/*          WRITE (MP,9000) ' KEEP(4N+1:4N+',KEEP(5*N+2),') =', */
/*     +      ' (Fathers of nodes in tree)', (KEEP(I),I=4*N+1,4*N+
K) */
/*          IF (K.LT.KEEP(5*N+2)) WRITE (MP,'(16X,A)') ' . . .' */
/*          WRITE (MP,9010) ' KEEP(5N+1:5N+2) =', */
/*     +      ' (Nos. of faulty entries and tree nodes)',KEEP(5*N+
1), */
/*     +      KEEP(5*N+2) */
/* 9010     FORMAT (A/A/ (12X,2I6)) */
	    k = min(10,*ne);
	    if (ldiag == 4) {
		k = *ne;
	    }
/*          WRITE (MP,9000) ' KEEP(5N+3:5N+2+',NE,') =',' (Map arr
ay)', */
/*     +      (KEEP(I),I=5*N+3,5*N+2+K) */
/*          IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .' */
	}
    }
} /* ma47ad_ */

/* Subroutine */ int ma47bd_(n, ne, jcn, a, la, iw, liw, keep, cntl, icntl, 
	iw1, rinfo, info)
integer *n, *ne, *jcn;
doublereal *a;
integer *la, *iw, *liw, *keep;
doublereal *cntl;
integer *icntl, *iw1;
doublereal *rinfo;
integer *info;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer node, mark, perm, lrow, i, k, ldiag;
    extern /* Subroutine */ int ma47od_(), ma47ud_();
    static integer sizes, ke, lp, mp, father, map;

    /* Parameter adjustments */
    --info;
    --rinfo;
    --iw1;
    --icntl;
    --cntl;
    --keep;
    --iw;
    --a;
    --jcn;

    /* Function Body */
    lp = icntl[1];
    mp = icntl[2];
    ldiag = icntl[3];
    if (*n < 1) {
	goto L30;
    }
    if ((*n + *n + *n) / 3 != *n) {
	goto L30;
    }
    if (*ne < 1) {
	goto L40;
    }
    if (*liw < *ne << 1) {
	goto L50;
    }
    if (*la < *ne << 1) {
	goto L70;
    }
    info[1] = 0;
    info[2] = 0;
    for (i = 15; i <= 24; ++i) {
	info[i] = 0;
/* L10: */
    }
    if (ldiag > 2 && mp > 0) {
/*        WRITE (MP,'(//A,I6,3(A,I7)/A,3I7/A,I7/A,1P2D12.4)') */
/*     +    ' Entering MA47BD with N =',N,'  NE =',NE,'  LA =',LA, */
/*     +    '  LIW =',LIW,' ICNTL(1:3)  =', (ICNTL(I),I=1,3), */
/*     +    ' ICNTL(5)    =',ICNTL(5),' CNTL(1:2)   =',CNTL(1),CNTL(2)
 */
	ke = min(6,*ne);
	if (ldiag == 4) {
	    ke = *ne;
	}
	if (ke > 0) {
/*          WRITE (MP,'(A/( 3(I6,'':'',1PD12.4) ) )') */
/*     +      ' Matrix entries:', (K,A(K),K=1,KE) */
/*          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .' */
/*          WRITE (MP,'(A/(10I6))') ' Column indices:', (JCN(K),K=
1,KE) */
/*          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .' */
	}
	k = min(10,*n);
	if (ldiag == 4) {
	    k = *n;
	}
/*        WRITE (MP,9000) ' KEEP(1:',N,') =',' (Permutation)', */
/*     +    (KEEP(I),I=1,K) */
/*        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .' */
/* 9000   FORMAT (A,I7,A/A/ (12X,10I6)) */
/*        WRITE (MP,9000) ' KEEP(N+1:N+',N,') =', */
/*     +    ' (No. of entries in permuted rows)', (KEEP(I),I=N+1,N+K) 
*/
/*        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .' */
/*        WRITE (MP,9000) ' KEEP(2N+1:2N+',N,') =', */
/*     +    ' (Tree nodes at which variables eliminated)', */
/*     +    (KEEP(I),I=2*N+1,2*N+K) */
/* Computing MIN */
	i__1 = k, i__2 = keep[*n * 5 + 2];
	k = min(i__1,i__2);
/*        IF (K.LT.KEEP(5*N+2)) WRITE (MP,'(16X,A)') ' . . .' */
/*        WRITE (MP,9000) ' KEEP(3N+1:3N+',KEEP(5*N+2),') =', */
/*     +    ' (Markowitz costs at tree nodes)', (KEEP(I),I=3*N+1,3*N+K
) */
/*        IF (K.LT.KEEP(5*N+2)) WRITE (MP,'(16X,A)') ' . . .' */
/*        WRITE (MP,9000) ' KEEP(4N+1:4N+',KEEP(5*N+2),') =', */
/*     +    ' (Fathers of nodes in tree)', (KEEP(I),I=4*N+1,4*N+K) */
/*        IF (K.LT.KEEP(5*N+2)) WRITE (MP,'(16X,A)') ' . . .' */
/*        WRITE (MP,9010) ' KEEP(5N+1:5N+2) =', */
/*     +    ' (Nos. of faulty entries and tree nodes)',KEEP(5*N+1), */
/*     +    KEEP(5*N+2) */
/* 9010   FORMAT (A/A/ (12X,2I6)) */
/*        WRITE (MP,9000) ' KEEP(5N+3:5N+2+',NE,') =',' (Map array)', 
*/
/*     +    (KEEP(I),I=5*N+3,5*N+2+KE) */
/*        IF (KE.LT.NE) WRITE (MP,'(16X,A)') ' . . .' */
    }
    perm = 1;
    lrow = perm + *n;
    node = lrow + *n;
    mark = node + *n;
    father = mark + *n;
    sizes = father + *n;
    map = sizes + 2;
/* DIR@ IVDEP */
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	a[keep[map + k - 1] + *la - *ne] = a[k];
	iw[keep[map + k - 1] + *liw - *ne] = jcn[k];
/* L20: */
    }
    i__1 = *ne - keep[sizes];
    ma47od_(n, &i__1, &a[1], la, &iw[1], liw, &keep[lrow], &keep[perm], &keep[
	    sizes + 1], &iw1[1], &keep[mark], &keep[node], &iw1[*n + 1], &
	    keep[father], &cntl[1], &icntl[5], &info[1], &rinfo[1]);
    if (info[1] == -3) {
	goto L60;
    }
    if (info[1] == -4) {
	goto L80;
    }
    goto L90;
L30:
    info[1] = -1;
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,9020) INFO(1) */
/* 9020 FORMAT (' **** Error return from MA47BD ****  INFO(1)=',I3) */
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) */
/*     +    WRITE (LP,'(A,I10)') ' N has value',N */
    goto L100;
L40:
    info[1] = -2;
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,9020) INFO(1) */
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP, */
/*     +    '(A,I10)') ' NE has value',NE */
    goto L100;
L50:
    info[1] = -3;
    info[2] = *ne << 1;
/*   60 IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,9020) INFO(1) */
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,'(A,I10,A,I10)') */
/*     +    ' LIW is too small. It must be increased from',LIW, */
/*     +    ' to at least',INFO(2) */
L60:
    goto L100;
L70:
    info[1] = -4;
    info[2] = *ne << 1;
/*   80 IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,9020) INFO(1) */
/*      IF (LP.GT.0 .AND. LDIAG.GT.0) WRITE (LP,'(A,I10,A,I10)') */
/*     +    ' LA is too small. It must be increased from',LA, */
/*     +    ' to at least',INFO(2) */
L80:
    goto L100;
L90:
    if (ldiag > 2 && mp > 0) {
/*        WRITE (MP,'(/A,I7,A,I7/A,2F9.0/A,2I8/A,5I8/(15X,5I8))') */
/*     +    ' Leaving MA47BD with  INFO(1) =',INFO(1),'  IERROR =', */
/*     +    INFO(2),' RINFO(3:4)  = ',RINFO(3),RINFO(4),' INFO(6:7)   
= ', */
/*     +     (INFO(I),I=6,7),' INFO(15:24) = ', (INFO(I),I=15,24) */
	ma47ud_(&a[1], la, &iw[1], liw, &icntl[1]);
    }
L100:
    return 0;
} /* ma47bd_ */

/* Subroutine */ int ma47cd_(n, a, la, iw, liw, w, rhs, iw1, icntl)
integer *n;
doublereal *a;
integer *la, *iw, *liw;
doublereal *w, *rhs;
integer *iw1, *icntl;
{
    static integer k, ldiag;
    extern /* Subroutine */ int ma47qd_(), ma47rd_(), ma47ud_();
    static integer mp;

    /* Parameter adjustments */
    --icntl;
    --iw1;
    --rhs;
    --w;
    --iw;
    --a;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (ldiag >= 3 && mp > 0) {
/*        WRITE (MP,'(//A,I6,3(A,I7))') ' Entering MA47CD with N =',N,
 */
/*     +    '  LA =',LA,'  LIW =',LIW */
	ma47ud_(&a[1], la, &iw[1], liw, &icntl[1]);
	k = min(10,*n);
	if (ldiag == 4) {
	    k = *n;
	}
/*        WRITE (MP,'(A, 1P,5D13.3/(4X, 1P,5D13.3))') ' RHS', */
/*     +    (RHS(I),I=1,K) */
/*        IF (K.LT.N) WRITE (MP,'(A)') '     . . .' */
    }
    if (iw[2] == 0) {
    } else {
	ma47qd_(n, &a[1], la, &iw[1], liw, &w[1], &rhs[1], &iw1[1], &icntl[1])
		;
	ma47rd_(n, &a[1], la, &iw[1], liw, &w[1], &rhs[1], &iw1[1], &icntl[1])
		;
    }
    if (ldiag >= 3 && mp > 0) {
/*        WRITE (MP,'(//A)') ' Leaving MA47CD with RHS:' */
/*        WRITE (MP,'(1P,5D13.3)') (RHS(I),I=1,K) */
/*        IF (K.LT.N) WRITE (MP,'(A)') '     . . .' */
    }
} /* ma47cd_ */

/* Subroutine */ int ma47fd_(n, ipe, flag_, iw, lw, iwfr, ncmpa)
integer *n, *ipe, *flag_, *iw, *lw, *iwfr, *ncmpa;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer lwfr, i, k, l, ir, len1, len2, len3;

    /* Parameter adjustments */
    --iw;
    --flag_;
    --ipe;

    /* Function Body */
    ++(*ncmpa);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	l = ipe[i];
	if (l <= 0 || flag_[i] == -1) {
	    goto L10;
	}
	ipe[i] = iw[l];
	iw[l] = -i;
L10:
	;
    }
    *iwfr = 1;
    lwfr = 1;
    i__1 = *n;
    for (ir = 1; ir <= i__1; ++ir) {
	i__2 = *lw;
	for (k = lwfr; k <= i__2; ++k) {
	    if (iw[k] < 0) {
		goto L30;
	    }
/* L20: */
	}
	goto L90;
L30:
	i = -iw[k];
	iw[*iwfr] = ipe[i];
	ipe[i] = *iwfr;
	++(*iwfr);
	if (flag_[i] <= -2) {
	    l = *iwfr - 1;
	    len1 = iw[k + 1];
	    iw[l + 1] = len1;
	    len3 = iw[k + 2];
	    iw[l + 2] = len3;
	    len2 = iw[l] - len1 - len3 - 2;
	    *iwfr = l + 3;
	    i__2 = k + 2 + len1;
	    for (lwfr = k + 3; lwfr <= i__2; ++lwfr) {
		if (flag_[iw[lwfr]] < 0) {
		    --iw[l + 1];
		    --iw[l];
		} else {
		    iw[*iwfr] = iw[lwfr];
		    ++(*iwfr);
		}
/* L40: */
	    }
	    k = lwfr;
	    i__2 = k - 1 + len2;
	    for (lwfr = k; lwfr <= i__2; ++lwfr) {
		if (flag_[iw[lwfr]] < 0) {
		    --iw[l];
		} else {
		    iw[*iwfr] = iw[lwfr];
		    ++(*iwfr);
		}
/* L50: */
	    }
	    k = lwfr;
	    i__2 = k - 1 + len3;
	    for (lwfr = k; lwfr <= i__2; ++lwfr) {
		if (flag_[iw[lwfr]] < 0) {
		    --iw[l + 2];
		    --iw[l];
		} else {
		    iw[*iwfr] = iw[lwfr];
		    ++(*iwfr);
		}
/* L60: */
	    }
	} else {
	    i__2 = k + iw[*iwfr - 1];
	    for (lwfr = k + 1; lwfr <= i__2; ++lwfr) {
		iw[*iwfr] = iw[lwfr];
		++(*iwfr);
/* L70: */
	    }
	}
/* L80: */
    }
L90:
    ;
} /* ma47fd_ */

/* Subroutine */ int ma47gd_(n, ne, irn, jcn, iw, ipe, count, nv, flag_, iwfr,
	 icntl, info)
integer *n, *ne, *irn, *jcn, *iw, *ipe, *count, *nv, *flag_, *iwfr, *icntl, *
	info;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j, k, l, ldiag, mp;

    /* Parameter adjustments */
    --info;
    --icntl;
    --flag_;
    --nv;
    --iw;
    --jcn;
    --irn;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (mp <= 0) {
	ldiag = 0;
    }
    info[1] = 0;
    count[0] = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	flag_[i] = 0;
	count[i] = 0;
	nv[i] = -1;
/* L10: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = irn[k];
	j = jcn[k];
	if (min(i,j) < 1 || max(i,j) > *n) {
	    irn[k] = 0;
	    jcn[k] = 0;
	    ++count[0];
	    info[1] = 1;
/*          IF (COUNT(0).LE.1 .AND. LDIAG.GT.1) WRITE (MP,'(2A,I2)
') */
/*     +        ' *** Warning message from subroutine MA47AD ***',
 */
/*     +        '  INFO(1) =',INFO(1) */
/*          IF (COUNT(0).LE.10 .AND. LDIAG.GT.1) WRITE (MP,'(3(I6,
A))') K, */
/*     +        'th entry (in row',I,' and column',J,') ignored' 
*/
	} else if (i != j) {
	    ++count[i];
	    ++count[j];
	} else {
	    ++count[i];
	    nv[i] = 1;
	}
/* L20: */
    }
    ipe[0] = count[0];
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ipe[i] = ipe[i - 1] + count[i] + 1;
/* L30: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = irn[k];
	j = jcn[k];
	iw[ipe[i]] = j;
	--ipe[i];
	if (i != j) {
	    iw[ipe[j]] = i;
	    --ipe[j];
	}
/* L40: */
    }
    info[3] = count[0];
    info[4] = 0;
    *iwfr = 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	l = ipe[i];
	ipe[i] = *iwfr;
	i__2 = l + count[i];
	for (k = l + 1; k <= i__2; ++k) {
	    j = iw[k];
	    if (flag_[j] != i) {
		flag_[j] = i;
		++(*iwfr);
		iw[*iwfr] = j;
	    } else {
		if (i <= j) {
		    ++info[4];
		}
	    }
/* L50: */
	}
	iw[ipe[i]] = *iwfr - ipe[i];
	++(*iwfr);
/* L60: */
    }
    if (info[4] > 0) {
	info[1] += 2;
/*        IF (LDIAG.GT.1) WRITE (MP,'(A/I6,A)') */
/*     +      ' *** Warning message from subroutine MA47AD ***',INFO(4
), */
/*     +      ' duplicate entries found.' */
    }
} /* ma47gd_ */

/* Subroutine */ int ma47hd_(n, ipe, iw, lw, iwfr, count, nv, next, last, ipr,
	 flag_, leaf, svars, icntl, info)
integer *n, *ipe, *iw, *lw, *iwfr, *count, *nv, *next, *last, *ipr, *flag_, *
	leaf, *svars, *icntl, *info;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;

    /* Builtin functions */
    double sqrt();
    integer i_sign();

    /* Local variables */
    static integer kind;
    static doublereal cmin;
    static integer nflg, minf, rank, ncmp, minr, ithr;
    static doublereal cost;
    static integer loop, part, list, pivt, i, k, ldiag;
    extern /* Subroutine */ int ma47fd_(), ma47td_(), ma47vd_(), ma47zd_();
    static logical bound;
    static integer shift, k1, k2, pivot[2], nvpiv, mrows, nrows, ie, ke, me, 
	    ip, jp, ir, is, js, kp, ks, ls, ml, mp, ms, np, nr[3], ns, thresh,
	     jp1, jp2, jp3, jp4, kp1, kp2, np0, nsvars, flg, iel, nel, nsc[3],
	     irl, nvc, piv[2], flg1;

    /* Parameter adjustments */
    --info;
    --icntl;
    --svars;
    --leaf;
    --flag_;
    --ipr;
    --last;
    --next;
    --nv;
    --count;
    --iw;
    --ipe;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (mp <= 0) {
	ldiag = 0;
    }
    mrows = *n;
    if (icntl[4] > 1) {
	mrows = icntl[4];
    }
    ma47td_(n, &ipe[1], &iw[1], lw, &nv[1], &next[1], &last[1], &leaf[1], &
	    flag_[1], &count[1], &svars[1]);
    if (ldiag > 4) {
/*        WRITE (MP,'(/A)') '    Row   NV   List' */
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    kp = ipe[i];
	    if (kp == 0) {
/*            WRITE (MP,'(I7,I5)') I,NV(I) */
	    } else {
/*            WRITE (MP,'(I7,I5,10I6:/ (12X,10I6))') I,NV(I), 
*/
/*     +        (IW(K),K=KP+1,KP+IW(KP)) */
	    }
/* L10: */
	}
    }
/* Computing MAX */
    r__1 = sqrt((real) (*n));
    thresh = dmax(r__1,(float)3.);
    rank = *n;
    minf = 1;
    minr = 1;
    ncmp = 0;
    nflg = *n * 3;
    nel = 0;
    i__1 = *n;
    for (is = 1; is <= i__1; ++is) {
	ipr[is] = 0;
	flag_[is] = nflg;
/* L20: */
    }
    i__1 = *n;
    for (is = 1; is <= i__1; ++is) {
	k = ipe[is];
	if (k == 0) {
	    flag_[is] = -1;
	} else {
	    ir = iw[k];
	    if (ir == 0) {
		rank += nv[is];
		nv[is] = -nv[is];
		ir = 1;
	    }
	    count[is] = ir;
	    ns = ipr[ir];
	    if (ns > 0) {
		last[ns] = is;
	    }
	    next[is] = ns;
	    ipr[ir] = is;
	    last[is] = 0;
	}
/* L30: */
    }
    i__1 = *n;
    for (ml = 1; ml <= i__1; ++ml) {
	if (nel >= *n) {
	    goto L680;
	}
	ir = minr;
	i__2 = *n;
	for (minr = ir; minr <= i__2; ++minr) {
	    if (ipr[minr] != 0) {
		goto L50;
	    }
/* L40: */
	}
L50:
	i__2 = *n;
	for (ithr = 1; ithr <= i__2; ++ithr) {
	    nrows = 0;
/* Computing 2nd power */
	    r__1 = (real) (*n);
	    cmin = r__1 * r__1;
	    i__3 = *n;
	    for (ir = minr; ir <= i__3; ++ir) {
		if (minf > ir) {
		    goto L70;
		}
		ms = ipr[ir];
		i__4 = *n;
		for (list = 1; list <= i__4; ++list) {
		    if (ms == 0) {
			minf = ir + 1;
			goto L70;
		    }
		    if (nv[ms] > 0) {
			if (ir <= thresh) {
			    goto L180;
			}
			goto L160;
		    }
		    ms = next[ms];
/* L60: */
		}
L70:
		ms = ipr[ir];
		i__4 = *n;
		for (list = 1; list <= i__4; ++list) {
		    if (ms == 0) {
			goto L130;
		    }
		    if (nflg <= 4) {
			ma47zd_(n, &flag_[1], &nflg);
		    }
		    ++nrows;
		    --nflg;
		    kp = ipe[ms];
		    kp1 = kp + 1;
		    kp2 = kp + iw[kp];
		    i__5 = kp2;
		    for (kp = kp1; kp <= i__5; ++kp) {
			part = (iw[kp] - 1) / *n;
			ke = iw[kp] - part * *n;
			if (flag_[ke] == -1) {
			    goto L100;
			}
			if (flag_[ke] <= -2) {
			    jp = ipe[ke];
			    if (part == 2) {
				jp1 = jp + 3 + iw[jp + 1];
				jp2 = jp + iw[jp];
			    } else {
				jp1 = jp + 3;
				jp2 = jp + iw[jp] - iw[jp + 2];
			    }
			} else {
			    jp1 = kp;
			    jp2 = kp2;
			}
			i__6 = jp2;
			for (jp = jp1; jp <= i__6; ++jp) {
			    is = iw[jp];
			    if (flag_[is] > nflg) {
				flag_[is] = nflg;
				if (nv[is] < 0) {
				    cost = (real) (count[is] - 1) * (real) (
					    ir - 1);
				} else {
				    cost = (real) (ir + count[is] - 3) * (
					    real) (ir - 1);
				}
				if (cost < cmin) {
				    cmin = cost;
				    pivot[0] = is;
				    pivot[1] = ms;
/* Computing 2nd power */
				    r__1 = (real) (ir - 1);
				    if (cmin <= r__1 * r__1) {
					goto L150;
				    }
				}
			    }
/* L90: */
			}
			if (jp2 == kp2) {
			    goto L110;
			}
L100:
			;
		    }
L110:
		    if (nrows >= mrows) {
			goto L150;
		    }
		    ms = next[ms];
/* L120: */
		}
L130:
/* Computing 2nd power */
		r__1 = (real) ir;
		if (cmin <= r__1 * r__1) {
		    goto L150;
		}
/* L140: */
	    }
L150:
/* Computing MAX */
	    i__3 = count[pivot[0]], i__4 = count[pivot[1]];
	    ir = max(i__3,i__4);
	    if (ir <= thresh) {
		goto L190;
	    }
L160:
	    i__3 = ir + *n / 10;
	    ma47vd_(&thresh, &i__3, n, &ipe[1], &iw[1], lw, &count[1], &nv[1],
		     &next[1], &last[1], &ipr[1], &flag_[1], &nflg);
/* L170: */
	}
L180:
	kind = 1;
	me = ms;
	pivot[0] = ms;
	pivot[1] = ms;
	pivt = 0;
	nvpiv = nv[ms];
/* Computing 2nd power */
	r__1 = (real) (ir - 1);
	cmin = r__1 * r__1;
	flag_[ms] = -1;
/*        IF (LDIAG.GT.4) WRITE (MP,'(A,I5,A,2I5,A,I5,A,F7.0)') */
/*     +      ' Pivot',NEL,':',PIVOT(1),NV(PIVOT(1)),' Row count:', */
/*     +      COUNT(PIVOT(1)),' M-cost=',CMIN */
	nel += nvpiv;
	goto L210;
L190:
	kind = 2;
	if (nv[pivot[0]] < 0) {
	    kind = 3;
	}
	flag_[pivot[0]] = -1;
	flag_[pivot[1]] = -1;
	piv[0] = pivot[0];
	piv[1] = pivot[1];
	nvpiv = (i__2 = nv[pivot[0]], abs(i__2));
	if (nvpiv == (i__2 = nv[pivot[1]], abs(i__2))) {
	    pivt = 0;
	} else {
	    if (nvpiv > (i__2 = nv[pivot[1]], abs(i__2))) {
		pivt = pivot[0];
		nvpiv = (i__2 = nv[pivot[1]], abs(i__2));
	    } else {
		pivt = pivot[1];
	    }
	    ms = leaf[pivt];
	    i__2 = nvpiv;
	    for (k = 2; k <= i__2; ++k) {
		ms = -next[ms];
/* L200: */
	    }
	    leaf[ms] = leaf[pivt];
	    leaf[pivt] = -next[ms];
	    flag_[pivt] = nflg;
	    next[ms] = 0;
	    i__3 = (i__2 = nv[pivt], abs(i__2)) - nvpiv;
	    nv[pivt] = i_sign(&i__3, &nv[pivt]);
	    nv[ms] = i_sign(&nvpiv, &nv[pivt]);
	    count[ms] = count[pivt];
	    if (pivt == pivot[0]) {
		piv[0] = ms;
	    } else {
		piv[1] = ms;
	    }
	}
	me = piv[0];
/*        IF (LDIAG.GT.4) WRITE (*,'(A,I5,A,3I5,A,2I5,A,F7.0)') */
/*     +      ' Pivot',NEL,':',PIV(1),PIV(2),NV(ME),' Row counts:', */
/*     +      COUNT(PIVOT(1)),COUNT(PIVOT(2)),' M-cost=',CMIN */
	nsc[1] = 0;
	nvc = 0;
	nel += nvpiv << 1;
L210:
	nsvars = 0;
	ir = 0;
	for (loop = min(kind,2); loop >= 1; --loop) {
	    ms = pivot[loop - 1];
	    nsc[0] = nsvars;
	    nr[0] = ir;
	    if (ms != pivt) {
		ns = next[ms];
		ls = last[ms];
		next[ms] = 0;
		if (ns > 0) {
		    last[ns] = ls;
		}
		if (ls > 0) {
		    next[ls] = ns;
		} else {
		    ipr[count[ms]] = ns;
		}
	    }
	    kp = ipe[ms];
	    kp1 = kp + 1;
	    kp2 = kp + iw[kp];
	    i__2 = kp2;
	    for (kp = kp1; kp <= i__2; ++kp) {
		part = (iw[kp] - 1) / *n;
		ke = iw[kp] - part * *n;
		if (flag_[ke] == -1) {
		    goto L270;
		}
		if (flag_[ke] <= -2) {
		    ie = ke;
		    i__3 = *n;
		    for (list = 1; list <= i__3; ++list) {
			if (next[ie] == 0) {
			    goto L230;
			}
			iel = ie;
			ie = -last[ie];
			if (ie == me) {
			    goto L240;
			}
			last[iel] = -me;
/* L220: */
		    }
L230:
		    next[ie] = -me;
		    last[ie] = -me;
L240:
		    jp = ipe[ke];
		    jp1 = jp + 3;
		    jp2 = jp + iw[jp];
		    if (part == 0) {
		    } else if (part == 2) {
			jp1 += iw[jp + 1];
		    } else {
			jp2 -= iw[jp + 2];
		    }
		} else {
		    jp1 = kp;
		    jp2 = kp2;
		}
		i__3 = jp2;
		for (jp = jp1; jp <= i__3; ++jp) {
		    is = iw[jp];
		    if (flag_[is] <= loop) {
		    } else if (flag_[is] == 2) {
			nvc += (i__4 = nv[is], abs(i__4));
			++nsc[1];
			flag_[is] = 0;
			count[is] -= nvpiv;
		    } else {
			ir += (i__4 = nv[is], abs(i__4));
			flag_[is] = loop;
			++nsvars;
			svars[nsvars] = is;
			ls = last[is];
			last[is] = 0;
			ns = next[is];
			next[is] = 0;
			if (ns > 0) {
			    last[ns] = ls;
			}
			if (ls > 0) {
			    next[ls] = ns;
			} else {
			    ipr[count[is]] = ns;
			}
			count[is] -= nvpiv;
		    }
/* L250: */
		}
		if (jp2 != kp2 && loop == 1) {
		    if (kind == 1) {
			if (part == 0) {
			    flag_[ke] = -1;
			}
		    } else {
			i__3 = jp2;
			for (jp = jp1; jp <= i__3; ++jp) {
			    if (iw[jp] == pivot[1]) {
				flag_[ke] = -1;
				goto L270;
			    }
/* L260: */
			}
		    }
		}
L270:
		;
	    }
/* L280: */
	}
	if (kind == 1) {
	    nsc[1] = nsvars;
	    nr[0] = ir;
	    nsc[2] = 0;
	    nr[2] = 0;
	} else {
	    next[piv[1]] = -me;
	    last[piv[1]] = -me;
	    count[piv[1]] = -count[piv[1]];
	    if (kind == 2) {
		nsc[2] = nsvars - nsc[0];
		nr[2] = ir - nr[0];
		nsc[1] = nsc[0];
		nr[1] = nr[0];
		nsc[0] = 0;
		nr[0] = ir;
	    } else {
		nsc[2] = nsvars - nsc[0];
		nr[2] = ir - nr[0] + nvc;
		nsc[0] -= nsc[1];
		nr[1] = nr[0];
		nr[0] = ir;
		k1 = 1;
		i__2 = nsc[0] + nsc[1];
		for (k = 1; k <= i__2; ++k) {
		    if (flag_[svars[k]] == 2) {
			ks = svars[k];
			svars[k] = svars[k1];
			svars[k1] = ks;
			++k1;
		    }
/* L310: */
		}
	    }
	}
	if (nsvars == 0) {
	    ipe[me] = 0;
	    goto L670;
	}
	if (nflg <= 4) {
	    ma47zd_(n, &flag_[1], &nflg);
	}
	--nflg;
	i__2 = nsvars;
	for (k = 1; k <= i__2; ++k) {
	    is = svars[k];
	    kp = ipe[is];
	    kp1 = kp + 1;
	    kp2 = iw[kp] + kp;
	    i__3 = kp2;
	    for (kp = kp1; kp <= i__3; ++kp) {
		part = (iw[kp] - 1) / *n;
		ke = iw[kp] - *n * part;
		if (flag_[ke] >= 0) {
		    goto L450;
		}
		if (flag_[ke] == -1) {
		    goto L440;
		}
		if (flag_[ke] == -nflg) {
		    goto L440;
		}
		flag_[ke] = -nflg;
		jp = ipe[ke];
		jp1 = jp + 3;
		jp2 = jp1 + iw[jp + 1];
		jp4 = jp + iw[jp];
		jp3 = jp4 - iw[jp + 2];
		if (kind == 1) {
		    i__4 = jp4;
		    for (jp = jp1; jp <= i__4; ++jp) {
			if (flag_[iw[jp]] > 2) {
			    goto L440;
			}
/* L340: */
		    }
		} else if (kind == 2) {
		    i__4 = jp3;
		    for (jp = jp2; jp <= i__4; ++jp) {
			if (flag_[iw[jp]] > 2) {
			    goto L440;
			}
			if (flag_[iw[jp]] == 1) {
			    goto L440;
			}
/* L350: */
		    }
		    flg1 = 0;
		    i__4 = jp4;
		    for (jp = jp3 + 1; jp <= i__4; ++jp) {
			flg = flag_[iw[jp]];
			if (flg > 2) {
			    goto L440;
			}
			if (flg <= 0) {
			    goto L360;
			}
			flg1 = flg;
L360:
			;
		    }
		    i__4 = jp2 - 1;
		    for (jp = jp1; jp <= i__4; ++jp) {
			flg = flag_[iw[jp]];
			if (flg > 2) {
			    goto L440;
			}
			if (flg <= 0) {
			    goto L370;
			}
			if (flg1 != 0) {
			    goto L440;
			}
L370:
			;
		    }
		} else {
		    i__4 = jp3;
		    for (jp = jp2; jp <= i__4; ++jp) {
			if (flag_[iw[jp]] > 0) {
			    goto L440;
			}
/* L380: */
		    }
		    flg1 = 0;
		    i__4 = jp4;
		    for (jp = jp3 + 1; jp <= i__4; ++jp) {
			flg = flag_[iw[jp]];
			if (flg > 2) {
			    goto L440;
			}
			if (flg <= 0) {
			    goto L390;
			}
			if (flg1 == 0) {
			    flg1 = flg;
			}
			if (flg != flg1) {
			    goto L440;
			}
L390:
			;
		    }
		    flg1 = 3 - flg1;
		    i__4 = jp2 - 1;
		    for (jp = jp1; jp <= i__4; ++jp) {
			flg = flag_[iw[jp]];
			if (flg > 2) {
			    goto L440;
			}
			if (flg <= 0) {
			    goto L400;
			}
			if (flg1 == 3) {
			    flg1 = flg;
			}
			if (flg != flg1) {
			    goto L440;
			}
L400:
			;
		    }
		}
		flag_[ke] = -1;
		ie = ke;
		i__4 = *n;
		for (list = 1; list <= i__4; ++list) {
		    if (next[ie] == 0) {
			goto L430;
		    }
		    iel = ie;
		    ie = -last[ie];
		    if (ie == me) {
			goto L440;
		    }
		    last[iel] = -me;
/* L420: */
		}
L430:
		next[ie] = -me;
		last[ie] = -me;
L440:
		;
	    }
L450:
	    ;
	}
	i__2 = kind;
	for (loop = 1; loop <= i__2; ++loop) {
	    if (loop == 1) {
		k1 = nsc[0] + 1;
		k2 = k1 + nsc[1] - 1;
	    } else if (loop == 2) {
		k1 = k2 + 1;
		k2 = nsvars;
		i__3 = k2;
		for (k = k1; k <= i__3; ++k) {
		    is = svars[k];
		    if (nv[is] != 0) {
			flag_[is] = nflg;
		    }
/* L460: */
		}
	    } else {
		i__3 = k2;
		for (k = k1; k <= i__3; ++k) {
		    is = svars[k];
		    if (nv[is] != 0) {
			flag_[is] = 1;
		    }
/* L470: */
		}
		k1 = 1;
		k2 = nsc[0];
		i__3 = k2;
		for (k = k1; k <= i__3; ++k) {
		    is = svars[k];
		    if (nv[is] != 0) {
			flag_[is] = nflg;
		    }
/* L480: */
		}
	    }
	    i__3 = k2;
	    for (k = k1; k <= i__3; ++k) {
		is = svars[k];
		bound = count[is] > thresh;
		if (loop == 1) {
		    nv[is] = (i__4 = nv[is], abs(i__4));
		}
		if (nflg <= 4) {
		    ma47zd_(n, &flag_[1], &nflg);
		}
		--nflg;
		ir = nr[loop - 1];
		kp = ipe[is];
		np = kp + 1;
		kp1 = kp + 1;
		kp2 = iw[kp] + kp;
		i__4 = kp2;
		for (kp = kp1; kp <= i__4; ++kp) {
		    part = (iw[kp] - 1) / *n;
		    ke = iw[kp] - *n * part;
		    if (flag_[ke] <= -2) {
			flag_[ke] = -nflg;
			if (bound) {
			    goto L510;
			}
			leaf[ke] = part;
			jp = ipe[ke];
			jp1 = jp + 3;
			jp2 = jp + iw[jp];
			if (part == 0) {
			} else if (part == 2) {
			    jp1 += iw[jp + 1];
			} else {
			    jp2 -= iw[jp + 2];
			}
			irl = ir;
			i__5 = jp2;
			for (jp = jp1; jp <= i__5; ++jp) {
			    js = iw[jp];
			    if (flag_[js] <= nflg) {
				goto L500;
			    }
			    ir += (i__6 = nv[js], abs(i__6));
			    flag_[js] = nflg;
L500:
			    ;
			}
			if (ir == irl && last[ke] == -me) {
			    goto L520;
			}
L510:
			iw[np] = iw[kp];
			++np;
		    } else if (flag_[ke] >= 0) {
			goto L530;
		    }
L520:
		    ;
		}
		np0 = np;
		goto L550;
L530:
		np0 = np;
		kp1 = kp;
		i__4 = kp2;
		for (kp = kp1; kp <= i__4; ++kp) {
		    ks = iw[kp];
		    if (flag_[ks] <= nflg) {
			goto L540;
		    }
		    ir += (i__5 = nv[ks], abs(i__5));
		    flag_[ks] = nflg;
		    iw[np] = ks;
		    ++np;
L540:
		    ;
		}
L550:
		if (bound) {
		    ir = count[is];
		}
		if (np > kp2) {
		    kp = ipe[is];
		    if (np + *iwfr - kp >= *lw) {
			i__4 = *iwfr - 1;
			ma47fd_(n, &ipe[1], &flag_[1], &iw[1], &i__4, iwfr, &
				ncmp);
			if (np + *iwfr - kp >= *lw) {
			    info[1] = -3;
			    info[2] = np + *iwfr - kp + 1 - *lw;
			}
			shift = ipe[is] - kp;
			kp += shift;
			kp2 += shift;
			np += shift;
			np0 += shift;
		    }
		    np = np + *iwfr - kp;
		    np0 = np0 + *iwfr - kp;
		    ipe[is] = *iwfr;
		    kp1 = kp;
		    i__4 = kp2;
		    for (kp = kp1; kp <= i__4; ++kp) {
			iw[*iwfr] = iw[kp];
			++(*iwfr);
/* L560: */
		    }
		    iw[*iwfr] = 0;
		    ++(*iwfr);
		}
		iw[np] = iw[np0];
		kp = ipe[is];
		iw[np0] = iw[kp + 1];
		iw[kp + 1] = me + (loop - 1) * *n;
		iw[kp] = np - kp;
		if (ir == 0) {
		    ir = -nv[is];
		    nv[is] = ir;
		    rank -= ir;
		    ir = 1;
		}
		if (cmin == 0.) {
		    goto L610;
		}
		if (ir > thresh) {
		    goto L610;
		}
		js = ipr[ir];
		i__4 = *n;
		for (list = 1; list <= i__4; ++list) {
		    if (js <= 0) {
			goto L610;
		    }
		    kp = ipe[js];
		    if (iw[kp + 1] != me + (loop - 1) * *n) {
			goto L610;
		    }
		    if (i_sign(&c__1, &nv[js]) != i_sign(&c__1, &nv[is])) {
			goto L580;
		    }
		    kp1 = kp;
		    i__5 = kp1 + iw[kp1];
		    for (kp = kp1 + 2; kp <= i__5; ++kp) {
			part = (iw[kp] - 1) / *n;
			ie = iw[kp] - part * *n;
			if ((i__6 = flag_[ie], abs(i__6)) > nflg) {
			    goto L580;
			}
			if (flag_[ie] == -nflg) {
			    if (part != leaf[ie]) {
				goto L580;
			    }
			}
/* L570: */
		    }
		    goto L600;
L580:
		    js = next[js];
/* L590: */
		}
L600:
		ipe[js] = 0;
		nv[is] += nv[js];
		nv[js] = 0;
		flag_[js] = -1;
		ns = next[js];
		ls = last[js];
		if (ns > 0) {
		    last[ns] = is;
		}
		if (ls > 0) {
		    next[ls] = is;
		}
		last[is] = ls;
		next[is] = ns;
		if (ipr[ir] == js) {
		    ipr[ir] = is;
		}
		count[is] = ir;
		next[js] = -leaf[is];
		leaf[is] = leaf[js];
		last[js] = -is;
		goto L620;
L610:
		ns = ipr[ir];
		if (ns > 0) {
		    last[ns] = is;
		}
		next[is] = ns;
		ipr[ir] = is;
		last[is] = 0;
		minr = min(minr,ir);
		if (nv[is] > 0) {
		    minf = min(minf,ir);
		}
		count[is] = ir;
L620:
		;
	    }
/* L630: */
	}
	if (*iwfr + nsvars + 3 >= *lw) {
	    i__2 = *iwfr - 1;
	    ma47fd_(n, &ipe[1], &flag_[1], &iw[1], &i__2, iwfr, &ncmp);
	    if (*iwfr + nsvars + 3 >= *lw) {
		info[1] = -3;
		info[2] = *iwfr + nsvars + 4 - *lw;
		return 0;
	    }
	}
	ip = *iwfr;
	*iwfr += 3;
	k2 = 0;
	for (loop = 1; loop <= 3; ++loop) {
	    k1 = k2 + 1;
	    k2 = k1 + nsc[loop - 1] - 1;
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		is = svars[k];
		if (nv[is] == 0) {
		    --nsc[loop - 1];
		} else {
		    flag_[is] = nflg;
		    iw[*iwfr] = is;
		    ++(*iwfr);
		}
/* L640: */
	    }
/* L650: */
	}
	iw[ip] = *iwfr - ip - 1;
	iw[ip + 1] = nsc[0];
	iw[ip + 2] = nsc[2];
	flag_[me] = -nflg;
	ipe[me] = ip;
L670:
	;
    }
L680:
    if (rank < *n) {
	info[1] += 4;
    }
    info[8] = *n - rank;
    info[12] = ncmp;
} /* ma47hd_ */

/* Subroutine */ int ma47id_(cntl, icntl)
doublereal *cntl;
integer *icntl;
{
    /* Parameter adjustments */
    --icntl;
    --cntl;

    /* Function Body */
    cntl[1] = .001;
    cntl[2] = 0.;
    icntl[1] = 6;
    icntl[2] = 6;
    icntl[3] = 1;
    icntl[4] = 0;
    icntl[5] = 5;
    icntl[6] = 5;
    icntl[7] = 4;
    return 0;
} /* ma47id_ */

/* Subroutine */ int ma47jd_(n, ne, irn, jcn, keep, iw, ipe, count, perm, 
	flag_, iwfr, icntl, info)
integer *n, *ne, *irn, *jcn, *keep, *iw, *ipe, *count, *perm, *flag_, *iwfr, *
	icntl, *info;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j, k, l, ldiag, mp;

    /* Parameter adjustments */
    --info;
    --icntl;
    --flag_;
    --iw;
    --keep;
    --jcn;
    --irn;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (mp <= 0) {
	ldiag = 0;
    }
    info[1] = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	perm[i] = 0;
	flag_[i] = 0;
	count[i] = 0;
/* L10: */
    }
    count[0] = 0;
    perm[0] = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	j = (i__2 = keep[i], abs(i__2));
	if (j < 1 || j > *n) {
	    goto L80;
	}
	if (perm[j] != 0) {
	    goto L80;
	}
	perm[j] = i;
/* L20: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = irn[k];
	j = jcn[k];
	if (min(i,j) < 1 || max(i,j) > *n) {
	    irn[k] = 0;
	    jcn[k] = 0;
	    ++count[0];
	    info[1] = 1;
/*          IF (COUNT(0).LE.1 .AND. LDIAG.GT.1) WRITE (MP,'(2A,I2)
') */
/*     +        ' *** Warning message from subroutine MA47AD ***',
 */
/*     +        '  INFO(1) =',INFO(1) */
/*          IF (COUNT(0).LE.10 .AND. LDIAG.GT.1) WRITE (MP,'(3(I6,
A))') K, */
/*     +        'th entry (in row',I,' and column',J,') ignored' 
*/
	} else if (perm[i] <= perm[j]) {
	    ++count[i];
	} else {
	    ++count[j];
	}
/* L30: */
    }
    ipe[0] = count[0];
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ipe[i] = ipe[i - 1] + count[i] + 1;
/* L40: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = irn[k];
	j = jcn[k];
	if (perm[i] <= perm[j]) {
	    iw[ipe[i]] = j;
	    --ipe[i];
	} else {
	    iw[ipe[j]] = i;
	    --ipe[j];
	}
/* L50: */
    }
    *iwfr = 1;
    info[4] = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	l = ipe[i];
	ipe[i] = *iwfr;
	i__2 = l + count[i];
	for (k = l + 1; k <= i__2; ++k) {
	    j = iw[k];
	    if (flag_[j] != i) {
		flag_[j] = i;
		++(*iwfr);
		iw[*iwfr] = j;
	    } else {
		++info[4];
	    }
/* L60: */
	}
	if (*iwfr > ipe[i]) {
	    iw[ipe[i]] = *iwfr - ipe[i];
	    ++(*iwfr);
	} else {
	    ipe[i] = 0;
	}
/* L70: */
    }
    if (info[4] > 0) {
	info[1] += 2;
/*        IF (LDIAG.GT.1) WRITE (MP,'(A/I6,A)') */
/*     +      ' *** Warning message from subroutine MA47AD ***', */
/*     +      INFO(4),' duplicate entries found.' */
    }
    info[3] = count[0];
    return 0;
L80:
    info[1] = -5;
    info[2] = i;
} /* ma47jd_ */

/* Subroutine */ int ma47kd_(n, ipe, iw, lw, iwfr, keep, count, next, last, 
	ipr, flag_, nexte, vars, icntl, info)
integer *n, *ipe, *iw, *lw, *iwfr, *keep, *count, *next, *last, *ipr, *flag_, 
	*nexte, *vars, *icntl, *info;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static logical diag;
    static integer kind, rank, ncmp, pive[2], loop, list, k, ldiag;
    extern /* Subroutine */ int ma47fd_();
    static integer nvars, k1, pivot[2];
    static logical a11, a21, a22;
    static integer ie, ke, me, ne, ip, jp, kp, is, js, ks, ml, mp, ms, jp1, 
	    jp2, jp3, jp4, kp1, kp2, nv1, nv2, iel, nel, nsc[3], nvc;

    /* Parameter adjustments */
    --info;
    --icntl;
    --vars;
    --nexte;
    --flag_;
    --ipr;
    --last;
    --next;
    --count;
    --keep;
    --iw;
    --ipe;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (mp <= 0) {
	ldiag = 0;
    }
    if (ldiag > 4) {
/*        WRITE (MP,'(/A)') '    Row   List' */
	i__1 = *n;
	for (is = 1; is <= i__1; ++is) {
	    kp = ipe[is];
	    if (kp == 0) {
/*            WRITE (MP,'(I7)') IS */
	    } else {
/*            WRITE (MP,'(I7,10I6:/ (7X,10I6))') IS, */
/*     +        (IW(K),K=KP+1,KP+IW(KP)) */
	    }
/* L1: */
	}
    }
    rank = *n;
    ncmp = 0;
    nel = 0;
    i__1 = *n;
    for (is = 1; is <= i__1; ++is) {
	ipr[is] = 0;
	flag_[is] = 4;
/* L10: */
    }
    i__1 = *n;
    for (ml = 1; ml <= i__1; ++ml) {
	if (nel >= *n) {
	    goto L220;
	}
	ms = keep[nel + 1];
	if (ms > 0) {
	    kind = 1;
	    pivot[0] = ms;
	    pivot[1] = ms;
	    ++nel;
	} else {
	    kind = 2;
	    ms = -ms;
	    if (nel + 2 > *n) {
		goto L230;
	    }
	    is = -keep[nel + 2];
	    if (is < 0) {
		goto L230;
	    }
	    pivot[0] = ms;
	    pivot[1] = is;
	    nvc = 0;
	    nel += 2;
	}
	nvars = 0;
	a21 = FALSE_;
	i__2 = kind;
	for (loop = 1; loop <= i__2; ++loop) {
	    diag = FALSE_;
	    ms = pivot[loop - 1];
	    ke = ipr[ms];
	    pive[loop - 1] = ke;
	    ipr[ms] = 1;
	    i__3 = *n + 1;
	    for (kp = 1; kp <= i__3; ++kp) {
		if (ke != 0) {
		    jp = ipe[ke];
		    jp1 = jp + 3;
		    jp4 = jp + iw[jp];
		    jp2 = jp1 + iw[jp + 1];
		    jp3 = jp4 - iw[jp + 2];
		    kp1 = jp1;
		    kp2 = jp4;
		    i__4 = jp2 - 1;
		    for (jp = jp1; jp <= i__4; ++jp) {
			if (iw[jp] == ms) {
			    goto L40;
			}
/* L20: */
		    }
		    i__4 = jp4;
		    for (jp = jp3 + 1; jp <= i__4; ++jp) {
			if (iw[jp] == ms) {
			    goto L50;
			}
/* L30: */
		    }
		    jp2 = jp4;
		    goto L60;
L40:
		    jp1 = jp2;
		    jp2 = jp4;
		    goto L60;
L50:
		    jp2 = jp3;
		} else {
		    jp = ipe[ms];
		    if (jp == 0) {
			goto L120;
		    }
		    jp1 = jp + 1;
		    jp2 = jp + iw[jp];
		}
L60:
		i__4 = jp2;
		for (jp = jp1; jp <= i__4; ++jp) {
		    is = iw[jp];
		    if (flag_[is] <= 3 - loop) {
		    } else if (flag_[is] == 2) {
			flag_[is] = 0;
			if (is == ms) {
			    diag = TRUE_;
			} else {
			    ++nvc;
			}
		    } else {
			flag_[is] = 3 - loop;
			if (is == ms) {
			    diag = TRUE_;
			} else if (is == pivot[1]) {
			    a21 = TRUE_;
			} else {
			    ++nvars;
			    vars[nvars] = is;
			}
		    }
/* L70: */
		}
		if (ke == 0) {
		    goto L120;
		}
		if (loop == 1) {
		    if (kind == 1) {
			if (jp1 == kp1 && jp2 == kp2) {
			    flag_[ke] = -2;
			}
		    } else {
			i__4 = jp2;
			for (jp = jp1; jp <= i__4; ++jp) {
			    if (iw[jp] == pivot[1]) {
				flag_[ke] = -2;
				goto L90;
			    }
/* L80: */
			}
L90:
			;
		    }
		}
		ke = nexte[ke];
/* L110: */
	    }
L120:
	    if (loop == 1) {
		a11 = diag;
		me = pivot[1];
		if (diag) {
		    me = pivot[0];
		}
		nv1 = nvars;
	    } else {
		a22 = diag;
	    }
	    ke = pive[loop - 1];
	    i__3 = *n + 1;
	    for (kp = 1; kp <= i__3; ++kp) {
		if (ke == 0) {
		    goto L180;
		}
		if (flag_[ke] <= -2) {
		    ie = ke;
		    i__4 = *n;
		    for (list = 1; list <= i__4; ++list) {
			if (next[ie] == 0) {
			    goto L140;
			}
			iel = ie;
			ie = -last[ie];
			if (ie == me) {
			    goto L150;
			}
			last[iel] = -me;
/* L130: */
		    }
L140:
		    next[ie] = -me;
		    last[ie] = -me;
		    if (flag_[ke] == -2 && loop == kind) {
			flag_[ke] = -1;
		    }
		}
L150:
		ne = nexte[ke];
		flag_[ms] = -1;
		if (flag_[ke] <= -2) {
		    jp = ipe[ke];
		    jp1 = jp + 3;
		    jp2 = jp + iw[jp];
		    js = *n + 1;
		    i__4 = jp2;
		    for (jp = jp1; jp <= i__4; ++jp) {
			is = iw[jp];
			if (flag_[is] >= 0) {
/* Computing MIN */
			    i__5 = js, i__6 = next[is];
			    js = min(i__5,i__6);
			}
/* L160: */
		    }
		    if (js <= *n) {
			js = (i__4 = keep[js], abs(i__4));
			nexte[ke] = ipr[js];
			ipr[js] = ke;
		    } else {
			flag_[ke] = -1;
		    }
		}
		ke = ne;
/* L170: */
	    }
L180:
	    ;
	}
	nsc[0] = 0;
	nsc[2] = 0;
	if (kind == 1) {
	    count[ms] = nvars + 1;
	    flag_[ms] = -1;
	} else {
	    nv2 = nvars - nv1 + nvc;
	    count[pivot[0]] = nv1;
	    count[pivot[1]] = nv2;
	    if (a11) {
		count[pivot[0]] = nv1 + 1;
	    }
	    if (a22) {
		count[pivot[1]] = nv2 + 1;
	    }
	    if (a21) {
		++count[pivot[0]];
		++count[pivot[1]];
	    }
	    flag_[pivot[0]] = -1;
	    flag_[pivot[1]] = -1;
	    kind = 4;
	    if (a11) {
		if (a21 && ! a22) {
		    kind = 2;
		    nsc[0] = nvars - nv2;
		    ipr[pivot[1]] = -1;
		}
	    } else {
		if (a21) {
		    if (a22) {
			kind = 2;
			nsc[2] = nv2 - nvc;
			ipr[pivot[0]] = -1;
		    } else {
			kind = 3;
			ipr[pivot[0]] = -1;
			ipr[pivot[1]] = -1;
			nsc[0] = nvars - nv2;
			nsc[2] = nv2 - nvc;
		    }
		}
	    }
	    if (! a22) {
		k1 = 1;
		i__2 = nv1;
		for (k = 1; k <= i__2; ++k) {
		    if (flag_[vars[k]] == 2) {
			ks = vars[k];
			vars[k] = vars[k1];
			vars[k1] = ks;
			++k1;
		    }
/* L190: */
		}
	    }
	    if (me == pivot[0]) {
		next[pivot[1]] = -me;
		last[pivot[1]] = -me;
		count[pivot[1]] = -count[pivot[1]];
	    } else {
		next[pivot[0]] = -me;
		last[pivot[0]] = -me;
		count[pivot[0]] = -count[pivot[0]];
	    }
	    if (kind == 4) {
		count[pivot[0]] = nvars + 2;
		count[pivot[1]] = nvars + 2;
		ipr[pivot[0]] = 0;
		ipr[pivot[1]] = 0;
		ipr[me] = 2;
	    }
	}
	next[me] = 0;
	if (nvars == 0) {
	    if (kind == 1) {
		if (! a11) {
		    --rank;
		}
	    } else {
		if (! a21) {
		    if (! a11) {
			--rank;
		    }
		    if (! a22) {
			--rank;
		    }
		}
	    }
	    ipe[me] = 0;
	    goto L210;
	}
	if (*iwfr + nvars + 3 >= *lw) {
	    i__2 = *iwfr - 1;
	    ma47fd_(n, &ipe[1], &flag_[1], &iw[1], &i__2, iwfr, &ncmp);
	    if (*iwfr + nvars + 3 >= *lw) {
		info[1] = -3;
		info[2] = *iwfr + nvars + 4 - *lw;
		return 0;
	    }
	}
	ip = *iwfr;
	*iwfr += 3;
	ms = *n + 1;
	i__2 = nvars;
	for (k = 1; k <= i__2; ++k) {
	    is = vars[k];
/* Computing MIN */
	    i__3 = ms, i__4 = next[is];
	    ms = min(i__3,i__4);
	    flag_[is] = 4;
	    iw[*iwfr] = is;
	    ++(*iwfr);
/* L200: */
	}
	iw[ip] = *iwfr - ip - 1;
	iw[ip + 1] = nsc[0];
	iw[ip + 2] = nsc[2];
	flag_[me] = -4;
	ipe[me] = ip;
	ms = (i__2 = keep[ms], abs(i__2));
	nexte[me] = ipr[ms];
	ipr[ms] = me;
L210:
	;
    }
L220:
    if (rank < *n) {
	info[1] += 4;
    }
    info[8] = *n - rank;
    info[12] = ncmp;
    return 0;
L230:
    info[2] = nel + 1;
    info[1] = -6;
} /* ma47kd_ */

/* Subroutine */ int ma47ld_(n, father, son, node, count, ne, na, mark, perm, 
	nodes, icntl)
integer *n, *father, *son, *node, *count, *ne, *na, *mark, *perm, *nodes, *
	icntl;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ison, i, j, k, l, ldiag, stage, nemin, level, nvpiv, 
	    count2, mp, nr, ifathr, ibrthr, nde, nrl, nst;

    /* Parameter adjustments */
    --icntl;
    --perm;
    --mark;
    --na;
    --ne;
    --count;
    --node;
    --son;
    --father;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    if (mp <= 0) {
	ldiag = 0;
    }
    nemin = icntl[6];
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	son[i] = 0;
	node[i] = 0;
/* L10: */
    }
    nrl = *n + 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ifathr = -father[i];
	na[i] = -ifathr;
	if (ifathr == 0) {
	    --nrl;
	    mark[nrl] = i;
	} else {
	    ibrthr = son[ifathr];
	    if (ibrthr > 0) {
		na[i] = ibrthr;
	    }
	    son[ifathr] = i;
	    if (count[i] < 0) {
		ne[ifathr] <<= 1;
		ne[i] = 0;
		node[ifathr] = -count[i];
	    }
	}
/* L20: */
    }
    nr = nrl;
    i = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (i <= 0) {
	    i = mark[nr];
	    ++nr;
	    level = *n;
	    nst = 0;
	    perm[level] = i;
	}
	i__2 = *n;
	for (l = 1; l <= i__2; ++l) {
	    if (son[i] <= 0) {
		goto L40;
	    }
	    ison = son[i];
	    son[i] = 0;
	    i = ison;
	    if (node[i] != 0) {
		++nst;
	    }
	    --level;
	    perm[level] = i;
/* L30: */
	}
L40:
	nvpiv = ne[i];
	if (nvpiv == 0) {
	    goto L60;
	}
	if (level == *n) {
	    goto L60;
	}
	ifathr = perm[level + 1];
	if (node[i] == 0 && node[ifathr] == 0) {
	    if (count[i] - nvpiv == count[ifathr]) {
		goto L50;
	    }
	    if (nst > 0) {
		goto L60;
	    }
	    if (nvpiv >= nemin) {
		goto L60;
	    }
	    if (ne[ifathr] >= nemin) {
		goto L60;
	    }
	} else {
	    if (node[i] == 0 || node[ifathr] == 0) {
		goto L60;
	    }
	    if (nvpiv > 0 && ne[ifathr] > 0) {
		if (node[i] - nvpiv / 2 != node[ifathr]) {
		    goto L60;
		}
		if (count[i] - nvpiv != count[ifathr]) {
		    goto L60;
		}
	    } else if (nvpiv < 0 && ne[ifathr] < 0) {
		nvpiv = -nvpiv / 2;
		if (node[i] - nvpiv != node[ifathr]) {
		    goto L60;
		}
		if (count[i] - nvpiv != count[ifathr]) {
		    goto L60;
		}
		if (count[i] == node[i]) {
		    goto L60;
		}
	    } else {
		goto L60;
	    }
	    node[ifathr] = node[i];
	    --nst;
	}
L50:
	ne[ifathr] += ne[i];
	ne[i] = 0;
/*        IF (LDIAG.GT.4) WRITE (MP,'(A,2I5)') ' Merging nodes',I,IFAT
HR */
	count[ifathr] += nvpiv;
	node[i] = -1;
L60:
	ibrthr = na[i];
	if (node[i] > 0) {
	    --nst;
	}
	if (ibrthr > 0) {
	    perm[level] = ibrthr;
	    i = ibrthr;
	    if (node[i] > 0) {
		++nst;
	    }
	} else {
	    ++level;
	    i = -ibrthr;
	}
/* L70: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	son[i] = 0;
/* L80: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ifathr = -father[i];
	na[i] = -ifathr;
	if (ifathr != 0) {
	    ibrthr = son[ifathr];
	    if (ibrthr > 0) {
		na[i] = ibrthr;
	    }
	    son[ifathr] = i;
	}
/* L90: */
    }
    i = 0;
    nr = nrl;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (i <= 0) {
	    i = mark[nr];
	    ++nr;
	    level = *n;
	    perm[*n] = i;
	}
	i__2 = *n;
	for (l = 1; l <= i__2; ++l) {
	    if (son[i] <= 0) {
		goto L110;
	    }
	    ison = son[i];
	    son[i] = 0;
	    i = ison;
	    father[i] = -perm[level];
	    if (node[i] >= 0) {
		--level;
		perm[level] = i;
	    }
/* L100: */
	}
L110:
	ibrthr = na[i];
	if (ibrthr > 0) {
	    if (node[i] < 0) {
		--level;
	    }
	    i = ibrthr;
	    perm[level] = i;
	    father[i] = -perm[level + 1];
	    if (node[i] < 0) {
		++level;
	    }
	} else {
	    if (node[i] >= 0) {
		++level;
	    }
	    i = -ibrthr;
	}
/* L120: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	son[i] = 0;
/* L130: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (ne[i] == 0 && count[i] >= 0) {
	    ifathr = -father[i];
	    ibrthr = son[ifathr];
	    if (ibrthr > 0) {
		father[i] = ibrthr;
	    }
	    son[ifathr] = i;
	}
/* L140: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (count[i] < 0) {
	    ifathr = -father[i];
	    ibrthr = son[ifathr];
	    if (ibrthr > 0) {
		father[i] = ibrthr;
	    }
	    son[ifathr] = i;
	}
/* L150: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (ne[i] != 0) {
	    ifathr = -father[i];
	    if (ifathr != 0) {
		ibrthr = son[ifathr];
		if (ibrthr > 0) {
		    father[i] = ibrthr;
		}
		son[ifathr] = i;
	    }
	}
/* L160: */
    }
    stage = 1;
    i = 0;
    nr = nrl;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (i <= 0) {
	    i = mark[nr];
	    ++nr;
	    level = *n;
	    na[*n] = 0;
	}
	l = level;
	for (level = l; level >= 1; --level) {
	    if (son[i] <= 0) {
		goto L180;
	    }
	    ison = son[i];
	    son[i] = 0;
	    i = ison;
	    na[level - 1] = 0;
/* L170: */
	}
L180:
	perm[k] = i;
	count2 = node[i];
	node[i] = stage;
	if (ne[i] == 0) {
	    goto L210;
	}
	if (level < *n) {
	    ++na[level + 1];
	}
	na[stage] = na[level];
	if (count2 == 0) {
/* Computing 2nd power */
	    i__2 = count[i] - 1;
	    mark[stage] = i__2 * i__2;
	} else if (ne[i] > 0) {
	    i__2 = k - ne[i] / 2;
	    for (j = k - ne[i] + 1; j <= i__2; ++j) {
		node[perm[j]] = -node[perm[j]];
/* L190: */
	    }
	    mark[stage] = (count[i] + count2 - 3) * (count2 - 1);
	} else {
	    i__2 = k;
	    for (j = k + ne[i] + 1; j <= i__2; ++j) {
		node[perm[j]] = -node[perm[j]];
/* L195: */
	    }
	    mark[stage] = (count[i] - 1) * (count2 - 1);
	}
	++stage;
L210:
	ibrthr = father[i];
	if (ibrthr > 0) {
	    na[level] = 0;
	    i = ibrthr;
	} else {
	    ++level;
	    ison = i;
	    i = -ibrthr;
	}
/* L220: */
    }
    *nodes = stage - 1;
    k = 0;
    l = 1;
    i__1 = *nodes;
    for (nde = 1; nde <= i__1; ++nde) {
	father[nde] = *nodes + 1;
	i__2 = na[nde];
	for (i = 1; i <= i__2; ++i) {
	    father[na[k]] = nde;
	    --k;
/* L230: */
	}
	++k;
	na[k] = nde;
/* L240: */
    }
} /* ma47ld_ */

/* Subroutine */ int ma47md_(n, ne, iw, liw, lrow, perm, nodes, elems, node, 
	ppos, father, nbloc, info, rinfo)
integer *n, *ne, *iw, *liw, *lrow, *perm, *nodes, *elems, *node, *ppos, *
	father, *nbloc, *info;
doublereal *rinfo;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer flag_, iell, kind, left, nell, ncmp, iass, iorg, apos, 
	    loop, istk, null, npiv, noxo, apos1, apos3, apos4, norg1, numm1, 
	    numm2, numm3, i, j, laell;
    extern /* Subroutine */ int ma47pd_();
    static integer liell, ielln, ntile, rlspa, iwnfs, j1, j2, n1, n2, iwpos, 
	    n3, notin1, jj, nstkac[2];
    static doublereal flopsb;
    static integer rstkac, intspa, maxfrt, nschur, nfront, numorg, iinput;
    static doublereal flopsx;
    static integer jp1, jp2, jp3, jp4;
    static logical nosure;
    static integer ptrirn, npotpv, ntotpv, elt, isw, nsc1, nsc2, num1, num2, 
	    num3;

    /* Parameter adjustments */
    --rinfo;
    --info;
    --father;
    --node;
    --elems;
    --perm;
    --lrow;
    --iw;

    /* Function Body */
    ncmp = 0;
    ntile = 0;
    noxo = 0;
    flopsb = 0.;
    flopsx = 0.;
    maxfrt = 0;
    nstkac[0] = *ne;
    nstkac[1] = *ne;
    rlspa = *ne;
    intspa = *ne;
    ptrirn = *liw - *ne + 1;
    istk = ptrirn - 1;
    iinput = 0;
    ntotpv = 0;
    npotpv = 0;
    apos = 1;
    iwpos = 8;
    i__1 = *n;
    for (i = 0; i <= i__1; ++i) {
	ppos[i] = *n + 1;
/* L10: */
    }
    i__1 = *nodes;
    for (i = 1; i <= i__1; ++i) {
	elems[i] = 0;
/* L20: */
    }
    i__1 = *nodes;
    for (iass = 1; iass <= i__1; ++iass) {
/* Computing MAX */
	i__2 = intspa, i__3 = iwpos + (*n - ntotpv << 1) + nstkac[1];
	intspa = max(i__2,i__3);
	if ((*n - ntotpv << 1) + 1 > istk) {
	    ma47pd_(&iw[1], &istk, &ptrirn, &iinput, &ncmp);
	    if ((*n - ntotpv << 1) + 1 > istk) {
		info[2] = nstkac[1] + 1 + (*n - ntotpv << 1);
		info[1] = -3;
		goto L390;
	    }
	}
	flag_ = istk - (*n - ntotpv);
	numorg = 0;
	i__2 = *n;
	for (i = npotpv + 1; i <= i__2; ++i) {
	    j = perm[i];
	    if ((i__3 = node[j], abs(i__3)) > iass) {
		goto L40;
	    }
	    ++numorg;
	    iw[numorg] = j;
	    ppos[j] = numorg;
	    iw[flag_ + ppos[j]] = 0;
/* L30: */
	}
L40:
	null = 0;
	i__2 = numorg;
	for (i = 1; i <= i__2; ++i) {
	    null += lrow[npotpv + i];
/* L50: */
	}
	if (null == 0 && elems[iass] == 0) {
	    npotpv += numorg;
	    goto L380;
	}
	kind = 1;
	if (numorg >= 2) {
	    if (node[perm[npotpv + numorg / 2]] < 0) {
		kind = 2;
		if (node[perm[npotpv + numorg]] < 0) {
		    kind = 3;
		}
	    }
	}
	iwnfs = numorg + 1;
	norg1 = numorg;
	if (kind > 1) {
	    norg1 = numorg / 2;
	}
	rstkac = 0;
	i__2 = min(2,kind);
	for (loop = 1; loop <= i__2; ++loop) {
	    num1 = iwnfs - numorg - 1;
	    num2 = 0;
	    if (null == 0) {
		goto L80;
	    }
	    j1 = ptrirn;
	    i__3 = loop * norg1;
	    for (iorg = (loop - 1) * norg1 + 1; iorg <= i__3; ++iorg) {
		j2 = j1 + lrow[npotpv + iorg] - 1;
		i__4 = j2;
		for (jj = j1; jj <= i__4; ++jj) {
		    j = iw[jj];
		    if (ppos[j] <= *n) {
			if (loop == 2) {
			    if (iw[flag_ + ppos[j]] == 1) {
				iw[flag_ + ppos[j]] = 3;
				++num2;
			    }
			}
		    } else {
			iw[iwnfs] = j;
			++iwnfs;
			ppos[j] = iwnfs - 1;
			iw[flag_ + ppos[j]] = loop;
		    }
/* L60: */
		}
		j1 = j2 + 1;
/* L70: */
	    }
	    nstkac[1] = nstkac[1] - j1 + ptrirn;
	    rstkac = rstkac + j1 - ptrirn;
	    iinput = iinput + j1 - ptrirn;
	    ptrirn = j1;
L80:
	    nell = elems[iass];
	    iell = istk + 1;
	    ppos[0] = 0;
	    i__3 = nell;
	    for (elt = 1; elt <= i__3; ++elt) {
		if (iw[iell] < 0) {
		    iell -= iw[iell];
		}
		jp1 = iell + 3;
		jp2 = jp1 + iw[iell + 1];
		jp3 = jp2 + iw[iell + 2];
		jp4 = iell + iw[iell] - 2;
		j1 = jp1;
		j2 = jp4;
		i__4 = jp3 - 1;
		for (jj = jp2; jj <= i__4; ++jj) {
		    j = iw[jj];
		    if (ppos[j] >= (loop - 1) * norg1 + 1 && ppos[j] <= loop *
			     norg1) {
			goto L130;
		    }
/* L90: */
		}
		notin1 = 0;
		i__4 = jp2 - 1;
		for (jj = jp1; jj <= i__4; ++jj) {
		    j = iw[jj];
		    if (ppos[j] >= (loop - 1) * norg1 + 1 && ppos[j] <= loop *
			     norg1) {
			goto L110;
		    }
/* L100: */
		}
		notin1 = 1;
		j2 = jp3 - 1;
L110:
		i__4 = jp4;
		for (jj = jp3; jj <= i__4; ++jj) {
		    j = iw[jj];
		    if (ppos[j] >= (loop - 1) * norg1 + 1 && ppos[j] <= loop *
			     norg1) {
			goto L130;
		    }
/* L120: */
		}
		if (notin1 == 1) {
		    goto L150;
		}
		j1 = jp2;
L130:
		i__4 = j2;
		for (jj = j1; jj <= i__4; ++jj) {
		    j = iw[jj];
		    if (j == 0) {
			goto L140;
		    }
		    if (ppos[j] <= *n) {
			null = 1;
			if (loop == 2) {
			    if (iw[flag_ + ppos[j]] == 1) {
				iw[flag_ + ppos[j]] = 3;
				++num2;
			    }
			}
		    } else {
			iw[iwnfs] = j;
			++iwnfs;
			ppos[j] = iwnfs - 1;
			iw[flag_ + ppos[j]] = loop;
		    }
L140:
		    ;
		}
L150:
		iell += iw[iell];
/* L160: */
	    }
/* L170: */
	}
	ppos[0] = *n + 1;
	if (null == 0) {
	    goto L190;
	}
	nfront = iwnfs - 1;
	maxfrt = max(maxfrt,nfront);
	if (kind == 1) {
	    num1 = 0;
	    num2 = nfront - numorg;
	    num3 = 0;
	} else {
	    if (kind == 3) {
		jj = numorg;
		j1 = jj;
		i__2 = num1;
		for (j = 1; j <= i__2; ++j) {
		    ++jj;
		    if (iw[flag_ + jj] == 1) {
			++j1;
			isw = iw[jj];
			iw[jj] = iw[j1];
			iw[j1] = isw;
			ppos[iw[jj]] = jj;
			ppos[iw[j1]] = j1;
		    }
/* L180: */
		}
		num3 = nfront - numorg - num1;
		num1 -= num2;
	    } else {
		num2 = num1;
		num1 = 0;
		num3 = nfront - numorg - num2;
	    }
	}
/* Computing MAX */
	i__2 = rlspa, i__3 = apos + numorg * (numorg + 1) / 2 + nfront * 
		numorg + nstkac[0];
	rlspa = max(i__2,i__3);
L190:
	if (null != 0) {
	    ntotpv += numorg;
	}
	npotpv += numorg;
	iell = istk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    left = 0;
	    jp1 = iell + 3;
	    jp4 = iell + iw[iell] - 2;
	    jp2 = jp1 + iw[iell + 1];
	    numm1 = iw[iell + 1];
	    numm2 = iw[iell + 2];
	    jp3 = jp2 + numm2;
	    numm3 = iw[iell] - numm1 - numm2 - 4;
	    i__3 = jp2 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L200;
		}
		if (ppos[j] <= numorg) {
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L200:
		;
	    }
	    i__3 = jp3 - 1;
	    for (jj = jp2; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L210;
		}
		if (ppos[j] <= numorg) {
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L210:
		;
	    }
	    if (left == 0) {
		goto L230;
	    }
	    i__3 = jp4;
	    for (jj = jp3; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L220;
		}
		if (ppos[j] <= numorg) {
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L220:
		;
	    }
L230:
	    ppos[0] = *n + 1;
	    ielln = iell + iw[iell];
	    if (left == 0) {
		--nell;
		--elems[iass];
		laell = numm2 * (numm2 + 1) / 2 + numm1 * (numm2 + numm3) + 
			numm2 * numm3 + 2;
		nstkac[0] -= laell;
		liell = iw[iell];
		nstkac[1] -= liell;
		j1 = iell;
		j2 = ielln - 1;
		if (ielln < ptrirn - iinput) {
		    if (iw[ielln] < 0) {
			liell -= iw[ielln];
			j2 -= iw[ielln];
		    }
		}
		if (iell - 1 == istk) {
		    istk += liell;
		} else {
		    if (iw[iell - 1] < 0) {
			liell -= iw[iell - 1];
			j1 = iell + iw[iell - 1];
		    }
		    iw[j1] = -liell;
		    iw[j2] = -liell;
		}
	    }
	    iell = ielln;
/* L250: */
	}
	nstkac[0] -= rstkac;
	if (null == 0) {
	    goto L380;
	}
	npiv = numorg;
	if (kind == 2) {
	    ntile += npiv / 2;
	}
	if (kind == 3) {
	    noxo += npiv / 2;
	}
	if (kind == 1) {
	    i__2 = npiv;
	    for (j = 1; j <= i__2; ++j) {
		flopsb += 1;
		j2 = npiv - j;
		i__3 = npiv - j;
		for (j1 = 1; j1 <= i__3; ++j1) {
		    flopsb = flopsb + 1 + (nfront - npiv + j2 << 1);
		    --j2;
/* L260: */
		}
		flopsb += num2;
/* L270: */
	    }
	} else {
	    i__2 = norg1;
	    for (j = 1; j <= i__2; ++j) {
		flopsb = flopsb + 3 + (norg1 - j) * ((nfront - npiv / 2 - j <<
			 1) + 1);
		j2 = norg1 - j;
		i__3 = norg1 - j;
		for (j1 = 1; j1 <= i__3; ++j1) {
		    flopsb += nfront - npiv + j2 + 1 << 2;
		    --j2;
/* L280: */
		}
		flopsb = flopsb + num1 + (num2 << 2) + num3;
/* L290: */
	    }
	}
	nosure = *nbloc >= nfront - numorg && kind == 1;
	if (npiv == nfront) {
	    goto L360;
	}
	if (kind > 1) {
	    nsc1 = num1 + num2;
	    nsc2 = num2 + num3;
	} else {
	    nsc1 = nfront - numorg;
	    nsc2 = nfront - numorg;
	}
	if (nsc1 * nsc2 == 0) {
	    goto L360;
	}
	apos3 = apos + numorg * (numorg + 1) / 2 + nfront * npiv;
	if (nosure) {
/* Computing MAX */
	    i__2 = rlspa, i__3 = apos3 + (nsc1 + 1) * nsc1 / 2 + 2 + nstkac[0]
		    ;
	    rlspa = max(i__2,i__3);
	} else {
	    apos4 = apos3 + npiv * nsc2;
	    if (kind > 1) {
/* Computing MAX */
		i__2 = nsc1 * nsc2 + 1, i__3 = num1 * (num2 + num3) + num2 * 
			num3 + num2 * (num2 + 1) / 2 + 2;
		nschur = max(i__2,i__3);
	    } else {
/* Computing MAX */
		i__2 = nsc1 * nsc2 + 1, i__3 = nsc1 * (nsc1 + 1) / 2 + 2;
		nschur = max(i__2,i__3);
	    }
/* Computing MAX */
	    i__2 = rlspa, i__3 = apos4 + nschur + nstkac[0];
	    rlspa = max(i__2,i__3);
	}
	if (kind > 1) {
	    flopsb += (npiv << 1) * num1 * (num2 + num3);
	    flopsb += npiv * ((num3 << 1) + num2 + 1) * num2;
	    flopsx = flopsx + npiv * num2 * num3 + *nbloc * (*nbloc - 1) * 
		    npiv * (num2 / *nbloc);
	} else {
	    flopsb += npiv * (nfront - numorg) * (nfront - numorg + 1);
	    if (npiv >= *nbloc) {
		flopsx += *nbloc * (*nbloc - 1) * npiv * (num2 / *nbloc);
	    }
	}
	if (kind > 1) {
	    numm1 = num1;
	    numm2 = num2;
	} else {
	    numm1 = 0;
	    numm2 = nfront - numorg;
	}
	nell = elems[iass];
	iell = istk + 1;
	ppos[0] = 0;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    ielln = iell + iw[iell];
	    jp1 = iell + 3;
	    jp2 = jp1 + iw[iell + 1];
	    jp3 = jp2 + iw[iell + 2];
	    jp4 = iell + iw[iell] - 2;
	    i__3 = jp3 - 1;
	    for (jj = jp2; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L300;
		}
		if (ppos[iw[jj]] - numorg <= numm1 || ppos[iw[jj]] - numorg > 
			numm1 + numm2) {
		    goto L330;
		}
L300:
		;
	    }
	    apos1 = 0;
	    i__3 = jp2 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L310;
		}
		if (ppos[iw[jj]] - numorg <= numm1) {
		    if (apos1 == 2) {
			goto L330;
		    }
		    apos1 = 1;
		}
		if (ppos[iw[jj]] - numorg > numm1 + numm2) {
		    if (ppos[iw[jj]] == *n + 1) {
			goto L330;
		    }
		    if (apos1 == 1) {
			goto L330;
		    }
		    apos1 = 2;
		}
L310:
		;
	    }
	    i__3 = jp4;
	    for (jj = jp3; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L320;
		}
		if (ppos[iw[jj]] - numorg <= numm1) {
		    if (apos1 == 1) {
			goto L330;
		    }
		}
		if (ppos[iw[jj]] - numorg > numm1 + numm2) {
		    if (ppos[iw[jj]] == *n + 1) {
			goto L330;
		    }
		    if (apos1 == 2) {
			goto L330;
		    }
		}
L320:
		;
	    }
	    --nell;
	    --elems[iass];
	    n1 = iw[iell + 1];
	    n2 = iw[iell + 2];
	    n3 = iw[iell] - 4 - n1 - n2;
	    laell = n2 * (n2 + 1) / 2 + n1 * (n2 + n3) + n2 * n3 + 2;
	    nstkac[0] -= laell;
	    liell = iw[iell];
	    nstkac[1] -= liell;
	    j1 = iell;
	    j2 = ielln - 1;
	    if (ielln < ptrirn - iinput) {
		if (iw[ielln] < 0) {
		    liell -= iw[ielln];
		    j2 -= iw[ielln];
		}
	    }
	    if (iell - 1 == istk) {
		istk += liell;
	    } else {
		if (iw[iell - 1] < 0) {
		    liell -= iw[iell - 1];
		    j1 = iell + iw[iell - 1];
		}
		iw[j1] = -liell;
		iw[j2] = -liell;
	    }
L330:
	    iell = ielln;
/* L340: */
	}
	ppos[0] = *n + 1;
	++elems[father[iass]];
	if (kind > 1) {
	    laell = num1 * (num2 + num3) + num2 * (num2 + 1) / 2 + num2 * 
		    num3 + 2;
	} else {
	    laell = (nsc1 + 1) * nsc1 / 2 + 2;
	}
	nstkac[0] += laell;
	liell = nfront - numorg + 4;
	nstkac[1] += liell;
/* Computing MAX */
	i__2 = intspa, i__3 = iwpos + (nfront << 1) + liell + nstkac[1];
	intspa = max(i__2,i__3);
	if ((nfront << 1) + 1 + liell > istk) {
	    ma47pd_(&iw[1], &istk, &ptrirn, &iinput, &ncmp);
	    if ((nfront << 1) + 1 + liell > istk) {
		info[2] = nstkac[1] + 1 + (nfront << 1) + liell;
		info[1] = -3;
		goto L390;
	    }
	}
	iw[istk] = liell;
	--istk;
	i__2 = nfront - numorg;
	for (i = 1; i <= i__2; ++i) {
	    iw[istk] = iw[nfront + 1 - i];
	    --istk;
/* L350: */
	}
	if (kind > 1) {
	    iw[istk] = num2;
	    iw[istk - 1] = num1;
	} else {
	    iw[istk] = nfront - numorg;
	    iw[istk - 1] = 0;
	}
	iw[istk - 2] = liell;
	istk += -3;
L360:
	i__2 = nfront;
	for (jj = npiv + 1; jj <= i__2; ++jj) {
	    j = (i__3 = iw[jj], abs(i__3));
	    ppos[j] = *n + 1;
/* L370: */
	}
	if (father[iass] <= *nodes) {
	    elems[father[iass]] += nell;
	}
/* ******************************** */
/* ******************************** */
	if (kind == 1) {
	    iwpos += -2;
	    apos = apos + npiv * (npiv + 1) / 2 + npiv * (nfront - npiv);
	} else {
	    if (kind == 2) {
		--iwpos;
		apos = apos + norg1 * (norg1 + 1) * 3 / 2 + norg1 * (num1 + (
			num2 << 1) + num3);
	    } else {
		apos = apos + norg1 * (norg1 + 2) + norg1 * (num1 + (num2 << 
			1) + num3);
	    }
	}
	iwpos = iwpos + nfront + 4;
L380:
	;
    }
    info[5] = *nodes;
/* Computing MAX */
    i__1 = *ne << 1;
    info[6] = max(i__1,rlspa);
/* Computing MAX */
    i__1 = *ne << 1;
    info[7] = max(i__1,intspa);
    info[8] = *n - ntotpv;
    info[9] = maxfrt;
    info[10] = apos - 1;
    info[11] = iwpos - 5;
    info[12] += ncmp;
    info[13] = ntile;
    info[14] = noxo;
    rinfo[1] = flopsb;
    rinfo[2] = flopsx;
L390:
    return 0;
} /* ma47md_ */

/* Subroutine */ int ma47nd_(n, ne, irn, jcn, map, lrow, perm, count, perm0)
integer *n, *ne, *irn, *jcn, *map, *lrow, *perm, *count, *perm0;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, j, k;

    /* Parameter adjustments */
    --perm;
    --lrow;
    --map;
    --jcn;
    --irn;

    /* Function Body */
    count[0] = 0;
    perm0[0] = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	count[i] = 0;
	perm0[perm[i]] = i;
/* L10: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = irn[k];
	j = jcn[k];
	if (perm0[i] <= perm0[j]) {
	    i = perm0[i];
	    ++count[i];
	} else {
	    jcn[k] = i;
	    irn[k] = j;
	    j = perm0[j];
	    ++count[j];
	}
/* L20: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	lrow[i] = count[i];
	count[i] = count[i - 1] + lrow[i];
/* L30: */
    }
    i__1 = *ne;
    for (k = 1; k <= i__1; ++k) {
	i = perm0[irn[k]];
	map[k] = count[i];
	--count[i];
/* L40: */
    }
    return 0;
} /* ma47nd_ */

/* Subroutine */ int ma47od_(n, ne, a, la, iw, liw, lrow, perm, nodes, elems, 
	mark, node, ppos, father, cntl, nbloc, info, rinfo)
integer *n, *ne;
doublereal *a;
integer *la, *iw, *liw, *lrow, *perm, *nodes, *elems, *mark, *node, *ppos, *
	father;
doublereal *cntl;
integer *nbloc, *info;
doublereal *rinfo;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer i_dnnt();

    /* Local variables */
    static integer ibeg, iend, neig, kblk;
    static doublereal amax;
    static integer iell, jcol, kind, left, nblk, ncol, iass, iorg, apos, astk,
	     jmax, nass, ipiv, jpiv, nell, ipos, istk, lpiv, npiv, ptra, irow,
	     krow;
    static doublereal rmax;
    static logical ltwo;
    static doublereal tmax;
    static integer zcol;
    static doublereal swop, zone, a1max;
    static integer j1max, apos1, apos2, apos3, apos4, numm1, numm2, i, j, k, 
	    l, laell;
    extern /* Subroutine */ int dgemm_(), ma47pd_(), ma47sd_();
    static integer iexch, liell, aposa, aposb, aposc, ielln, nnell, aposh, 
	    aposi, aposj;
    extern /* Subroutine */ int ma47yd_(), ma47wd_(), ma47xd_();
    static integer rlspa, noxob, iwnfs, lnpiv, nnull, i1, j1, j2, k1, k2, 
	    iswop, iwpos;
    static doublereal pivot;
    static integer npivd2;
    static doublereal amult1, amult2;
    static integer notin1, cntnz1, cntnz2, jj, offdag, pospv1, pospv2, kr;
    extern integer idamax_();
    static integer aposbb, ncmpbi;
    static doublereal uu;
    static integer ntileb, ncmpbr, atrash, nfullb, nirbdu;
    static doublereal flopsb;
    static integer nrlbdu, nstkac[2];
    static doublereal detpiv;
    static integer ja1, ja2, intspa, ksweep, ainput, ma1, maxfrt, nschur;
    static doublereal maxpiv;
    static integer nfront, ldummy, iinput, numorg, poselt;
    static doublereal flopsx;
    static integer jp1, jp2, jp3, jp4;
    static logical nosure;
    static integer nstruc, ptrelt, ptrirn, pivsiz, npotpv, numpvt, ntotpv, 
	    inc, jjj, jay, elt, nfs;
    static doublereal tol;
    static integer nsc1, nsc2, num1, num2, num3;

/* ?? We could refine this count and have NCOLST, NCOLP etc. */
    /* Parameter adjustments */
    --rinfo;
    --info;
    --cntl;
    --father;
    --node;
    --mark;
    --elems;
    --perm;
    --lrow;
    --iw;
    --a;

    /* Function Body */
    tol = cntl[2];
    nblk = 0;
    ntileb = 0;
    noxob = 0;
    nfullb = 0;
    ncmpbi = 0;
    ncmpbr = 0;
    flopsb = 0.;
    flopsx = 0.;
    neig = 0;
    maxfrt = 0;
    uu = min(cntl[1],.5);
    uu = max(uu,0.);
    i__1 = *n;
    for (i = 0; i <= i__1; ++i) {
	ppos[i] = *n + 1;
/* L10: */
    }
    i__1 = *nodes;
    for (i = 1; i <= i__1; ++i) {
	elems[i] = 0;
/* L20: */
    }
    iwpos = 8;
    aposbb = 1;
    nstkac[0] = *ne;
    nstkac[1] = *ne;
    intspa = *ne;
    rlspa = *ne;
    ptrirn = *liw - *ne + 1;
    ptra = *la - *ne + 1;
    istk = ptrirn - 1;
    astk = ptra - 1;
    ainput = 0;
    iinput = 0;
    ntotpv = 0;
    npotpv = 0;
    i__1 = *nodes;
    for (iass = 1; iass <= i__1; ++iass) {
	nnull = 0;
/* Computing MAX */
	i__2 = intspa, i__3 = iwpos + (*n - ntotpv << 1) + nstkac[1];
	intspa = max(i__2,i__3);
	if (iwpos + (*n - ntotpv << 1) > istk) {
	    ma47pd_(&iw[1], &istk, &ptrirn, &iinput, &ncmpbi);
	    if (iwpos + (*n - ntotpv << 1) > istk) {
		info[2] = intspa;
		info[1] = -3;
		goto L2170;
	    }
	}
	numorg = 0;
	i__2 = *n;
	for (i = npotpv + 1; i <= i__2; ++i) {
	    j = perm[i];
	    if ((i__3 = node[j], abs(i__3)) > iass) {
		goto L40;
	    }
	    iw[iwpos + numorg] = j;
	    ++numorg;
	    ppos[j] = numorg;
/* L30: */
	}
L40:
	kind = 1;
	if (numorg >= 2) {
	    if (node[perm[npotpv + numorg / 2]] < 0) {
		kind = 2;
		numpvt = numorg / 2;
		if (node[perm[npotpv + numorg]] < 0) {
		    kind = 3;
		}
	    }
	}
	nass = numorg;
	nell = elems[iass];
	nnell = nell;
	iell = istk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    if (iw[iell + 1] != -1) {
		goto L60;
	    }
	    i__3 = iell + 2 + iw[iell + 2];
	    for (jj = iell + 3; jj <= i__3; ++jj) {
		j = iw[jj];
		iw[iwpos + nass] = j;
		++nass;
		ppos[j] = nass;
/* L50: */
	    }
L60:
	    iell += iw[iell];
/* L70: */
	}
	iwnfs = iwpos + nass;
	j1 = ptrirn;
	i__2 = numorg;
	for (iorg = 1; iorg <= i__2; ++iorg) {
	    j2 = j1 + lrow[npotpv + iorg] - 1;
	    i__3 = j2;
	    for (jj = j1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= *n) {
		    goto L80;
		}
		iw[iwnfs] = j;
		++iwnfs;
		ppos[j] = iwnfs - iwpos;
L80:
		;
	    }
	    j1 = j2 + 1;
/* L90: */
	}
	iell = istk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    if (iw[iell + 1] == -1) {
		goto L160;
	    }
	    jp1 = iell + 3;
	    jp2 = jp1 + iw[iell + 1];
	    jp3 = jp2 + iw[iell + 2];
	    jp4 = iell + iw[iell] - 2;
	    j1 = jp1;
	    j2 = jp4;
	    i__3 = jp3 - 1;
	    for (jj = jp2; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= nass) {
		    goto L140;
		}
/* L100: */
	    }
	    notin1 = 0;
	    i__3 = jp2 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= nass) {
		    goto L120;
		}
/* L110: */
	    }
	    notin1 = 1;
	    j2 = jp3 - 1;
L120:
	    i__3 = jp4;
	    for (jj = jp3; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= nass) {
		    goto L140;
		}
/* L130: */
	    }
	    if (notin1 == 1) {
		goto L160;
	    }
	    j1 = jp2;
L140:
	    ppos[0] = 0;
	    i__3 = j2;
	    for (jj = j1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= *n) {
		    goto L150;
		}
		iw[iwnfs] = j;
		++iwnfs;
		ppos[j] = iwnfs - iwpos;
L150:
		;
	    }
	    ppos[0] = *n + 1;
L160:
	    iell += iw[iell];
/* L170: */
	}
	ncol = iwnfs - iwpos;
	iell = istk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    if (iw[iell + 1] != -1) {
		goto L190;
	    }
	    ppos[0] = 0;
	    i__3 = iell + iw[iell] - 2;
	    for (jj = iell + 3 + iw[iell + 2]; jj <= i__3; ++jj) {
		j = iw[jj];
		if (ppos[j] <= *n) {
		    goto L180;
		}
		iw[iwnfs] = j;
		++iwnfs;
		ppos[j] = iwnfs - iwpos;
L180:
		;
	    }
	    ppos[0] = *n + 1;
L190:
	    iell += iw[iell];
/* L200: */
	}
	nfront = iwnfs - iwpos;
	maxfrt = max(maxfrt,nfront);
	if (info[1] != -4) {
	    apos = aposbb + nass * (nass + 1) / 2;
	} else {
	    apos = 1;
	}
/* Computing MAX */
	i__2 = rlspa, i__3 = info[2] + apos + nfront * nass + nstkac[0];
	rlspa = max(i__2,i__3);
	if (apos + nfront * nass > astk) {
	    ma47sd_(&a[1], &astk, &ptra, &ainput, &ncmpbr);
	    if (apos + nfront * nass > astk) {
		info[2] = info[2] + apos - 1;
		apos = 1;
		aposbb = 1;
		info[1] = -4;
		if (nfront * nass > astk) {
		    info[2] = rlspa;
		    goto L2170;
		}
	    }
	}
	atrash = apos + nfront * nass;
	i__2 = atrash;
	for (jj = apos; jj <= i__2; ++jj) {
	    a[jj] = 0.;
/* L210: */
	}
	j1 = ptrirn;
	i__2 = numorg;
	for (iorg = 1; iorg <= i__2; ++iorg) {
	    j = perm[npotpv + iorg];
	    aposi = apos + (ppos[j] - 1) * nfront - 1;
	    j2 = j1 + lrow[npotpv + iorg] - 1;
	    i__3 = j2;
	    for (jj = j1; jj <= i__3; ++jj) {
		apos2 = aposi + ppos[iw[jj]];
		a[apos2] += a[ptra];
		++ptra;
/* L220: */
	    }
	    ainput = ainput + j2 - j1 + 1;
	    nstkac[0] = nstkac[0] - j2 + j1 - 1;
	    j1 = j2 + 1;
/* L230: */
	}
	iinput = iinput + j1 - ptrirn;
	nstkac[1] = nstkac[1] - j1 + ptrirn;
	ptrirn = j1;
	npotpv += numorg;
	nfs = numorg;
	iell = istk + 1;
	ptrelt = astk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (i_dnnt(&a[ptrelt]) < 0) {
		ptrelt -= i_dnnt(&a[ptrelt]);
	    }
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    left = 0;
	    poselt = ptrelt + 1;
	    jp1 = iell + 3;
	    jp4 = iell + iw[iell] - 2;
	    if (iw[iell + 1] != -1) {
		goto L260;
	    }
	    liell = jp4 - jp1 + 1;
	    i1 = 1;
	    i__3 = jp4;
	    for (jj = jp1; jj <= i__3; ++jj) {
		j = iw[jj];
		apos2 = apos + nfs * nfront - 1 + ppos[j];
		apos1 = poselt;
/* Computing MIN */
		i__5 = i1, i__6 = iw[iell + 2];
		i__4 = min(i__5,i__6);
		for (i = 1; i <= i__4; ++i) {
		    a[apos2] = a[apos1];
		    apos1 += liell;
		    apos2 += nfront;
/* L240: */
		}
		++i1;
		++poselt;
/* L250: */
	    }
	    nfs += iw[iell + 2];
	    goto L370;
L260:
	    jp2 = jp1 + iw[iell + 1];
	    num1 = iw[iell + 1];
	    num2 = iw[iell + 2];
	    jp3 = jp2 + num2;
	    num3 = iw[iell] - num1 - num2 - 4;
	    i__3 = jp2 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L280;
		}
		if (ppos[j] <= nass) {
		    poselt = ptrelt + 1 + (jj - jp1) * (num2 + num3);
		    aposi = apos + (ppos[j] - 1) * nfront - 1;
		    ppos[0] = atrash - aposi;
		    i__4 = jp4;
		    for (jjj = jp2; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			++poselt;
/* L270: */
		    }
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L280:
		;
	    }
	    i__3 = jp3 - 1;
	    for (jj = jp2; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L320;
		}
		if (ppos[j] <= nass) {
		    poselt = ptrelt + 1 + jj - jp2;
		    inc = num2 + num3;
		    aposi = apos + (ppos[j] - 1) * nfront - 1;
		    ppos[0] = atrash - aposi;
		    i__4 = jp2 - 1;
		    for (jjj = jp1; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			poselt += inc;
/* L290: */
		    }
		    i__4 = jj - 1;
		    for (jjj = jp2; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			--inc;
			poselt += inc;
/* L300: */
		    }
		    i__4 = jp4;
		    for (jjj = jj; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			++poselt;
/* L310: */
		    }
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L320:
		;
	    }
	    if (left == 0) {
		goto L360;
	    }
	    i__3 = jp4;
	    for (jj = jp3; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    goto L350;
		}
		if (ppos[j] <= nass) {
		    poselt = ptrelt + 1 + num2 + jj - jp3;
		    aposi = apos + (ppos[j] - 1) * nfront - 1;
		    ppos[0] = atrash - aposi;
		    inc = num2 + num3;
		    i__4 = jp2 - 1;
		    for (jjj = jp1; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			poselt += inc;
/* L330: */
		    }
		    i__4 = jp3 - 1;
		    for (jjj = jp2; jjj <= i__4; ++jjj) {
			apos2 = aposi + ppos[iw[jjj]];
			a[apos2] += a[poselt];
			--inc;
			poselt += inc;
/* L340: */
		    }
		    iw[jj] = 0;
		} else {
		    ++left;
		}
L350:
		;
	    }
L360:
	    ppos[0] = *n + 1;
L370:
	    ielln = iell + iw[iell];
	    poselt = ptrelt + i_dnnt(&a[ptrelt]);
	    if (left == 0) {
		--nell;
		nnell = nell;
		--elems[iass];
		liell = iw[iell];
		nstkac[1] -= liell;
		j1 = iell;
		j2 = ielln - 1;
		if (ielln < ptrirn - iinput) {
		    if (iw[ielln] < 0) {
			liell -= iw[ielln];
			j2 -= iw[ielln];
		    }
		}
		if (iell - 1 == istk) {
		    istk += liell;
		} else {
		    if (iw[iell - 1] < 0) {
			liell -= iw[iell - 1];
			j1 = iell + iw[iell - 1];
		    }
		    iw[j1] = -liell;
		    iw[j2] = -liell;
		}
		laell = i_dnnt(&a[ptrelt]);
		nstkac[0] -= laell;
		ja1 = ptrelt;
		ja2 = poselt - 1;
		if (poselt < ptra - ainput) {
		    if (i_dnnt(&a[poselt]) < 0) {
			laell -= i_dnnt(&a[poselt]);
			ja2 -= i_dnnt(&a[poselt]);
		    }
		}
		if (ptrelt - 1 == astk) {
		    astk += laell;
		} else {
		    if (i_dnnt(&a[ptrelt - 1]) < 0) {
			laell -= i_dnnt(&a[ptrelt - 1]);
			ja1 = ptrelt + i_dnnt(&a[ptrelt - 1]);
		    }
		    a[ja1] = (doublereal) (-laell);
		    a[ja2] = (doublereal) (-laell);
		}
	    }
	    iell = ielln;
	    ptrelt = poselt;
/* L380: */
	}
	i__2 = numorg;
	for (i = 1; i <= i__2; ++i) {
	    apos2 = apos + (nfront + 1) * (i - 1);
	    i__3 = numorg - i;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2 + j * nfront] += a[apos2 + j];
/* L390: */
	    }
	    i__3 = nass - i;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2 + j] = a[apos2 + j * nfront];
/* L400: */
	    }
/* L410: */
	}
	npiv = 0;
	nstruc = 0;
	if (kind == 1) {
	    goto L760;
	}
	i__2 = kind - 1;
	for (k = 1; k <= i__2; ++k) {
	    poselt = apos + (k - 1) * (numorg / 2 * (nfront + 1)) - 1;
	    i__3 = numorg / 2;
	    for (i = 1; i <= i__3; ++i) {
		i__4 = numorg / 2;
		for (j = i; j <= i__4; ++j) {
		    if (a[poselt + j] != 0.) {
			goto L760;
		    }
/* L420: */
		}
		poselt += nfront;
/* L430: */
	    }
/* L440: */
	}
	cntnz1 = ncol - nass;
	cntnz2 = ncol - nass;
	i__2 = ncol;
	for (j = nass + 1; j <= i__2; ++j) {
	    ja1 = apos + j - 1;
	    i__3 = numpvt;
	    for (i = 1; i <= i__3; ++i) {
		if (a[ja1] != 0.) {
		    goto L460;
		}
		ja1 += nfront;
/* L450: */
	    }
	    --cntnz1;
L460:
	    if (kind == 3) {
		ja2 = apos + nfront * numpvt + j - 1;
		i__3 = numpvt;
		for (i = 1; i <= i__3; ++i) {
		    if (a[ja2] != 0.) {
			goto L480;
		    }
		    ja2 += nfront;
/* L470: */
		}
		--cntnz2;
	    }
L480:
	    ;
	}
	ma1 = (cntnz1 + numpvt - 1) * (cntnz2 + numpvt - 1);
	if (father[iass] != *nodes + 1) {
	}
	i__2 = (numorg - npiv) / 2;
	for (ksweep = 1; ksweep <= i__2; ++ksweep) {
	    nnull = 0;
	    lnpiv = npiv;
	    npivd2 = npiv / 2;
	    i__3 = numpvt;
	    for (ipiv = npivd2 + 1; ipiv <= i__3; ++ipiv) {
		aposi = apos + numpvt + npivd2 + (ipiv - 1) * nfront;
		i__4 = numpvt - npivd2;
		jmax = idamax_(&i__4, &a[aposi], &c__1);
		offdag = aposi + jmax - 1;
		jmax = jmax + numpvt + npivd2;
		pivot = a[offdag];
		a[offdag] = 0.;
		i__4 = ncol - (numpvt + npivd2);
		jj = idamax_(&i__4, &a[aposi], &c__1);
		rmax = (d__1 = a[aposi + jj - 1], abs(d__1));
/* Computing MAX */
		d__1 = rmax, d__2 = abs(pivot);
		if (max(d__1,d__2) <= tol) {
		    ++nnull;
		    goto L640;
		}
		if (abs(pivot) <= tol) {
		    goto L640;
		}
		pospv1 = apos + (ipiv - 1) * (nfront + 1);
		pospv2 = apos + (jmax - 1) * (nfront + 1);
		if (ma1 == 0) {
		    detpiv = -pivot * pivot;
		    if ((d__1 = detpiv / pivot, abs(d__1)) <= tol) {
			goto L640;
		    }
		    a[offdag] = pivot;
		    goto L530;
		}
		aposj = apos + (jmax - 1) * nfront;
		amult2 = a[pospv2];
		a[pospv2] = 0.;
		i__4 = numpvt - npivd2;
		j1max = npivd2 + idamax_(&i__4, &a[apos + npivd2 * nfront + 
			jmax - 1], &nfront);
		tmax = (d__1 = a[apos + jmax - 1 + (j1max - 1) * nfront], abs(
			d__1));
		if (jmax - 1 - (numpvt + npivd2) > 0) {
		    i__4 = jmax - 1 - (numpvt + npivd2);
		    j1max = numpvt + npivd2 + idamax_(&i__4, &a[apos + (
			    numpvt + npivd2) * nfront + jmax - 1], &nfront);
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[apos + jmax - 1 + (j1max - 
			    1) * nfront], abs(d__1));
		    tmax = max(d__2,d__3);
		}
		i__4 = ncol - jmax + 1;
		j1max = jmax - 1 + idamax_(&i__4, &a[pospv2], &c__1);
/* Computing MAX */
		d__2 = tmax, d__3 = (d__1 = a[aposj + j1max - 1], abs(d__1));
		tmax = max(d__2,d__3);
		a[pospv2] = amult2;
		a[offdag] = pivot;
		detpiv = -pivot * pivot;
		if ((d__1 = detpiv / pivot, abs(d__1)) <= tol) {
		    goto L640;
		}
		if (((d__1 = a[pospv2], abs(d__1)) * rmax + abs(pivot) * tmax)
			 * uu > abs(detpiv)) {
		    goto L490;
		}
		if (abs(pivot) * rmax * uu > abs(detpiv)) {
		    goto L490;
		}
		goto L530;
L490:
		if (kind == 3) {
		    goto L640;
		}
/* Computing MAX */
		d__2 = abs(pivot);
		if ((d__1 = a[pospv2], abs(d__1)) < uu * max(d__2,tmax)) {
		    goto L640;
		}
		tmax = 0.;
		amult2 = -(a[offdag] / a[pospv2]);
		ja1 = apos + npivd2 * nfront + (jmax - 1);
		i__4 = numpvt;
		for (j = npivd2 + 1; j <= i__4; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja1] * amult2, abs(d__1));
		    tmax = max(d__2,d__3);
		    ja1 += nfront;
/* L500: */
		}
		ja1 = apos + (numpvt + npivd2) * nfront + (jmax - 1);
		ja2 = apos + (ipiv - 1) * nfront + numpvt + npivd2;
		i__4 = jmax - 1;
		for (j = numpvt + npivd2 + 1; j <= i__4; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja2] + a[ja1] * amult2, abs(
			    d__1));
		    tmax = max(d__2,d__3);
		    ja1 += nfront;
		    ++ja2;
/* L510: */
		}
		ja1 = apos + (jmax - 1) * nfront + jmax;
		ja2 = apos + (ipiv - 1) * nfront + jmax;
		i__4 = ncol;
		for (j = jmax + 1; j <= i__4; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja2] + a[ja1] * amult2, abs(
			    d__1));
		    tmax = max(d__2,d__3);
		    ++ja1;
		    ++ja2;
/* L520: */
		}
		if ((d__1 = detpiv / a[pospv2], abs(d__1)) < uu * tmax) {
		    goto L640;
		}
L530:
		if (ipiv == npivd2 + 1) {
		    goto L550;
		}
		ja1 = apos + (ipiv - 1) * nfront + numpvt;
		j1 = apos + npivd2 * nfront + numpvt;
		i__4 = ncol - numpvt;
		for (jj = 1; jj <= i__4; ++jj) {
		    swop = a[ja1];
		    a[ja1] = a[j1];
		    a[j1] = swop;
		    ++ja1;
		    ++j1;
/* L540: */
		}
		ipos = iwpos + npivd2;
		iexch = iwpos + ipiv - 1;
		iswop = iw[ipos];
		iw[ipos] = iw[iexch];
		iw[iexch] = iswop;
L550:
		if (jmax == numpvt + npivd2 + 1) {
		    goto L590;
		}
		ja1 = apos + jmax - 1;
		j1 = apos + numpvt + npivd2;
		i__4 = numpvt + npivd2;
		for (jj = 1; jj <= i__4; ++jj) {
		    swop = a[ja1];
		    a[ja1] = a[j1];
		    a[j1] = swop;
		    ja1 += nfront;
		    j1 += nfront;
/* L560: */
		}
		ja1 = apos + (numpvt + npivd2 + 1) * nfront + jmax - 1;
		j1 = apos + (numpvt + npivd2) * (nfront + 1) + 1;
		i__4 = jmax - 1;
		for (jj = numpvt + npivd2 + 2; jj <= i__4; ++jj) {
		    swop = a[ja1];
		    a[ja1] = a[j1];
		    a[j1] = swop;
		    ja1 += nfront;
		    ++j1;
/* L570: */
		}
		i__4 = ncol - jmax;
		for (jj = 1; jj <= i__4; ++jj) {
		    ++ja1;
		    ++j1;
		    swop = a[ja1];
		    a[ja1] = a[j1];
		    a[j1] = swop;
/* L580: */
		}
		swop = a[apos + (numpvt + npivd2) * (nfront + 1)];
		a[apos + (numpvt + npivd2) * (nfront + 1)] = a[apos + (jmax - 
			1) * (nfront + 1)];
		a[apos + (jmax - 1) * (nfront + 1)] = swop;
		ipos = iwpos + numpvt + npivd2;
		iexch = iwpos + jmax - 1;
		iswop = iw[ipos];
		iw[ipos] = iw[iexch];
		iw[iexch] = iswop;
L590:
		pospv1 = apos + npivd2 * (nfront + 1);
		pospv2 = apos + (numpvt + npivd2) * (nfront + 1);
		offdag = pospv1 + numpvt;
		flopsb += (float)3.;
		++neig;
/* ?? We could decide not to swop these over. */
/* ?? Also could save if pivot was oxo. */
		a[pospv1] = a[pospv2] / detpiv;
		a[pospv2] = 0.;
		a[offdag] = -a[offdag] / detpiv;
		j1 = offdag + nfront;
		j2 = offdag + (numpvt - npivd2 - 1) * nfront;
		ibeg = j1 + 1;
		iend = apos + (npivd2 + 1) * nfront + ncol - 1;
		i__4 = j2;
		i__5 = nfront;
		for (jj = j1; i__5 < 0 ? jj >= i__4 : jj <= i__4; jj += i__5) 
			{
		    amult1 = -(a[jj] * a[offdag]);
		    flopsb = flopsb + (iend - ibeg + 1 << 1) + 1;
		    k1 = offdag + 1;
/* DIR$            IVDEP */
		    i__6 = iend;
		    for (irow = ibeg; irow <= i__6; ++irow) {
			a[irow] += amult1 * a[k1];
			++k1;
/* L600: */
		    }
		    a[jj] = -amult1;
		    ibeg += nfront;
		    iend += nfront;
/* L610: */
		}
		j1 = pospv2 + 1;
		j2 = pospv2 + nass - numpvt - npivd2 - 1;
		ibeg = j1 + nfront;
		iend = apos + (numpvt + npivd2 + 2) * nfront - 1;
		i__5 = j2;
		for (jj = j1; jj <= i__5; ++jj) {
		    amult1 = -(a[jj - numpvt * nfront] * a[pospv1] + a[jj] * 
			    a[offdag]);
		    amult2 = -(a[jj - numpvt * nfront] * a[offdag]);
		    flopsb += iend - ibeg + 2 << 2;
		    k1 = jj - nfront * numpvt;
		    k2 = jj;
/* DIR$            IVDEP */
		    i__4 = iend;
		    for (irow = ibeg; irow <= i__4; ++irow) {
			a[irow] = a[irow] + amult1 * a[k1] + amult2 * a[k2];
			++k1;
			++k2;
/* L620: */
		    }
		    if (jj <= pospv2 + numpvt - npivd2 - 1) {
			a[jj - numpvt * nfront] = -amult2;
			a[jj] = amult1;
		    } else {
			a[jj - numpvt * nfront] = amult1;
			a[jj] = amult2;
		    }
		    ibeg = ibeg + nfront + 1;
		    iend += nfront;
/* L630: */
		}
		iw[iwpos + numpvt + npivd2] = -iw[iwpos + numpvt + npivd2];
		iw[iwpos + npivd2] = -iw[iwpos + npivd2];
		npiv += 2;
		++npivd2;
		ntotpv += 2;
		if (kind == 3) {
		    ++noxob;
		}
		if (kind == 2) {
		    ++ntileb;
		}
		nstruc += 2;
L640:
		;
	    }
	    if (lnpiv == npiv) {
		goto L660;
	    }
/* L650: */
	}
L660:
	if (npiv == nass) {
	    goto L930;
	}
	if (npiv == numorg || npiv == 0) {
	    goto L750;
	}
	i__2 = npivd2;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + numpvt + npivd2 + (i - 1) * nfront;
	    i__3 = numpvt - npivd2;
	    for (j = 1; j <= i__3; ++j) {
		a[ja1] = -a[ja1];
		++ja1;
/* L670: */
	    }
/* L680: */
	}
	i__2 = numpvt - npivd2;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + npivd2 * nfront + numpvt + (i - 1) * nfront;
	    i__3 = npivd2;
	    for (j = 1; j <= i__3; ++j) {
		a[ja1] = -a[ja1];
		++ja1;
/* L690: */
	    }
/* L700: */
	}
	i__2 = numpvt - npivd2;
	ma47yd_(&a[apos + npivd2], &nfront, &npivd2, &i__2, &a[apos + numpvt *
		 nfront], &nfront, &c__0, &c__0);
	ma47yd_(&a[apos + numpvt], &nfront, &npivd2, &npivd2, &a[apos + 
		npivd2], &nfront, &c__0, &c__0);
	i__2 = numpvt - npivd2;
	ma47yd_(&a[apos + numpvt * nfront], &nfront, &npivd2, &i__2, &a[apos 
		+ npiv], &nfront, &c__0, &c__0);
	i__2 = numpvt - npivd2;
	ma47yd_(&a[apos + npivd2 * nfront + numpvt], &nfront, &i__2, &npivd2, 
		&a[apos + numpvt * nfront], &nfront, &c__0, &c__0);
	ma47yd_(&a[apos + numpvt * (nfront + 1)], &nfront, &npivd2, &npivd2, &
		a[apos + npivd2 * (nfront + 1)], &nfront, &c__1, &c__0);
	i__2 = numpvt - npivd2;
	i__3 = numpvt - npivd2;
	ma47xd_(&a[apos + npiv * (nfront + 1)], &nfront, &i__2, &i__3);
	i__2 = numpvt - npivd2;
	ma47yd_(&a[apos + numpvt * nfront], &nfront, &i__2, &npivd2, &a[apos 
		+ npivd2 * nfront + npiv], &nfront, &c__0, &c__1);
	atrash = (numpvt - npivd2) * (nfront - numpvt - npivd2);
/* Computing MAX */
	i__2 = rlspa, i__3 = info[2] + apos + nfront * nass + atrash + nstkac[
		0];
	rlspa = max(i__2,i__3);
	if (apos + nass * nfront + atrash > astk) {
	    ma47sd_(&a[1], &astk, &ptra, &ainput, &ncmpbr);
	    if (apos + nass * nfront + atrash > astk) {
		info[2] = info[2] + apos - 1;
		i__2 = nass * nfront;
		for (i = 1; i <= i__2; ++i) {
		    a[i] = a[apos + i - 1];
/* L710: */
		}
		apos = 1;
		aposbb = 1;
		info[1] = -4;
		if (nfront * nass + atrash > astk) {
		    info[2] = rlspa;
		    goto L2170;
		}
	    }
	}
	astk -= atrash;
	i__2 = numpvt - npivd2;
	i__3 = nfront - numpvt - npivd2;
	i__5 = nfront - numpvt - npivd2;
	ma47yd_(&a[apos + npivd2 * nfront + npivd2 + numpvt], &nfront, &i__2, 
		&i__3, &a[astk], &i__5, &c__0, &c__0);
	i__2 = nfront - numpvt - npivd2;
	ma47yd_(&a[apos + numpvt * nfront + npivd2 + numpvt], &nfront, &
		npivd2, &i__2, &a[apos + npivd2 * nfront + npivd2 + numpvt], &
		nfront, &c__0, &c__0);
	i__2 = nfront - numpvt - npivd2;
	i__3 = numpvt - npivd2;
	i__5 = nfront - numpvt - npivd2;
	ma47yd_(&a[astk], &i__2, &i__3, &i__5, &a[apos + npiv * nfront + 
		npivd2 + numpvt], &nfront, &c__0, &c__0);
	astk += atrash;
	i__2 = numpvt - npivd2 << 1;
	ma47yd_(&a[apos + npiv], &nfront, &npivd2, &i__2, &a[apos + npiv * 
		nfront], &nfront, &c__0, &c__1);
	i__2 = numpvt - npivd2 << 1;
	ma47yd_(&a[apos + npivd2 * nfront + npiv], &nfront, &npivd2, &i__2, &
		a[apos + npiv], &nfront, &c__0, &c__0);
	i__2 = numpvt - npivd2 << 1;
	ma47yd_(&a[apos + npiv * nfront], &nfront, &i__2, &npivd2, &a[apos + 
		npivd2 * nfront + npiv], &nfront, &c__0, &c__1);
	j1 = iwpos + *n - ntotpv + npiv;
	i__2 = numpvt - 1;
	for (jj = npivd2; jj <= i__2; ++jj) {
	    iw[j1] = iw[iwpos + jj];
	    ++j1;
/* L720: */
	}
	j1 = iwpos + npivd2;
	i__2 = numpvt + npivd2 - 1;
	for (jj = numpvt; jj <= i__2; ++jj) {
	    iw[j1] = iw[iwpos + jj];
	    ++j1;
/* L730: */
	}
	i__2 = numpvt - npivd2;
	for (jj = 1; jj <= i__2; ++jj) {
	    iw[j1] = iw[iwpos + *n - ntotpv + npiv + jj - 1];
	    ++j1;
/* L740: */
	}
L750:
	ncol = nfront;
L760:
/* Computing 2nd power */
	i__2 = ncol - npiv - 1;
	ma1 = i__2 * i__2;
	if (father[iass] != *nodes + 1) {
	}
	i__2 = nass - npiv;
	for (ksweep = 1; ksweep <= i__2; ++ksweep) {
	    nnull = 0;
	    lnpiv = npiv;
	    jpiv = 1;
	    i__3 = nass;
	    for (ipiv = npiv + 1; ipiv <= i__3; ++ipiv) {
		--jpiv;
		if (jpiv == 1) {
		    goto L910;
		}
		if (ipiv > numorg) {
		    ncol = nfront;
		}
		aposi = apos + (ipiv - 1) * nfront;
		pospv1 = aposi + ipiv - 1;
		pivot = a[pospv1];
		if (ma1 == 0) {
		    if (abs(pivot) > tol) {
			pivsiz = 1;
			jmax = ipiv;
			goto L810;
		    }
		}
		a[pospv1] = 0.;
		i__5 = ipiv - npiv;
		j1max = npiv + idamax_(&i__5, &a[apos + npiv * nfront + ipiv 
			- 1], &nfront);
		a1max = (d__1 = a[apos + ipiv - 1 + (j1max - 1) * nfront], 
			abs(d__1));
		i__5 = nass - ipiv + 1;
		jmax = ipiv - 1 + idamax_(&i__5, &a[pospv1], &c__1);
		amax = (d__1 = a[aposi + jmax - 1], abs(d__1));
		if (a1max > amax) {
		    amax = a1max;
		    jmax = j1max;
		}
		if (nass == nfront) {
		    rmax = 0.;
		} else {
		    i__5 = nfront - nass;
		    jj = idamax_(&i__5, &a[aposi + nass], &c__1);
		    rmax = (d__1 = a[aposi + nass + jj - 1], abs(d__1));
		}
/* Computing MAX */
		d__1 = max(amax,rmax), d__2 = abs(pivot);
		if (max(d__1,d__2) <= tol) {
		    ++nnull;
		    goto L910;
		}
/* Computing MAX */
		d__1 = amax, d__2 = abs(pivot);
		if (max(d__1,d__2) <= tol) {
		    goto L910;
		}
		pivsiz = 0;
		if (abs(pivot) > uu * max(rmax,amax)) {
		    pivsiz = 1;
		    a[pospv1] = pivot;
		    goto L810;
		}
		if (npiv + 1 == nass) {
		    a[pospv1] = pivot;
		    goto L910;
		}
		if (rmax < amax) {
		    amult1 = a[aposi + jmax - 1];
		    a[aposi + jmax - 1] = 0.;
		    i__5 = ipiv - npiv;
		    j1max = npiv + idamax_(&i__5, &a[apos + npiv * nfront + 
			    ipiv - 1], &nfront);
		    a1max = (d__1 = a[apos + ipiv - 1 + (j1max - 1) * nfront],
			     abs(d__1));
		    i__5 = nass - ipiv + 1;
		    j1max = ipiv - 1 + idamax_(&i__5, &a[pospv1], &c__1);
/* Computing MAX */
		    d__2 = max(rmax,a1max), d__3 = (d__1 = a[aposi + j1max - 
			    1], abs(d__1));
		    rmax = max(d__2,d__3);
		    a[aposi + jmax - 1] = amult1;
		}
		aposj = apos + (jmax - 1) * nfront;
		if (jmax > numorg) {
		    ncol = nfront;
		}
		pospv2 = aposj + jmax - 1;
		if (ipiv > jmax) {
		    offdag = aposj + ipiv - 1;
		} else {
		    offdag = aposi + jmax - 1;
		}
		amult1 = a[offdag];
		a[offdag] = 0.;
		amult2 = a[pospv2];
		a[pospv2] = 0.;
		i__5 = jmax - npiv;
		j1max = npiv + idamax_(&i__5, &a[apos + npiv * nfront + jmax 
			- 1], &nfront);
		tmax = (d__1 = a[apos + jmax - 1 + (j1max - 1) * nfront], abs(
			d__1));
		i__5 = ncol - jmax + 1;
		j1max = jmax - 1 + idamax_(&i__5, &a[pospv2], &c__1);
/* Computing MAX */
		d__2 = tmax, d__3 = (d__1 = a[aposj + j1max - 1], abs(d__1));
		tmax = max(d__2,d__3);
		a[offdag] = amult1;
		a[pospv2] = amult2;
		a[pospv1] = pivot;
		detpiv = a[pospv1] * a[pospv2] - amax * amax;
/* Computing MAX */
		d__3 = (d__1 = a[pospv1], abs(d__1)), d__4 = (d__2 = a[pospv2]
			, abs(d__2));
		maxpiv = max(d__3,d__4);
		if (maxpiv == 0.) {
		    maxpiv = 1.;
		}
		if (abs(detpiv) / maxpiv <= tol) {
		    goto L910;
		}
		pivsiz = 2;
		if (((d__1 = a[pospv2], abs(d__1)) * rmax + amax * tmax) * uu 
			> abs(detpiv)) {
		    goto L770;
		}
		if (((d__1 = a[pospv1], abs(d__1)) * tmax + amax * rmax) * uu 
			> abs(detpiv)) {
		    goto L770;
		}
		goto L810;
L770:
		if ((d__1 = a[pospv2], abs(d__1)) < uu * max(amax,tmax)) {
		    goto L910;
		}
		if ((d__1 = a[pospv2], abs(d__1)) <= tol) {
		    goto L910;
		}
		tmax = 0.;
		amult2 = -(a[offdag] / a[pospv2]);
		ja1 = apos + (jmax - 1) + npiv * nfront;
		ja2 = apos + (ipiv - 1) + npiv * nfront;
		i__5 = min(ipiv,jmax) - 1;
		for (j = npiv + 1; j <= i__5; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja2] + a[ja1] * amult2, abs(
			    d__1));
		    tmax = max(d__2,d__3);
		    ja1 += nfront;
		    ja2 += nfront;
/* L780: */
		}
		if (jmax > ipiv) {
		    j1 = nfront;
		    j2 = 1;
		    ja1 = apos + ipiv * nfront + jmax - 1;
		    ja2 = apos + (ipiv - 1) * nfront + ipiv;
		} else {
		    j1 = 1;
		    j2 = nfront;
		    ja1 = apos + (jmax - 1) * nfront + jmax;
		    ja2 = apos + jmax * nfront + ipiv - 1;
		}
		i__5 = max(ipiv,jmax) - 1;
		for (j = min(ipiv,jmax) + 1; j <= i__5; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja2] + a[ja1] * amult2, abs(
			    d__1));
		    tmax = max(d__2,d__3);
		    ja1 += j1;
		    ja2 += j2;
/* L790: */
		}
		ja1 = apos + (jmax - 1) * nfront + jmax;
		ja2 = apos + (ipiv - 1) * nfront + jmax;
		i__5 = ncol;
		for (j = max(ipiv,jmax) + 1; j <= i__5; ++j) {
/* Computing MAX */
		    d__2 = tmax, d__3 = (d__1 = a[ja2] + a[ja1] * amult2, abs(
			    d__1));
		    tmax = max(d__2,d__3);
		    ++ja1;
		    ++ja2;
/* L800: */
		}
		if ((d__1 = detpiv / a[pospv2], abs(d__1)) < uu * tmax) {
		    goto L910;
		}
L810:
		lpiv = ipiv;
		if (pivsiz == 2) {
		    lpiv = min(ipiv,jmax);
		}
		i__5 = npiv + pivsiz - 1;
		for (krow = npiv; krow <= i__5; ++krow) {
		    if (lpiv == krow + 1) {
			goto L850;
		    }
		    ja1 = apos + (lpiv - 1);
		    j1 = apos + krow;
		    i__4 = krow;
		    for (jj = 1; jj <= i__4; ++jj) {
			swop = a[ja1];
			a[ja1] = a[j1];
			a[j1] = swop;
			ja1 += nfront;
			j1 += nfront;
/* L820: */
		    }
		    ja1 += nfront;
		    ++j1;
		    i__4 = lpiv - krow - 2;
		    for (jj = 1; jj <= i__4; ++jj) {
			swop = a[ja1];
			a[ja1] = a[j1];
			a[j1] = swop;
			ja1 += nfront;
			++j1;
/* L830: */
		    }
		    swop = a[apos + krow * (nfront + 1)];
		    a[apos + krow * (nfront + 1)] = a[ja1];
		    a[ja1] = swop;
		    i__4 = ncol - lpiv;
		    for (jj = 1; jj <= i__4; ++jj) {
			++ja1;
			++j1;
			swop = a[ja1];
			a[ja1] = a[j1];
			a[j1] = swop;
/* L840: */
		    }
		    ipos = iwpos + krow;
		    iexch = iwpos + lpiv - 1;
		    iswop = iw[ipos];
		    iw[ipos] = iw[iexch];
		    iw[iexch] = iswop;
L850:
		    lpiv = max(ipiv,jmax);
/* L860: */
		}
		pospv1 = apos + npiv * (nfront + 1);
		pospv2 = pospv1 + nfront + 1;
		if (pivsiz == 1) {
		    flopsb += 1.;
		    a[pospv1] = 1. / a[pospv1];
		    if (a[pospv1] < 0.) {
			++neig;
		    }
		    j1 = pospv1 + 1;
		    j2 = pospv1 + nass - (npiv + 1);
		    ibeg = pospv1 + nfront + 1;
		    iend = apos + (npiv + 1) * nfront + ncol - 1;
		    i__5 = j2;
		    for (jj = j1; jj <= i__5; ++jj) {
			amult1 = -a[jj] * a[pospv1];
			jcol = jj;
			flopsb = flopsb + (iend - ibeg + 1 << 1) + 1;
/* DIR$            IVDEP */
			i__4 = iend;
			for (irow = ibeg; irow <= i__4; ++irow) {
			    a[irow] += amult1 * a[jcol];
			    ++jcol;
/* L870: */
			}
			a[jj] = amult1;
			ibeg = ibeg + nfront + 1;
			iend += nfront;
/* L880: */
		    }
		    ++npiv;
		    ++ntotpv;
		    jpiv = 1;
		} else {
		    offdag = pospv1 + 1;
		    flopsb += (float)6.;
		    swop = a[pospv2];
		    if (detpiv < 0.) {
			++neig;
		    } else {
			if (swop < 0.) {
			    neig += 2;
			}
		    }
		    a[pospv2] = a[pospv1] / detpiv;
		    a[pospv1] = swop / detpiv;
		    a[offdag] = -a[offdag] / detpiv;
		    j1 = pospv1 + 2;
		    j2 = pospv1 + nass - (npiv + 1);
		    ibeg = pospv2 + nfront + 1;
		    iend = apos + (npiv + 2) * nfront + ncol - 1;
		    i__5 = j2;
		    for (jj = j1; jj <= i__5; ++jj) {
			k1 = jj;
			k2 = jj + nfront;
			amult1 = -(a[pospv1] * a[k1] + a[pospv1 + 1] * a[k2]);
			amult2 = -(a[pospv1 + 1] * a[k1] + a[pospv2] * a[k2]);
			flopsb = flopsb + (iend - ibeg + 1 << 2) + 6;
/* DIR$            IVDEP */
			i__4 = iend;
			for (irow = ibeg; irow <= i__4; ++irow) {
			    a[irow] = a[irow] + amult1 * a[k1] + amult2 * a[
				    k2];
			    ++k1;
			    ++k2;
/* L890: */
			}
			a[jj] = amult1;
			a[jj + nfront] = amult2;
			ibeg = ibeg + nfront + 1;
			iend += nfront;
/* L900: */
		    }
		    ipos = iwpos + npiv;
		    iw[ipos] = -iw[ipos];
		    iw[ipos + 1] = -iw[ipos + 1];
		    npiv += 2;
		    ntotpv += 2;
		    ++nfullb;
		    jpiv = 2;
		}
L910:
		;
	    }
	    if (lnpiv == npiv) {
		goto L930;
	    }
/* L920: */
	}
L930:
	if (npiv == 0) {
	    goto L1730;
	}
	nosure = *nbloc >= ncol - nass && nstruc == 0;
	zone = 0.;
	if (nstruc > 0 && npiv != nstruc) {
	    zone = 1.;
	}
	num1 = 0;
	num2 = 0;
	num3 = 0;
	if (npiv == nfront) {
	    goto L1830;
	}
	if (nstruc > 0) {
	    j2 = ncol;
	    i__2 = ncol;
	    for (j1 = nass + 1; j1 <= i__2; ++j1) {
		ja1 = apos + j1 - 1;
		i__3 = nstruc / 2;
		for (k = 1; k <= i__3; ++k) {
		    if (a[ja1] != 0.) {
			goto L990;
		    }
		    ja1 += nfront;
/* L940: */
		}
		i__3 = j1 + 1;
		for (jj = j2; jj >= i__3; --jj) {
		    ja2 = apos + jj - 1;
		    i__5 = nstruc / 2;
		    for (k = 1; k <= i__5; ++k) {
			if (a[ja2] != 0.) {
			    j2 = jj;
			    goto L970;
			}
			ja2 += nfront;
/* L950: */
		    }
/* L960: */
		}
		j2 = j1 - 1;
		goto L1010;
L970:
		iswop = iw[iwpos - 1 + j1];
		iw[iwpos - 1 + j1] = iw[iwpos - 1 + j2];
		iw[iwpos - 1 + j2] = iswop;
		ppos[iswop] = j2;
		ppos[iw[iwpos - 1 + j1]] = j1;
		ja1 = apos + j1 - 1;
		ja2 = apos + j2 - 1;
		i__3 = nass;
		for (ldummy = 1; ldummy <= i__3; ++ldummy) {
		    swop = a[ja1];
		    a[ja1] = a[ja2];
		    a[ja2] = swop;
		    ja1 += nfront;
		    ja2 += nfront;
/* L980: */
		}
		--j2;
L990:
		if (j1 >= j2) {
		    goto L1010;
		}
/* L1000: */
	    }
L1010:
	    num1 = 0;
	    num2 = j2 - nass;
	    num3 = ncol - j2;
	    if (kind == 3) {
		i__2 = nass + num2;
		for (j1 = nass + 1; j1 <= i__2; ++j1) {
		    ja1 = apos + nstruc / 2 * nfront + j1 - 1;
		    i__3 = nstruc / 2;
		    for (k = 1; k <= i__3; ++k) {
			if (a[ja1] != 0.) {
			    goto L1030;
			}
			ja1 += nfront;
/* L1020: */
		    }
		    goto L1080;
L1030:
		    i__3 = j1 + 1;
		    for (jj = j2; jj >= i__3; --jj) {
			ja2 = apos + nstruc / 2 * nfront + jj - 1;
			i__5 = nstruc / 2;
			for (k = 1; k <= i__5; ++k) {
			    if (a[ja2] != 0.) {
				goto L1050;
			    }
			    ja2 += nfront;
/* L1040: */
			}
			j2 = jj;
			goto L1060;
L1050:
			;
		    }
		    j2 = j1 - 1;
		    goto L1100;
L1060:
		    iswop = iw[iwpos - 1 + j1];
		    iw[iwpos - 1 + j1] = iw[iwpos - 1 + j2];
		    iw[iwpos - 1 + j2] = iswop;
		    ppos[iswop] = j2;
		    ppos[iw[iwpos - 1 + j1]] = j1;
		    ja1 = apos + j1 - 1;
		    ja2 = apos + j2 - 1;
		    i__3 = nass;
		    for (ldummy = 1; ldummy <= i__3; ++ldummy) {
			swop = a[ja1];
			a[ja1] = a[ja2];
			a[ja2] = swop;
			ja1 += nfront;
			ja2 += nfront;
/* L1070: */
		    }
		    --j2;
L1080:
		    if (j1 >= j2) {
			goto L1100;
		    }
/* L1090: */
		}
L1100:
		num1 = j2 - nass;
		num2 -= num1;
	    }
	}
	zcol = 0;
	apos1 = apos + ncol - 1;
	apos2 = apos + ncol;
	i__2 = ncol - (num1 + num2 + nass);
	for (j = 1; j <= i__2; ++j) {
	    i__3 = npiv;
	    for (ja1 = 1; ja1 <= i__3; ++ja1) {
		if (a[apos1 + (ja1 - 1) * nfront] != 0.) {
		    goto L1130;
		}
/* L1110: */
	    }
	    ++zcol;
	    --apos2;
	    if (apos1 == apos2) {
		goto L1130;
	    }
	    i__3 = nass;
	    for (ja1 = 1; ja1 <= i__3; ++ja1) {
		swop = a[apos1 + (ja1 - 1) * nfront];
		a[apos1 + (ja1 - 1) * nfront] = a[apos2 + (ja1 - 1) * nfront];
		a[apos2 + (ja1 - 1) * nfront] = swop;
/* L1120: */
	    }
	    iswop = iw[iwpos + apos2 - apos];
	    iw[iwpos + apos2 - apos] = iw[iwpos + apos1 - apos];
	    iw[iwpos + apos1 - apos] = iswop;
	    ppos[iswop] = apos1 - apos + 1;
	    ppos[iw[iwpos + apos2 - apos]] = apos2 - apos + 1;
L1130:
	    --apos1;
/* L1140: */
	}
	ncol -= zcol;
	num3 -= zcol;
	if (ncol == nass) {
	    goto L1730;
	}
	apos3 = apos + nfront * nass;
	if (npiv == nstruc) {
	    nsc1 = num1 + num2;
	    nsc2 = num2 + num3;
	    if (nsc1 * nsc2 == 0) {
		aposi = apos + nass;
		offdag = apos + nstruc / 2;
		flopsb += (num1 + num3) * nstruc / 2;
		i__2 = nstruc / 2;
		for (k = 1; k <= i__2; ++k) {
		    i__3 = aposi + num1 - 1;
		    for (jj = aposi; jj <= i__3; ++jj) {
			a[jj + nfront * nstruc / 2] = -a[offdag] * a[jj];
			a[jj] = 0.;
/* L1150: */
		    }
		    i__3 = aposi + num1 + num3 - 1;
		    for (jj = aposi + num1; jj <= i__3; ++jj) {
			a[jj] = -(a[offdag] * a[jj + nfront * nstruc / 2]);
			a[jj + nfront * nstruc / 2] = 0.;
/* L1160: */
		    }
		    aposi += nfront;
		    offdag = offdag + nfront + 1;
/* L1170: */
		}
		goto L1730;
	    }
	} else {
	    nsc1 = ncol - nass;
	    nsc2 = ncol - nass;
	}
	if (nosure) {

	    apos4 = astk - (nsc1 + 1) * nsc1 / 2 - 1;
/* Computing MAX */
	    i__2 = rlspa, i__3 = info[2] + apos3 + (nsc1 + 1) * nsc1 / 2 + 2 
		    + nstkac[0];
	    rlspa = max(i__2,i__3);
	    if (apos3 > apos4) {
		ma47sd_(&a[1], &astk, &ptra, &ainput, &ncmpbr);

		apos4 = astk - (nsc1 + 1) * nsc1 / 2 - 1;
		if (apos3 > apos4) {
		    info[2] = info[2] + apos - 1;
		    i__2 = nass * nfront;
		    for (i = 1; i <= i__2; ++i) {
			a[i] = a[apos + i - 1];
/* L1180: */
		    }
		    aposbb = 1;
		    apos = 1;
		    apos3 = apos + nfront * nass;
		    info[1] = -4;
		    if (apos3 > apos4) {
			info[2] = rlspa;
			goto L2170;
		    }
		}
	    }
	    a[astk] = (doublereal) ((nsc1 + 1) * nsc1 / 2 + 2);
	    a[apos4] = (doublereal) ((nsc1 + 1) * nsc1 / 2 + 2);
	    nstkac[0] = nstkac[0] + (nsc1 + 1) * nsc1 / 2 + 2;
	} else {
	    apos4 = apos3 + npiv * nsc2;
	    if (npiv == nstruc) {
/* Computing MAX */
		i__2 = nsc1 * nsc2 + 1, i__3 = num1 * (num2 + num3) + num2 * 
			num3 + num2 * (num2 + 1) / 2 + 2;
		nschur = max(i__2,i__3);
	    } else {
/* Computing MAX */
		i__2 = nsc1 * nsc2 + 1, i__3 = nsc1 * (nsc1 + 1) / 2 + 2;
		nschur = max(i__2,i__3);
	    }
/* Computing MAX */
	    i__2 = rlspa, i__3 = info[2] + apos4 + nschur + nstkac[0];
	    rlspa = max(i__2,i__3);
	    if (apos4 + nschur - 1 > astk) {
		ma47sd_(&a[1], &astk, &ptra, &ainput, &ncmpbr);
		if (apos4 + nschur - 1 > astk) {
		    info[2] = info[2] + apos - 1;
		    i__2 = nass * nfront;
		    for (i = 1; i <= i__2; ++i) {
			a[i] = a[apos + i - 1];
/* L1190: */
		    }
		    aposbb = 1;
		    apos = 1;
		    apos3 = apos + nfront * nass;
		    info[1] = -4;
		    apos4 = apos3 + npiv * nsc2;
		    if (apos4 + nschur - 1 > astk) {
			info[2] = rlspa;
			goto L2170;
		    }
		}
	    }
	    if (zone == 1.) {
		i__2 = apos4 + nsc1 * nsc2 - 1;
		for (k = apos4; k <= i__2; ++k) {
		    a[k] = 0.;
/* L1200: */
		}
	    }
	}
	if (nstruc > 0) {
	    pospv1 = apos;
	    aposi = apos;
	    offdag = pospv1 + nstruc / 2;
	    flopsb += (num1 + (num2 << 2) + num3) * nstruc / 2;
	    i__2 = nstruc / 2;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = aposi + nass + num1 - 1;
		for (jj = aposi + nass; jj <= i__3; ++jj) {
		    a[jj + nfront * nstruc / 2] = -a[offdag] * a[jj];
		    a[jj] = 0.;
/* L1210: */
		}
		apos2 = apos3 + (k - 1) * (num2 + num3);
		i__3 = aposi + nass + num1 + num2 - 1;
		for (jj = aposi + nass + num1; jj <= i__3; ++jj) {
		    a[apos2] = a[jj];
		    a[apos2 + nstruc / 2 * (num2 + num3)] = a[jj + nfront * 
			    nstruc / 2];
		    a[jj] = -(a[pospv1] * a[jj] + a[offdag] * a[jj + nfront * 
			    nstruc / 2]);
		    a[jj + nfront * nstruc / 2] = -a[offdag] * a[apos2];
		    ++apos2;
/* L1220: */
		}
		apos2 = apos3 + (k - 1) * (num2 + num3) + num2;
		i__3 = aposi + ncol - 1;
		for (jj = aposi + nass + num1 + num2; jj <= i__3; ++jj) {
		    a[apos2] = 0.;
		    a[apos2 + nstruc / 2 * (num2 + num3)] = a[jj + nfront * 
			    nstruc / 2];
		    a[jj] = -(a[offdag] * a[jj + nfront * nstruc / 2]);
		    a[jj + nfront * nstruc / 2] = 0.;
		    ++apos2;
/* L1230: */
		}
		pospv1 = pospv1 + nfront + 1;
		aposi += nfront;
		offdag = offdag + nfront + 1;
/* L1240: */
	    }
	    if (npiv == nstruc) {
		apos1 = apos4;
	    } else {
		apos1 = apos4 + num1;
	    }
	    flopsb += (nstruc << 1) * num1 * (num2 + num3);
	    if (num1 >= *nbloc) {
		i__2 = num2 + num3;
		i__3 = num2 + num3;
		dgemm_("N", "T", &i__2, &num1, &nstruc, &c_b371, &a[apos3], &
			i__3, &a[apos + nass], &nfront, &c_b372, &a[apos1], &
			nsc2, 1L, 1L);
	    } else {
		aposi = apos;
		i__2 = nstruc / 2;
		for (k = 1; k <= i__2; ++k) {
		    inc = nsc2;
		    ja1 = apos1;
		    j1 = aposi + nass + nfront * nstruc / 2;
		    i__3 = num1;
		    for (jj = 1; jj <= i__3; ++jj) {
			amult2 = a[j1];
			++j1;
			k2 = apos3 + (k + nstruc / 2 - 1) * (num2 + num3);
			if (k == 1) {
			    i__5 = ja1 + num2 + num3 - 1;
			    for (irow = ja1; irow <= i__5; ++irow) {
				a[irow] = amult2 * a[k2];
				++k2;
/* L1250: */
			    }
			} else {
			    i__5 = ja1 + num2 + num3 - 1;
			    for (irow = ja1; irow <= i__5; ++irow) {
				a[irow] += amult2 * a[k2];
				++k2;
/* L1260: */
			    }
			}
			ja1 += inc;
/* L1270: */
		    }
		    aposi += nfront;
/* L1280: */
		}
	    }
	    flopsb += nstruc * ((num3 << 1) + num2 + 1) * num2;
	    flopsx += nstruc * num2 * num3;
	    apos1 += num1 * nsc2;
	    ja1 = apos1;
	    kblk = num2 / *nbloc;
	    i__2 = kblk;
	    for (k = 1; k <= i__2; ++k) {
		ja1 = (k - 1) * *nbloc;
		flopsx += *nbloc * (*nbloc - 1) * nstruc;
		i__3 = num2 + num3 - ja1;
		i__5 = num2 + num3;
		dgemm_("N", "T", &i__3, nbloc, &nstruc, &c_b371, &a[apos3 + 
			ja1], &i__5, &a[apos + nass + num1 + ja1], &nfront, &
			c_b372, &a[apos1 + ja1 * (nsc2 + 1)], &nsc2, 1L, 1L);
/* L1290: */
	    }
	    aposi = apos;
	    apos1 += kblk * *nbloc * (nsc2 + 1);
	    i__2 = nstruc / 2;
	    for (k = 1; k <= i__2; ++k) {
		inc = nsc2;
		ja1 = apos1;
		ja2 = ja1 + num2 - 1 - kblk * *nbloc;
		j1 = aposi + nass + num1 + kblk * *nbloc;
		i__3 = num2 - kblk * *nbloc;
		for (jj = 1; jj <= i__3; ++jj) {
		    amult1 = a[j1];
		    amult2 = a[j1 + nfront * nstruc / 2];
		    ++j1;
		    k1 = apos3 + kblk * *nbloc + (k - 1) * (num2 + num3) + jj 
			    - 1;
		    k2 = k1 + nstruc / 2 * (num2 + num3);
		    if (k == 1) {
			i__5 = ja2;
			for (irow = ja1; irow <= i__5; ++irow) {
			    a[irow] = amult1 * a[k1] + amult2 * a[k2];
			    ++k1;
			    ++k2;
/* L1300: */
			}
		    } else {
			i__5 = ja2;
			for (irow = ja1; irow <= i__5; ++irow) {
			    a[irow] = a[irow] + amult1 * a[k1] + amult2 * a[
				    k2];
			    ++k1;
			    ++k2;
/* L1310: */
			}
		    }
		    if (k == 1) {
			i__5 = ja2 + num3;
			for (irow = ja2 + 1; irow <= i__5; ++irow) {
			    a[irow] = amult2 * a[k2];
			    ++k2;
/* L1320: */
			}
		    } else {
			i__5 = ja2 + num3;
			for (irow = ja2 + 1; irow <= i__5; ++irow) {
			    a[irow] += amult2 * a[k2];
			    ++k2;
/* L1330: */
			}
		    }
		    ja1 = ja1 + inc + 1;
		    ja2 += inc;
/* L1340: */
		}
		aposi += nfront;
/* L1350: */
	    }
	}
	if (nstruc == npiv) {
	    goto L1560;
	}
	j1 = iwpos + nstruc;
	ltwo = FALSE_;
	pospv1 = apos + nfront * nstruc + nstruc;
	if (nosure) {
	    i = nstruc + 2;
	    aposi = apos + nstruc * nfront + nass;
	    j2 = apos + nfront * nstruc + ncol - 1;
	    aposc = apos4 + 1;
	    if (iw[j1] > 0) {
		flopsb = flopsb + (ncol - nass) + (ncol - nass) * (ncol - 
			nass + 1);
		i__2 = j2;
		for (jj = aposi; jj <= i__2; ++jj) {
		    amult1 = -a[jj] * a[pospv1];
		    i__3 = j2;
		    for (jjj = jj; jjj <= i__3; ++jjj) {
			a[aposc] = amult1 * a[jjj];
			++aposc;
/* L1360: */
		    }
		    a[jj] = amult1;
/* L1370: */
		}
		++j1;
	    } else {
		pospv2 = pospv1 + nfront + 1;
		offdag = pospv1 + 1;
		flopsb = flopsb + (ncol - nass) * 6 + (ncol - nass << 1) * (
			ncol - nass + 1);
		i__2 = j2;
		for (jj = aposi; jj <= i__2; ++jj) {
		    amult1 = -(a[pospv1] * a[jj] + a[offdag] * a[jj + nfront])
			    ;
		    amult2 = -a[pospv2] * a[jj + nfront] - a[offdag] * a[jj];
		    i__3 = j2;
		    for (jjj = jj; jjj <= i__3; ++jjj) {
			a[aposc] = amult1 * a[jjj] + amult2 * a[jjj + nfront];
			++aposc;
/* L1380: */
		    }
		    a[jj] = amult1;
		    a[jj + nfront] = amult2;
/* L1390: */
		}
		j1 += 2;
		pospv1 = pospv2;
		++i;
	    }
	    pospv1 = pospv1 + nfront + 1;
	    i__2 = npiv;
	    for (i1 = i; i1 <= i__2; ++i1) {
		if (ltwo) {
		    goto L1440;
		}
		aposi = apos + (i1 - 1) * nfront + nass;
		j2 = apos + nfront * (i1 - 1) + ncol - 1;
		aposc = apos4 + 1;
		if (iw[j1] > 0) {
		    flopsb = flopsb + (ncol - nass) + (ncol - nass) * (ncol - 
			    nass + 1);
		    i__3 = j2;
		    for (jj = aposi; jj <= i__3; ++jj) {
			amult1 = -a[jj] * a[pospv1];
			i__5 = j2;
			for (jjj = jj; jjj <= i__5; ++jjj) {
			    a[aposc] += amult1 * a[jjj];
			    ++aposc;
/* L1400: */
			}
			a[jj] = amult1;
/* L1410: */
		    }
		    ++j1;
		} else {
		    pospv2 = pospv1 + nfront + 1;
		    offdag = pospv1 + 1;
		    flopsb = flopsb + (ncol - nass) * 6 + (ncol - nass << 1) *
			     (ncol - nass + 1);
		    i__3 = j2;
		    for (jj = aposi; jj <= i__3; ++jj) {
			amult1 = -(a[pospv1] * a[jj] + a[offdag] * a[jj + 
				nfront]);
			amult2 = -a[pospv2] * a[jj + nfront] - a[offdag] * a[
				jj];
			i__5 = j2;
			for (jjj = jj; jjj <= i__5; ++jjj) {
			    a[aposc] = a[aposc] + amult1 * a[jjj] + amult2 * 
				    a[jjj + nfront];
			    ++aposc;
/* L1420: */
			}
			a[jj] = amult1;
			a[jj + nfront] = amult2;
/* L1430: */
		    }
		    j1 += 2;
		    pospv1 = pospv2;
		    ltwo = TRUE_;
		    goto L1450;
		}
L1440:
		ltwo = FALSE_;
		pospv1 = pospv1 + nfront + 1;
L1450:
		;
	    }
	} else {
	    i__2 = npiv;
	    for (i = nstruc + 1; i <= i__2; ++i) {
		if (ltwo) {
		    goto L1480;
		}
		aposi = apos + (i - 1) * nfront + nass;
		poselt = apos3 + (i - nstruc - 1) * (ncol - nass);
		if (iw[j1] > 0) {
		    flopsb += ncol - nass;
		    i__3 = apos + nfront * (i - 1) + ncol - 1;
		    for (jj = aposi; jj <= i__3; ++jj) {
			a[poselt] = a[jj];
			a[jj] = -a[jj] * a[pospv1];
			++poselt;
/* L1460: */
		    }
		    ++j1;
		} else {
		    pospv2 = pospv1 + nfront + 1;
		    offdag = pospv1 + 1;
		    flopsb += (ncol - nass) * 6;
		    i__3 = apos + nfront * (i - 1) + ncol - 1;
		    for (jj = aposi; jj <= i__3; ++jj) {
			a[poselt] = a[jj];
			a[poselt + ncol - nass] = a[jj + nfront];
			a[jj] = -(a[pospv1] * a[jj] + a[offdag] * a[jj + 
				nfront]);
			a[jj + nfront] = -a[pospv2] * a[jj + nfront] - a[
				offdag] * a[poselt];
			++poselt;
/* L1470: */
		    }
		    j1 += 2;
		    pospv1 = pospv2;
		    ltwo = TRUE_;
		    goto L1490;
		}
L1480:
		ltwo = FALSE_;
		pospv1 = pospv1 + nfront + 1;
L1490:
		;
	    }
/* Computing 2nd power */
	    i__2 = ncol - nass;
	    flopsb = flopsb + npiv * (i__2 * i__2) + npiv * (ncol - nass);
	    kblk = (ncol - nass) / *nbloc;
	    if (npiv - nstruc < *nbloc) {
		kblk = 0;
	    }
	    aposh = apos + nstruc * nfront;
	    l = ncol - nass;
	    i__2 = kblk;
	    for (kr = 1; kr <= i__2; ++kr) {
		flopsx += *nbloc * (*nbloc - 1) * (npiv - nstruc);
		i__3 = l - (kr - 1) * *nbloc;
		i__5 = npiv - nstruc;
		dgemm_("N", "T", &i__3, nbloc, &i__5, &c_b371, &a[aposh + 
			nass + *nbloc * (kr - 1)], &nfront, &a[apos3 + *nbloc 
			* (kr - 1)], &l, &zone, &a[apos4 + *nbloc * (l + 1) * 
			(kr - 1)], &l, 1L, 1L);
/* L1500: */
	    }
	    i__2 = l;
	    for (i = kblk * *nbloc + 1; i <= i__2; ++i) {
		aposa = apos + nstruc * nfront + nass + i - 1;
		aposb = apos3 - 1;
		aposc = apos4 + (i - 1) * l - 1;
		if (zone == 1.) {
		    i__3 = l;
		    for (j = i; j <= i__3; ++j) {
			a[aposc + j] += a[aposa] * a[aposb + j];
/* L1510: */
		    }
		} else {
		    i__3 = l;
		    for (j = i; j <= i__3; ++j) {
			a[aposc + j] = a[aposa] * a[aposb + j];
/* L1520: */
		    }
		}
		i__3 = npiv - nstruc;
		for (k = 2; k <= i__3; ++k) {
		    aposa += nfront;
		    aposb += l;
		    i__5 = l;
		    for (j = i; j <= i__5; ++j) {
			a[aposc + j] += a[aposa] * a[aposb + j];
/* L1530: */
		    }
/* L1540: */
		}
/* L1550: */
	    }
	}
L1560:
	if (ncol < nfront) {
	    i__2 = iwpos + nfront - 1;
	    for (jj = iwpos + ncol; jj <= i__2; ++jj) {
		j = iw[jj];
		ppos[j] = *n + 1;
/* L1570: */
	    }
	}
	if (npiv == nstruc) {
	    numm1 = num1;
	    numm2 = num2;
	} else {
	    numm1 = 0;
	    numm2 = ncol - nass;
	}
	nell = elems[iass];
	iell = istk + 1;
	ptrelt = astk + 1;
	i__2 = nell;
	for (elt = 1; elt <= i__2; ++elt) {
	    if (i_dnnt(&a[ptrelt]) < 0) {
		ptrelt -= i_dnnt(&a[ptrelt]);
	    }
	    if (iw[iell] < 0) {
		iell -= iw[iell];
	    }
	    ielln = iell + iw[iell];
	    poselt = ptrelt + i_dnnt(&a[ptrelt]);
	    jp1 = iell + 3;
	    jp2 = jp1 + iw[iell + 1];
	    jp3 = jp2 + iw[iell + 2];
	    jp4 = iell + iw[iell] - 2;
	    ppos[0] = 0;
	    i__3 = jp3 - 1;
	    for (jj = jp2; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L1580;
		}
		if (ppos[iw[jj]] - nass <= numm1 || ppos[iw[jj]] - nass > 
			numm1 + numm2) {
		    goto L1640;
		}
L1580:
		;
	    }
	    apos1 = 0;
	    i__3 = jp2 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L1590;
		}
		if (ppos[iw[jj]] - nass <= numm1) {
		    if (apos1 == 2) {
			goto L1640;
		    }
		    apos1 = 1;
		}
		if (ppos[iw[jj]] - nass > numm1 + numm2) {
		    if (ppos[iw[jj]] == *n + 1) {
			goto L1640;
		    }
		    if (apos1 == 1) {
			goto L1640;
		    }
		    apos1 = 2;
		}
L1590:
		;
	    }
	    i__3 = jp4;
	    for (jj = jp3; jj <= i__3; ++jj) {
		if (iw[jj] == 0) {
		    goto L1600;
		}
		if (ppos[iw[jj]] - nass <= numm1) {
		    if (apos1 == 1) {
			goto L1640;
		    }
		}
		if (ppos[iw[jj]] - nass > numm1 + numm2) {
		    if (ppos[iw[jj]] == *n + 1) {
			goto L1640;
		    }
		    if (apos1 == 2) {
			goto L1640;
		    }
		}
L1600:
		;
	    }
	    apos1 = ptrelt + 1;
	    ja2 = jp2;
	    i__3 = jp3 - 1;
	    for (jj = jp1; jj <= i__3; ++jj) {
		j = iw[jj];
		if (j == 0) {
		    if (jj <= jp2) {
			apos1 = apos1 + jp4 - jp2 + 1;
		    } else {
			apos1 = apos1 + jp4 - jj + 1;
		    }
		    goto L1630;
		}
		if (jj > jp2) {
		    ja2 = jj;
		}
		if (nosure) {
		    i__5 = jp4;
		    for (jjj = ja2; jjj <= i__5; ++jjj) {
			jay = iw[jjj];
			if (ppos[jay] >= ppos[j]) {
			    apos2 = apos4 + 1 + (ppos[j] - nass - 1) * ((nsc2 
				    << 1) - ppos[j] + nass + 2) / 2 + ppos[
				    jay] - ppos[j];
			    a[apos2] += a[apos1];
			} else if (jay > 0) {
			    apos2 = apos4 + 1 + (ppos[jay] - nass - 1) * ((
				    nsc2 << 1) - ppos[jay] + nass + 2) / 2 + 
				    ppos[j] - ppos[jay];
			    a[apos2] += a[apos1];
			}
			++apos1;
/* L1610: */
		    }
		} else {
		    i__5 = jp4;
		    for (jjj = ja2; jjj <= i__5; ++jjj) {
			jay = iw[jjj];
			if (ppos[jay] >= ppos[j]) {
			    apos2 = apos4 + (ppos[j] - nass - 1) * nsc2 + 
				    ppos[jay] - nass - 1;
			    if (nstruc == npiv) {
				apos2 -= num1;
			    }
			    a[apos2] += a[apos1];
			} else if (jay > 0) {
			    apos2 = apos4 + (ppos[jay] - nass - 1) * nsc2 + 
				    ppos[j] - nass - 1;
			    if (nstruc == npiv) {
				apos2 -= num1;
			    }
			    a[apos2] += a[apos1];
			}
			++apos1;
/* L1620: */
		    }
		}
L1630:
		;
	    }
	    ppos[0] = *n + 1;
	    --nell;
	    --nnell;
	    --elems[iass];
	    liell = iw[iell];
	    nstkac[1] -= liell;
	    j1 = iell;
	    j2 = ielln - 1;
	    if (ielln < ptrirn - iinput) {
		if (iw[ielln] < 0) {
		    liell -= iw[ielln];
		    j2 -= iw[ielln];
		}
	    }
	    if (iell - 1 == istk) {
		istk += liell;
	    } else {
		if (iw[iell - 1] < 0) {
		    liell -= iw[iell - 1];
		    j1 = iell + iw[iell - 1];
		}
		iw[j1] = -liell;
		iw[j2] = -liell;
	    }
	    laell = i_dnnt(&a[ptrelt]);
	    nstkac[0] -= laell;
	    ja1 = ptrelt;
	    ja2 = poselt - 1;
	    if (poselt < ptra - ainput) {
		if (i_dnnt(&a[poselt]) < 0) {
		    laell -= i_dnnt(&a[poselt]);
		    ja2 -= i_dnnt(&a[poselt]);
		}
	    }
	    if (ptrelt - 1 == astk && ! nosure) {
		astk += laell;
	    } else {
		if (i_dnnt(&a[ptrelt - 1]) < 0) {
		    laell -= i_dnnt(&a[ptrelt - 1]);
		    ja1 = ptrelt + i_dnnt(&a[ptrelt - 1]);
		}
		a[ja1] = (doublereal) (-laell);
		a[ja2] = (doublereal) (-laell);
	    }
L1640:
	    iell = ielln;
	    ptrelt = poselt;
/* L1650: */
	}
	ppos[0] = *n + 1;
	++elems[father[iass]];
	if (nosure) {
	    astk = apos4 - 1;
	    laell = (nsc1 + 1) * nsc1 / 2 + 2;
	} else {
	    if (nstruc == npiv) {
		laell = num1 * (num2 + num3) + num2 * (num2 + 1) / 2 + num2 * 
			num3 + 2;
		ja1 = apos4 + (num1 + num2) * (num2 + num3) - 1;
	    } else {
		laell = (nsc1 + 1) * nsc1 / 2 + 2;
		ja1 = apos4 + nsc1 * nsc1 - 1;
	    }
	    a[astk] = (doublereal) laell;
	    nstkac[0] += laell;
	    --astk;
	    if (nstruc == npiv) {
		i__2 = num2;
		for (i = 1; i <= i__2; ++i) {
		    i__3 = ja1 - num3 - i + 1;
		    for (jj = ja1; jj >= i__3; --jj) {
			a[astk] = a[jj];
			--astk;
/* L1660: */
		    }
		    ja1 -= num2 + num3;
/* L1670: */
		}
		i__2 = num1;
		for (i = 1; i <= i__2; ++i) {
		    i__3 = num2 + num3;
		    for (jj = 1; jj <= i__3; ++jj) {
			a[astk] = a[ja1];
			--astk;
			--ja1;
/* L1680: */
		    }
/* L1690: */
		}
	    } else {
		i__2 = nsc1;
		for (i = 1; i <= i__2; ++i) {
		    i__3 = ja1 - i + 1;
		    for (jj = ja1; jj >= i__3; --jj) {
			a[astk] = a[jj];
			--astk;
/* L1700: */
		    }
		    ja1 -= nsc1;
/* L1710: */
		}
	    }
	    a[astk] = (doublereal) laell;
	    --astk;
	}
	liell = ncol - nass + 4;
	nstkac[1] += liell;
/* Computing MAX */
	i__2 = intspa, i__3 = iwpos + (nfront << 1) + liell + nstkac[1];
	intspa = max(i__2,i__3);
	if (iwpos + (nfront << 1) + liell > istk) {
	    ma47pd_(&iw[1], &istk, &ptrirn, &iinput, &ncmpbi);
	    if (iwpos + (nfront << 1) + liell > istk) {
		info[2] = intspa;
		info[1] = -3;
		goto L2170;
	    }
	}
	iw[istk] = liell;
	--istk;
	++nnell;
	i__2 = ncol - nass;
	for (i = 1; i <= i__2; ++i) {
	    iw[istk] = iw[iwpos + ncol - i];
	    --istk;
/* L1720: */
	}
	if (nstruc == npiv) {
	    iw[istk] = num2;
	    iw[istk - 1] = num1;
	} else {
	    iw[istk] = ncol - nass;
	    iw[istk - 1] = 0;
	}
	iw[istk - 2] = liell;
	istk += -3;
L1730:
	if (npiv + nnull < nass) {
	    zcol = 0;
	    apos1 = apos + npiv * nfront + nass - 1;
	    i__2 = nfront;
	    for (j = nass + 1; j <= i__2; ++j) {
		++apos1;
		i__3 = nass - npiv;
		for (ja1 = 1; ja1 <= i__3; ++ja1) {
		    if (a[apos1 + (ja1 - 1) * nfront] != 0.) {
			goto L1750;
		    }
/* L1740: */
		}
		iw[iwpos + j - 1] = -iw[iwpos + j - 1];
		++zcol;
L1750:
		;
	    }
	    ++elems[father[iass]];
	    laell = (nass - npiv) * (nfront - npiv - zcol) + 2;
/* Computing MAX */
	    i__2 = rlspa, i__3 = info[2] + apos + nass * nfront + laell + 
		    nstkac[0];
	    rlspa = max(i__2,i__3);
	    if (apos + nass * nfront + laell > astk) {
		ma47sd_(&a[1], &astk, &ptra, &ainput, &ncmpbr);
		if (apos + nass * nfront + laell > astk) {
		    info[2] = info[2] + apos - 1;
		    i__2 = nass * nfront;
		    for (i = 1; i <= i__2; ++i) {
			a[i] = a[apos + i - 1];
/* L1760: */
		    }
		    aposbb = 1;
		    apos = 1;
		    info[1] = -4;
		    if (nass * nfront + laell > astk) {
			info[2] = rlspa;
			goto L2170;
		    }
		}
	    }
	    a[astk] = (doublereal) laell;
	    nstkac[0] += laell;
	    astk -= laell;
	    a[astk + 1] = (doublereal) laell;
	    ja1 = astk + 2;
	    ja2 = apos + npiv * nfront + npiv;
	    if (zcol == 0) {
		i__2 = nass - npiv;
		for (i = 1; i <= i__2; ++i) {
		    i__3 = ja2 + nfront - npiv - 1;
		    for (jj = ja2; jj <= i__3; ++jj) {
			a[ja1] = a[jj];
			++ja1;
/* L1770: */
		    }
		    ja2 += nfront;
/* L1780: */
		}
	    } else {
		i__2 = nass - npiv;
		for (i = 1; i <= i__2; ++i) {
		    i__3 = ja2 + nfront - npiv - 1;
		    for (jj = ja2; jj <= i__3; ++jj) {
			if (iw[iwpos + npiv + jj - ja2] > 0) {
			    a[ja1] = a[jj];
			    ++ja1;
			}
/* L1790: */
		    }
		    ja2 += nfront;
/* L1800: */
		}
	    }
	    liell = nfront - npiv - zcol + 4;
/* Computing MAX */
	    i__2 = intspa, i__3 = iwpos + (nfront << 1) + 2 + liell + nstkac[
		    1];
	    intspa = max(i__2,i__3);
	    if (iwpos + (nfront << 1) + 2 + liell > istk) {
		ma47pd_(&iw[1], &istk, &ptrirn, &iinput, &ncmpbi);
		if (iwpos + (nfront << 1) + 2 + liell > istk) {
		    info[2] = intspa;
		    info[1] = -3;
		    goto L2170;
		}
	    }
	    ++nnell;
	    iw[istk] = liell;
	    nstkac[1] += liell;
	    istk -= liell;
	    iw[istk + 1] = liell;
	    iw[istk + 2] = -1;
	    iw[istk + 3] = nass - npiv;
	    j1 = istk + 4;
	    if (zcol == 0) {
		i__2 = iwpos + nfront - 1;
		for (jj = iwpos + npiv; jj <= i__2; ++jj) {
		    iw[j1] = iw[jj];
		    ++j1;
/* L1810: */
		}
	    } else {
		i__2 = iwpos + nfront - 1;
		for (jj = iwpos + npiv; jj <= i__2; ++jj) {
		    if (iw[jj] <= 0) {
			iw[jj] = -iw[jj];
		    } else {
			iw[j1] = iw[jj];
			++j1;
		    }
/* L1820: */
		}
	    }
	}
L1830:
	i__2 = iwpos + nfront - 1;
	for (jj = iwpos + npiv; jj <= i__2; ++jj) {
	    j = (i__3 = iw[jj], abs(i__3));
	    ppos[j] = *n + 1;
/* L1840: */
	}
	if (father[iass] <= *nodes) {
	    elems[father[iass]] += nell;
	}
/* ******************************** */
/* ******************************** */
	if (npiv == 0) {
	    goto L2160;
	}
	++nblk;
	zcol = 0;
	apos1 = apos + npiv;
	i__2 = nass - npiv;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = npiv;
	    for (ja1 = 1; ja1 <= i__3; ++ja1) {
		if (a[apos1 + (ja1 - 1) * nfront] != 0.) {
		    goto L1870;
		}
/* L1850: */
	    }
	    ++zcol;
	    apos2 = apos + nass - zcol;
	    i__3 = npiv;
	    for (ja1 = 1; ja1 <= i__3; ++ja1) {
		swop = a[apos1 + (ja1 - 1) * nfront];
		a[apos1 + (ja1 - 1) * nfront] = a[apos2 + (ja1 - 1) * nfront];
		a[apos2 + (ja1 - 1) * nfront] = swop;
/* L1860: */
	    }
	    iswop = iw[iwpos + npiv + j - zcol];
	    iw[iwpos + npiv + j - zcol] = iw[iwpos + nass - zcol];
	    iw[iwpos + nass - zcol] = iswop;
	    goto L1880;
L1870:
	    ++apos1;
L1880:
	    ;
	}
	if (zcol > 0) {
	    apos1 = apos + nass - zcol;
	    apos2 = apos + nass;
	    i__2 = ncol - nass;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = npiv;
		for (ja1 = 1; ja1 <= i__3; ++ja1) {
		    a[apos1 + (ja1 - 1) * nfront] = a[apos2 + (ja1 - 1) * 
			    nfront];
/* L1890: */
		}
		iw[iwpos + nass - zcol + j - 1] = iw[iwpos + nass + j - 1];
		++apos1;
		++apos2;
/* L1900: */
	    }
	}
	ncol -= zcol;
	nass -= zcol;
	if (nstruc > 0) {
	    iw[iwpos - 4] = -ncol;
	    iw[iwpos - 3] = nstruc;
	    if (kind == 3) {
		iw[iwpos - 3] = -nstruc;
	    }
	    iw[iwpos - 2] = num3;
	    if (kind == 2) {
		i__2 = ncol;
		for (i = 1; i <= i__2; ++i) {
		    iw[iwpos + i - 2] = iw[iwpos + i - 1];
/* L1910: */
		}
		iwpos = iwpos + ncol + 3;
		num1 = 0;
	    } else {
		iw[iwpos - 1] = num1;
		iwpos = iwpos + ncol + 4;
	    }
	    if (npiv != nstruc) {
		++nblk;
		j1 = iwpos - ncol - 4 + nstruc;
		j2 = iwpos - 2;
		i__2 = ncol;
		for (i = nstruc + 1; i <= i__2; ++i) {
		    iw[j2] = iw[j1];
		    ++j2;
		    ++j1;
/* L1920: */
		}
	    }
	    if (nstruc < nass) {
		if (num1 != 0) {
		    j1 = iwpos - 4 - ncol + nstruc;
		    j2 = j1;
		    i__2 = nass - nstruc;
		    for (i = 1; i <= i__2; ++i) {
			ppos[i] = iw[j2];
			++j2;
/* L1930: */
		    }
		    j2 = j1 + nass - nstruc;
		    i__2 = num1;
		    for (i = 1; i <= i__2; ++i) {
			iw[j1] = iw[j2];
			++j1;
			++j2;
/* L1940: */
		    }
		    i__2 = nass - nstruc;
		    for (i = 1; i <= i__2; ++i) {
			iw[j1] = ppos[i];
			ppos[i] = *n + 1;
			++j1;
/* L1950: */
		    }
		}
	    }
	    if (npiv == nstruc) {
		goto L1970;
	    }
	    iwpos += -2;
	} else {
	    j1 = iwpos;
	    j2 = iwpos - 2;
	    i__2 = ncol;
	    for (i = 1; i <= i__2; ++i) {
		iw[j2] = iw[j1];
		++j2;
		++j1;
/* L1960: */
	    }
	    iwpos += -2;
	}
	iw[iwpos - 2] = ncol - nstruc;
	iw[iwpos - 1] = npiv - nstruc;
	iwpos = iwpos + ncol - nstruc + 4;
L1970:
	if (info[1] == -4) {
	    info[2] += (npiv - nstruc) * ((ncol << 1) - npiv + nstruc + 1) / 
		    2;
	    goto L2160;
	}
	apos2 = aposbb;
	i__2 = nstruc / 2;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + (i - 1) * nfront + nstruc / 2;
	    i__3 = nstruc / 2;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L1980: */
	    }
/* L1990: */
	}
	if (kind == 2) {
	    i__2 = nstruc / 2;
	    for (i = 1; i <= i__2; ++i) {
		ja1 = apos + (nstruc / 2 + i - 1) * (nfront + 1) + 1;
		i__3 = nstruc / 2;
		for (j = i + 1; j <= i__3; ++j) {
		    a[apos2] = a[ja1];
		    ++apos2;
		    ++ja1;
/* L2000: */
		}
/* L2010: */
	    }
	    ja1 = apos;
	    i__2 = nstruc / 2;
	    for (i = 1; i <= i__2; ++i) {
		a[apos2] = a[ja1];
		ja1 = ja1 + nfront + 1;
		++apos2;
/* L2020: */
	    }
	    i__2 = nstruc / 2;
	    for (i = 1; i <= i__2; ++i) {
		a[apos2] = 0.;
		++apos2;
/* L2030: */
	    }
	} else {
	    i__2 = nstruc;
	    for (i = 1; i <= i__2; ++i) {
		a[apos2] = 0.;
		++apos2;
/* L2040: */
	    }
	}
	i__2 = nstruc / 2;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + (i - 1) * nfront + nstruc;
	    i__3 = nass;
	    for (j = nstruc + 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2050: */
	    }
	    ja1 += num1;
	    i__3 = ncol - nass - num1;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2060: */
	    }
/* L2070: */
	}
	i__2 = nstruc;
	for (i = nstruc / 2 + 1; i <= i__2; ++i) {
	    ja1 = apos + (i - 1) * nfront + nass;
	    i__3 = num1;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2080: */
	    }
	    ja1 = apos + (i - 1) * nfront + nstruc;
	    i__3 = nass - nstruc;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2090: */
	    }
	    ja1 = apos + (i - 1) * nfront + nass + num1;
	    i__3 = ncol - nass - num1 - num3;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2100: */
	    }
/* L2110: */
	}
	i__2 = npiv - nstruc;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + (nstruc + i - 1) * (nfront + 1);
	    i__3 = npiv - nstruc;
	    for (j = i; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2120: */
	    }
/* L2130: */
	}
	i__2 = npiv - nstruc;
	for (i = 1; i <= i__2; ++i) {
	    ja1 = apos + (nstruc + i - 1) * nfront + npiv;
	    i__3 = ncol - npiv;
	    for (j = 1; j <= i__3; ++j) {
		a[apos2] = a[ja1];
		++apos2;
		++ja1;
/* L2140: */
	    }
/* L2150: */
	}
	aposbb = apos2;
L2160:
	;
    }
    if (info[1] == -4) {
	info[2] = rlspa;
	goto L2170;
    }
    nrlbdu = aposbb - 1;
    nirbdu = iwpos - 5;
    iw[1] = nrlbdu + 1;
    iw[2] = nrlbdu + nfullb;
    iw[3] = nblk;
    ma47wd_(&a[1], la, &iw[1], liw, &nrlbdu);
/* Computing MAX */
    i__1 = *ne << 1;
    info[6] = max(i__1,rlspa);
/* Computing MAX */
    i__1 = *ne << 1;
    info[7] = max(i__1,intspa);
    info[15] = maxfrt;
    info[16] = nrlbdu;
    info[17] = nirbdu;
    info[18] = ncmpbr;
    info[19] = ncmpbi;
    info[20] = ntileb;
    info[21] = noxob;
    info[22] = nfullb;
    info[23] = neig;
    info[24] = *n - ntotpv;
    rinfo[3] = flopsb;
    rinfo[4] = flopsx;
L2170:
    return 0;
} /* ma47od_ */

/* Subroutine */ int ma47pd_(iw, bottom, top, move, ncmp)
integer *iw, *bottom, *top, *move, *ncmp;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer size, k, jj, idummy;

    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ++(*ncmp);
    k = *top - *move - 1;
    i__1 = *top;
    for (idummy = 1; idummy <= i__1; ++idummy) {
	if (k <= *bottom) {
	    goto L30;
	}
	size = iw[k];
	if (size < 0) {
	    *move -= size;
	    k += size;
	} else {
	    if (*move > 0) {
		i__2 = k - size + 1;
		for (jj = k; jj >= i__2; --jj) {
		    iw[jj + *move] = iw[jj];
/* L10: */
		}
	    }
	    k -= size;
	}
/* L20: */
    }
L30:
    *bottom += *move;
    *move = 0;
    return 0;
} /* ma47pd_ */

/* Subroutine */ int ma47qd_(n, a, la, iw, liw, w, rhs, iw1, icntl)
integer *n;
doublereal *a;
integer *la, *iw, *liw;
doublereal *w, *rhs;
integer *iw1, *icntl;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer iblk, kind, apos, irhs, ipiv, i, j, k;
    extern /* Subroutine */ int dgemv_();
    static integer ncols, aposm, aposo;
    extern /* Subroutine */ int dtpmv_();
    static integer j1, j2;
    extern /* Subroutine */ int dtpsv_();
    static integer iwpos;
    extern /* Subroutine */ int dtrsv_();
    static integer n2, nrows;
    static doublereal w1;
    static integer num1, num3;

    /* Parameter adjustments */
    --icntl;
    --iw1;
    --rhs;
    --w;
    --iw;
    --a;

    /* Function Body */
    apos = 1;
    iwpos = 4;
    i__1 = iw[3];
    for (iblk = 1; iblk <= i__1; ++iblk) {
	iw1[iblk] = iwpos;
	ncols = iw[iwpos];
	nrows = iw[iwpos + 1];
	iwpos += 2;
	num1 = 0;
	num3 = 0;
	if (ncols > 0) {
	    kind = 1;
	    n2 = nrows;
	} else {
	    ncols = -ncols;
	    num3 = iw[iwpos];
	    ++iwpos;
	    if (nrows > 0) {
		kind = 2;
	    } else {
		nrows = -nrows;
		num1 = iw[iwpos];
		++iwpos;
		kind = 3;
	    }
	    n2 = nrows / 2;
	}
	if (n2 > icntl[7]) {
	    i__2 = ncols;
	    for (i = 1; i <= i__2; ++i) {
		w[i] = rhs[i__3 = iw[iwpos + i - 1], abs(i__3)];
/* L10: */
	    }
	    if (kind == 1) {
		dtpsv_("L", "N", "U", &nrows, &a[apos], &w[1], &c__1, 1L, 1L, 
			1L);
		apos += nrows * (nrows + 1) / 2;
		if (ncols > nrows) {
		    i__2 = ncols - nrows;
		    i__3 = ncols - nrows;
		    dgemv_("N", &i__2, &nrows, &c_b371, &a[apos], &i__3, &w[1]
			    , &c__1, &c_b371, &w[nrows + 1], &c__1, 1L);
		}
		apos += nrows * (ncols - nrows);
		i__2 = ncols;
		for (i = 1; i <= i__2; ++i) {
		    rhs[i__3 = iw[iwpos + i - 1], abs(i__3)] = w[i];
/* L35: */
		}
	    } else if (kind == 2) {
		dtrsv_("U", "T", "U", &n2, &a[apos], &n2, &w[1], &c__1, 1L, 
			1L, 1L);
		aposm = apos + n2 * n2;
		aposo = aposm + n2 * (n2 + 1) / 2 + n2;
		k = ncols - nrows;
		if (k > num1) {
		    i__2 = k - num1;
		    i__3 = k - num1;
		    dgemv_("N", &i__2, &n2, &c_b371, &a[aposo], &i__3, &w[1], 
			    &c__1, &c_b371, &w[nrows + 1 + num1], &c__1, 1L);
		}
		i__2 = n2;
		for (i = 1; i <= i__2; ++i) {
		    rhs[i__3 = iw[iwpos + i - 1], abs(i__3)] = w[i];
/* L40: */
		}
		i__2 = n2 - 1;
		dtpmv_("L", "N", "N", &i__2, &a[aposm], &w[1], &c__1, 1L, 1L, 
			1L);
		i__2 = nrows;
		for (i = n2 + 2; i <= i__2; ++i) {
		    w[i] += w[i - n2 - 1];
/* L50: */
		}
		dtrsv_("L", "N", "U", &n2, &a[apos], &n2, &w[n2 + 1], &c__1, 
			1L, 1L, 1L);
		apos = aposo + n2 * (k - num1);
		if (k > num3) {
		    i__2 = k - num3;
		    i__3 = k - num3;
		    dgemv_("N", &i__2, &n2, &c_b371, &a[apos], &i__3, &w[n2 + 
			    1], &c__1, &c_b371, &w[nrows + 1], &c__1, 1L);
		}
		apos += n2 * (k - num3);
		i__2 = ncols;
		for (i = n2 + 1; i <= i__2; ++i) {
		    rhs[i__3 = iw[iwpos + i - 1], abs(i__3)] = w[i];
/* L60: */
		}
	    } else {
		dtrsv_("U", "T", "U", &n2, &a[apos], &n2, &w[1], &c__1, 1L, 
			1L, 1L);
		dtrsv_("L", "N", "U", &n2, &a[apos], &n2, &w[n2 + 1], &c__1, 
			1L, 1L, 1L);
		apos = apos + n2 * n2 + (n2 << 1);
		k = ncols - nrows;
		if (k > num1) {
		    i__2 = k - num1;
		    i__3 = k - num1;
		    dgemv_("N", &i__2, &n2, &c_b371, &a[apos], &i__3, &w[1], &
			    c__1, &c_b371, &w[nrows + 1 + num1], &c__1, 1L);
		}
		apos += n2 * (k - num1);
		if (k > num3) {
		    i__2 = k - num3;
		    i__3 = k - num3;
		    dgemv_("N", &i__2, &n2, &c_b371, &a[apos], &i__3, &w[n2 + 
			    1], &c__1, &c_b371, &w[nrows + 1], &c__1, 1L);
		}
		apos += n2 * (k - num3);
		i__2 = ncols;
		for (i = 1; i <= i__2; ++i) {
		    rhs[i__3 = iw[iwpos + i - 1], abs(i__3)] = w[i];
/* L90: */
		}
	    }
	} else {
	    j1 = iwpos;
	    j2 = iwpos + nrows - 1;
	    if (kind == 1) {
		i__2 = nrows;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    ++apos;
		    w1 = rhs[i__3 = iw[j1], abs(i__3)];
		    ++j1;
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			rhs[irhs] -= a[apos] * w1;
			++apos;
/* L100: */
		    }
/* L130: */
		}
		j2 = iwpos + ncols - 1;
		i__2 = nrows;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    w1 = rhs[i__3 = iw[iwpos + ipiv - 1], abs(i__3)];
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			rhs[irhs] += w1 * a[apos];
			++apos;
/* L133: */
		    }
/* L136: */
		}
	    } else {
		j2 = iwpos + n2 - 1;
		aposo = apos + n2 * n2;
		i__2 = n2 - 1;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    k = apos + (ipiv - 1) * (n2 + 1);
		    w1 = rhs[-iw[j1]];
		    ++j1;
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			k += n2;
			rhs[irhs] -= w1 * a[k];
/* L140: */
		    }
/* L145: */
		}
		if (kind == 2) {
		    k = aposo;
		    aposo = aposo + n2 * (n2 + 1) / 2 + n2;
		    j1 = iwpos;
		    i__2 = n2 - 1;
		    for (ipiv = 1; ipiv <= i__2; ++ipiv) {
			w1 = rhs[-iw[j1]];
			++j1;
			i__3 = j2;
			for (j = j1; j <= i__3; ++j) {
			    irhs = (i__4 = iw[j + n2], abs(i__4));
			    rhs[irhs] += w1 * a[k];
			    ++k;
/* L150: */
			}
/* L155: */
		    }
		} else {
		    aposo += n2 << 1;
		}
		j1 = iwpos + n2;
		j2 += n2;
		i__2 = n2 - 1;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    k = apos + (ipiv - 1) * (n2 + 1);
		    w1 = rhs[-iw[j1]];
		    ++j1;
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			++k;
			rhs[irhs] -= w1 * a[k];
/* L160: */
		    }
/* L165: */
		}
		apos = aposo;
		j1 = iwpos + nrows;
		j2 = iwpos + ncols - 1;
		i__2 = n2;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    w1 = rhs[i__3 = iw[iwpos + ipiv - 1], abs(i__3)];
		    i__3 = j2;
		    for (j = j1 + num1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			rhs[irhs] += w1 * a[apos];
			++apos;
/* L190: */
		    }
/* L195: */
		}
		i__2 = n2;
		for (ipiv = 1; ipiv <= i__2; ++ipiv) {
		    w1 = rhs[i__3 = iw[iwpos + ipiv - 1 + n2], abs(i__3)];
		    i__3 = j2 - num3;
		    for (j = j1; j <= i__3; ++j) {
			irhs = (i__4 = iw[j], abs(i__4));
			rhs[irhs] += w1 * a[apos];
			++apos;
/* L200: */
		    }
/* L210: */
		}
	    }
	}
	iwpos += ncols;
/* L270: */
    }
} /* ma47qd_ */

/* Subroutine */ int ma47rd_(n, a, la, iw, liw, w, rhs, iw1, icntl)
integer *n;
doublereal *a;
integer *la, *iw, *liw;
doublereal *w, *rhs;
integer *iw1, *icntl;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer iblk, kind, apos, irhs, ipiv, jpiv, lrow, apos2, irhs1, 
	    irhs2, i, j, k;
    extern /* Subroutine */ int dgemv_();
    static integer aposd, apose, aposf, ncols, aposm;
    extern /* Subroutine */ int dtpmv_();
    static integer j1, j2;
    extern /* Subroutine */ int dtpsv_();
    static integer iwpos;
    extern /* Subroutine */ int dtrsv_();
    static integer n2, nrows;
    static doublereal w1;
    static integer num1, num3;

    /* Parameter adjustments */
    --icntl;
    --iw1;
    --rhs;
    --w;
    --iw;
    --a;

    /* Function Body */
    apos = iw[1];
    apos2 = iw[2];
    for (iblk = iw[3]; iblk >= 1; --iblk) {
	iwpos = iw1[iblk];
	ncols = (i__1 = iw[iwpos], abs(i__1));
	nrows = (i__1 = iw[iwpos + 1], abs(i__1));
	num1 = 0;
	num3 = 0;
	if (iw[iwpos] > 0) {
	    kind = 1;
	    apos -= nrows * (ncols - nrows);
	    n2 = nrows;
	} else {
	    n2 = nrows / 2;
	    apos -= nrows * (ncols + ncols - nrows + 1) / 2;
	    num3 = iw[iwpos + 2];
	    apos = apos + num3 * n2 + n2 * (n2 + 1) / 2;
	    ++iwpos;
	    if (iw[iwpos] > 0) {
		apos -= n2;
		kind = 2;
	    } else {
		num1 = iw[iwpos + 2];
		apos = apos + num1 * n2 + n2 * (n2 + 1) / 2;
		apos -= n2 << 1;
		++iwpos;
		kind = 3;
	    }
	}
	iwpos += 2;
	if (n2 > icntl[7]) {
	    i__1 = ncols;
	    for (i = nrows + 1; i <= i__1; ++i) {
		w[i] = rhs[i__2 = iw[iwpos + i - 1], abs(i__2)];
/* L5: */
	    }
	    if (kind == 1) {
		for (ipiv = nrows; ipiv >= 1; --ipiv) {
		    irhs = (i__1 = iw[iwpos + ipiv - 1], abs(i__1));
		    apos -= nrows + 1 - ipiv;
		    w[ipiv] = rhs[irhs] * a[apos];
/* L10: */
		}
		jpiv = -1;
		for (ipiv = nrows; ipiv >= 1; --ipiv) {
		    irhs = iw[iwpos + ipiv - 1];
		    if (irhs < 0) {
			irhs1 = -iw[iwpos + ipiv - 1 + jpiv];
			w[ipiv] = rhs[irhs1] * a[apos2] + w[ipiv];
			if (jpiv == 1) {
			    --apos2;
			}
			jpiv = -jpiv;
		    }
/* L20: */
		}
		k = ncols - nrows;
		if (k > 0) {
		    dgemv_("T", &k, &nrows, &c_b371, &a[apos + nrows * (nrows 
			    + 1) / 2], &k, &w[nrows + 1], &c__1, &c_b371, &w[
			    1], &c__1, 1L);
		}
		dtpsv_("L", "T", "U", &nrows, &a[apos], &w[1], &c__1, 1L, 1L, 
			1L);
		i__1 = nrows;
		for (i = 1; i <= i__1; ++i) {
		    rhs[i__2 = iw[iwpos + i - 1], abs(i__2)] = w[i];
/* L60: */
		}
	    } else {
		aposd = apos;
		aposm = apos + n2 * n2;
		apose = aposm + n2 * (n2 - 1) / 2;
		if (kind == 3) {
		    apose = aposm;
		}
		aposf = apose + n2;
		i__1 = n2;
		for (ipiv = 1; ipiv <= i__1; ++ipiv) {
		    irhs1 = -iw[iwpos + ipiv - 1];
		    irhs2 = -iw[iwpos + ipiv - 1 + n2];
		    w[ipiv] = rhs[irhs1] * a[apose] + rhs[irhs2] * a[aposd];
		    w[ipiv + n2] = rhs[irhs1] * a[aposd] + rhs[irhs2] * a[
			    aposf];
		    aposd = aposd + n2 + 1;
		    ++apose;
		    ++aposf;
/* L70: */
		}
		aposd = apos;
		apos = aposf;
		k = ncols - nrows;
		if (k > num1) {
		    i__1 = k - num1;
		    i__2 = k - num1;
		    dgemv_("T", &i__1, &n2, &c_b371, &a[apos], &i__2, &w[
			    nrows + num1 + 1], &c__1, &c_b371, &w[1], &c__1, 
			    1L);
		}
		if (k > num3) {
		    i__1 = k - num3;
		    i__2 = k - num3;
		    dgemv_("T", &i__1, &n2, &c_b371, &a[apos + n2 * (k - num1)
			    ], &i__2, &w[nrows + 1], &c__1, &c_b371, &w[n2 + 
			    1], &c__1, 1L);
		}
		if (kind == 2) {
		    apos = aposd;
		    dtrsv_("L", "T", "U", &n2, &a[aposd], &n2, &w[n2 + 1], &
			    c__1, 1L, 1L, 1L);
		    i__1 = nrows;
		    for (i = n2 + 1; i <= i__1; ++i) {
			rhs[i__2 = iw[iwpos + i - 1], abs(i__2)] = w[i];
/* L90: */
		    }
		    i__1 = n2 - 1;
		    dtpmv_("L", "T", "N", &i__1, &a[aposm], &w[n2 + 2], &c__1,
			     1L, 1L, 1L);
		    i__1 = n2 - 1;
		    for (i = 1; i <= i__1; ++i) {
			w[i] += w[i + n2 + 1];
/* L100: */
		    }
		    dtrsv_("U", "N", "U", &n2, &a[aposd], &n2, &w[1], &c__1, 
			    1L, 1L, 1L);
		    i__1 = n2;
		    for (i = 1; i <= i__1; ++i) {
			rhs[i__2 = iw[iwpos + i - 1], abs(i__2)] = w[i];
/* L110: */
		    }
		} else {
		    apos = aposd;
		    dtrsv_("L", "T", "U", &n2, &a[apos], &n2, &w[n2 + 1], &
			    c__1, 1L, 1L, 1L);
		    dtrsv_("U", "N", "U", &n2, &a[apos], &n2, &w[1], &c__1, 
			    1L, 1L, 1L);
		    i__1 = nrows;
		    for (i = 1; i <= i__1; ++i) {
			rhs[i__2 = iw[iwpos + i - 1], abs(i__2)] = w[i];
/* L180: */
		    }
		}
	    }
	} else {
	    j1 = iwpos;
	    j2 = iwpos + ncols - 1;
	    k = apos;
	    if (kind == 1) {
		jpiv = -1;
		for (ipiv = nrows; ipiv >= 1; --ipiv) {
		    irhs = iw[iwpos + ipiv - 1];
		    lrow = nrows + 1 - ipiv;
		    if (irhs > 0) {
			apos -= lrow;
			rhs[irhs] *= a[apos];
		    } else {
			if (jpiv == -1) {
			    irhs1 = -iw[iwpos + ipiv - 2];
			    irhs2 = -irhs;
			    apos = apos - lrow - lrow - 1;
			    w1 = rhs[irhs1] * a[apos] + rhs[irhs2] * a[apos2];
			    rhs[irhs2] = rhs[irhs1] * a[apos2] + rhs[irhs2] * 
				    a[apos + lrow + 1];
			    rhs[irhs1] = w1;
			    --apos2;
			}
			jpiv = -jpiv;
		    }
/* L210: */
		}
		apos += nrows * (nrows + 1) / 2;
		k = apos;
		j1 = iwpos + nrows;
		i__1 = nrows;
		for (ipiv = 1; ipiv <= i__1; ++ipiv) {
		    irhs = (i__2 = iw[iwpos + ipiv - 1], abs(i__2));
		    w1 = rhs[irhs];
		    i__2 = j2;
		    for (j = j1; j <= i__2; ++j) {
			w1 += a[k] * rhs[i__3 = iw[j], abs(i__3)];
			++k;
/* L215: */
		    }
		    rhs[irhs] = w1;
/* L220: */
		}
		j2 = iwpos + nrows - 1;
		i__1 = nrows;
		for (ipiv = 1; ipiv <= i__1; ++ipiv) {
		    irhs = (i__2 = iw[j1 - 1], abs(i__2));
		    apos -= ipiv;
		    w1 = rhs[irhs];
		    k = apos + 1;
		    i__2 = j2;
		    for (j = j1; j <= i__2; ++j) {
			w1 -= a[k] * rhs[i__3 = iw[j], abs(i__3)];
			++k;
/* L230: */
		    }
		    rhs[irhs] = w1;
		    --j1;
/* L260: */
		}
	    } else {
		aposd = apos;
		aposm = apos + n2 * n2;
		apose = aposm + n2 * (n2 - 1) / 2;
		if (kind == 3) {
		    apose = aposm;
		}
		aposf = apose + n2;
		for (ipiv = n2; ipiv >= 1; --ipiv) {
		    irhs1 = -iw[j1];
		    irhs2 = -iw[j1 + n2];
		    w1 = rhs[irhs1] * a[apose] + rhs[irhs2] * a[aposd];
		    rhs[irhs2] = rhs[irhs1] * a[aposd] + rhs[irhs2] * a[aposf]
			    ;
		    rhs[irhs1] = w1;
		    aposd = aposd + n2 + 1;
		    ++apose;
		    ++aposf;
		    ++j1;
/* L270: */
		}
		k = aposf;
		j1 = iwpos + nrows;
		i__1 = n2;
		for (ipiv = 1; ipiv <= i__1; ++ipiv) {
		    irhs1 = (i__2 = iw[iwpos + ipiv - 1], abs(i__2));
		    w1 = rhs[irhs1];
		    i__2 = j2;
		    for (j = j1 + num1; j <= i__2; ++j) {
			irhs = (i__3 = iw[j], abs(i__3));
			w1 += rhs[irhs] * a[k];
			++k;
/* L290: */
		    }
		    rhs[irhs1] = w1;
/* L295: */
		}
		i__1 = n2;
		for (ipiv = 1; ipiv <= i__1; ++ipiv) {
		    irhs2 = (i__2 = iw[iwpos + ipiv - 1 + n2], abs(i__2));
		    w1 = rhs[irhs2];
		    i__2 = j2 - num3;
		    for (j = j1; j <= i__2; ++j) {
			irhs = (i__3 = iw[j], abs(i__3));
			w1 += rhs[irhs] * a[k];
			++k;
/* L300: */
		    }
		    rhs[irhs2] = w1;
/* L310: */
		}
		j1 = iwpos + nrows - 1;
		j2 = j1;
		for (ipiv = n2 - 1; ipiv >= 1; --ipiv) {
		    irhs1 = -iw[j1 - 1];
		    k = apos + (ipiv - 1) * (n2 + 1);
		    i__1 = j2;
		    for (j = j1; j <= i__1; ++j) {
			++k;
			rhs[irhs1] -= rhs[i__2 = iw[j], abs(i__2)] * a[k];
/* L320: */
		    }
		    --j1;
/* L330: */
		}
		if (kind == 2) {
		    k = aposm;
		    i__1 = n2 - 1;
		    for (ipiv = 1; ipiv <= i__1; ++ipiv) {
			++j1;
			irhs2 = -iw[j1 - 1 - n2];
			i__2 = j2;
			for (j = j1; j <= i__2; ++j) {
			    rhs[irhs2] += rhs[i__3 = iw[j], abs(i__3)] * a[k];
			    ++k;
/* L340: */
			}
/* L350: */
		    }
		}
		j1 = j2 - n2;
		j2 = j1;
		for (ipiv = n2 - 1; ipiv >= 1; --ipiv) {
		    irhs2 = -iw[j1 - 1];
		    k = apos + (ipiv - 1) * (n2 + 1);
		    i__1 = j2;
		    for (j = j1; j <= i__1; ++j) {
			k += n2;
			rhs[irhs2] -= rhs[i__2 = iw[j], abs(i__2)] * a[k];
/* L360: */
		    }
		    --j1;
/* L370: */
		}
	    }
	}
/* L380: */
    }
} /* ma47rd_ */

/* Subroutine */ int ma47sd_(a, bottom, top, move, ncmp)
doublereal *a;
integer *bottom, *top, *move, *ncmp;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer size, k, jj, idummy;

    /* Parameter adjustments */
    --a;

    /* Function Body */
    ++(*ncmp);
    k = *top - *move - 1;
    i__1 = *top;
    for (idummy = 1; idummy <= i__1; ++idummy) {
	if (k <= *bottom) {
	    goto L30;
	}
	size = (integer) a[k];
	if (size < 0) {
	    *move -= size;
	    k += size;
	} else {
	    if (*move > 0) {
		i__2 = k - size + 1;
		for (jj = k; jj >= i__2; --jj) {
		    a[jj + *move] = a[jj];
/* L10: */
		}
	    }
	    k -= size;
	}
/* L20: */
    }
L30:
    *bottom += *move;
    *move = 0;
    return 0;
} /* ma47sd_ */

/* Subroutine */ int ma47td_(n, ipe, iw, lw, nv, next, last, leaf, flag_, var,
	 svar)
integer *n, *ipe, *iw, *lw, *nv, *next, *last, *leaf, *flag_, *var, *svar;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_sign();

    /* Local variables */
    static integer free, i, j, k, kk, is, js, ls, ns;

    /* Parameter adjustments */
    --svar;
    --var;
    --flag_;
    --leaf;
    --last;
    --next;
    --nv;
    --iw;
    --ipe;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	svar[i] = 1;
	last[i] = i - 1;
	next[i] = i + 1;
	flag_[i] = 0;
	var[i] = i + 1;
/* L10: */
    }
    last[1] = *n;
    next[*n] = 1;
    free = 2;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	kk = ipe[j];
	i__2 = kk + iw[kk];
	for (k = kk + 1; k <= i__2; ++k) {
	    i = iw[k];
	    is = svar[i];
	    if (flag_[is] != j) {
		flag_[is] = j;
		if (next[i] == i) {
		} else {
		    js = free;
		    free = var[js];
		    var[is] = i;
		    svar[i] = js;
		    ns = next[i];
		    ls = last[i];
		    next[ls] = ns;
		    last[ns] = ls;
		    next[i] = i;
		    last[i] = i;
		}
	    } else {
		if (next[i] == i) {
		    ns = var[is];
		    var[is] = free;
		    free = is;
		} else {
		    ns = next[i];
		    ls = last[i];
		    next[ls] = ns;
		    last[ns] = ls;
		    ns = var[is];
		}
		ls = last[ns];
		next[ls] = i;
		next[i] = ns;
		last[ns] = i;
		last[i] = ls;
		svar[i] = svar[ns];
	    }
/* L20: */
	}
/* L30: */
    }
    i__1 = *n;
    for (is = 1; is <= i__1; ++is) {
	leaf[is] = is;
	if (last[is] == is) {
	    goto L60;
	}
	if (last[is] < 0) {
	    goto L60;
	}
	ls = last[is];
	leaf[is] = ls;
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    i = ls;
	    if (i == is) {
		goto L50;
	    }
	    ipe[i] = 0;
	    ls = last[i];
	    next[i] = -ls;
	    last[i] = -ls;
	    nv[i] = 0;
/* L40: */
	}
L50:
	nv[is] = i_sign(&k, &nv[is]);
L60:
	;
    }
} /* ma47td_ */

/* Subroutine */ int ma47ud_(a, la, iw, liw, icntl)
doublereal *a;
integer *la, *iw, *liw, *icntl;
{

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer iblk, kind, nblk, apos, jpiv, irow, apos2, j, ldiag, iline,
	     ncols, aposm, j1, j2, n2, iwpos, nrows, mp, len, num1, num3;

    /* Parameter adjustments */
    --icntl;
    --iw;
    --a;

    /* Function Body */
    mp = icntl[2];
    ldiag = icntl[3];
    apos2 = iw[1];
    nblk = iw[3];
    if (ldiag <= 2) {
	nblk = 0;
    }
    if (ldiag == 3) {
	nblk = min(1,nblk);
    }
    len = 12;
    if (ldiag == 5) {
	len = 1;
    }
    if (len == 12) {
	if (nblk == iw[3]) {
/*          WRITE (MP,'(A)') */
/*     +      ' For each block, the following information is provi
ded:' */
	} else {
/*          WRITE (MP,'(A,A)') ' For the first block only,', */
/*     +      ' the following information is provided:' */
	}
    }
/*      IF (LEN.EQ.12) WRITE (MP,'(A)') */
/*     +    '   1. Block number, number of rows, number of columns', */
/*     +    '   2. Whether the block is full, tile, or oxo', */
/*     +    '   3. List of indices for the pivot, each negated if part of'
 */
/*     +    ,'      a 2x2 pivot', */
/*     +    '   4. The factorized block pivot (see below)', */
/*     +    '   5. List of indices for the non-pivot columns', */
/*     +    '   6. The non-pivot part, with structural zeros printed as', 
*/
/*     +    '      0000000000',' A factorized full pivot has the form', */
/*     +    '               T','        L  D  L ', */
/*     +    '                          T', */
/*     +    ' and is printed as D and L  packed together.' */
/*      IF (LEN.EQ.12) WRITE (MP,'(A)') */
/*     +    ' A factorized tile or oxo pivot has the form', */
/*     +    '                              T', */
/*     +    '       ( L     ) ( F   D ) ( L   M )', */
/*     +    '       (  T  T ) (       ) (       )', */
/*     +    '       ( M  U  ) ( D   E ) (     U )', */
/*     +    ' where L is unit lower triangular; D, E, and F are diagonal;'
 */
/*     +    , */
/*     +   ' M is zero for an oxo pivot and is upper triangular with zero'
 */
/*     +    ,' diagonal for a tile pivot; and U is upper triangular with',
 */
/*     +    ' unit diagonal. It is printed as', */
/*     +    '    a. L,D,U  packed together,', */
/*     +    '    b. M (tile case only) in packed form,','    c. F and E.' 
*/
    iwpos = 4;
    apos = 1;
    i__1 = nblk;
    for (iblk = 1; iblk <= i__1; ++iblk) {
	ncols = iw[iwpos];
	nrows = iw[iwpos + 1];
	iwpos += 2;
	if (ncols > 0) {
	    kind = 1;
	} else {
	    ncols = -ncols;
	    num3 = iw[iwpos];
	    ++iwpos;
	    if (nrows > 0) {
		kind = 2;
	    } else {
		nrows = -nrows;
		num1 = iw[iwpos];
		++iwpos;
		kind = 3;
	    }
	}
/*        WRITE (MP,'(4(A,I6))') ' Block',IBLK,' with',NROWS,' rows an
d', */
/*     +    NCOLS,' columns' */
/*        IF (KIND.EQ.1) WRITE (MP,'(A)') ' Full pivot' */
/*        IF (KIND.EQ.2) WRITE (MP,'(A)') ' Tile pivot' */
/*        IF (KIND.EQ.3) WRITE (MP,'(A)') ' Oxo pivot' */
/*        IF (LEN.EQ.12) WRITE (MP,'(6I12)') (IW(K),K=IWPOS,IWPOS+NROW
S-1) */
/*        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (SGN(IW(K)),K=IWPOS, */
/*     +      IWPOS+NROWS-1) */
	if (kind == 1) {
	    jpiv = 0;
	    i__2 = nrows;
	    for (irow = 1; irow <= i__2; ++irow) {
		if (jpiv == 1) {
		    jpiv = 0;
		} else {
		    if (iw[iwpos + irow - 1] < 0) {
			jpiv = 1;
		    }
		}
		iline = 1;
		i__3 = irow - 1;
		for (j = 1; j <= i__3; ++j) {
/*              WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' ' 
*/
		    iline += len;
		    if (iline > 72) {
/*                WRITE (MP,'(A)') LINE */
			iline = 1;
		    }
/* L10: */
		}
		i__3 = nrows;
		for (j = irow; j <= i__3; ++j) {
/*              IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11)
, */
/*     +            '(1P,D12.4)') A(APOS) */
/*              IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A
)') NZ(A(APOS)) */
		    ++apos;
		    if (j == irow + 1) {
			if (jpiv == 1) {
/*                  IF (LEN.EQ.12) WRITE (LINE(ILI
NE:ILINE+11), */
/*     +                '(1P,D12.4)') A(APOS2) */
/*                  IF (LEN.EQ.1) WRITE (LINE(ILIN
E:ILINE), */
/*     +                '(A)') NZ(A(APOS2)) */
			    ++apos2;
			}
		    }
		    iline += len;
		    if (iline > 72) {
/*                WRITE (MP,'(A)') LINE */
			iline = 1;
		    }
/* L20: */
		}
		if (iline > 1) {
/*              LINE(ILINE:) = ' ' */
/*              WRITE (MP,'(A)') LINE */
		}
/* L30: */
	    }
	} else {
	    n2 = nrows / 2;
	    i__2 = n2;
	    for (irow = 1; irow <= i__2; ++irow) {
		iline = 1;
		i__3 = n2;
		for (j = 1; j <= i__3; ++j) {
/*              IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11)
, */
/*     +            '(1P,D12.4)') A(APOS) */
/*              IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A
)') NZ(A(APOS)) */
		    ++apos;
		    iline += len;
		    if (iline > 72) {
/*                WRITE (MP,'(A)') LINE */
			iline = 1;
		    }
/* L40: */
		}
		if (iline > 1) {
/*              LINE(ILINE:) = ' ' */
/*              WRITE (MP,'(A)') LINE */
		}
/* L50: */
	    }
	    if (kind == 2) {
		i__2 = n2;
		for (irow = 2; irow <= i__2; ++irow) {
		    iline = 1;
		    aposm = apos + irow - 2;
		    i__3 = irow - 1;
		    for (j = 1; j <= i__3; ++j) {
/*                IF (LEN.EQ.12) WRITE (LINE(ILINE:ILI
NE+11), */
/*     +              '(1P,D12.4)') A(APOSM) */
/*                IF (LEN.EQ.1) WRITE (LINE(ILINE:ILIN
E), */
/*     +              '(A)') NZ(A(APOSM)) */
			aposm = aposm + n2 - j - 1;
			iline += len;
			if (iline > 72) {
/*                  WRITE (MP,'(A)') LINE */
			    iline = 1;
			}
/* L60: */
		    }
		    i__3 = n2;
		    for (j = irow + 1; j <= i__3; ++j) {
/*                WRITE (LINE(ILINE:ILINE+LEN-1),'(A)'
) ' ' */
			iline += len;
			if (iline > 72) {
/*                  WRITE (MP,'(A)') LINE */
			    iline = 1;
			}
/* L70: */
		    }
		    if (iline > 1) {
/*                LINE(ILINE:) = ' ' */
/*                WRITE (MP,'(A)') LINE */
		    }
/* L80: */
		}
		apos += n2 * (n2 - 1) / 2;
	    }
	    iline = 1;
	    apos += n2;
	    i__2 = nrows;
	    for (j = 1; j <= i__2; ++j) {
/*            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11), */
/*     +          '(1P,D12.4)') A(APOS) */
/*            IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A)') NZ
(A(APOS)) */
		++apos;
		iline += len;
		if (iline > 72) {
/*              WRITE (MP,'(A)') LINE */
		    iline = 1;
		}
		if (j == n2) {
		    apos -= nrows;
		}
/* L90: */
	    }
	    apos += n2;
	    if (iline > 1) {
/*            LINE(ILINE:) = ' ' */
/*            WRITE (MP,'(A)') LINE */
	    }
	}
	iwpos += nrows;
/*        IF (LEN.EQ.12) WRITE (MP,'(6I12)') (IW(K),K=IWPOS, */
/*     +      IWPOS+NCOLS-NROWS-1) */
/*        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (SGN(IW(K)),K=IWPOS, */
/*     +      IWPOS+NCOLS-NROWS-1) */
	iwpos = iwpos + ncols - nrows;
	i__2 = nrows;
	for (irow = 1; irow <= i__2; ++irow) {
	    j1 = nrows;
	    j2 = ncols;
	    if (kind > 1) {
		if (irow > n2) {
		    j2 = ncols - num3;
		} else if (kind == 3) {
		    j1 = nrows + num1;
		}
	    }
	    iline = 1;
	    i__3 = j1;
	    for (j = nrows + 1; j <= i__3; ++j) {
/*            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11), */
/*     +          '(A)') '  0000000000' */
/*            IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A)') '0
' */
		iline += len;
		if (iline > 72) {
/*              WRITE (MP,'(A)') LINE */
		    iline = 1;
		}
/* L100: */
	    }
	    i__3 = j2;
	    for (j = j1 + 1; j <= i__3; ++j) {
/*            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11), */
/*     +          '(1P,D12.4)') A(APOS) */
/*            IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A)') NZ
(A(APOS)) */
		++apos;
		iline += len;
		if (iline > 72) {
/*              WRITE (MP,'(A)') LINE */
		    iline = 1;
		}
/* L110: */
	    }
	    i__3 = ncols;
	    for (j = j2 + 1; j <= i__3; ++j) {
/*            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11), */
/*     +          '(A)') '  0000000000' */
/*            IF (LEN.EQ.1) WRITE (LINE(ILINE:ILINE),'(A)') '0
' */
		iline += len;
		if (iline > 72) {
/*              WRITE (MP,'(A)') LINE */
		    iline = 1;
		}
/* L120: */
	    }
	    if (iline > 1) {
/*            LINE(ILINE:) = ' ' */
/*            WRITE (MP,'(A)') LINE */
	    }
/* L280: */
	}
/* L300: */
    }
} /* ma47ud_ */

/* Subroutine */ int ma47vd_(thresh, newthr, n, ipe, iw, lw, count, nv, next, 
	last, ipr, flag_, nflg)
integer *thresh, *newthr, *n, *ipe, *iw, *lw, *count, *nv, *next, *last, *ipr,
	 *flag_, *nflg;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer part, k;
    extern /* Subroutine */ int ma47zd_();
    static integer ke, jp, ir, is, kp, ls, ms, ns, jp1, jp2, kp2;

    /* Parameter adjustments */
    --flag_;
    --ipr;
    --last;
    --next;
    --nv;
    --count;
    --iw;
    --ipe;

    /* Function Body */
    i__1 = *n;
    for (ms = 1; ms <= i__1; ++ms) {
	if (flag_[ms] < 0) {
	    goto L50;
	}
	if (count[ms] <= *thresh) {
	    goto L50;
	}
	if (count[ms] > *newthr) {
	    goto L50;
	}
	ir = 0;
	if (*nflg <= 4) {
	    ma47zd_(n, &flag_[1], nflg);
	}
	--(*nflg);
	k = ipe[ms];
	kp2 = k + iw[k];
	i__2 = kp2;
	for (kp = k + 1; kp <= i__2; ++kp) {
	    part = (iw[kp] - 1) / *n;
	    ke = iw[kp] - part * *n;
	    if (flag_[ke] == -1) {
		goto L30;
	    }
	    if (flag_[ke] <= -2) {
		jp = ipe[ke];
		jp1 = jp + 3;
		jp2 = jp + iw[jp];
		if (part == 0) {
		} else if (part == 2) {
		    jp1 += iw[jp + 1];
		} else {
		    jp2 -= iw[jp + 2];
		}
	    } else {
		jp1 = kp;
		jp2 = kp2;
	    }
	    i__3 = jp2;
	    for (jp = jp1; jp <= i__3; ++jp) {
		is = iw[jp];
		if (flag_[is] > *nflg) {
		    flag_[is] = *nflg;
		    ir += (i__4 = nv[is], abs(i__4));
		}
/* L20: */
	    }
	    if (jp2 == kp2 || ir > *newthr) {
		goto L40;
	    }
L30:
	    ;
	}
L40:
	if (ir != count[ms]) {
	    ns = next[ms];
	    ls = last[ms];
	    next[ms] = 0;
	    if (ns > 0) {
		last[ns] = ls;
	    }
	    if (ls > 0) {
		next[ls] = ns;
	    } else {
		ipr[count[ms]] = ns;
	    }
	    ns = ipr[ir];
	    if (ns > 0) {
		last[ns] = ms;
	    }
	    next[ms] = ns;
	    ipr[ir] = ms;
	    last[ms] = 0;
	    count[ms] = ir;
	}
L50:
	;
    }
    *thresh = *newthr;
} /* ma47vd_ */

/* Subroutine */ int ma47wd_(a, la, iw, liw, nrlbdu)
doublereal *a;
integer *la, *iw, *liw, *nrlbdu;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer iblk, apos, jpiv, irow, j, ncols, n2, iwpos, nrows, num1, 
	    num3;

    /* Parameter adjustments */
    --iw;
    --a;

    /* Function Body */
    apos = 1;
    iwpos = 6;
    i__1 = iw[3];
    for (iblk = 1; iblk <= i__1; ++iblk) {
	ncols = iw[iwpos - 2];
	nrows = iw[iwpos - 1];
	n2 = abs(nrows) / 2;
	if (ncols > 0) {
	    jpiv = 1;
	    i__2 = nrows;
	    for (irow = 1; irow <= i__2; ++irow) {
		--jpiv;
		if (jpiv == 1) {
		    goto L10;
		}
		if (iw[iwpos + irow - 1] < 0) {
		    jpiv = 2;
		    ++(*nrlbdu);
		    a[*nrlbdu] = a[apos + 1];
		    a[apos + 1] = 0.;
		}
L10:
		i__3 = apos + nrows - irow;
		for (j = apos + 1; j <= i__3; ++j) {
		    a[j] = -a[j];
/* L20: */
		}
		apos = apos + nrows - irow + 1;
/* L30: */
	    }
	    apos += nrows * (ncols - nrows);
	} else {
	    ncols = -ncols;
	    num3 = iw[iwpos];
	    apos = apos - num3 * n2 - n2 * (n2 + 1) / 2 + n2;
	    ++iwpos;
	    if (nrows < 0) {
		nrows = -nrows;
		num1 = iw[iwpos];
		apos = apos - num1 * n2 - n2 * (n2 + 1) / 2 + n2;
		++iwpos;
	    }
	    apos += nrows * (ncols + ncols - nrows + 1) / 2;
	}
	iwpos = iwpos + abs(ncols) + 2;
/* L40: */
    }
} /* ma47wd_ */

/* Subroutine */ int ma47xd_(a, lda, m, n)
doublereal *a;
integer *lda, *m, *n;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] = 0.;
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* ma47xd_ */

/* Subroutine */ int ma47yd_(a, lda, m, n, b, ldb, triang, trans)
doublereal *a;
integer *lda, *m, *n;
doublereal *b;
integer *ldb;
logical *triang, *trans;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;

    /* Parameter adjustments */
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (! (*triang)) {
	if (! (*trans)) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i = 1; i <= i__2; ++i) {
		    b[i + j * b_dim1] = a[i + j * a_dim1];
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i = 1; i <= i__2; ++i) {
		    b[j + i * b_dim1] = a[i + j * a_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
    } else {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i = j; i <= i__2; ++i) {
		b[i + j * b_dim1] = a[i + j * a_dim1];
/* L50: */
	    }
/* L60: */
	}
    }
    return 0;
} /* ma47yd_ */

/* Subroutine */ int ma47zd_(n, flag_, nflg)
integer *n, *flag_, *nflg;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, n3;

    /* Parameter adjustments */
    --flag_;

    /* Function Body */
    n3 = *n * 3;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (flag_[i] > 2) {
	    flag_[i] = n3;
	}
	if (flag_[i] <= -2) {
	    flag_[i] = -n3;
	}
/* L80: */
    }
    *nflg = n3;
} /* ma47zd_ */

