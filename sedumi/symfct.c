/* symfct.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "mex.h"
typedef mwSignedIndex integer;                  /* removed "long" */

#if !defined(max)
#define  max(A, B)   ((A) > (B) ? (A) : (B))
#endif

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP */
/*       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION. */

/*   NOTE: */
/*       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E., */
/*       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES). */

/*   INPUT PARAMETERS: */
/*       NEQNS       -   NUMBER OF EQUATIONS. */
/*       NNZA        -   LENGTH OF ADJACENCY STRUCTURE. */
/*       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                       TO THE ADJACENCY STRUCTURE. */
/*       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                       THE ADJACENCY STRUCTURE. */
/*       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   UPDATED PARAMETERS: */
/* (* JFS Sept 1, 1998: moved "PERM" and "INVP" to updated section *) */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE EQUIVALENT */
/*                           ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                       OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                       INCLUDING THE DIAGONAL ENTRY. */
/*       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                       THE DIAGONAL ENTRIES. */
/*       NSUB        -   NUMBER OF SUBSCRIPTS. */
/*       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                       SUPERNODE MEMBERSHIP. */
/*       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                       SUPERNODE PARTITIONING. */
/*       IFLAG(*)    -   ERROR FLAG. */
/*                          0: SUCCESSFUL SF INITIALIZATION. */
/*                         -1: INSUFFICENT WORKING STORAGE */
/*                             [IWORK(*)]. */

/*   WORK PARAMETERS: */
/*       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3. */

/*   FIRST CREATED ON    NOVEMBER 14, 1994. */
/*   LAST UPDATED ON     January 12, 1995. */

/* *********************************************************************** */

/* Subroutine */ int sfinit_(neqns, nnza, xadj, adjncy, perm, invp, colcnt, 
	nnzl, nsub, nsuper, snode, xsuper, iwsiz, iwork, iflag)
integer *neqns, *nnza, *xadj, *adjncy, *perm, *invp, *colcnt, *nnzl, *nsub, *
	nsuper, *snode, *xsuper, *iwsiz, *iwork, *iflag;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int fsup1_(), fsup2_();
    static integer i__;
    extern /* Subroutine */ int fcnthn_(), chordr_(), etordr_();


/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* ***********************************************************************
 */

/*       -------------------------------------------------------- */
/*       RETURN IF THERE IS INSUFFICIENT INTEGER WORKING STORAGE. */
/*       -------------------------------------------------------- */
    /* Parameter adjustments */
    --iwork;
    --xsuper;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --xadj;
    --adjncy;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < *neqns * 7 + 3) {
	*iflag = -1;
	return 0;
    }
/* (* JFS: */
/* ------------------------------------------------------------ */
/* Handle the case of diagonal matrices separately, to avoid */
/* an unsolved BUG in FCNTHN. */
/* This patch is due to Jos F. Sturm, Sept 1, 1998. */
/* ------------------------------------------------------------ */
    if (xadj[*neqns + 1] - 1 == 0) {
	i__1 = *neqns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    colcnt[i__] = 1;
	    snode[i__] = i__;
	    xsuper[i__] = i__;
/* L10: */
	}
	xsuper[*neqns + 1] = *neqns + 1;
	*nnzl = *neqns;
	*nsub = *neqns;
	*nsuper = *neqns;
	*iflag = 0;
	return 0;
    }
/* ------------------------------------------------------------ */
/* end of patch JFS *) */

/*       ------------------------------------------ */
/*       COMPUTE ELIMINATION TREE AND POSTORDERING. */
/*       ------------------------------------------ */
    etordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], &
	    iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*neqns * 3 + 
	    1]);

/*       --------------------------------------------- */
/*       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS. */
/*       --------------------------------------------- */
    fcnthn_(neqns, nnza, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], 
	    &snode[1], &colcnt[1], nnzl, &iwork[*neqns + 1], &iwork[(*neqns <<
	     1) + 1], &xsuper[1], &iwork[*neqns * 3 + 1], &iwork[(*neqns << 2)
	     + 2], &iwork[*neqns * 5 + 3], &iwork[*neqns * 6 + 4]);

/*       --------------------------------------------------------- */
/*       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM */
/*       NUMBER OF NONZEROS IN ITS COLUMN OF L. */
/*       --------------------------------------------------------- */
    chordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[1], &
	    iwork[1], &iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*
	    neqns * 3 + 1]);

/*       ---------------- */
/*       FIND SUPERNODES. */
/*       ---------------- */
    fsup1_(neqns, &iwork[1], &colcnt[1], nsub, nsuper, &snode[1]);
    fsup2_(neqns, nsuper, &iwork[1], &snode[1], &xsuper[1]);

    return 0;
} /* sfinit_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **********     ETORDR ..... ELIMINATION TREE REORDERING     *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE AN EQUIVALENT REORDERING BASED ON THE STRUCTURE OF */
/*       THE ELIMINATION TREE.  A POSTORDERING OF THE GIVEN ELIMINATION */
/*       TREE IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE EQUIVALENT */
/*                           ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BETREE, ETPOST, ETREE , INVINV. */

/* *********************************************************************** */

/* Subroutine */ int etordr_(neqns, xadj, adjncy, perm, invp, parent, fson, 
	brothr, invpos)
integer *neqns, *xadj, *adjncy, *perm, *invp, *parent, *fson, *brothr, *
	invpos;
{
    extern /* Subroutine */ int etree_(), betree_(), invinv_(), etpost_();


/* ***********************************************************************
 */



/* ***********************************************************************
 */

/*       ----------------------------- */
/*       COMPUTE THE ELIMINATION TREE. */
/*       ----------------------------- */
    /* Parameter adjustments */
    --invpos;
    --brothr;
    --fson;
    --parent;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    etree_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &parent[1], &
	    invpos[1]);

/*       -------------------------------------------------------- */
/*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE. */
/*       -------------------------------------------------------- */
    betree_(neqns, &parent[1], &fson[1], &brothr[1]);

/*       ------------------------------- */
/*       POSTORDER THE ELIMINATION TREE. */
/*       ------------------------------- */
    etpost_(neqns, &fson[1], &brothr[1], &invpos[1], &parent[1], &perm[1]);

/*       -------------------------------------------------------- */
/*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
/*       -------------------------------------------------------- */
    invinv_(neqns, &invp[1], &invpos[1], &perm[1]);

    return 0;
} /* etordr_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************     ETREE ..... ELIMINATION TREE     ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND */
/*       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */
/*       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */

/*   WORKING PARAMETERS: */
/*       ANCSTR          -   THE ANCESTOR VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int etree_(neqns, xadj, adjncy, perm, invp, parent, ancstr)
integer *neqns, *xadj, *adjncy, *perm, *invp, *parent, *ancstr;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer node, next, i__, j, jstop, jstrt, nbr;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --ancstr;
    --parent;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	parent[i__] = 0;
	ancstr[i__] = 0;
	node = perm[i__];

	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	if (jstrt <= jstop) {
	    i__2 = jstop;
	    for (j = jstrt; j <= i__2; ++j) {
		nbr = adjncy[j];
		nbr = invp[nbr];
		if (nbr < i__) {
/*                       --------------------------------
----------- */
/*                       FOR EACH NBR, FIND THE ROOT OF IT
S CURRENT */
/*                       ELIMINATION TREE.  PERFORM PATH C
OMPRESSION */
/*                       AS THE SUBTREE IS TRAVERSED. */
/*                       --------------------------------
----------- */
L100:
		    if (ancstr[nbr] == i__) {
			goto L300;
		    }
		    if (ancstr[nbr] > 0) {
			next = ancstr[nbr];
			ancstr[nbr] = i__;
			nbr = next;
			goto L100;
		    }
/*                       --------------------------------
------------ */
/*                       NOW, NBR IS THE ROOT OF THE SUBTR
EE.  MAKE I */
/*                       THE PARENT NODE OF THIS ROOT. */
/*                       --------------------------------
------------ */
		    parent[nbr] = i__;
		    ancstr[nbr] = i__;
		}
L300:
		;
	    }
	}
/* L400: */
    }

    return 0;
} /* etree_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BETREE ..... BINARY TREE REPRESENTATION OF ETREE     ******* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE BINARY TREE REPRESENTATION OF THE ELIMINATION */
/*       TREE GIVEN BY THE PARENT VECTOR.  THE RETURNED REPRESENTATION */
/*       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT */
/*       OF THE BINARY TREE IS ALWAYS NEQNS. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */
/*                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF */
/*                           THE ROOTS. */

/*   OUTPUT PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int betree_(neqns, parent, fson, brothr)
integer *neqns, *parent, *fson, *brothr;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, ndpar, lroot;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --brothr;
    --fson;
    --parent;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (node = 1; node <= i__1; ++node) {
	fson[node] = 0;
	brothr[node] = 0;
/* L100: */
    }
    lroot = *neqns;
/*       ------------------------------------------------------------ */
/*       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING. */
/*       ------------------------------------------------------------ */
    if (*neqns <= 1) {
	return 0;
    }
    for (node = *neqns - 1; node >= 1; --node) {
	ndpar = parent[node];
	if (ndpar <= 0 || ndpar == node) {
/*               -------------------------------------------------
 */
/*               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
 */
/*               SET NODE TO BE ONE OF THE ROOTS OF THE TREES. */
/*               -------------------------------------------------
 */
	    brothr[lroot] = node;
	    lroot = node;
	} else {
/*               ------------------------------------------- */
/*               OTHERWISE, BECOMES FIRST SON OF ITS PARENT. */
/*               ------------------------------------------- */
	    brothr[node] = fson[ndpar];
	    fson[ndpar] = node;
	}
/* L300: */
    }
    brothr[lroot] = 0;

    return 0;
} /* betree_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***************     ETPOST ..... ETREE POSTORDERING     *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (SEPT 17, 1986) */

/*   PURPOSE: */
/*       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF */
/*       THE ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE */
/*       CORRESPONDING PARENT VECTOR IS ALSO MODIFIED TO REFLECT */
/*       THE REORDERING. */

/*   INPUT PARAMETERS: */
/*       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT */
/*                           IS NEQNS). */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHR VECTOR. */

/*   UPDATED PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR. */

/*   OUTPUT PARAMETERS: */
/*       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING. */

/*   WORKING PARAMETERS: */
/*       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE */
/*                           TREE. */

/* *********************************************************************** */

/* Subroutine */ int etpost_(root, fson, brothr, invpos, parent, stack)
integer *root, *fson, *brothr, *invpos, *parent, *stack;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, itop, ndpar, nunode, num;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --stack;
    --parent;
    --invpos;
    --brothr;
    --fson;

    /* Function Body */
    num = 0;
    itop = 0;
    node = *root;
/*       ------------------------------------------------------------- */
/*       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES */
/*       ALONG THE TRAVERSAL INTO THE STACK. */
/*       ------------------------------------------------------------- */
L100:
    ++itop;
    stack[itop] = node;
    node = fson[node];
    if (node > 0) {
	goto L100;
    }
/*           ---------------------------------------------------------- */
/*           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT. */
/*           ---------------------------------------------------------- */
L200:
    if (itop <= 0) {
	goto L300;
    }
    node = stack[itop];
    --itop;
    ++num;
    invpos[node] = num;
/*               ---------------------------------------------------- */
/*               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE. */
/*               ---------------------------------------------------- */
    node = brothr[node];
    if (node <= 0) {
	goto L200;
    }
    goto L100;

L300:
/*       ------------------------------------------------------------ */
/*       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR */
/*       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR. */
/*       ------------------------------------------------------------ */
    i__1 = num;
    for (node = 1; node <= i__1; ++node) {
	nunode = invpos[node];
	ndpar = parent[node];
	if (ndpar > 0) {
	    ndpar = invpos[ndpar];
	}
	brothr[nunode] = ndpar;
/* L400: */
    }

    i__1 = num;
    for (nunode = 1; nunode <= i__1; ++nunode) {
	parent[nunode] = brothr[nunode];
/* L500: */
    }

    return 0;
} /* etpost_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***********     INVINV ..... CONCATENATION OF TWO INVP     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO PERFORM THE MAPPING OF */
/*           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP */
/*       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION */
/*       VECTOR PERM IS ALSO COMPUTED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR. */

/*   UPDATED PARAMETERS: */
/*       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON */
/*                           OUTPUT, IT CONTAINS THE NEW INVERSE */
/*                           PERMUTATION. */

/*   OUTPUT PARAMETER: */
/*       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS */
/*                           INVP2). */

/* *********************************************************************** */

/* Subroutine */ int invinv_(neqns, invp, invp2, perm)
integer *neqns, *invp, *invp2, *perm;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, i__, interm;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --perm;
    --invp2;
    --invp;

    /* Function Body */
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	interm = invp[i__];
	invp[i__] = invp2[interm];
/* L100: */
    }

    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = invp[i__];
	perm[node] = i__;
/* L200: */
    }

    return 0;
} /* invinv_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **********     CHORDR ..... CHILD REORDERING                *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       REARRANGE THE CHILDREN OF EACH VERTEX SO THAT THE LAST ONE */
/*       MAXIMIZES (AMONG THE CHILDREN) THE NUMBER OF NONZEROS IN THE */
/*       CORRESPONDING COLUMN OF L.  ALSO DETERMINE AN NEW POSTORDERING */
/*       BASED ON THE STRUCTURE OF THE MODIFIED ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE NEW */
/*                           POSTORDERING. */
/*       COLCNT          -   COLUMN COUNTS IN L UNDER INITIAL ORDERING; */
/*                           MODIFIED TO REFLECT THE NEW ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BTREE2, EPOST2, INVINV. */

/* *********************************************************************** */

/* Subroutine */ int chordr_(neqns, xadj, adjncy, perm, invp, colcnt, parent, 
	fson, brothr, invpos)
integer *neqns, *xadj, *adjncy, *perm, *invp, *colcnt, *parent, *fson, *
	brothr, *invpos;
{
    extern /* Subroutine */ int btree2_(), epost2_(), invinv_();


/* ***********************************************************************
 */



/* ***********************************************************************
 */

/*       ---------------------------------------------------------- */
/*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE, */
/*       SO THAT EACH "LAST CHILD" MAXIMIZES AMONG ITS SIBLINGS THE */
/*       NUMBER OF NONZEROS IN THE CORRESPONDING COLUMNS OF L. */
/*       ---------------------------------------------------------- */
    /* Parameter adjustments */
    --invpos;
    --brothr;
    --fson;
    --parent;
    --colcnt;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    btree2_(neqns, &parent[1], &colcnt[1], &fson[1], &brothr[1], &invpos[1]);

/*       ---------------------------------------------------- */
/*       POSTORDER THE ELIMINATION TREE (USING THE NEW BINARY */
/*       REPRESENTATION. */
/*       ---------------------------------------------------- */
    epost2_(neqns, &fson[1], &brothr[1], &invpos[1], &parent[1], &colcnt[1], &
	    perm[1]);

/*       -------------------------------------------------------- */
/*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
/*       -------------------------------------------------------- */
    invinv_(neqns, &invp[1], &invpos[1], &perm[1]);

    return 0;
} /* chordr_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BTREE2 ..... BINARY TREE REPRESENTATION OF ETREE     ******* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       TO DETERMINE A BINARY TREE REPRESENTATION OF THE ELIMINATION */
/*       TREE, FOR WHICH EVERY "LAST CHILD" HAS THE MAXIMUM POSSIBLE */
/*       COLUMN NONZERO COUNT IN THE FACTOR.  THE RETURNED REPRESENTATION */
/*       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT OF */
/*       THE BINARY TREE IS ALWAYS NEQNS. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */
/*                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF */
/*                           THE ROOTS. */
/*       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR. */

/*   OUTPUT PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */

/*   WORKING PARAMETERS: */
/*       LSON            -   LAST SON VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int btree2_(neqns, parent, colcnt, fson, brothr, lson)
integer *neqns, *parent, *colcnt, *fson, *brothr, *lson;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, ndpar, lroot, ndlson;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --lson;
    --brothr;
    --fson;
    --colcnt;
    --parent;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (node = 1; node <= i__1; ++node) {
	fson[node] = 0;
	brothr[node] = 0;
	lson[node] = 0;
/* L100: */
    }
    lroot = *neqns;
/*       ------------------------------------------------------------ */
/*       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING. */
/*       ------------------------------------------------------------ */
    if (*neqns <= 1) {
	return 0;
    }
    for (node = *neqns - 1; node >= 1; --node) {
	ndpar = parent[node];
	if (ndpar <= 0 || ndpar == node) {
/*               -------------------------------------------------
 */
/*               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
 */
/*               SET NODE TO BE ONE OF THE ROOTS OF THE TREES. */
/*               -------------------------------------------------
 */
	    brothr[lroot] = node;
	    lroot = node;
	} else {
/*               ------------------------------------------- */
/*               OTHERWISE, BECOMES FIRST SON OF ITS PARENT. */
/*               ------------------------------------------- */
	    ndlson = lson[ndpar];
	    if (ndlson != 0) {
		if (colcnt[node] >= colcnt[ndlson]) {
		    brothr[node] = fson[ndpar];
		    fson[ndpar] = node;
		} else {
		    brothr[ndlson] = node;
		    lson[ndpar] = node;
		}
	    } else {
		fson[ndpar] = node;
		lson[ndpar] = node;
	    }
	}
/* L300: */
    }
    brothr[lroot] = 0;

    return 0;
} /* btree2_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***************     EPOST2 ..... ETREE POSTORDERING #2  *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF THE */
/*       ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE */
/*       CORRESPONDING PARENT AND COLCNT VECTORS ARE ALSO MODIFIED TO */
/*       REFLECT THE REORDERING. */

/*   INPUT PARAMETERS: */
/*       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT */
/*                           IS NEQNS). */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHR VECTOR. */

/*   UPDATED PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR. */
/*       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR. */

/*   OUTPUT PARAMETERS: */
/*       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING. */

/*   WORKING PARAMETERS: */
/*       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE */
/*                           TREE. */

/* *********************************************************************** */

/* Subroutine */ int epost2_(root, fson, brothr, invpos, parent, colcnt, 
	stack)
integer *root, *fson, *brothr, *invpos, *parent, *colcnt, *stack;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, itop, ndpar, nunode, num;


/* ***********************************************************************
 */



/* ***********************************************************************
 */


/* ***********************************************************************
 */

    /* Parameter adjustments */
    --stack;
    --colcnt;
    --parent;
    --invpos;
    --brothr;
    --fson;

    /* Function Body */
    num = 0;
    itop = 0;
    node = *root;
/*       ------------------------------------------------------------- */
/*       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES */
/*       ALONG THE TRAVERSAL INTO THE STACK. */
/*       ------------------------------------------------------------- */
L100:
    ++itop;
    stack[itop] = node;
    node = fson[node];
    if (node > 0) {
	goto L100;
    }
/*           ---------------------------------------------------------- */
/*           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT. */
/*           ---------------------------------------------------------- */
L200:
    if (itop <= 0) {
	goto L300;
    }
    node = stack[itop];
    --itop;
    ++num;
    invpos[node] = num;
/*               ---------------------------------------------------- */
/*               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE. */
/*               ---------------------------------------------------- */
    node = brothr[node];
    if (node <= 0) {
	goto L200;
    }
    goto L100;

L300:
/*       ------------------------------------------------------------ */
/*       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR */
/*       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR. */
/*       ------------------------------------------------------------ */
    i__1 = num;
    for (node = 1; node <= i__1; ++node) {
	nunode = invpos[node];
	ndpar = parent[node];
	if (ndpar > 0) {
	    ndpar = invpos[ndpar];
	}
	brothr[nunode] = ndpar;
/* L400: */
    }

    i__1 = num;
    for (nunode = 1; nunode <= i__1; ++nunode) {
	parent[nunode] = brothr[nunode];
/* L500: */
    }

/*       ---------------------------------------------- */
/*       PERMUTE COLCNT(*) TO REFLECT THE NEW ORDERING. */
/*       ---------------------------------------------- */
    i__1 = num;
    for (node = 1; node <= i__1; ++node) {
	nunode = invpos[node];
	stack[nunode] = colcnt[node];
/* L600: */
    }

    i__1 = num;
    for (node = 1; node <= i__1; ++node) {
	colcnt[node] = stack[node];
/* L700: */
    }

    return 0;
} /* epost2_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************     FCNTHN  ..... FIND NONZERO COUNTS    *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN */
/*       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM. */

/*       TECHNIQUES: */
/*       1) SUPERNODE DETECTION. */
/*       2) PATH HALVING. */
/*       3) NO UNION BY RANK. */

/*   NOTES: */
/*       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH ROW OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                           THE DIAGONAL ENTRIES. */

/*   WORK PARAMETERS: */
/*       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE */
/*                           DISJOINT SETS (I.E., SUBTREES). */
/*       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS LEAF OF EACH ROW SUBTREE. */
/*       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL */
/*                           (DISTANCE FROM THE ROOT). */
/*       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS */
/*                           USED TO COMPUTE COLUMN COUNTS. */
/*       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT. */
/*       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           NUMBER OF CHILDREN. */
/*       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE. */

/*   FIRST CREATED ON    APRIL 12, 1990. */
/*   LAST UPDATED ON     JANUARY 12, 1995. */

/* (*JFS Sept 1, 1998: there is a BUG in fcnthn: if adjlen = 0, i.e. */
/*    the matrix is purely diagonal, then "segment violation" *) */
/* *********************************************************************** */

/* Subroutine */ int fcnthn_(neqns, adjlen, xadj, adjncy, perm, invp, etpar, 
	rowcnt, colcnt, nlnz, set, prvlf, level, weight, fdesc, nchild, 
	prvnbr)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *etpar, *rowcnt, *
	colcnt, *nlnz, *set, *prvlf, *level, *weight, *fdesc, *nchild, *
	prvnbr;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer temp, xsup, last1, last2, j, k, lflag, pleaf, hinbr, jstop,
	     jstrt, ifdesc, oldnbr, parent, lownbr, lca;


/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

/*       -------------------------------------------------- */
/*       COMPUTE LEVEL(*), FDESC(*), NCHILD(*). */
/*       INITIALIZE ROWCNT(*), COLCNT(*), */
/*                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*). */
/*       -------------------------------------------------- */
    /* Parameter adjustments */
    --prvnbr;
    --prvlf;
    --set;
    --colcnt;
    --rowcnt;
    --etpar;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    level[0] = 0;
    for (k = *neqns; k >= 1; --k) {
	rowcnt[k] = 1;
	colcnt[k] = 0;
	set[k] = k;
	prvlf[k] = 0;
	level[k] = level[etpar[k]] + 1;
	weight[k] = 1;
	fdesc[k] = k;
	nchild[k] = 0;
	prvnbr[k] = 0;
/* L100: */
    }
    nchild[0] = 0;
    fdesc[0] = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	parent = etpar[k];
	weight[parent] = 0;
	++nchild[parent];
	ifdesc = fdesc[k];
	if (ifdesc < fdesc[parent]) {
	    fdesc[parent] = ifdesc;
	}
/* L200: */
    }
/*       ------------------------------------ */
/*       FOR EACH ``LOW NEIGHBOR'' LOWNBR ... */
/*       ------------------------------------ */
    i__1 = *neqns;
    for (lownbr = 1; lownbr <= i__1; ++lownbr) {
	lflag = 0;
	ifdesc = fdesc[lownbr];
	oldnbr = perm[lownbr];
	jstrt = xadj[oldnbr];
	jstop = xadj[oldnbr + 1] - 1;
/*           ----------------------------------------------- */
/*           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ... */
/*           ----------------------------------------------- */
	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    hinbr = invp[adjncy[j]];
	    if (hinbr > lownbr) {
		if (ifdesc > prvnbr[hinbr]) {
/*                       ------------------------- */
/*                       INCREMENT WEIGHT(LOWNBR). */
/*                       ------------------------- */
		    ++weight[lownbr];
		    pleaf = prvlf[hinbr];
/*                       --------------------------------
--------- */
/*                       IF HINBR HAS NO PREVIOUS ``LOW NE
IGHBOR'' */
/*                       THEN ... */
/*                       --------------------------------
--------- */
		    if (pleaf == 0) {
/*                           ------------------------
----------------- */
/*                           ... ACCUMULATE LOWNBR-->H
INBR PATH LENGTH */
/*                               IN ROWCNT(HINBR). */
/*                           ------------------------
----------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				hinbr];
		    } else {
/*                           ------------------------
----------------- */
/*                           ... OTHERWISE, LCA <-- FI
ND(PLEAF), WHICH */
/*                               IS THE LEAST COMMON A
NCESTOR OF PLEAF */
/*                               AND LOWNBR. */
/*                               (PATH HALVING.) */
/*                           ------------------------
----------------- */
			last1 = pleaf;
			last2 = set[last1];
			lca = set[last2];
L300:
			if (lca != last2) {
			    set[last1] = lca;
			    last1 = lca;
			    last2 = set[last1];
			    lca = set[last2];
			    goto L300;
			}
/*                           ------------------------
------------- */
/*                           ACCUMULATE PLEAF-->LCA PA
TH LENGTH IN */
/*                           ROWCNT(HINBR). */
/*                           DECREMENT WEIGHT(LCA). */
/*                           ------------------------
------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				lca];
			--weight[lca];
		    }
/*                       --------------------------------
-------------- */
/*                       LOWNBR NOW BECOMES ``PREVIOUS LEA
F'' OF HINBR. */
/*                       --------------------------------
-------------- */
		    prvlf[hinbr] = lownbr;
		    lflag = 1;
		}
/*                   ----------------------------------------
---------- */
/*                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' 
OF HINBR. */
/*                   ----------------------------------------
---------- */
		prvnbr[hinbr] = lownbr;
	    }
/* L500: */
	}
/*           ---------------------------------------------------- */
/*           DECREMENT WEIGHT ( PARENT(LOWNBR) ). */
/*           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP). */
/*           ---------------------------------------------------- */
	parent = etpar[lownbr];
	--weight[parent];
	if (lflag == 1 || nchild[lownbr] >= 2) {
	    xsup = lownbr;
	}
	set[xsup] = parent;
/* L600: */
    }
/*       --------------------------------------------------------- */
/*       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS. */
/*       --------------------------------------------------------- */
    *nlnz = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	temp = colcnt[k] + weight[k];
	colcnt[k] = temp;
	*nlnz += temp;
	parent = etpar[k];
	if (parent != 0) {
	    colcnt[parent] += temp;
	}
/* L700: */
    }

    return 0;
} /* fcnthn_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************    FSUP1 ..... FIND SUPERNODES #1    ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE FIRST OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT RETURNS ONLY THE NUMBER OF */
/*       SUPERNODES NSUPER AND THE SUPERNODE MEMBERSHIP VECTOR SNODE(*), */
/*       WHICH IS OF LENGTH NEQNS.  THE VECTORS OF LENGTH NSUPER ARE */
/*       COMPUTED SUBSEQUENTLY BY THE COMPANION ROUTINE FSUP2. */

/*   METHOD AND ASSUMPTIONS: */
/*       THIS ROUTINE USES THE ELIMINATION TREE AND THE FACTOR COLUMN */
/*       COUNTS TO COMPUTE THE SUPERNODE PARTITION; IT ALSO ASSUMES A */
/*       POSTORDERING OF THE ELIMINATION TREE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           FACTOR COLUMN COUNTS: I.E., THE NUMBER OF */
/*                           NONZERO ENTRIES IN EACH COLUMN OF L */
/*                           (INCLUDING THE DIAGONAL ENTRY). */

/*   OUTPUT PARAMETERS: */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMBER 11, 1994. */

/* *********************************************************************** */

/* Subroutine */ int fsup1_(neqns, etpar, colcnt, nofsub, nsuper, snode)
integer *neqns, *etpar, *colcnt, *nofsub, *nsuper, *snode;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kcol;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

/*       -------------------------------------------- */
/*       COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION. */
/*       -------------------------------------------- */
    /* Parameter adjustments */
    --snode;
    --colcnt;
    --etpar;

    /* Function Body */
    *nsuper = 1;
    snode[1] = 1;
    *nofsub = colcnt[1];
    i__1 = *neqns;
    for (kcol = 2; kcol <= i__1; ++kcol) {
	if (etpar[kcol - 1] == kcol) {
	    if (colcnt[kcol - 1] == colcnt[kcol] + 1) {
		snode[kcol] = *nsuper;
		goto L300;
	    }
	}
	++(*nsuper);
	snode[kcol] = *nsuper;
	*nofsub += colcnt[kcol];
L300:
	;
    }

    return 0;
} /* fsup1_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****************    FSUP2  ..... FIND SUPERNODES #2   ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO */
/*       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE */
/*       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE */
/*       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS. */


/*   ASSUMPTIONS: */
/*       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT */
/*       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   OUTPUT PARAMETERS: */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           SUPERNODE PARTITIONING. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMEBER 22, 1994. */

/* *********************************************************************** */

/* Subroutine */ int fsup2_(neqns, nsuper, etpar, snode, xsuper)
integer *neqns, *nsuper, *etpar, *snode, *xsuper;
{
    static integer kcol, ksup, lstsup;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

/*       ------------------------------------------------- */
/*       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*). */
/*       ------------------------------------------------- */
    /* Parameter adjustments */
    --xsuper;
    --snode;
    --etpar;

    /* Function Body */
    lstsup = *nsuper + 1;
    for (kcol = *neqns; kcol >= 1; --kcol) {
	ksup = snode[kcol];
	if (ksup != lstsup) {
	    xsuper[lstsup] = kcol + 1;
	}
	lstsup = ksup;
/* L100: */
    }
    xsuper[1] = 1;

    return 0;
} /* fsup2_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFCT ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE CALLS SYMFC2 WHICH PERFORMS SUPERNODAL SYMBOLIC */
/*       FACTORIZATION ON A REORDERED LINEAR SYSTEM. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */
/*       (I) IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                              -1 - INSUFFICIENT INTEGER WORKING SPACE. */
/*                              -2 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) IWORK       -   WORKING ARRAY OF LENGTH NSUPER+2*NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int symfct_(neqns, adjlen, xadj, adjncy, perm, invp, colcnt, 
	nsuper, xsuper, snode, nofsub, xlindx, lindx, xlnz, iwsiz, iwork, 
	flag__)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *colcnt, *nsuper, *
	xsuper, *snode, *nofsub, *xlindx, *lindx, *xlnz, *iwsiz, *iwork, *
	flag__;
{
    extern /* Subroutine */ int symfc2_();


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    --xlnz;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --xadj;
    --adjncy;
    --iwork;
    --xlindx;
    --xsuper;
    --lindx;

    /* Function Body */
    *flag__ = 0;
    if (*iwsiz < *nsuper + (*neqns << 1) + 1) {
	*flag__ = -1;
	return 0;
    }
    symfc2_(neqns, adjlen, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[
	    1], nsuper, &xsuper[1], &snode[1], nofsub, &xlindx[1], &lindx[1], 
	    &xlnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[*nsuper + *neqns 
	    + 2], flag__);
    return 0;
} /* symfct_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFC2 ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A */
/*       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS */
/*       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE */
/*       ASSOCIATED WITH THE FACTOR MATRIX L. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                               1 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE */
/*                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST. */
/*       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                           CURRENT LINKED LIST OF MERGED INDICES (THE */
/*                           "REACH" SET). */
/*       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES */
/*                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S */
/*                           INDEX SET. */

/* *********************************************************************** */

/* Subroutine */ int symfc2_(neqns, adjlen, xadj, adjncy, perm, invp, colcnt, 
	nsuper, xsuper, snode, nofsub, xlindx, lindx, xlnz, mrglnk, rchlnk, 
	marker, flag__)
integer *neqns, *adjlen, *xadj, *adjncy, *perm, *invp, *colcnt, *nsuper, *
	xsuper, *snode, *nofsub, *xlindx, *lindx, *xlnz, *mrglnk, *rchlnk, *
	marker, *flag__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer head, node, tail, pcol, newi, jptr, kptr, jsup, ksup, psup,
	     i__, nzbeg, nzend, width, nexti, point, jnzbeg, knzbeg, length, 
	    jnzend, jwidth, fstcol, knzend, lstcol, knz;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    --marker;
    --xlnz;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --xadj;
    --adjncy;
    --mrglnk;
    --xlindx;
    --xsuper;
    --lindx;

    /* Function Body */
    *flag__ = 0;
    if (*neqns <= 0) {
	return 0;
    }

/*       --------------------------------------------------- */
/*       INITIALIZATIONS ... */
/*           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX. */
/*           TAIL   : END OF LIST INDICATOR */
/*                    (IN RCHLNK(*), NOT MRGLNK(*)). */
/*           MRGLNK : CREATE EMPTY LISTS. */
/*           MARKER : "UNMARK" THE INDICES. */
/*       --------------------------------------------------- */
    nzend = 0;
    head = 0;
    tail = *neqns + 1;
    point = 1;
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	marker[i__] = 0;
	xlnz[i__] = point;
	point += colcnt[i__];
/* L50: */
    }
    xlnz[*neqns + 1] = point;
    point = 1;
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
	mrglnk[ksup] = 0;
	fstcol = xsuper[ksup];
	xlindx[ksup] = point;
	point += colcnt[fstcol];
/* L100: */
    }
    xlindx[*nsuper + 1] = point;

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {

/*           ---------------------------------------------------------
 */
/*           INITIALIZATIONS ... */
/*               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP. */
/*               LSTCOL : LAST COLUMN OF SUPERNODE KSUP. */
/*               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL.
 */
/*               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL. */
/*           ---------------------------------------------------------
 */
	fstcol = xsuper[ksup];
	lstcol = xsuper[ksup + 1] - 1;
	width = lstcol - fstcol + 1;
	length = colcnt[fstcol];
	knz = 0;
	rchlnk[head] = tail;
	jsup = mrglnk[ksup];

/*           ------------------------------------------------- */
/*           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ... */
/*           ------------------------------------------------- */
	if (jsup > 0) {
/*               --------------------------------------------- */
/*               COPY THE INDICES OF THE FIRST CHILD JSUP INTO */
/*               THE LINKED LIST, AND MARK EACH WITH THE VALUE */
/*               KSUP. */
/*               --------------------------------------------- */
	    jwidth = xsuper[jsup + 1] - xsuper[jsup];
	    jnzbeg = xlindx[jsup] + jwidth;
	    jnzend = xlindx[jsup + 1] - 1;
	    i__2 = jnzbeg;
	    for (jptr = jnzend; jptr >= i__2; --jptr) {
		newi = lindx[jptr];
		++knz;
		marker[newi] = ksup;
		rchlnk[newi] = rchlnk[head];
		rchlnk[head] = newi;
/* L200: */
	    }
/*               ------------------------------------------ */
/*               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ... */
/*               ------------------------------------------ */
	    jsup = mrglnk[jsup];
L300:
	    if (jsup != 0 && knz < length) {
/*                   ---------------------------------------- 
*/
/*                   MERGE THE INDICES OF JSUP INTO THE LIST, 
*/
/*                   AND MARK NEW INDICES WITH VALUE KSUP. */
/*                   ---------------------------------------- 
*/
		jwidth = xsuper[jsup + 1] - xsuper[jsup];
		jnzbeg = xlindx[jsup] + jwidth;
		jnzend = xlindx[jsup + 1] - 1;
		nexti = head;
		i__2 = jnzend;
		for (jptr = jnzbeg; jptr <= i__2; ++jptr) {
		    newi = lindx[jptr];
L400:
		    i__ = nexti;
		    nexti = rchlnk[i__];
		    if (newi > nexti) {
			goto L400;
		    }
		    if (newi < nexti) {
			++knz;
			rchlnk[i__] = newi;
			rchlnk[newi] = nexti;
			marker[newi] = ksup;
			nexti = newi;
		    }
/* L500: */
		}
		jsup = mrglnk[jsup];
		goto L300;
	    }
	}
/*           --------------------------------------------------- */
/*           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET. */
/*           "SORT" ITS STRUCTURE INTO THE LINKED LIST, */
/*           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE */
/*           LIST. */
/*           --------------------------------------------------- */
	if (knz < length) {
	    node = perm[fstcol];
	    knzbeg = xadj[node];
	    knzend = xadj[node + 1] - 1;
	    i__2 = knzend;
	    for (kptr = knzbeg; kptr <= i__2; ++kptr) {
		newi = adjncy[kptr];
		newi = invp[newi];
		if (newi > fstcol && marker[newi] != ksup) {
/*                       -------------------------------- 
*/
/*                       POSITION AND INSERT NEWI IN LIST 
*/
/*                       AND MARK IT WITH KCOL. */
/*                       -------------------------------- 
*/
		    nexti = head;
L600:
		    i__ = nexti;
		    nexti = rchlnk[i__];
		    if (newi > nexti) {
			goto L600;
		    }
		    ++knz;
		    rchlnk[i__] = newi;
		    rchlnk[newi] = nexti;
		    marker[newi] = ksup;
		}
/* L700: */
	    }
	}
/*           --------------------------------------------------------
---- */
/*           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LI
ST. */
/*           --------------------------------------------------------
---- */
	if (rchlnk[head] != fstcol) {
	    rchlnk[fstcol] = rchlnk[head];
	    rchlnk[head] = fstcol;
	    ++knz;
	}

/*           -------------------------------------------- */
/*           COPY INDICES FROM LINKED LIST INTO LINDX(*). */
/*           -------------------------------------------- */
	nzbeg = nzend + 1;
	nzend += knz;
	if (nzend + 1 != xlindx[ksup + 1]) {
	    goto L8000;
	}
	i__ = head;
	i__2 = nzend;
	for (kptr = nzbeg; kptr <= i__2; ++kptr) {
	    i__ = rchlnk[i__];
	    lindx[kptr] = i__;
/* L800: */
	}

/*           --------------------------------------------------- */
/*           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S */
/*           "MERGE" LIST. */
/*           --------------------------------------------------- */
	if (length > width) {
	    pcol = lindx[xlindx[ksup] + width];
	    psup = snode[pcol];
	    mrglnk[ksup] = mrglnk[psup];
	    mrglnk[psup] = ksup;
	}

/* L1000: */
    }

    return 0;

/*       ----------------------------------------------- */
/*       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED. */
/*       ----------------------------------------------- */
L8000:
    *flag__ = -2;
    return 0;

} /* symfc2_ */

#ifdef DO_BFINIT
/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING */
/*       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */
/*       CACHSZ          -   CACHE SIZE (IN KBYTES). */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */

/* *********************************************************************** */

/* Subroutine */ int bfinit_(neqns, nsuper, xsuper, snode, xlindx, lindx, 
	cachsz, tmpsiz, split)
integer *neqns, *nsuper, *xsuper, *snode, *xlindx, *lindx, *cachsz, *tmpsiz, *
	split;
{
    extern /* Subroutine */ int fnsplt_(), fntsiz_();


/* ***********************************************************************
 */


/* ***********************************************************************
 */

/*       --------------------------------------------------- */
/*       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT. */
/*       --------------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --lindx;
    --xlindx;
    --snode;
    --xsuper;

    /* Function Body */
    fntsiz_(nsuper, &xsuper[1], &snode[1], &xlindx[1], &lindx[1], tmpsiz);

/*       ------------------------------- */
/*       PARTITION SUPERNODES FOR CACHE. */
/*       ------------------------------- */
    fnsplt_(neqns, nsuper, &xsuper[1], &xlindx[1], cachsz, &split[1]);

    return 0;
} /* bfinit_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE */
/*       REQUIRED BY BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */

/* *********************************************************************** */

/* Subroutine */ int fntsiz_(nsuper, xsuper, snode, xlindx, lindx, tmpsiz)
integer *nsuper, *xsuper, *snode, *xlindx, *lindx, *tmpsiz;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iend, clen, ksup, i__, bound, ncols, width, tsize, ibegin, 
	    length, cursup, nxtsup;


/* ***********************************************************************
 */



/* ***********************************************************************
 */

/*       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE. 
*/
/*       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY
 */
/*       TIGHT. */

/*       ---------------------------------------- */
/*       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR */
/*       NEEDED BY BLKFCT. */
/*       ---------------------------------------- */
    /* Parameter adjustments */
    --lindx;
    --xlindx;
    --snode;
    --xsuper;

    /* Function Body */
    *tmpsiz = 0;
    for (ksup = *nsuper; ksup >= 1; --ksup) {
	ncols = xsuper[ksup + 1] - xsuper[ksup];
	ibegin = xlindx[ksup] + ncols;
	iend = xlindx[ksup + 1] - 1;
	length = iend - ibegin + 1;
	bound = length * (length + 1) / 2;
	if (bound > *tmpsiz) {
	    cursup = snode[lindx[ibegin]];
	    clen = xlindx[cursup + 1] - xlindx[cursup];
	    width = 0;
	    i__1 = iend;
	    for (i__ = ibegin; i__ <= i__1; ++i__) {
		nxtsup = snode[lindx[i__]];
		if (nxtsup == cursup) {
		    ++width;
		    if (i__ == iend) {
			if (clen > length) {
			    tsize = length * width - (width - 1) * width / 2;
			    *tmpsiz = max(tsize,*tmpsiz);
			}
		    }
		} else {
		    if (clen > length) {
			tsize = length * width - (width - 1) * width / 2;
			*tmpsiz = max(tsize,*tmpsiz);
		    }
		    length -= width;
		    bound = length * (length + 1) / 2;
		    if (bound <= *tmpsiz) {
			goto L500;
		    }
		    width = 1;
		    cursup = nxtsup;
		    clen = xlindx[cursup + 1] - xlindx[cursup];
		}
/* L400: */
	    }
	}
L500:
	;
    }

    return 0;
} /* fntsiz_ */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     ***** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES */
/*       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE */
/*       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           POINTERS IN THE SUPERNODE INDICES. */
/*       CACHSZ          -   CACHE SIZE IN KILO BYTES. */
/*                           IF THERE IS NO CACHE, SET CACHSZ = 0. */

/*   OUTPUT PARAMETERS: */
/*       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE */
/*                           FINE PARTITIONING. */

/* *********************************************************************** */

/* Subroutine */ int fnsplt_(neqns, nsuper, xsuper, xlindx, cachsz, split)
integer *neqns, *nsuper, *xsuper, *xlindx, *cachsz, *split;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kcol, used, ksup, cache, ncols, width, height, curcol, 
	    fstcol, lstcol, nxtblk;


/* ***********************************************************************
 */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ******************************************************************* */

/*       -------------------------------------------- */
/*       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE. */
/*       -------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --xlindx;
    --xsuper;

    /* Function Body */
    if (*cachsz <= 0) {
	cache = 2000000000;
    } else {
	cache = (float) (*cachsz) * (float)1024. / (float)8. * (float).9;
    }

/*       --------------- */
/*       INITIALIZATION. */
/*       --------------- */
    i__1 = *neqns;
    for (kcol = 1; kcol <= i__1; ++kcol) {
	split[kcol] = 0;
/* L100: */
    }

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
/*           ----------------------- */
/*           ... GET SUPERNODE INFO. */
/*           ----------------------- */
	height = xlindx[ksup + 1] - xlindx[ksup];
	fstcol = xsuper[ksup];
	lstcol = xsuper[ksup + 1] - 1;
	width = lstcol - fstcol + 1;
	nxtblk = fstcol;
/*           -------------------------------------- */
/*           ... UNTIL ALL COLUMNS OF THE SUPERNODE */
/*               HAVE BEEN PROCESSED ... */
/*           -------------------------------------- */
	curcol = fstcol - 1;
L200:
/*               ------------------------------------------- */
/*               ... PLACE THE FIRST COLUMN(S) IN THE CACHE. */
/*               ------------------------------------------- */
	++curcol;
	if (curcol < lstcol) {
	    ++curcol;
	    ncols = 2;
	    used = height * 3 - 1;
	    height += -2;
	} else {
	    ncols = 1;
	    used = height << 1;
	    --height;
	}

/*               -------------------------------------- */
/*               ... WHILE THE CACHE IS NOT FILLED AND */
/*                   THERE ARE COLUMNS OF THE SUPERNODE */
/*                   REMAINING TO BE PROCESSED ... */
/*               -------------------------------------- */
L300:
	if (used + height < cache && curcol < lstcol) {
/*                   -------------------------------- */
/*                   ... ADD ANOTHER COLUMN TO CACHE. */
/*                   -------------------------------- */
	    ++curcol;
	    ++ncols;
	    used += height;
	    --height;
	    goto L300;
	}
/*               ------------------------------------- */
/*               ... RECORD THE NUMBER OF COLUMNS THAT */
/*                   FILLED THE CACHE. */
/*               ------------------------------------- */
	split[nxtblk] = ncols;
	++nxtblk;
/*               -------------------------- */
/*               ... GO PROCESS NEXT BLOCK. */
/*               -------------------------- */
	if (curcol < lstcol) {
	    goto L200;
	}
/* L1000: */
    }

    return 0;
} /* fnsplt_ */
#endif

