.. _licensing2:

=======
License
=======

| CVX: A system for disciplined convex programming   
| © 2012 CVX Research, Inc., Austin, TX.
| http://cvxr.com
| info@cvxr.com

Thank you for using CVX!

The files contained in the CVX distribution come from several different sources and
are covered under a variety of licenses. The files owned by CVX Research, Inc. are
covered under one of two licenses: the *CVX Standard License* and the *CVX Professional
License*. The CVX Standard License is effectively the GNU General Public License, Version
2 (GPLv2), but with additional terms that govern modifications that connect CVX to 
additional solvers. The added terms are discussed in "Solver Interfaces" below.

CVX Professional License
------------------------

The standard CVX distribution includes several files in Matlab *p-code* format, which
contain encrypted binary versions of Matlab bytecode. As their name implies, they are 
recognized by their ``.p`` suffix. Currently the following files are distributed this way:

::

    shims/cvx_mosek.p
    shims/cvx_gurobi.p
    cvx_license.p
    
You may redistribute these files only as part of the complete, unmodified CVX package as
distributed by CVX Research, Inc. itself.

CVX Standard License
--------------------

If you wish, you may distribute a version of CVX with all of the p-code files *removed*.
The resulting package retains full functionality with the exception of its ability to
connect to commercial solvers. This modified package is covered by the CVX Standard
License. Under this license, you are free to redistribute and/or modify the files under
the terms of the GPLv2, plus the additional terms discussed in "Solver Interfaces" below. 

You must include the files ``LICENSE.txt`` and ``GPL.txt`` in unmodified form when 
redistributing this software. If you did not receive a copy of either of these files
with your distribution, please contact us.

Solver Interfaces
-----------------

CVX relies upon other software packages, called *solvers*, to perform many of its 
underlying calculations. Currently CVX supports free solvers SeDuMi and SDPT3, and 
commercial solvers Gurobi and MOSEK. The resulting nexus of free and commercial
software presents a licensing challenge. Our vision is guided by three goals:

- to ensure that CVX remains free to use by *all* users with any compatible *free* solver.
- to generate revenue by selling interfaces to *commercial* solvers to *commercial* customers.
- to provide the academic community with the full commercial capability at no charge.

The terms we lay out here are intended to support these trifold goals.  

We invite our users to create new interfaces between CVX and other *free* solvers. By 
"free", we mean that the solver must be made available at no charge for to *all* users,
including commercial users, without restriction. Please contact us if you are interested
in creating such an interface; we can offer assistance. If you do create one, please 
consider submitting it to us for inclusion in the standard CVX distribution. But you are 
under no obligation to do this. Instead, you can ship the interface code with the solver
itself; or you can construct a modified version of CVX with your interface included.

We do not permit the creation and distribution of new interfaces between CVX and 
*non-free* solvers---even if those solvers are made available to academic users at no 
charge. If you are a vendor or developer of a commercial solver, and would like to develop
or offer a CVX interface to your users, please contact us at info@cvxr.com. We welcome
the opportunity to support a wider variety of commercial solvers with CVX, and are
willing to devote engineering resources to make those connections.

If you are a user of a particular commercial solver and would like to see it supported by 
CVX, please contact your solver vendor---but please contact us at info@cvxr.com as well. 
If there is sufficient demand, and it proves technically and financially feasible, we will
reach out to the solver vendor to work on an implementation.

Bundled solvers
----------------

The solvers SDPT3 and SeDuMi are distributed with CVX in the ``sdpt3/`` and ``sedumi/``
subdirectories, respectively. Neither of these packages is owned by CVX Research, Inc.
Both are included with permission of the authors, and licensed under the terms of
the GPLv2. Please consult the plain-text documentation contained in each of these
directories for more information about copying, citation, and so forth.

Example library
---------------

The contents of the example library, which is distributed with CVX in the ``examples/``
subdirectory, is *public domain*. You are free to use them in any way you wish; but when 
you do, we request that you give appropriate credit to the authors. A number of people 
have contributed to the examples in this library, including Lieven Vandenberghe, 
Joëlle Skaf, Argyris Zymnis, Almir Mutapcic, Michael Grant, and Stephen Boyd.

No Warranty
-----------

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.