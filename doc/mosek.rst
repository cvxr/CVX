.. _mosek:

====================
Using MOSEK with CVX
====================

About MOSEK
-----------

`MOSEK ApS`_ is widely considered the leader in commercial software
for nonlinear convex optimization. The company was established in 1997, and is
led by founding CEO `Erling Andersen <http://www.linkedin.com/in/edandersen>`_ and
a technical advisory board chaired by
Stanford Professor `Yinyu Ye <http://www.stanford.edu/~yyye/>`_. Both are internationally
recognized for their contributions to the field of convex optimization, and remain active 
in research and publication. With its support for integer variables, the semidefinite cone,
and (with version 9.0) the exponential cone, the MOSEK solver has native support for a
wider variety of CVX models than any other solver.

Using MOSEK with CVX requires a valid license:

* *Academic users*: request an license from the `MOSEK Academic Licensing`_ page.
* *Commercial users* must purchase one of our CVX Professional licenses:

  * A *bundled CVX + MOSEK license* allows MOSEK to be used exclusively within
    CVX. This is the most cost-effective approach for users who do not intend
    to use MOSEK outside of CVX and/or MATLAB.
  * A *bring-your-own-solver (BYOS)* license allows CVX to be paired with a
    separate MOSEK license, enabling the same installation to be used within
    CVX and separate from it.

  Please contact `CVX Sales`_ for more information about either option, and
  `MOSEK ApS Sales`_ for pricing information for standalone MOSEK licenses.

Using MOSEK with CVX
--------------------

1. Download the the appropriate CVX bundle from the `CVX download page`_
   and following the regular installation instructions at :ref:`install`.
   The standard bundles include a CVX-specific version of the MOSEK version 9.1.

2. Obtain the licenses for MOSEK and/or CVX, as needed:
   
   * A MOSEK license should be installed in the location ``mosek/mosek.lic``
     in your home directory.
   * A CVX license should be saved in a convenient location for Step 4. You
     will need to be able to supply its full path to the ``cvx_setup`` command.

3. If you need a full installation of MOSEK—either because you wish to use a
   different version than is bundled with CVX, or because you wish to use MOSEK
   outside of CVX—obtain an appropriate installer from the `MOSEK download page`_ 
   and follow their instructions. Confirm that it can be successfully run from
   the MATLAB command line *before* proceeding with Step 2.

4. Re-run ``cvx_setup`` so that the new MOSEK and/or CVX licenses can be detected.
   If a CVX Professional license was obtained, supply the path to this file as the
   argument to the ``cvx_setup`` command, as discussed in :ref:`licinstall`.

If successful, the output of step 4 should show that MOSEK is among the list
of available solvers. If you installed both a standalone and bundled version of MOSEK,
they should both be available after setup.

Selecting MOSEK as your default solver
--------------------------------------

Even if MOSEK is successfully added to your solver list, it will not automatically
be selected as your default solver. To change this, type the following two commands
on the MATLAB command line:

::

    cvx_solver mosek
    cvx_save_prefs

The first command changes the active solver to MOSEK, but only for the current session.
The second line saves that change to CVX's preference file, so that MOSEK will be 
selected as the active solver every time you start MATLAB.

If multiple versions of MOSEK were found on the MATLAB path, then CVX will append a
numeral to the end of the solver name, allowing you to switch between them; e.g.,

::

    cvx_solver mosek
    cvx_solver mosek_2
    cvx_solver mosek_3

and so forth.
    
Obtaining support for CVX and MOSEK
------------------------------------

If you encounter problems using CVX and MOSEK, please contact 
`CVX Support`_ first instead of MOSEK ApS.
If we can reproduce your problem, we will determine whether or not it is an
issue that is unique to CVX or needs to be forwarded to MOSEK ApS for further
analysis.

.. _MOSEK ApS: https://mosek.com/
.. _CVX Support: http://support.cvxr.com/
.. _CVX download page: http://cvxr.com/cvx/download
.. _CVX Sales: mailto:sales@cvxr.com   
.. _MOSEK Academic Licensing: https://www.mosek.com/products/academic-licenses/
.. _MOSEK download page: https://www.mosek.com/downloads/
.. _MOSEK ApS Sales: mailto:sales@mosek.com
