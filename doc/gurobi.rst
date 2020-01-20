.. _gurobi:

=====================
Using Gurobi with CVX
=====================

About Gurobi
------------

`Gurobi Optimization <http://www.gurobi.com>`_ was founded in 2008 by some of the most
experienced and respected members of the optimization community. The Gurobi solver quickly
became an industry performance leader in linear, quadratic, and mixed-integer programming.
Gurobi is a fantastic solver for use with CVX, particularly with the integer and binary
variable capability added in CVX 2.0.

Using Gurobi with CVX requires a valid license:

* *Academic users*: information about obtaining a license can be found on the
  `Gurobi Academic Program`_ page.
* *Commercial users* must purchase one of our CVX Professional licenses:

  * A *bundled CVX + Gurobi license* allows Gurobi to be used exclusively within
    CVX. This is the most cost-effective approach for users who do not intend
    to use Gurobi outside of CVX and/or MATLAB.
  * A *bring-your-own-solver (BYOS)* license allows CVX to be paired with a
    separate Gurobi license, enabling the same installation to be used within
    CVX and separate from it.

  Please contact `CVX Sales`_ for more information about either option, and
  `Gurobi Sales`_ for pricing information for standalone Gurobi licenses.

.. _gurobilic:

Using Gurobi with CVX
---------------------

1. Download the the appropriate CVX bundle from the `CVX download page`_
   and following the regular installation instructions at :ref:`install`.
   The standard bundles include a CVX-specific version of the Gurobi version 9.0.

2. Obtain the licenses for Gurobi and/or CVX, as needed:
   
   * A Gurobi license code, which is composed of 32 hexidecimal digits in the format
     ``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``. If you purchase a commercial CVX+Gurobi
     package, you will receive this code in an email from CVX Research. If you are an
     academic user, you will receive it directly from the `Gurobi academic license request`_ page.
   * A CVX license should be saved in a convenient location for Step 4. You
     will need to be able to supply its full path to the ``cvx_setup`` command.

3. If you need a full installation of Gurobi—because you wish to use a
   different version than is bundled with CVX, or because you wish to use Gurobi
   outside of CVX—obtain an appropriate installer from the `Gurobi Download Center`_ 
   and follow their instructions. Confirm that it can be successfully run from
   the MATLAB command line *before* proceeding with Step 2.

4. Next, retrieve your Gurobi license key by running the command ``cvx_grbgetkey``
   *{code}*, where *{code}* is the 32-digit Gurobi key. This is a convenience
   wrapper around Gurobi's own ``grbgetkey`` script. The command will look 
   something like this::

    cvx_grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
    
   *Important note for academic users:* this step must be run from a computer
   connected to your university network (a VPN is usually sufficient).
   Please consult the Gurobi documentation for more details.

5. Re-run ``cvx_setup`` so that the new Gurobi and/or CVX licenses can be detected.
   If a CVX Professional license was obtained, supply the path to this file as the
   argument to the ``cvx_setup`` command, as discussed in :ref:`licinstall`.

If successful, the output of step 4 should show that Gurobi is among the list
of available solvers. If you installed both a standalone Gurobi bundled version of Gurobi,
they should both be available after setup.

Selecting Gurobi as your default solver
--------------------------------------

Even if Gurobi is successfully added to your solver list, it will not automatically
be selected as your default solver. To change this, type the following two commands
on the MATLAB command line:

::

    cvx_solver gurobi
    cvx_save_prefs

The first command changes the active solver to Gurobi, but only for the current session.
The second line saves that change to CVX's preference file, so that Gurobi will be 
selected as the active solver every time you start MATLAB.

If multiple versions of Gurobi were found on the MATLAB path, then CVX will append a
numeral to the end of the solver name, allowing you to switch between them; e.g.,

::

    cvx_solver gurobi
    cvx_solver gurobi_2
    cvx_solver gurobi_3

and so forth.
    
Obtaining support for CVX and Gurobi
------------------------------------

If you encounter problems using CVX and Gurobi, please contact 
`CVX Support`_ first instead of Gurobi.
If we can reproduce your problem, we will determine whether or not it is an
issue that is unique to CVX or needs to be forwarded to Gurobi for further
analysis.

.. _Gurobi Optimization LLC: https://gurobi.com/
.. _CVX Support: http://support.cvxr.com/
.. _CVX download page: http://cvxr.com/cvx/download
.. _CVX Sales: mailto:sales@cvxr.com   
.. _Gurobi Academic Program: https://www.gurobi.com/academia/academic-program-and-licenses/
.. _Gurobi academic license request: https://www.gurobi.com/downloads/end-user-license-agreement-academic/
.. _Gurobi Download Center: https://www.gurobi.com/downloads/
.. _Gurobi Sales: mailto:sales@gurobi.com
