.. index::
	single: Installation
	single: CVX; installing

.. _install:

============
Installation
============
	
Supported platforms
-------------------

.. index::
	single: Platforms
	single: Matlab; versions
	single: Windows
	single: Linux
	single: Mac
	
CVX is supported on 32-bit and 64-bit versions of Linux, Mac OSX, and Windows. For 32-bit
platforms, MATLAB version 7.5 (R2007b) or later is required; for 64-bit platforms, MATLAB
version 7.8 (R2009a) or later is required. There are some important platform-specific
cautions, however:

- Gurobi support requires Matlab 7.7 (R2008b) or later.

- 32-bit Linux: the Gurobi solver is not available for this platform, as Gurobi is phasing
  out support for 32-bit Linux altogether.

- Older versions of Mac OS X (e.g. 10.5) ship with Java 1.5. The standard version of
  CVX works properly on this platform, but CVX Professional support requires Java 1.6.
  To restore this support, upgrade your operating system or Java installation.
  
As of version 2.0, support for versions 7.4 (R2007a) or older has been discontinued.
If you need to use CVX with these older versions of Matlab, please use CVX 1.22 or 
earlier, which will remain available indefinitely on the CVX Research web site. However,
this version is no longer supported, and will not receive bug fixes or improvements. 
We strongly encourage you to update your Matlab installation to the latest version
possible.

Installation instructions
-------------------------

.. index:: cvx_setup

.. note ::

	If you wish to use CVX with Gurobi or MOSEK, they must be installed and accessible
	from MATLAB *before* running ``cvx_setup``. See :ref:`below <extsolv>` for more details.

1. Retrieve the latest version of CVX from `the web site <http://cvxr.com/cvx/download>`_.
   You can download the package as either a ``.zip`` file or a ``.tar.gz`` file.
   
2. Unpack the file anywhere you like; a directory called ``cvx`` will be
   created. There are two important exceptions: 
   
   - *Do not* place CVX in Matlab's own ``toolbox`` directory.
   - *Do not* unpack a new version of CVX on top of an old one. We recommend moving the
     old version out of the way, but do not delete it until you are sure the new 
     version is working as you expect.

3. Start Matlab.

4. Change directories to the top of the CVX distribution, and run  the ``cvx_setup``
   command. For example, if you installed CVX into ``C:\Matlab\personal\cvx`` on
   Windows, type these commands:

   ::

       cd C:\Matlab\personal\cvx
       cvx_setup

   at the Matlab command prompt. If you installed CVX into
   ``~/MATLAB/cvx`` on Linux or a Mac, type these commands:
   
   ::

       cd ~/MATLAB/cvx
       cvx_setup
       
   The ``cvx_setup`` function performs a variety of tasks to verify that your 
   installation is correct, sets your Matlab search path so it can find all of the CVX 
   program files, and runs a simple test problem to verify the installation.       
       
5. In some cases---usually on Linux---the ``cvx_setup`` command may instruct you to 
   create or modify a ``startup.m`` file that allows you to use CVX without having
   to type ``cvx_setup`` every time you re-start Matlab.

.. index:: License; installing

.. _licinstall:

Installing a CVX Professional license
--------------------------------------

If you acquire a license key for CVX Professional, the only change required to the above
steps is to include the name of the license file as an input to the ``cvx_setup`` command.
For example, if you saved your license file to ``~/licenses/cvx_license.mat`` on a Mac,
this would be the modified command:

::

       cd ~/MATLAB/cvx
       cvx_setup ~/licenses/cvx_license.mat
       
If you have previously run ``cvx_setup`` without a license, or you need to replace your
current license with a new one, simply run ``cvx_setup`` again with the filename.
Once the license has been accepted and installed, you are free to move your license 
file anywhere you wish for safekeeping---CVX saves a copy in its preferences.

.. index::
	single: SeDuMi
	single: Solvers; SeDuMi
	single: SDPT3
	single: Solvers; SDPT3
	single: MOSEK
	single: Solvers; MOSEK
	single: Gurobi
	single: Solvers; Gurobi
	single: Solvers; included
	single: Solvers
	
.. _extsolv:

Solvers included with CVX
-------------------------

All versions of CVX include copies of the solvers
`SeDuMi <http://sedumi.ie.lehigh.edu/>`_
and 
`SDPT3 <http://www.math.nus.edu.sg/~mattohkc/sdpt3.html>`_
in the directories :file:`cvx/sedumi` and :file:`cvx/sdpt3`, respectively. When you
run `cvx_setup`, CVX will automatically add these solvers to its solver list.

If you have downloaded a CVX Professional Solver Bundle, then the solvers 
`Gurobi <http://gurobi.com>`_
and/or 
`MOSEK <http://mosek.com>`_ will be included with CVX as well. Use of these
solvers requires a CVX Professional license. You may also use your existing
copies of these solvers with CVX as well. We have created special sections of
this users' guide for each solver:

* Gurobi: :ref:`gurobi`
* MOSEK: :ref:`mosek`

For more general information on the solvers supported by CVX, an how to select a
solver for your particular problem, see the :ref:`Solvers <solvers>` section.
