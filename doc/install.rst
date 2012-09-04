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

CVX is supported on the following platforms:

* Microsoft Windows, Linux, and Mac OSX
* MATLAB versions 7.5 (R2007b) or later
* 32-bit and 64-bit platforms

As of version 2.0, support for versions of Matlab more than five years old---specifically,
Matlab 7.4 (R2007a) or older---has been discontinued. If you need to use CVX with these 
older versions of Matlab, please use CVX 1.22 or earlier, which will remain available
indefinitely on the CVX Research web site.

It is theoretically possible to use CVX on other platforms, but we are
unable to provide support of any kind. In particular, you will have to compile 
the MEX files for both CVX and the solvers.

Installation instructions
-------------------------

.. index:: cvx_setup

1. Retrieve the latest version of CVX from `the web site <http://cvxr.com/cvx/download>`_.
   You can download the package as either a ``.zip`` file or a ``.tar.gz`` file.

2. Unpack the file anywhere you like; a directory called ``cvx`` will be
   created. There is one important exception: *do not* place CVX in Matlab's
   own ``toolbox`` directory.

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
	single: Solvers; included
	single: Solvers

About the included solvers
---------------------------

The CVX distribution includes copies of the solvers 
`SeDuMi <http://sedumi.ie.lehigh.edu/>`_
and 
`SDPT3 <http://www.math.nus.edu.sg/~mattohkc/sdpt3.html>`_
in the directories :file:`cvx/sedumi` and :file:`cvx/sdpt3`, respectively. We have
designed CVX to use its own copy of these solvers, because we can better support the 
specific version that we have chosen. Indeed, CVX has generated quite a few bug reports
for these solvers! However, you are free to keep your alternate copy in your
MATLAB path. When you are not constructing a CVX model, MATLAB will rely on your
copy of the solver instead.
