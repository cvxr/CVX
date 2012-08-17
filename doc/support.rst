.. _support:

=======
Support
=======

The user base for CVX has grown to such an extent that email-based
support is no longer tenable. Therefore, we have created several avenues
for obtaining support.

For help on how to *use* CVX, this document is your first line of support.
Please make sure that you have attempted to find an answer to your question
here before you pursue other avenues. We have placed this document
online and made it searchable in order to help you find the answers to 
the questions you may have.

The CVX user community
-----------------------

If your answers cannot be found here, consider posting your question
to the `CVX Forum <http://ask.cvxr.com>`_. This is a community forum
that allows our users to submit questions and to answer other people's questions.
The forum uses the open-source `Askbot <http://www.askbot.com>`_ system, and its format
should be familiar to anyone who participates in `OR-Exchange <http://www.or-exchange.com>`_,
`Stack Overflow <http://stackoverflow.com>`_, or any one of the `Stack Exchange <http://stackexchange.com>`_
family of sites. 

We highly encourage our expert users who enjoy helping others to participate in
this forum. We hope that it will not only serve as a resource for diagnosing problems
and issues, but a clearinghouse for advanced usage tips and tricks.

Bug reports
-----------

If you believe you have found a *bug* in CVX or in one of the underlying solvers, 
then we encourage you to  submit a bug report directly to
`CVX Research Support (support@cvxr.com) <mailto:support@cvxr.com>`_. Please include the following in your
bug report so we can fully reproduce the problem:

1. the CVX model and supporting data that caused the error. 
2. a copy of any error messages that it produced
3. the CVX version number and build number
4. the version number of Matlab that you are running
5. the name and version of the operating system you are using

The easiest way to supply items 3-5 is to type ``cvx_version`` at the command
prompt and copy its output into your email message.

Please note that we do not own all of Matlab's toolboxes. We cannot debug a model that
employs functions from a toolbox we do not own.

What *is* a bug?
~~~~~~~~~~~~~~~~

Certain issues are unambiguously bugs, and you should feel free to report them 
immediately. In particular, CVX often attempts to catch unexpected errors in key
places---including ``cvx_setup``, ``cvx_end``, etc. It will report those errors and
instruct you to report them to us.

If your model produces a MATLAB error that CVX did not itself generate, and you cannot
readily tie it to a syntax error in your model, please report that as well.

That said, because disciplined convex programming is a new concept for many, we often 
receive support requests for problems that are in fact *not* bugs. 

A particularly common
class of support requests are due to
``Disciplined convex programming error`` messages. This message indicates that the model fails
to adhere to the rules in the :ref:`DCP rulesert <dcp>`. In nearly all cases,
the underlying issue is in one of two categories:

1. The problem is *not convex* (and not an :ref:`MIDCP <what-is-midcp>`). In this case,
   no amount of manipulation of the expressions will cause CVX to solve the problem. 
   In some rare cases, it
   is possible to transform a problem that is non-convex into one that is convex (for
   example, :ref:`geometric programs <gp-mode>`). This has not been the case for any
   problem submitted by our users---so far.
   
2. The problem *is* convex (or an MIDCP), but requires the use of functions that do not
   exist in the CVX library. We have attempted to supply all of the commonly used 
   functions that the underlying solvers can support; so if you cannot easily rewrite 
   your problem using the functions supplied, it may not be possible. If you think this
   is a possibility, you may wish to see if the wizards on the
   `CVX Forum <http://ask.cvxr.com>`_ have suggestions for you.
   
In rare cases, users have discovered that certain models were rejected 
even though they satisfied the :ref:`DCP ruleset <dcp>`.
We have not received a bug report of this type in quite a long time, however, so we
suspect our users have helped us iron out any of these issues.

CVX Professional support
-------------------------

Paid CVX Professional users will receive support through a trouble-ticket support system,
so that they can have confidence that their issues are being addressed promptly. This
infrastructure is still under development; we will update this section and the Web site
with more information once it has been completed.
