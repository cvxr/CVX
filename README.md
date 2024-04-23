## CVX: A system for disciplined convex programming

#### [Click here](https://github.com/cvxr/cvx/releases/latest) to download a bundle of this repository, including pre-compiled MEX files.

### Important Update

This repository is in the process of becoming the *official* source
for CVX. We will be updating the [original site](https://cvxr.com/cvx)
and its [download page](https://cvxr.com/cvx/download) with links back
to this site. When this work is complete, the CVX bundles hosted here
will provide pre-compiled Matlab MEX files for Windows, Linux, macOS
Intel, and macOS Apple Silicon. Furtherore, open versions of the Mosek
and Gurobi solver shims will be available. Please stay tuned, here and
on the [web site](https://cvxr.com/cvx), for further developments.

### Introduction

CVX is a Matlab package for convex optimizaton.
To learn more about what CVX is and how to use it, please visit our
[web site](https://cvxr.com/cvx), read the
[users' guide](https://cvxr.com/cvx/doc), and browse the
[example library](https://cvxr.com/cvx/examples). 

The best way to obtain CVX is to visit the
[Releases Page](https://github.com/cvxr/CVX/releases/latest/) of
this repository, and download _either_ `cvx.tgz` or `cvx.zip`.
These archives contain:

- The full CVX code base
- Full copies of the [SeDuMi](https://github.com/sqlp/sedumi) and
  SDPT3 [SDPT3](https://github.com/sqlp/sdpt3) solvers
- Additional shims for Gurobi, Mosek, and GLPK. The solvers (and
  any needed licenses) must be be supplied by the user.
- Pre-compiled MEX files for Windows, macOS, and Linux
- HTML and PDF versions of the documentation

### About this repository

For now, this is a _filtered_ mirror of the main branch of our
internal development repository. This internal repository includes
code not ready or available for redistribution. Some of that code will
was built strictly to support our dual-source approach, and is
therefore no longer needed. There remains some additional code there,
however, that would likely be valuable to users, and we intend to
bring that here over time.

### Support

There are four primary mechanisms for obtaining support for CVX:

- The [user guide](https://cvxr.com/cvx/doc). PDF and HTML versions
  of this guide are included in `cvx.zip` and `cvx.tgz` bundles.
  The [documentation source](https://github.com/cvxr/CVX/tree/master/doc)
  is available on the repository as well.

- The [example library](https://cvxr.com/cvx/examples/). Many user
  problems can be solved as slight modifications to one of these
  examples. The example code is also included in the `cvx.zip`
  and `cvx.tgz` bundles, in the `examples/` subdirectory, and
  [in the repository](https://github.com/cvxr/CVX/tree/master/examples).
  
- The [CVX Forum](https://ask.cvxr.com/) is a Discourse-based
  server that is focused completely on CVX usage.
  
- The [Computational Science Stack Exchange](https://scicomp.stackexchange.com/)
  is a great community-driven Q&A site for a variety of
  computational science topics, including convex optimization.
  This would be a perfect choice for questions that are not
  necessarily specific to CVX.

Easily the most important page on the CVX Forum is the FAQ:

> [***Why isn't CVX accepting my model? READ THIS FIRST!***](https://ask.cvxr.com/t/why-isnt-cvx-accepting-my-model-read-this-first/570)
 
*Everyone* who attempts to use CVX should read that page! It should
save much frustration.

Supporting CVX users is a challenging exercise, because there are
multiple *categories* of issues that occur when using the software;
including, but not limited to:

1. Generic challenges involving convex optimization, especially the
   challenge of determining whether or not a given model is convex.
   See the [FAQ](https://ask.cvxr.com/t/why-isnt-cvx-accepting-my-model-read-this-first/570).
2. Difficulties caused by the strictness imposed by the disciplined
   convex programming framwork. Once again, see the
   [FAQ](https://ask.cvxr.com/t/why-isnt-cvx-accepting-my-model-read-this-first/570).
3. Numerical issues caused by model scaling issues, or even bugs or
   limitations in the underling solver software. The
   [forum](https://ask.cvxr.com/) can sometimes offer assistance here.
4. _Actual bugs_ caused by _unintended_ behavior of the software.

Our experience is that the _vast_ majority of user issues fall into
one of the first three categories; a small fraction fall into
category 4. In other words, _most user issues are not bugs_.

That does not mean CVX is bug-free! If you have truly eliminated
issue categories 1-3 above from contention, including a thorough
understanding of the
[FAQ](https://ask.cvxr.com/t/why-isnt-cvx-accepting-my-model-read-this-first/570),
you may wish to submit a bug request to this repository's 
[issue tracker](https://github.com/cvxr/CVX/issues).

Unfortunately, urgent attention to bugs cannot be promised, yet, but
as the contributor community is cultivated, they will begin to be
able to tackle issues and perhaps even offer new features. Furthermore,
if we determine an issue is _not_ a bug, it will be closed, often
without a response. For this reason we truly commend users to the
[CVX Forum](https://ask.cvxr.com/) before considering filing an issue.

_The original authors of this package are
not available to assist directly._ Support emails sent to the authors
will go unanswered.

### Citing CVX

Are you using CVX in research work to be published? If so, please include explicit
mention of our work in your publication. We have provided example language, citation
entries, and even BiBTeX citation code on the
[Citing CVX](https://cvxr.com/cvx/doc/citing.html) page of the documentation.

### License

Most of the files in this repository are governed by the terms of our 
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)-based 
[CVX Standard License](https://cvxr.com/cvx/doc/license.html). Please
see the files 
[LICENSE.txt](https://github.com/cvxr/CVX/blob/master/README.txt) and 
[GPL.txt](https://github.com/cvxr/CVX/blob/master/GPL.txt), 
included with the distribution, for more details about the license.

The contents of the example library in the `examples/` subdirectory
are *public domain*. We do ask that if you use any of this content in
your own work, that you acknowledge the source and any specific authors
cited therein.

Thank you for your interest in CVX!    
Michael Grant and Stephen Boyd
