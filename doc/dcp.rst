.. _dcp:

===============
The DCP ruleset
===============

*Disciplined convex programming* (DCP) requires all models obey a set of 
rules, or conventions, that govern how expressions and functions can appear 
in objectives and constraints. These rules, which we call the *DCP ruleset*, 
are drawn from basic principles of convex analysis and are relatively easy 
to learn. But they are not exhaustive, which means that it is possible to 
construct expressions and models that are known to be convex, but still 
violate the rules.

To illustrate the difference between *convexity* and *disciplined convexity*,
consider the function :math:`f(x)=\sqrt{x^2+1}`. It is simple to prove that 
this function is convex; its second derivative 
:math:`f^{(2)}(x)=(x^2+1)^{-3/2}` is positive for all :math:`x`. The ruleset
is silent on derivatives, however. And if you attempt to express :math:`f`
in the obvious manner in CVX, you will see this error:

::

    >> sqrt(x^2+1)
    Error using cvx/sqrt (line 61)
    Disciplined convex programming error:
        Illegal operation: sqrt( {nonnegative convex} ).

The problem is that the :ref:`composition rules <compositions>` forbid the 
application of a concave function (like ``sqrt``) to a convex expression 
(like ``x^2+1``), since that *usually* produces a nonconvex result. So CVX
rejects the expression on this basis.

Fortunately, most situations like this can be resolved simply by rewriting 
the expression. In this case, :math:`f` can also be written as 
:math:`f(x)=\|[x~1]\|_2`; a form which CVX accepts:

::

    >> norm([x 1])
    ans =     
        cvx nonnegative convex expression (scalar)

This expression is acceptable because ``norm`` is among the
:ref:`functions supported by CVX <funcref>`, and it is being used in a 
manner compliant with the :ref:`composition rules <compositions>`. So  
:math:`f(x)=\sqrt{x^2+1}` can indeed be used in CVX models, as long as it is
expressed in a compliant manner. 

At first, the DCP ruleset  may seem arbitrary or restrictive, but it serves a
very important purpose. Each rule corresponds to a specific step that CVX
takes to convert models to solvable form. When CVX rejects a model for a rule
violation, then, it is doing so *because it does not know how to solve it*.
Put another way, by complying with the rules, you are not only proving
that your model is convex, you are also giving CVX *a detailed recipe for
solving it.* (We thank you for the help!)

In truth, it is never the *rules* that ultimately prevent a model from being
represented in CVX. Rather, it is the finite size of the function library,
which is in turn limited by the capabilities of the underlying numerical 
solvers. Still, the greater your mastery of the DCP ruleset, the more 
productive you will be.

Top-level rules
---------------

Objectives
~~~~~~~~~~

Acceptable objective expressions come in one of two forms:

- Convex minimization: ``minimize(`` *expr* ``)``, 
  where *expr* is convex (or real affine).
- Concave maximization: ``maximize(`` *expr* ``)``, 
  where *expr* is concave (or real affine).

A model does not need to have an objective; such a problem is called a
*feasibility problem*. CVX will attempt to find a single point that  satisfies
all of the constraints.

Constraints
~~~~~~~~~~~

Acceptable constraints come in one of four forms:

- A *less-than inequality constraint*, using ``<=``, where the left
  side is convex and the right side is concave.
- A *greater-than inequality constraint*, using ``>=``, where the left
  side is concave and the right side is convex.
- An *equality constraint*, using ``==``, where both the left and
  right-hand sides are affine.
- A *set membership constraint*, using ``<In>``, involving affine
  expressions. (See :ref:`sets` for more details.)

*Non*-equality constraints, constructed using ``~=``, are *never allowed*.
(Such constraints are not convex.)

CVX treats strict ``<`` ``>`` inequalities identically to non-strict  ``<=``
``>=`` inequalities, so to avoid confusion the use of strict  inequalities is
*strongly discouraged*. For more information,  see :ref:`strict` below.

Inequality constraints must be real. Equality constraints, on the
other hand, may be complex. Complex equality constraints are equivalent
to two real equality constraints, one for the real part and one for
the imaginary part. An equality constraint with a real side and a complex 
side, then, has the effect of constraining the imaginary part
of the complex side to be zero.

.. _expressions:

Expression rules
----------------

Each scalar expression and subexpression is analyzed to determin
its *curvature* and *sign*. Vectors, matrices, and arrays are analyzed on an 
elementwise basis.

Curvature
~~~~~~~~~

CVX considers four types of *curvature*: curvature:  *constant*, *affine*,
*convex*, and *concave*. The reader should already be familiar with these
definitions. But for review, a function
:math:`f:\mathbf{R}^n\rightarrow\mathbf{R}` defined on all
:math:`\mathbf{R}^n`, the categories have the following meanings:

.. math::

  \begin{array}{l@{\quad}ll}
    \text{constant} & f(\alpha x + (1-\alpha)y) =  f(x) 
    & \forall x,y\in\mathbf{R}^n,~\alpha\in\mathbf{R} \\
    \text{affine}   & f(\alpha x + (1-\alpha)y) = \alpha f(x) + (1-\alpha) f(y)
    & \forall x,y\in\mathbf{R}^n,~\alpha\in\mathbf{R} \\
    \text{convex}   & f(\alpha x + (1-\alpha)y) \leq \alpha f(x) + (1-\alpha) f(y) 
    & \forall x,y\in\mathbf{R}^n,~\alpha\in[0,1] \\
    \text{concave}  & f(\alpha x + (1-\alpha)y) \geq \alpha f(x) + (1-\alpha) f(y) 
    & \forall x,y\in\mathbf{R}^n,~\alpha\in[0,1]
  \end{array}

There is, of course, significant overlap in these  categories: constant
expressions are also affine, and (real) affine expressions are both convex and
concave. Convex and concave expressions are real by definition, but constants
and affine expressions can be complex.

CVX does *not* determine convexity using the above definitions. Instead,
curvature is determined recursively applying the following rules.  While this
list may seem long, it is for the most part  an enumeration of basic rules of
convex analysis for combining convex, concave, and affine forms: sums,
multiplication by scalars, and so forth.

-  A valid constant expression is

   -  any well-formed expression that immediately evaluates to a finite
      value.

-  A valid affine expression is

   -  a valid constant expression;
   -  a declared variable;
   -  the sum or difference of affine expressions;
   -  the product of an affine expression and a constant.
   -  a valid affine function expression---see :ref:`compositions`;

-  A valid convex expression is

   -  a valid constant or affine expression;
   -  the sum of two or more convex expressions;
   -  the difference between a convex expression and a concave
      expression;
   -  the product of a convex expression and a nonnegative constant;
   -  the product of a concave expression and a nonpositive constant;
   -  the negation of a concave expression;
   -  a valid convex function expression---see :ref:`compositions`;
   -  an affine scalar raised to a constant power :math:`p\geq 1`,
      :math:`p\neq3,5,7,9,...`;
   -  a convex scalar quadratic form---see :ref:`quadforms`.

-  A valid concave expression is

   -  a valid constant or affine expression;
   -  the sum of two or more concave expressions;
   -  the difference between a concave expression and a convex expression;
   -  the product of a concave expression and a nonnegative constant;
   -  the product of a convex expression and a nonpositive constant;
   -  the negation of a convex expression;
   -  a valid concave function expression---see :ref:`compositions`;
   -  a concave scalar raised to a power :math:`p\in(0,1)`;
   -  a concave scalar quadratic form---see :ref:`quadforms`.

We note that the set of rules listed above is redundant; there are much
smaller, equivalent sets of rules. For matrix and array expressions, these
rules are applied on an elementwise basis.

Of particular note is that these expression rules generally forbid *products*
between nonconstant expressions, with the exception of scalar quadratic forms.
For example, the expression ``x*sqrt(x)`` happens to be a convex function of
``x``, but its convexity cannot be verified using the CVX ruleset, and so is
rejected. (It can be expressed as ``pow_p(x,3/2)``, however.)  We call this
the *no-product rule*:

- The product or ratio of two non-constant (affine, convex, concave)
  expressions is forbidden.

Adherence to the no-product rule will go a long way to insuring that you
construct valid expressions. There is one notable exception to this rule,
however: see :ref:`quadforms` below. But quadratic forms are, strictly
speaking, an unnecessary convenience, since CVX includes a ``quad_form``
function that provides the same functionality.

.. _sign:

Sign
~~~~

CVX also keeps track of the *sign* of an expression as well. Expressions are
classified as *positive*, *negative*, and *unknown sign*. In a slight abuse of
notation, nonnegative expressions are also treated as positive, and
nonpositive expressions are also treated as negative. It should be noted that
CVX does *not* perform any sort of advanced interval analysis to determine if
an expression is positive or negative. As with curvature, it draws its
conclusions by applying a simple set of rules:

- A "positive" expression is
  
  - a positive constant (or zero);
  - a variable *declared* `nonnegative` (see :ref:`variables`);
  - a diagonal element of a variable declared `semidefinite`  
    (see :ref:`variables`);
  - a call to any function specifically labeled as *positive* 
    (see :ref:`functions` below);
  - a negative expression multiplied by a negative constant;
  - a positive expression multiplied by a positive constant;
  - the sum of positive expressions.

- A "negative" expression is 

  - a negative constant (or zero);
  - a call to any function specifically labeled as *negative* 
    (see :ref:`functions` below);
  - a negative expression multiplied by a positive constant;
  - a positive expression multiplied by a negative constant;
  - the sum of negative expressions.

That's it! Any expression whose sign cannot be determined from these rules is
classified as having *unknown sign*. For example, the expression ``x - 1`` has
unknown sign---even if a constraint in the model ensures that ``x >= 1``.
These rules provide just enough information to CVX to give the user more
flexibility in how it combines functions together; more on this in 
:ref:`sign-monotonicity` below.

Function expressions
--------------------

Now let us consider how CVX classifies an expression of the form
:math:`f(\arg_1,\arg_2,\dots,\arg_n)`, where :math:`f` is a function from
CVX's function library, and each argument :math:`arg_k` is an otherwise well-
posed scalar  CVX expression. In the case where a MATLAB function accepts
vector, matrix, or array arguments, everything we discuss here is applied in
an elementwise fashion. For instance, the  expression `norm(x)`, where `x` is
a vector of length :math:`n`, can be thought of as a function expression
involving :math:`n` separate scalar arguments.

Function classification
~~~~~~~~~~~~~~~~~~~~~~~

In order to proceed, we must first understand the properties of the function
:math:`f` itself. As with basic expressions, CVX categorizes functions
according to their *curvature* and *sign*. They also obtain two more
attributes as well: *monotonicity* and *domain*. For functions with only one
argument, the categorization is straightforward. Some examples are given in
the table below.

.. tabularcolumns:: CCCCCCC

================== ==================== =========== ================ ========== ===========================
 Function           Meaning              Curvature   Monotonicity     Sign       Domain
================== ==================== =========== ================ ========== ===========================
 ``sum( x )``       :math:`\sum_i x_i`   affine      increasing       unknown    :math:`\mathbb{R}`
 ``abs( x )``       :math:`|x|`          convex      sign-dependent   positive   :math:`\mathbb{R}`
 ``inv_pos( x )``       :math:`1/x`          convex      decreasing       unknown    :math:`\{x\,|\,x>0\}`
 ``sqrt( x )``      :math:`\sqrt x`      concave     increasing       positive   :math:`\{x\,|\,x\geq 0\}`
 ``log( x )``   :math:`\log x`       concave     increasing       unknown    :math:`\{x\,|\,x>0\}`
 ``entr( x )``      :math:`-x\log x`     concave     non-monotonic    unknown    :math:`\{x\,|\,x\geq 0\}`
================== ==================== =========== ================ ========== ===========================

Domain
======

The *domain* of a function is simply the set of points over which a function
is well-defined. For a convex or concave function, this set is always convex.
The domain serves as an *implicit constraint* on the function's input. For
instance, if we form  ``sqrt(x+1)`` in a CVX specification, a new constraint
``x+1>=0`` is automatically assumed. There is no need to add such a constraint
separately. Monotonicity is also considered with respect to the function's
domain; so, for instance, ``sqrt(x)`` is considered increasing, since that is
indeed the case for all nonnegative inputs.

CVX does *not* consider a function to be convex or concave if it is so only
over a portion of its domain, even if the argument is constrained to lie in
one of these portions. For example, consider the function :math:`1/x`. This
function is convex for :math:`x>0`, and concave for :math:`x<0`. But you can
never write ``1/x`` in CVX (unless ``x`` is constant), even if you have
imposed a constraint such as ``x>=1``, which restricts ``x`` to lie in the
convex portion of function. You *can*, however, use the CVX function
``inv_pos(x)``, listed above, which is defined to have the domain
`:math:\mathbb{R}_{++}`. CVX recognizes this function as convex and
decreasing.

Monotonicity
============

CVX considers two types of monotonicity: *increasing* and *decreasing*. In a
slight abuse of notation, we classify nondecreasing functions as increasing, 
and nonincreasing functions as decreasing. These categories have the following
meanings:

.. math::

  \begin{array}{l@{\quad}l}
    \text{increasing} & x \geq y ~~\Longrightarrow~~ f(x) \geq f(y) \\
    \text{decreasing} & x \geq y ~~\Longrightarrow~~ f(x) \leq f(y)
  \end{array}   

A function that is neither increasing or decreasing is called *nonmonotonic*.
In more recent versiojns of CVX, we also consider *sign-dependent* 
monotonicity. For example, the functions ``square(x)`` and ``abs(x)`` 
are decreasing for negative :math:`x` and increasing for positive :math:`x`. 
Previous versions of CVX classifed these functions as nonmonotonic, which 
affected their use in compositions; more on this in
:ref:`sign-monotonicity` below.

For functions with multiple arguments, curvature is always considered 
*jointly*, but monotonicity can be considered on an *argument-by-argument* 
basis. For example, the function ``quad_over_lin(x,y)``

.. math:: 

	f_{\text{quad\_over\_lin}}(x,y) = \begin{cases} |x|^2/y & y > 0 \\
                                    +\infty & y\leq 0  \end{cases}

is jointly convex in both :math:`x` and :math:`y` and decreasing
in :math:`y`, and exhibits sign-dependent monotonicity in `x`.

Some functions are convex, concave, or affine only for a *subset* of its 
arguments. For example, the function ``norm(x,p)`` where ``p \geq 1`` is 
convex only in its first argument. Whenever this function is used in a CVX
specification, then, the remaining arguments must be constant, or CVX will 
issue an error message. Such arguments correspond to a function's parameters
in mathematical terminology; *e.g.*,

.. math:: 

	f_p(x):\mathbf{R}^n\rightarrow\mathbf{R}, \quad f_p(x) \triangleq \|x\|_p

So it seems fitting that we should refer to such arguments as *parameters* in
this context as well. Henceforth, whenever we speak of a CVX function as being 
convex, concave, or affine, we will assume that its parameters are known and 
have been given appropriate, constant values.

.. _compositions:

Composition rules
~~~~~~~~~~~~~~~~~

Armed with relevant information about :math:`f` and the classification of the
arguments :math:`\arg_k` according to the rules in :ref:`expressions`, we may 
proceed to classify the full expression. We call the rules that govern these function expressions the *composition rules*.

Perhaps the most basic composition rule in convex anaysis is  that convexity
is closed under composition with an affine mapping. For example, function
``square(x)``---which, as its name implies, computes :math:`f(x)=x^2`---is
convex for real arguments `x`.  So if ``x`` is a real variable of dimension
:math:`n`, ``a`` is a  constant :math:`n`-vector, and ``b`` is a constant, the
expression

::

    square( a' * x + b )

is accepted by CVX, which knows that it is convex. 

The affine composition rule is just one one case in a more
sophisiticated composition ruleset. Here is the complete set:

- The function expression :math:`f(\arg_1,\arg_2,\dots,\arg_n)` is affine
  if :math:`f` is affine and every expression
  :math:`\arg_k` is affine.

- The function expression :math:`f(\arg_1,\arg_2,\dots,\arg_n)` is convex
  if :math:`f` is convex (or affine), and if one of the following is true
  for *every* expression :math:`\arg_k`:

  - :math:`\arg_k` is affine.
  - :math:`\arg_k` is convex, *and* the function is increasing 
    in argument :math:`k`.
  - :math:`\arg_k` is concave, *and* the function is decreasing
    in argument :math:`k`.

- The function expression :math:`f(\arg_1,\arg_2,\dots,\arg_n)` is concave
  if :math:`f` is concave (or affine), and if one of the following is true
  for *every* expression :math:`\arg_k`:

  - :math:`\arg_k` is affine.
  - :math:`\arg_k` is concave, *and* the function is increasing 
    in argument :math:`k`.
  - :math:`\arg_k` is convex, *and* the function is decreasing 
    in argument :math:`k`.

For more background on these composition rules, see 
`Convex Optimization <http://www.stanford.edu/~boyd/cvxbook>`_, Section 3.2.4.

Let us examine some examples. The maximum function is convex and increasing 
in every argument, so it can accept any convex expressions as arguments. For
example, if ``x`` is a vector variable, then

::

    max( abs( x ) )

obeys the "convex/increasing/convex" composition rule, and is therefore 
accepted by CVX, and classified as convex. As another example, consider the 
sum function, which is both convex and concave (since it is affine), and 
increasing in each argument. Therefore the expressions

::

    sum( square( x ) )
    sum( sqrt( x ) )

are recognized as valid in CVX, and classified as convex and concave, 
respectively. The first one follows from the "convex/increasing/convex" rule, 
while the second follows from the "concave/increasing/concave" rule.

Most people who know basic convex analysis like to think of these simpler 
examples in terms of more specific rules: a maximum of convex functions is 
convex, and a sum of convex (concave) functions is convex (concave). But as
you can see, these rules are just *special cases* of the this general 
composition ruleset. In fact, with the exception of scalar quadratic 
expressions, the entire DCP ruleset can be thought of as special cases 
of these rules.

For a more complex example, suppose ``x`` is a vector variable, and ``A``, 
``b``, and ``f`` are constants with appropriate dimensions. CVX recognizes
the expression

::

    sqrt(f'*x) + min(4,1.3-norm(A*x-b))

as concave. Consider the term ``sqrt(f'*x)``. CVX recognizes that ``sqrt`` is
concave and ``f'*x`` is affine, so it concludes that ``sqrt(f'*x)`` is
concave. Now consider the second term ``min(4,1.3-norm(A*x-b))``. CVX
recognizes that ``min`` is concave and increasing, so it can accept concave
arguments. CVX recognizes that ``1.3-norm(A*x-b)`` is concave, since it is the
difference of a constant and a convex function. So CVX concludes that the
second term is also concave. The whole expression is then recognized as
concave, since it is the sum of two concave functions.

For a negative example, we can return to the original expression presented in
the beginnnig of this chapter, ``sqrt( x^2 + 1 )``. Assuming that ``x`` is a
scalar variable, this is the composition of a concave, increasing function
``sqrt`` and a convex expression ``x^2+1``. According to the composition
rules, ``sqrt`` can accept a *concave* argument, not a convex argument, so CVX
rejects it. On the other hand, ``norm([x 1])`` is the composition of a convex
function ``norm`` and an affine expression ``[x 1]``, so CVX can indeed accept
that.

.. _sign-monotonicity:

Sign-dependent monotonicity 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Monotonicity is clearly a critical aspect of the rules for nonlinear 
compositions. Previous versions of CVX enforced these rules in a way that 
occasionally produced some unfortunate consequences. For  example, consider 
the expression

::

    square( square( x ) + 1 )

where ``x`` is a scalar variable. This expression is in fact convex, since
:math:`(x^2+1)^2 = x^4+2x^2+1` is convex. However, previous versions of CVX
used to *reject* this expression, because ``square`` is nonmontonic; and so 
it may not accept a convex argument according to the strictest reading of the
composition rules above. Indeed, the  square of a convex function is not, in
general, convex: for example, :math:`(x^2-1)^2 = x^4-2x^2+1` is not convex.

In practice, this explanation may proved unsatisfying. After all, even though
``square`` is nonmonotonic over the entire real line, the expression
``square(x)+1`` has a range of :math:`[1,+\infty)`. And *over that interval*,
``square`` is increasing. Therefore, one could justifiably claim that the
composition rules are satisfied it this case.

The latest versions of CVX implement a simple but effective approach for
extending the composition rules to cover such cases:  *sign-dependent
monotonicity*. To accomplish this, functions that are positive or negative
over their entire domain are noted as such, so this information can be used 
in the sign analysis described in :ref:`sign` above. Furthermore, each 
functions monotonicity is considered *with respect to the sign of its 
input*. So, for example, ``square`` is increasing for positive inputs, and 
decreasing for negative inputs.

Under this new regime, we can now see how ``square(square(x)+1)`` can be
accepted by CVX. First, CVX knows that ``square`` is nonnegative; and as the
sum of two nonnegative terms, it draws the same conclusion about
``square(x)+1``. Because of this, CVX can conclude that the outer instance to
``square`` is increasing. CVX determines that this expression is the
composition of a convex, increasing function and a convex argument, and it is
accepted by the ruleset.

Clearly, sign-dependent monotonicity, and the simple rule-based sign analysis
performed in CVX, is limited. For example, `entr( x )` defined above is
increasing for :math:`x\geq 1/e` and decreasing for :math:`x\leq 1/e`, but 
CVX does not consider that. But our experience with implementations found in
`CVXPY <https://github.com/cvxgrp/cvxpy)>`_, the `Stanford DCP expression
analyzer <http://dcp.stanford.edu/>`_, and our internal version of CVX 
suggest that this covers nearly all of the cases CVX users are likely to 
encounter.

.. _quadforms:

Scalar quadratic forms
----------------------

In its pure form, the DCP ruleset forbids even the use of simple quadratic 
expressions such as ``x * x`` (assuming ``x`` is a scalar variable). For 
practical reasons, we have chosen to make an exception to the ruleset to 
allow for the recognition of certain specific quadratic forms that map 
directly to certain convex quadratic functions (or their concave negatives)
in the CVX atom library:

=====================   =============================
``x .* x``              ``square( x )`` (real ``x``)
``conj( x ) .* x``      ``square_abs( x )``                
``y' * y``              ``sum_square_abs( y )``            
``(A*x-b)'*Q*(Ax-b)``   ``quad_form( A*x - b, Q )`` 
=====================   =============================

CVX detects the quadratic expressions such as those on the left above, and 
determines whether or not they are convex or concave; and if so, translates 
them to an equivalent function call, such as those on the right above.

CVX examines each *single* product of affine expressions, and each *single* 
squaring of an affine expression, checking for convexity; it will not check,
for example, sums of products of affine expressions. For example, given 
scalar variables ``x`` and ``y``, the expression

::

    x ^ 2 + 2 * x * y + y ^2

will cause an error in CVX, because the second of the three terms
``2 * x * y``, is neither convex nor concave. But the equivalent
expressions

::

    ( x + y ) ^ 2
    ( x + y ) * ( x + y )

will be accepted. 

CVX actually completes the square when it comes across a scalar quadratic
form, so the form need not be symmetric. For example, if ``z`` is a vector 
variable, ``a``, ``b`` are constants, and ``Q`` is positive definite, then

::

    ( z + a )' * Q * ( z + b )

will be recognized as convex. Once a quadratic form has been verified by  
CVX, it can be freely used in any way that a normal convex or concave  
expression can be, as described in :ref:`expressions`.

Quadratic forms should actually be used *less frequently* in disciplined
convex programming than in a more traditional mathematical programming
framework, where a quadratic form is often a smooth substitute for a
nonsmooth form that one truly wishes to use. In CVX, such
substitutions are rarely necessary, because of its support for nonsmooth
functions. For example, the constraint

::

    sum( ( A * x - b ) .^ 2 ) <= 1

is equivalently represented using the Euclidean norm:

::

    norm( A * x - b ) <= 1

With modern solvers, the second form is more naturally represented using
a second-order cone constraint---so the second form may actually be more
efficient. In fact, in our experience, the non-squared form will often
be handled more accurately. So we strongly encourage you to re-evaluate
the use of quadratic forms in your models, in light of the new
capabilities afforded by disciplined convex programming.

.. _strict:

Strict inequalities
-------------------

As mentioned in :ref:`constraints`, strict inequalities ``<``, ``>`` are 
interpreted in an identical fashion to nonstrict inequalities ``>=``, 
``<=``. It is important to note that CVX cannot guarantee that an inequality 
will be strictly satisfied at the solution it computes. This is not simply a 
choice we have made in CVX; it is a natural consequence of both the 
underlying mathematics and the design of convex optimization solvers. For 
that reason, we *strongly* discourage the use of strict inequalities in CVX, 
and a future version may remove them altogether.

When a strict inequality is essential to your model, you may need to take 
additional steps to ensure compliance. In some cases, this can be 
accomplished through *normalization*. For instance, consider a set of 
homogeneous equations and inequalities:

.. math::

  A x = 0, \quad C x \preceq 0, \quad x \succ 0
  
Except for the strict inequality, :math:`x=0` would be an acceptable 
solution; indeed the need to avoid the origin is the very reason for the 
strict inequality. However, note that if a given :math:`x` satisfies these 
constraints, then so does :math:`\alpha x` for all :math:`\alpha>0`. By 
eliminating this degree of freedom with normalization, we can eliminate the 
strict inequality; for instance:

.. math::

  A x = 0, \quad C x \preceq 0, \quad x \succ 0, \quad \mathbf{1}^T x = 1
  
If normalization is not a valid approach for your model, you may simply need 
to convert the strict inequality into a non-strict one by adding a small 
offset; *e.g.*, convert ``x > 0`` to, say, ``x >= 1e-4``. Note that the 
bound needs to be large enough so that the underlying solver considers it 
numerically significant.

Finally, note that for some functions like ``log(x)`` and ``inv_pos(x)``, 
which have domains defined by strict inequalities, the domain restriction is 
handled *by the function itself*. You do not need to add an explicit 
constraint ``x > 0`` to your model to guarantee that the solution is 
positive.

Log convexity
-------------

Given our strong emphasis on adherence to the DCP ruleset, experienced users 
of CVX may be surprised to accidentally stumble upon certain expressions 
involving ``log`` and ``exp`` that violate the ruleset *but are accepted 
anyway*; for example, ``log(exp(x)+1)``. It turns out that this is an 
artifact of CVX's support for :ref:`geometric programming <gp-mode>`; and 
since it also requires the use of CVX's experimental :ref:`successive 
approximation approach <successive>`, it is unsupported. Nevertheless, 
advanced users may be interested in reading more about these "hidden" rules 
in the :ref:`Advanced topics <log-convexity>` chapter.

