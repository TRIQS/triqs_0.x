Setting the parameters
======================

Here's a step-by-step guide that should show you how to prepare a CTQMC run.
Take the time to go through this little guide as it can help you avoid doing
simple mistakes. At first you need to create an instance of the CTQMC solver
class. This is done with::

    from pytriqs.applications.impurity_solvers.operators import *
    from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver

    # Create a solver instance
    S = Solver( Beta = invT, parameters )

The first parameter of the solver is the inverse temperature. In order to set
all the other necessay parameters, you need to go through the following steps:

Step 1 - the structure of the Green's function
----------------------------------------------

You need to figure out is the structure of the Green's function for the problem
you are considering. As emphasized earlier, you should always try to take
advantage of a possible block structure of the Green's function.  In a spin
system, the Green's function can often be cut into *up* and *down* spint
sectors.

When the structure is clear you can set the parameter ``GFstruct`` which is
a list of tuples describing the Green's function blocks. Each of these tuples
has two elements: a string giving the name of the block and a list of the indices
of the block.

Examples
........

* For a single-band Hubbard model with a local Coulomb interaction, the Green's function
  can be cut in two up/down blocks of size 1. We would have::

    GFstruct = [ ('up', [1]), ('down', [1]) ]

* For a two-band Hubbard model with a hybridization between the bands, the Green's function
  can be cut in two up/down blocks, but there are off-diagonal orbital elements. We have::

    GFstruct = [ ('up', [1,2]), ('down', [1,2]) ]

Step 2 - the Hamiltonian
------------------------

The next step is to describe the local Hamiltonian. This is the Hamiltonian
acting on the effective impurity sites/orbitals. It is very important to
**include only the quartic terms and not the quadratic terms** in this
Hamiltonian. The latter terms will be computed from the kownledge of the
non-interacting Green's function ``S.G0`` as explained below (see Step 7). The
Hamiltonian is given in the parameter ``H_Local`` with the use of
:ref:`operators`. The arguments in the parenthesis of the ``C()``, ``Cdag()``
and ``N()`` operators must be compatible with the structure of the Green's
function ``GFstruct``.

Examples
........

* For a single-band Hubbard model with a local Coulom interaction::

    H_Local = U * N('up',1) * N('down',1)

* Two-orbital Hubbard model, no inter-orbital interaction, but a hybridization
  between the bands (this term will not appear in the local Hamiltonian!)::

    H_Local = U * (N('up',1) * N('down',1) + N('up',2) * N('down',2))

Step 3 - the quantum numbers
----------------------------

In order to improve the speed of the solver, it is important to use as many
quantum numbers as possible. These quantum numbers are the eigenvalues of a set
of symmetry *operators* that commute with the Hamiltonian. In other words, the
local eigenstates will be organized so that they are both eigenvectors of the
Hamiltonian and of the quantum number operators.  Note that the symmetries
corresponding to the quantum numbers must be *abelian* symmetries.  The quantum
number operators are given as a dictionary in the parameter
``Quantum_Numbers``.  The key is a string that describes the quantum number and
the value is the corresponding operator (see :ref:`operators`).


Examples
........

* For a single-band Hubbard model with a local Coulom interaction::

    Quantum_Numbers = { 'Number of up spins' : N('up',1), 'Number of down spins' : N('down',1) }

* If we consider a two-band Hubbard model with a fully SU(2) Hund's exchange::

    Ntot = N('up',1) + N('down',1) + N('up',2) + N('down',2)
    Sz = N('up',1) + N('down',2) - N('down',1) - N('down',2)
    Quantum_Numbers = { 'Total N' : Ntot, 'Total spin along z' : Sz }

Step 4 - can I use the segment picture?
---------------------------------------

In the special case where *all* the operators ``N(a,alpha)`` are good quantum
number operators (i.e. they all commute with the Hamiltonian), one can improve
the algorithm and use the so-called *segment picture*. This case is special
because then all construction and destruction operators are just a number
between the local eigenstates. If you are in this situation you should turn the
parameter ``Use_Segment_Picture`` to ``True``. For example, the segment picture
can be used for all Hamiltonians that are purely density-density. On the
contrary the segment picture cannot be used for a two-band Hubbard model with a
full SU(2) Hund's exchange.

Step 5 - the Monte Carlo parameters
-----------------------------------

There are three parameters that control how many steps and measures are done
during the Monte Carlo sampling. ``N_Cycles`` is the total number of measurment
cycles done. A cycle has ``Length_Cycle`` Monte Carlo moves in it.  The
measurments start after ``N_Warmup_Cycles`` cycles have been made and there is
a measurment at the end of every cycle. At the end of the run, there will be
``N_Cycles`` measurments and a total of (``N_Cycles`` + ``N_Warmup_Cycles``) x
``Length_Cycle`` moves.

When the solver is spread on a parallel machine, each core will do ``N_Cycles``
measurments cycles and ``N_Warmup_Cycles`` warmup cycles. Therefore the same
input run on a larger number of cores will yield a larger statistics.

Step 6 - how many Legendre coefficients?
----------------------------------------

The CTQMC algorithm computes the Green's function on the imaginary-time
interval :math:`[0,\beta]`. In order to gain memory and to reduce
high-frequency noise, the Green's function is expanded on a basis of
``N_Legendre_Coeffs`` Legendre polynomials.  The question is, how many of these
polynomials should one use? Our recommendation is to do a first *test* run
with a large number of coefficients, say 80. When the run is over, one
can inspect the Legendre Green's function and decide how many coefficients
should be kept. This will be detailed below.


Step 7 - prepare the non-interacting Green's function
-----------------------------------------------------

The last step before starting the solver is to prepare the non-interacting
Green's function of the problem. From the knowledge of this Green's function,
the solver can extract the hybridization function used in the algorithm and the
quadratic terms of the local Hamiltonian. The non-interacting Green's function
must be initialized in the member ``G0`` of the solver instance. For example,
one would write::

  for spin, g0 in S.G0:
    g0 <<= inverse( iOmega_n - e_f - V**2 * Wilson(D) )

to initialize the Green's function of an impurity imbedded in a flat
conduction bath.

Step 8 - we're ready to go!
---------------------------

Everything is ready at this stage and you just need to call the ``Solve()``
member of the solver::

  S.Solve()

When you call the solver, the local Hamiltonian (with the quadratic terms) is
shown. Be careful to check that this is indeed the Hamiltonian that you expect!
At the end of the run, the solver has computed the following objects:

  * The interacting Green's function of the problem on the Matsubara frequency
    axis. This is in the class member ``G``.

  * The interacting Legendre Green's function of the problem. This is put
    in the member ``G_Legendre``. This output is useful to decide how many
    Legendre coefficients should be used. 

  * The self-energy on the Matsubara frequency axis. This is in the member ``Sigma``.


Final Step - analyze the output
-------------------------------

The final step is to decide on the ideal number of Legendre coefficients to keep
for the following runs. If you have saved the Legendre Green's function in an
archive, you can then plot it:

.. plot:: tutorials/applications/solvers/legendre_plot.py
  :include-source:
  :scale: 70

From this plot you see that for coefficient :math:`l > 30`, the value of the
coefficient is of the order of the statistical noise. There is therefore no
information in the coefficients :math:`l > 30` and one can set
``N_Legendre_Coeffs = 30`` for the following runs. Of course, if you will use a
larger statistics or a larger number of cores, you might have to reajust this
value.
