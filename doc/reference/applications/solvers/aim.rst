.. index:: Anderson impurity model

.. _aim:


An example: the Anderson impurity model
========================================

To illustrate how the CTQMC solver works in practice, we show the example of a
one-orbital Anderson impurity embedded in a flat (Wilson) conduction bath. The
local Hamiltonian of this problem is simply:

.. math::

  \mathcal{H}_\mathrm{loc} = U n_\uparrow n_\downarrow,

and the non-interacting Green's function is:

.. math::

  G^{-1}_{0,\sigma} (i \omega_n) = i \omega_n - \epsilon_f - V^2 \Gamma_\sigma(i \omega_n).

In this example, there is a Coulomb interaction :math:`U` on the impurity level
which is at an energy :math:`e_f`. The bath Green's function is :math:`\Gamma(i
\omega_n)` and it has a flat density of states over the interval
:math:`[-1,1]`.  Finally, :math:`V` is the hybridization amplitude between the
impurity and the bath. Let us solve this problem with the CTQMC solver. Here is
the python :download:`script <aim.py>`:

.. literalinclude:: aim.py

Running this script on a single processor takes about 5 minutes and generates
an HDF5 archive file called :file:`solution.h5`. This file contains the
Green's function found by the solver as well as the computed average for the
occupation of the impurity level (the average of the operator :math:`n_\uparrow
+ n_\downarrow`). Let us plot the Green's function and print the value of the
average:

.. plot:: reference/applications/solvers/aim_plot.py
   :include-source:
   :scale: 70

As expected the result shows a particle-hole symmetric impurity Green's
function (the real part vanishes up to the statistical noise).

Let us now go through the script in some more details.

.. literalinclude:: aim.py
  :lines: 1-3

These lines import the classes to manipulate Green's functions, fermionic
operators, and the CTQMC solver.

.. literalinclude:: aim.py
  :lines: 5-6

This just sets the parameters of the problem.

.. literalinclude:: aim.py
  :lines: 8-25

This is the construction of the solver object. The class is described
in more detail in the section: :ref:`ctqmc_ref`. Basically, the constructor
of the solver needs a set of keyword arguments, essentially:

  - ``Beta``:      The inverse temperature
  - ``GFstruct``:  A tuple of lists that describe the block structure of the Green's function
  - ``H_Local``:      The local Hamiltonian written with operators. See the section: :ref:`operators`
  - ``Quantum_Numbers``:   A dict of operators describing the quantum numbers
  - ``N_Cycles``: The number of Monte Carlo cycles
  - ``Length_Cycle``: The number Monte Carlo moves in a cycle
  - ``Random_Generator_Name``: The name of the random number generator
  - ``Use_Segment_Picture``: Use the segment picture
  - ``Measured_Operators``: A dictionary of operator that will be averaged
  - ``Global_Moves``: A tuple of lists that describe the global moves to be used

After the solver is constructed it needs to know what the non-interacting Green's function
of the impurity is. From this information, the solver will deduce the hybridization function
which is used in the algorithm. The non-interacting Green's function must be put in the
class member ``S.G0``:

.. literalinclude:: aim.py
  :lines: 27-29

At this stage, everything is ready for the solver and we just run it calling its member
function ``Solve``:

.. literalinclude:: aim.py
  :lines: 31-32

When the solver has finished, it puts the result for the interacting Green's function in its
member ``S.G``. If the user asked for an operator to be averaged, the average will be put
in ``S.Measured_Operators_Results[operatorname]``. The last lines of the script save the Green's
function and the average of :math:`n_\uparrow+n_\downarrow` in the HDF archive.
