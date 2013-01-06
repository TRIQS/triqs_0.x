.. _solver_concept:


The solver concept
===================

A solver class must obey a specific concept so that inputs don't depend on
which specific solver you are using. The solver class must contain the members
`G` and `G0` corresponding to the interacting and non-interacting Green's
function, and a member function `Solve()` that uses `G0` to compute `G`.
Here's the concept explicitly:

.. class:: GenericSolver

  The generic concept of a solver class

  .. attribute:: G

    The interacting Green's function. This is a :ref:`full Green's function
    <fullgreen>` with each block being a :ref:`GfImFreq`. G is set after
    the Solve() method has been called.

  .. attribute:: G0

    The non-interacting Green's function. This is a :ref:`full Green's function
    <fullgreen>` with each block being a :ref:`GfImFreq`. G0 is needed
    by the solver and should be set before Solve() is called.

  .. function:: Solve()

    This function computes the interacting Green's function of the problem
    using G0 as an input. The solution is stored in G.

