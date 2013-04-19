.. index:: Quantum impurity solvers & DMFT

.. module:: pytriqs.applications.impurity_solvers

.. _solvers:

Quantum impurity solvers
========================

TRIQS has a set of extension modules that allow to solve the generic problem of
a quantum "impurity" embedded in a conduction bath. The "impurity" can be any
set of orbitals, on one or several atoms. To be more specific, the Hamiltonian
of the problem has the form:

.. math::

  \mathcal{H} = \sum_{k,\alpha} \epsilon_{k,\alpha} c^\dagger_{k,\alpha} c_{k,\alpha} + \sum_{k,\alpha}
                V_{k,\alpha} ( c^\dagger_{k,\alpha} d_{\alpha} + h.c. ) +
                \mu \sum_\alpha d^\dagger_\alpha d_\alpha +
                \sum{} U_{1234} d^\dagger_1 d^\dagger_2 d_3 d_4.

Here the operators :math:`c^\dagger` construct a fermion in the bath, while
the operators :math:`d^\dagger` construct a fermion on the impurity.
In this problem, the hybridization function :math:`\Delta` between the bath
and the impurity is given by:

.. math::

  \Delta_{\alpha,\beta} (i \omega_n) = \sum_k \frac{V_{k,\alpha} V^*_{k,\beta}}{i \omega_n + \epsilon_{k,\alpha}},

so that the non-interacting Green's function of the impurity is:

.. math::

  \hat{G}^{-1}_0 (i \omega_n) = i \omega_n + \mu - \hat{\Delta}(i \omega_n).

With the knowledge of :math:`G_0` and the matrix :math:`U_{1234}`, the quantum
impurity solvers find the interacting Green's function :math:`G` of the
problem. TRIQS has two such solvers available and they are described below.
Both these solvers have a similar concept, so that it is very simple to switch
from one solver to the other in an input file. This generic concept is
described below as well.

.. toctree::
   :maxdepth: 1

   ctqmc_hyb
   hubbardI
   concept
   dmft/dmft
   ipt/ipt
