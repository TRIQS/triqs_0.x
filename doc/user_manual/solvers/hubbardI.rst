.. _hubbardI:


Hubbard I solver
================

Contrary to the hybridization-expansion algorithm, the Hubbard I solver does
not solve the quantum impurity problem exactly. It provides an approximation in
which the hybridization between the impurity and the electronic bath is
neglected.  Hence, solving the impurity problem is reduced to finding the
"atomic" Green's function of an effective atomic Hamiltonian with a correlated
impurity subject to a Coulomb repulsion. The advantage of this solver is that
it is very fast and can provide a reasonable approximation of the Green's
function. It is mainly used in the LDA+DMFT framework.

There, the one-particle terms of this Hamiltonian are given by
:math:`\hat{\epsilon} = -\mu+\langle H^{ff} \rangle - \Sigma_{DC}`, where
:math:`\langle H^{ff} \rangle` is the `LDA` Hamiltonian for the impurity shell
averaged over the Brillouin zone and :math:`\Sigma_{DC}` is the double counting
correction to the self-energy. 

Note that the Hubbard-I approximation is only reasonable for strongly-localized
systems, e.g. local moment rare-earth compounds.  It is very fast and can be
employed both on the Matsubara grid and on the real axis. The TRIQS interface
to the Hubbard-I solver is described in :ref:`HubITRIQS`.  A detailed example
of LDA+DMFT calculations of :math:`\gamma`-Ce within Hubbard-I approximation is
presented in :ref:`LDADMFTtutorial`. 

