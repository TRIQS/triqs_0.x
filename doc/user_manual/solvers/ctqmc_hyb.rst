.. _ctqmc_hyb:


Continuous-time quantum Monte Carlo algorithm
==============================================

Introduction
------------------------

The continuous-time quantum Monte Carlo (CTQMC) algorithm that is implemented
in TRIQS is based on a hybridization expansion of the partition function
[#ctqmc1]_, [#ctqmc2]_.  The principle of the algorithm is to sample
stochastically the diagrams of this expansion with the correct Monte Carlo
weights and to compute the Green's function. A Monte Carlo configuration
:math:`\mathcal{C}` is a set of fermionic operators (in interaction
representation) at different imaginary times:

.. math::

  \mathcal{C} = d^\dagger_{\alpha_1}(\tau_1) d_{\alpha'_1}(\tau'_1) d^\dagger_{\alpha_2}(\tau_2)
                d^\dagger_{\alpha_3}(\tau_3) \ldots d_{\alpha}(\tau_N)

The algorithm samples new configurations by inserting/removing pairs of
operators, or by moving operators in the configuration. Note that it is a
finite-temperature algorithm so that the :math:`\tau \in [0,\beta]`, where
:math:`\beta` is the inverse temperature. The Monte Carlo weight of a
configuration is essentially the product of the trace :math:`\mathrm{Tr} \,
\mathcal{C}` and the determinant of a matrix whose elements are the
hybridization functions :math:`\Delta_{\alpha_i \alpha_j'} (\tau_i - \tau_j')`.

The main inputs of the solver are the hybridization functions
:math:`\Delta(i\omega_n)` and the local Hamiltonian
:math:`\mathcal{H}_\mathrm{loc}` on the impurity. The solver then computes the
Green's function on the imaginary-time interval :math:`[0,\beta]`.  This is done on a
basis of Legendre polynomials, as described in Ref. [#legendre]_.  Note that
our implementation of the algorithm uses a *matrix* representation [#ctqmc3]_ of the
operators :math:`d^\dagger_{\alpha}`. This allows to use any local Hamiltonian
:math:`\mathcal{H}_\mathrm{loc}` in the algorithm.

There is a special case, in which the operators are just numbers. This
happens when :math:`\mathcal{H}_\mathrm{loc}` commutes with every number
operator of the problem :math:`[ \mathcal{H}_\mathrm{loc}, n_\alpha]=0, \,
\forall \alpha`. In this situation, only a restricted number of configurations
can exist (those where two construction (destruction) operators don't lie
next to each other). This is called the *segment picture* and the code can
be optimized. There is an input flag to specify if the segment picture should
be used.

In the following, we illustrate how to use the solver with an example. Then we
describe the Operator class, used to define the local Hamiltonian, give
the full reference of the Solver class and finally show how to use
different random number generators.


.. toctree::
   :maxdepth: 1

   aim
   operators
   reference
   settingparameters
   random


Journal references
------------------------

.. [#ctqmc1] P. Werner, A. Comanac, L. de’ Medici, M. Troyer, and
             A. J. Millis, Phys. Rev. Lett. 97, 076405 (2006).
.. [#ctqmc2] E. Gull, A. J. Millis, A. I. Lichtenstein, A. N. Rubtsov,
             M. Troyer, and P.Werner, Rev. Mod. Phys. 83, 349 (2011).
.. [#legendre] L. Boehnke, H. Hafermann, M. Ferrero, F. Lechermann, and O. Parcollet,
               Phys. Rev. B 84, 075145 (2011).
.. [#ctqmc3] P. Werner and A. J. Millis,
             Phys. Rev. B 74, 155107 (2006).


 


