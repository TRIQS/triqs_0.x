The Move concept 
===================

* **Purpose**  : Defines the basic MonteCarlo move.
* **Definition** : 

  ===========================================  =============================================================================================
  Elements                                     Comment
  ===========================================  =============================================================================================
  * mc_weight_type                             - The type of the MC weight : complex or real
  * mc_weight_type Try()                       - First part of the Move.
                                               - Returns the probability to accept the move.  If :
                                                  - the move is :math:`x\rightarrow x'`, proposed with proba :math:`T_{x\rightarrow x'}` 
                                                  - the probability of the :math:`x` config is denoted :math:`p_x`
                                               - then Try should return the Metropolis ratio:
                                                    .. math::
                                                           \dfrac{p_{x'} T_{x'\rightarrow x}}{p_x T_{x\rightarrow x'}} 
   
                                                 with the sign. (The sign will be separated from the absolute value and kept by the MonteCarlo class).
                                                 In other words: Try() = move_sign * move_rate with abs(move_sign) = 1 
  * mc_weight_type Accept()                    - Called iif the Move is accepted.
                                               - Returns a number r such that :math:`|r| =1`
                                               - There is a possibility to update the sign
                                                 here as well, so that the final sign is: move_sign * Accept()
                                                 CAREFUL that we need the ratio new_sign / old_sign here just like
                                                 we needed a ratio in Try()
  * mc_weight_type Reject()                    - Called iif the Move is rejected (for cleaning).
  ===========================================  =============================================================================================


