Introduction & Motivations
**************************************

* **Purpose** : 

  The purpose of this little class is too facilitate the writing and maintenance
  of Metropolis MonteCarlo algorithms, with several **advantages** :

   * Easy to add new moves and measures.
   * For complex MC algorithms, with 6 or 7 different moves, it helps to clearly separate the code for each move.
   * Parallelism is automatic.
   * The random generator is a simple parameter, it can be chosen dynamically. 

* **Principle**

  The `mc_generic` class is a generic version of the algorithms, with `moves` and `measures`.
  The user  : 
    
    - writes move classes, modelling the Move concept.
    - writes measure classes, modelling the Measure concepts.
    - register them in a mc_generic object MC
    - call `MC.run(...)`  ... and that is (almost) it  !

