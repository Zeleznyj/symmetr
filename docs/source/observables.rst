Observables
============

The key thing that matters for symmetry is how observables transform under symmetry operations. When dealing with relativistic symmetry
there are 4 possibilities, depending on transformation under time-reversal (T) and spatial inversion (P). The following table gives example
of observables for each combination.

.. list-table:: 
   :header-rows: 1
   :stub-columns: 1

   * - 
     - P-even
     - P-odd
   * - T-even
     - t
     - r
   * - T-odd
     - L
     - v

Here 't' is torque, 'r' is position, 'L' is the orbital momentum and 'v' is velocity.

In absence of spin-orbit coupling, the real and spin space are decoupled and the spin then transforms differently. Then a fifth type of
transformation is necessary, which is denoted by 's'. Using these 5 observable types any tensor can be in principle specified. 

The type of the observable used does not matter as long as it transform the same. The only exception to this is that the program in some
cases applies the Onsager relations and symmetrizes indices corresponding to same observables. For example, using ``j E.E``, the last two
indices of the response tensor will be symmetrized, however, if ``j E.gT`` is used the last two indices are not symmetrized, even though
'gT' and 'E' transform the same. This is described in detail :doc:`here <syminds>`.





