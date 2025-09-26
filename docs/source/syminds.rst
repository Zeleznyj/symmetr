Symmetric indices and Onsager relations
========================================

Symmetric and anti-symmetric indices
------------------------------------

The program automatically symmetrizes indices of the tensor in two cases:

- When op1 or op2 contains same observables. So for example, if ``j E.E`` is specified then the last two indices are symmetrizes. Similarly,
  if ``s.s E`` was specified, then the first two indices are symmetrized. This does not apply for the case when op1 and op2 are the same,
  as this is described by the :ref:`Onsager relations <Onsager_relations>`, discussed in the following. 
  
  When such symmetrization is applied, the program inputs a message with the indices that are symmetrized. This can be overriden using the
  flags ``--sym-inds``, ``--asym-inds``, ``--T-sym-inds`` and ``--T-asym-inds``. When one of these flags is used the program does not apply
  symmetrization of indices in this case.
- In the case of the expansion in the order parameter, the last indices of the tensor correspond to the order parameter and thus the tensor
  must be symmetric in this case. This cannot be overriden.

Specify the symmetric or anti-symmetric indices like this:

.. code::

   --sym-inds 1,2:3,4

This means that the tensor will be symmetric under intechanging indices 1 and 2 as well as 3 and 4.

.. warning::

   The indices are numbered starting from 1!

This works similarly for ``--asym-inds`` except the tensor will be anti-symmetric.

The flags ``--T-sym-inds`` and ``--T-asym-inds`` specify indices which are symmetric for the T-even component of the tensor and
anti-symmetric for the T-odd component or vice versa. Specifically, the ``--T-sym-inds`` will make the tensor symmetric for the first
component and anti-symmetric for the second, and vice versa for the ``--T-asym-inds``.

.. _Onsager_relations:
Onsager relations
-----------------

The Onsager reciprocal relations relate an out-of-equilibrum process with an inverse process where the force and flow are reversed. Thus for
example, the consequence of Onsager relations is that the Peltier effect (thermal gradient caused by a voltage) and the Seebeck effect
(electric current caused by thermal gradient) are equivalent. This is not implemented in the code, except in the case when the force and the
flow are of the same type in which case the reverse process is the same and the Onsger relations give then constrains on the response tensor
itself.

For example, in the case of conductivity tensor, the consequence of Onsager relations is that the T-even component of the tensor is
symmetric and the T-odd is anti-symmetric. The same applies for the thermal conduction tensor. At this point these are the only two cases
where the Onsager relations are automatically applied.

In other cases this can be specified manualy using the flags ``--T-sym-inds`` and ``--T-asym-inds``. The ``--T-sym-inds`` flag specifies
indices which are symmetric for the first part of the response tensor (T-even component in the case of conductivity) and anti-symmetric for the second part (T-odd component for conductivity tensor). The ``--T-asym-inds`` flag is the opposite.


