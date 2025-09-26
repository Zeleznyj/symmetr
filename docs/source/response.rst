Response mode
==============


This is specified using the **res** keyword.

Example:

.. code-block:: none

    symmetry res j E -f input.in

Will determine the symmetry of the conductivity tensor for system specified in 'input.in'.

The response tensor is specified using arguments op1, and op2, which determine the response observable and the perturbing field, respectively.

Currently the following are implemented:

- 'j' or 'v' denote electric current
- 'jq' denotes heat current
- 's' denotes spin 
- this can also be combined, so 's.v' is a spin current
- 'L' denotes orbital momentum
- 't' denotes torque. Note that torque acting on the magnetic moment is perpendicular to the magnetic moment, however this constraint is not
  implemented in the code!
- 'E' denotes electric field (polarization), 'V' denotes voltage
- 'B' denotes magnetic field
- 'gT' denotes thermal gradient
- 'x' denotes strain 
- higher-order response is specified as 'E.E'
- when '0' is used for op2 it specifies equilibrium properties

Note that all operators, except for the '0' can be used in op1 and op2.

Further information about the observable types can be found :doc:`here <observables>`.

The program automatically determines which indices of the tensor are symmetric and in some cases apllies the Onsager relations. In such case
the program prints a message about the symmetrization or application of the Onsager relations. However, this is not neccesarily applied in
every case and thus this can also be specified manually. For this the flags ``--sym-inds``, ``--asym-inds``, ``--T-sym-inds``,
``--T-asym-inds`` are utilized, as described :doc:`here <syminds>`.


The code outputs always two response tensor, corresponding to the time-reversal even and time-reversal odd parts. Note that the order of the two can change depending on the response tensor! The order is always such that the second tensor is the one that corresponds to the intrinsc transformation of the observables under time-reversal, whereas the first has additional minus sign for symmetries with time-reversal. Such response exists in general due to the irreversible nature of dissipative processes. 

Example output for conductivity of the bcc Fe with magnetization along *z*:

.. code-block:: none

    even part of the response tensor:
    ⎡x₀₀   0    0 ⎤
    ⎢             ⎥
    ⎢ 0   x₀₀   0 ⎥
    ⎢             ⎥
    ⎣ 0    0   x₂₂⎦
    odd part of the response tensor:
    ⎡ 0   -x₁₀  0⎤
    ⎢            ⎥
    ⎢x₁₀   0    0⎥
    ⎢            ⎥
    ⎣ 0    0    0⎦
    
The :math:`x_{ij}` variables are indepedent components of the response tensor. These can take any value. Note that the components in the even and odd parts are independent!

Other flags relevant for the response mode:

- ``--noso`` - this turns on the non-relativistic symmetry. By default the symmetry is relativistic.
- ``--exp n`` - Here *n* is an integer. This does expansion in the order paramater in the n-th order. Only works for collinear systems!!!

   For example, the first-order expansion for conductivity in bcc Fe (corresponding to lowest order anomalous Hall effect) is determined by:

   .. code:: none

       symmetr res j E -f Fe.in --exp 1

   and the output is:

   .. code::

        ⎡   0      -m₂⋅x₂₁₀  m₁⋅x₂₁₀ ⎤
        ⎢                            ⎥
        ⎢m₂⋅x₂₁₀      0      -m₀⋅x₂₁₀⎥
        ⎢                            ⎥
        ⎣-m₁⋅x₂₁₀  m₀⋅x₂₁₀      0    ⎦
       
   Note that in this case only one tensor is returned, since each term in the expansion corresponds to either T-even or T-odd component.

- ``-e`` - this will return symmetry of all magnetic configurations related by symmetry. This corresponds to different magnetic domains.
- ``-p a`` - this will turn on projection on a specific atom, specified by integer *a*. What this means specifically is that when an atom is chosen, only symmetries that keep this atom invariant will be considered. This is useful for the spin-orbit torque. For example, to detemrine the effective magnetic field corresponding to the spin-orbit torque acting on site 1:

   .. code:: none

      symmetr res B E -f input.in -p 1

   .. warning::

      The sites are numbered starting from 1!


- ``-p2 a2`` - this can only be used together with -p a. If a symmetry is found that connects atoms *a* and *a2* it will als return the projected tensor transformed to atom 2. Note in this case the independent variables :math:`x_{ij}` for the tensors project on *a* and *a2* are the same.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   observables
   syminds

