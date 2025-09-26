Quickstart
===========

For example, to determine the symmetry of the conductivity tensor effect of bcc Fe you need an input file that specifies its crystal structure:

.. code:: none

   Fe
   0.001
   2
   3 3 3 90 90 90
   2
   I
   1
   Fe
   magnetic
   0 0 0 0 0 1

Save this as Fe.in. The input file format is described in detail :doc:`here <input>`.

.. code:: none

   symmetr res j E -f Fe.in

This specifies that we want the response tensor corresponding to the current induced by electric field, i.e. the conductivity tensor.

This returns:

.. code:: none

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

Here the even part refers to the T-even part of the conductivity tensor, which is the normal conductivity. The T-odd part refers to the anomalous Hall effect.

The :math:`x_{ij}` variables are independent components of the tensor. 

If we are instead interested in the current due to a thermal gradient, we specify:

.. code:: none

   symmetr res j gT -f Fe.in

The output in this case is the same.

The observable types are described in detail :doc:`here <observables>`.
