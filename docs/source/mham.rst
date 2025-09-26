Magnetic interactions
=======================

This mode is used to determine the symmetry of semiclassical magnetic Hamiltonians, that is Hamiltonians of the type:

.. math::

   H = \sum_{nm} H^{(2),nm}_{ij} M^n_i M^m_j + \sum_{nmpq} H^{(4),nmpq}_{ijkl} M^n_i M^m_j M^n_k M^m_l + \dots 

Here :math:`M^n_i` is the i-th component of the magnetic moment on site :math:`n`. :math:`H^{2}` and :math:`H^{(n)}` are the expansion coefficients that correspond to the magnetic interactions.

.. note::
   Typically it is assumed that the magnitude of the magnetic moments is constant. This is, however, not in any way enforced or taken into account by the code. In practice this means that the expansion is done in terms of variables that are not independent of each other. A consequence of this is that when transformed to independent variables such as spherical coordinates of each magnetic moment, n-th order term in magnetic moment will correspond to n-th order term in spherical coordinates, but may also contain lower order terms. This lower order terms correspond to what is directly obtained by symmetry for lower order terms, so this is not a problem.

To utilize this mode, it is necessary to specify the specific expansion term by specifying the sites it correspond to using the :code:`--sites` keyword. For example:

.. code:: none

   symmetr mham -f input.in --sites 1,2

will determine the interaction between sites 1 and 2, i.e. the term :math:`H^{(2),1 2}`.

The interaction can also be specified for sites in different unit cells. In general, the symmetry of the interaction does not depend on the cell the atoms are in, with one key exception. When considering interaction for the same atom, such as :code:`--sites 1,1`, the resulting tensor has to be symmetric under the interchanging of indices :math:`i,j`. However, when the sites correspon to the same atom but in different unit cells, the tensor does not have to be symmetric. To distinguish the two cases, the atom number can be specified with an apostrophe ': such as :code:`--sites "1,1'"`. Note that in this case you will most likely wrap the whole argument in quotation marks. It is also possible to specify several apostrophes to distinguish several different sites. 

Example
---------


For example, for bcc Fe, this is returned by the code when specifying sites 1,1:

.. code:: none

   ⎡x₀₀   0    0 ⎤
   ⎢             ⎥
   ⎢ 0   x₀₀   0 ⎥
   ⎢             ⎥
   ⎣ 0    0   x₀₀⎦

      2           2           2    
   M₁ₓ ⋅x₀₀ + M_1y ⋅x₀₀ + M_1z ⋅x₀₀

Here the matrix is only printed in the code of 2-order terms and is directly the :math:`H^{(2)}` matrix. The expression below the matrix is always printed and correspond directly to the full expansion term with the magnetic moments directly included. As for the response mode, the :math:`x_{ij}` variables are independent variables that can take any value.

In this case the interaction correspond to the 2-nd order anisotropy. However, because M is normalized, in this case the expression becomes M-independent. This is a correct result: since this is cubic crystal, the second order anisotropy vanishes. If we instead specify sites 1,1', the result becomes:

.. code:: none

   ⎡x₀₀   0    0 ⎤
   ⎢             ⎥
   ⎢ 0   x₀₀   0 ⎥
   ⎢             ⎥
   ⎣ 0    0   x₀₀⎦

   M_1'x⋅M₁ₓ⋅x₀₀ + M_1'y⋅M_1y⋅x₀₀ + M_1'z⋅M_1z⋅x₀₀

In this case the matrix is the same, however, since now the the magnetic moments are different this expression does not become constant. It corresponds to the isotropic Heisenberg exchange:

.. math::

   H^{(2)} = x_{00} \mathbf{M}^1 \cdot \mathbf{M}^2,

where :math:`x_{00}` is the exchange interaction parameter.



