Tensors
=======

The tensors are represented using the Tensor class, which is a custom class for storing symbolic tensors, which uses the sympy library under the hood.

 To print the tensor use the ``pprint`` method of the class, for example, for tensor ``X``:

.. code:: python

   X.pprint()

will print the tensor.

To access the individual component of the tensor use normal python indexing, with indices starting from 0, such as:

.. code:: python

   X[0,1]

The components of the tensor are `sympy`_ expression and thus sympy functionalities can be used. For example, to multiply a conductivity tensor by electric field to get the resultant current, you can use:

.. code:: python

   import sympy as sp

   X = res[0]

   j = sp.zeros(1,3)
   Ex,Ey,Ez = sp.symbols('Ex Ey Ez')
   E = sp.Matrix([Ex,Ey,Ez])
   for i in range(3):
       for k in range(3):
           j[i] += X[i,k] * E[k]


.. _sympy: https://www.sympy.org
