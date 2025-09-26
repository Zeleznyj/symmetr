Symmetries
===========

You can also use the program to obtain the symmetries of the system. In the case of relativistic symmetry, this just parses the symmetries from findsym output. For non-relativistic symmetries, the output is, however, non-trivial since these are determined by the program.

The functions to obtain the symmetries are in the :mod:`symmetr.symT` module. These functions take as an input the :class:`symmetr.input.options` object, which is generated using the :func:`symmetr.input.parse` function.

The following functions can be useful:

- :func:`symmetr.symT.get_syms`: this function returns the relativistic (magnetic space group) symmetries
- :func:`symmetr.symT.get_syms_nonmag`: this returns the relativistic (magnetic space group) symmetries of the non-magnetic system
- :func:`symmetr.symT.get_syms_noso`: this returns the non-relativistic symmetries (spin space group)

Example:

.. code:: python

   from symmetr.input import parse
   from symmetr.symT import get_syms

   opt = parse('res j E -f input.in')
   syms = get_syms(opt)

Note that it is always necessary to specify valid program input for the parse function, which means response type has to be specified, even though the type of response selected has no influence on the symmetries themselves.

Symmetry representation
~~~~~~~~~~~~~~~~~~~~~~~~

The symmetries are stored using the :class:`symmetr.symmetry.Symmetry` class. This stores the symmetry in the following form:

- "R" is the real space transformation matrix
- "Rs" is the spin space transformation matrix
- "has_T" a boolean, if `True` the symmetry also contains time-reversal (i.e. is anti-unitary), if `False` the symmetry does not contain time-reversal (i.e. is unitary)
- "permutations" the atom permutations, i.e. how does the symmetry transform the atomic positions

Note that the translation part of the symmetry is not stored since it is not used anywhere in the code, the only important is the atomic permutations.

You can print the symmetry using the python `print` function, which will return for example:

.. code:: none

   R: Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
   Rs: Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
   has_T: False
   permutations: {1: 1}

Coordinate system
~~~~~~~~~~~~~~~~~~~

The symmetries are returned in the conventional coordinate system for the corresponding magnetic space group. This is the coordinate system used by findsym. This is true even in the case of the non-relativistic symmetries!

To transform to the cartesian coordinate system that is used by default for symmetrizing the tensors you can use the :func:`symmetr.symT.get_T` function, which also takes the :class:`symmetr.input.options` object as an input and returns the transformation matrix, which can be directly used to convert the symmetries to the proper coordinate system:

.. code:: python

   from symmetr.symT import get_T
   T = get_T(opt)
   for sym in syms:
      sym.convert(T,in_place=True)
