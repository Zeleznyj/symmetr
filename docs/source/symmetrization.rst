Symmetrization of tensors
=========================

You can also use the code for any tensor and list of symmetries, provided that you define how the tensor transforms under the symmetries. For this, the function :func:`symmetr.symmetrize.symmetr` function is used.

This takes as an input:

- "syms": a list of symmetries, which can be in principle anything and is not restricted to the :class:`symmetr.symmetry.Symmetry` class
- "X": the symbolic tensor represented by the :class:`symmetr.tensors.Tensor` class
- "trans_func": function that transforms the tensor X by symmetry. It must work like this:
   X_trans = trans_func(X,sym,params)
- "params": parameters to be sent to the trans_func.

