Usage
=====

For a quick introduction to the code see :doc:`here <quickstart>`.

Input
------

The input to the program can be either:

- crystal structure, inclding the magnetic moments. This is specified in a separate file, which is described :doc:`here </input>`.

  To use this option use

.. code-block:: none

    -f findsym.in

where *findsym.in* is the name of the input file.

- symmetry group, however, this can be used only for limited number of functionalities as many features of the code will not work with this option.


.. code-block:: none

    -g "P4mm"

This has to be one of the magnetic space groups specified either by its name or number as given in the `list <https://bitbucket.org/zeleznyj/linear-response-symmetry/src/master/findsym/mag_groups.txt>`__ of magnetic space groups distributed with findsym. 

- id of a material in the MAGNDATA_ database of magnetic materials. This is specified using the `--magndata` keyword:

.. code:: none

   --magndata "0.222"

.. _MAGNDATA: https://www.cryst.ehu.es/magndata/

This will download the structure from MAGNDATA and store it by default in a file named "magndata_{id}.in", where {id} is the id of the material. This file can directly be used as an input.

When using the MAGNDATA database in publication, cite the corresponding papers ([J. Appl. Cryst. (2016). 49, 1750-1776], [J. Appl. Cryst. (2016). 49, 1941-1956])

One of the -f, -g or --magndata parameters must always be present.

The program has two different modes: :doc:`response mode <response>` used for response tensors (though also equilibrium tensors can also be determined) and a mode for determining the :doc:`symmetry of magnetic interactions <mham>`.

Coordinate system
------------------

The code uses by default a Cartesian coordinate system with a fixed orientation with respect to the lattice. When using Findsym input, this coordinate system is defined with respect to the lattice as defined in the Findsym input. In case the lattice vectors are explicitly given in the input, the coordinate system is simply the Cartesian coordinate system with respect to which the lattice vectors are defined. In case the type 2 lattice input is used, in which only the lattice vector lengths and angles are given, the Cartesian coordinate system is chosen such that the :math:`\hat{z}` axis is oriented along the :math:`c` lattice vector and the :math:`y` axis lies in the :math:`b-c` plane.

It is possible to specify a different coordinate system using the flag :code:`-b` or :code:`--basis`. The options are:

- :code:`abc`: this uses the conventional Bravais lattice corresponding to the magnetic space group.
- :code:`i`: Uses the Bravais lattice given in the Findsym input.
- :code:`abc_c`: Uses the convential Bravais lattice as in :code:`abc`, but normalizes the lattice vectors.

Note that this is included mainly for testing purposes and is rarely needed.


Other flags
------------

:code:`--noso`
~~~~~~~~~~~~~~

This switches on the non-relativistic mode, which utilizes the spin groups.

:code:`--generators`
~~~~~~~~~~~~~~~~~~~~~~

The program is quite fast for small tensors, however, the computation time grows roughly exponentially with the rank of the tensor. To speed
up the computeation it is possible to utilize the flag :code:`--generators`. This will search for the generators of the symmetry group and
utilize only those for the symmetrization since the generators are sufficient to determine the full symmetry. This provides a significant
speedup for larger tensor, however, the scaling remains the same.

:code:`--latex`
~~~~~~~~~~~~~~~

Prints the matrices in a latex format.

:code:`--num-prec`
~~~~~~~~~~~~~~~~~~

This specifies numerical precision. The numerical coefficients of the tensors are rounded to this precision. By default it is :math:`10^{-3}`.




   
.. warning::

   Beside the features and flags discussed here there are many other command line parameters implemented in the code. You can see them using the '-h' flag. Many of these are included only for development or testing purposes and some may be obsolete. Use at your own risk.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   quickstart
   input
   response
   mham

