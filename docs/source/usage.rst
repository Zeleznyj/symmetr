Usage
=====

For a quick introduction to the code see :doc:`here <quickstart>`.

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

   
.. warning::

   Beside the features and flags discussed here there are many other command line parameters implemented in the code. You can see them using the '-h' flag. Many of these are included only for development or testing purposes and some may be obsolete. Use at your own risk.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   quickstart
   input
   response
   mham

