Introduction
=============

Symmetr is a python package for determining symmetry properties of crystals. The program takes as an input the crystal structure or the symmetry group and generates the symmetry restricted of various tensors. The main focus is on magnetic systems, however, the code can be used for non-magnetic symmetry. The code is used from the terminal or from within python. Many of the features of the code only work on Linux, although using the code from other OSs.

The capabilities include:

* Support for collinear and non-collinear magnetism
* Symmetry of response tensors such as the anomalous Hall effect, spin-Hall effect or spin-orbit torques, including projection on sublattices
* Higher-order response
* Symmetry of magnetic interactions
* Relativistic symmetry (based on magnetic groups) and non-relativistic symmetry (based on spin-groups)
* Expansions in the order parameter for collinear systems
* Python API
* Determination of symmetry related magnetic configurations

.. note::

   When using the code in scientific publications please cite it as:

   J. Železný, Symmetr, bitbucket.org/zeleznyj/linear-response-symmetry/

   You should also cite the findsym package as that is included with the code. Use the following references:

   H. T. Stokes, D. M. Hatch, and B. J. Campbell, FINDSYM, ISOTROPY Software Suite, iso.byu.edu.
   H. T. Stokes and D. M. Hatch, "Program for Identifying the Space Group Symmetry of a Crystal", J. Appl. Cryst. 38, 237-238 (2005).

   If you are using MAGNDATA you should also cite it: [J. Appl. Cryst. (2016). 49, 1750-1776], [J. Appl. Cryst. (2016). 49, 1941-1956].



.. warning::

   The features of the program described in this documentation have been tested extensively, nevertheless bugs can occur or mistakes can happen when using the program incorrectly. Use caution and if you need help with using the code, contact me at jakub.zelezny@gmail.com.

