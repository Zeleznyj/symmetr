Python package for determining symmetry properties of crystals. See the [documentation](https://symmetr.readthedocs.io/) or the [manuscript](https://arxiv.org/abs/2602.21034). The source code can be found at the [GitHub repository](https://github.com/zeleznyj/symmetr).

The code internally uses [findsym](https://stokes.byu.edu/iso/isolinux.php) to determine the symmetry of a given crystal. Since findsym only works on Linux, most features of Symmetr will also only work on Linux. However, it is possible to run the code on any platform through [Docker](https://symmetr.readthedocs.io/en/latest/docker.html). The Windows Subsystem for Linux could likely be used for running the code, but this has not been tested. 

When using the code in scientific publications please cite it as:

- J. Železný, Symmetr: a Python package for determining symmetry properties of crystals, arXiv:2602.21034

You should also cite the findsym package as that is included with the code. Use the following references:

- H. T. Stokes, D. M. Hatch, and B. J. Campbell, FINDSYM, ISOTROPY Software Suite, iso.byu.edu
- H. T. Stokes, D. M. Hatch, “Program for Identifying the Space Group Symmetry of a Crystal”, J. Appl. Cryst. 38, 237-238 (2005).

If you are using MAGNDATA you should also cite it:
- J. Appl. Cryst. (2016). 49, 1750-1776
- J. Appl. Cryst. (2016). 49, 1941-1956


