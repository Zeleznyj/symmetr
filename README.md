# Symmetr

Symmetr is a program for determining various symmetry properties of
magnetic and nonmagnetic crystals.

Full documentation for the program is available
[here](https://symmetr.readthedocs.io). The physics behind the code and
the code structure is described in detail in the
[manuscript](https://arxiv.org/abs/2602.21034).

## Installation

The code only works on Linux and requires Python3.

To install the code use pip

```bash
  pip install symmetr
```

or download the repository and then run

```bash
pip install .
```

from the repository directory.

> [!NOTE]
> When using the code in scientific publications please cite it as:
>
> J. Železný, Symmetr: a Python package for determining symmetry
>     properties of crystals, arXiv:2602.21034
>
> You should also cite the findsym package as that is included with the
> code. Use the following references:
>
>H. T. Stokes, D. M. Hatch, and B. J. Campbell, FINDSYM, ISOTROPY Software Suite, iso.byu.edu. H. T. Stokes
> D. M. Hatch, “Program for Identifying the Space Group Symmetry of a Crystal”, J. Appl. Cryst. 38, 237-238 (2005).
>
> If you are using MAGNDATA you should also cite it: \[J. Appl. Cryst.
> (2016). 49, 1750-1776\], \[J. Appl. Cryst. (2016). 49, 1941-1956\].

## Examples

Examples of using Symmetr are available in the examples directory.

## Tests

Tests are included in the tests directory. To run them use pytest

``` bash
cd tests
pytest
```

