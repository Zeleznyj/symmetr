### Program for determining symmetry of linear transport coefficients ###

## How to install ##

Download the findsym code from the [isotropy suite](http://stokes.byu.edu/iso/isolinux.php). This is used for determining the symmetry. See the findsym.txt file for info on the input format for findsym.

Install python and sympy. Sympy is a symbolic library for python. Sympy is in the repositories of most distributions or you can download the latest version from [github](https://github.com/sympy/sympy). If installing sympy from github, you will probably also need to install [mpmath](https://github.com/fredrik-johansson/mpmath#1-download--installation)

Download the latest version of the code using git:

```
#!bash

git clone https://zeleznyj@bitbucket.org/zeleznyj/linear-response-symmetry.git
```

This will create a folder called linear-response-symmetry.

## How to use ##

Create input file for findsym and run findsym:

```
#!bash

findsym < findsym.in > findsym.out
```

Main executable is called main.py. It takes as input the type of operators, use v v for conductivity tensor, s v for non-equilibrium spin-polarization:

```
#!bash

main.py v v < findsym.out
main.py s v < findsym.out

```

You can also do:

```
#!bash

findsym < findsym.in | main.py v v
```


Non-equilibrium spin-polarization can also be calculated projected on a particular atom, run with s v atom_number. Atom numbers start from 0. For example for atom with number 0 run:

```
#!bash

main.py s v 0 < findsym.out

```

The code outputs symmetrized linear response matrix in three different coordinate systems:

1. Conventional unit cell - this is the unit cell used in crystallographic tables. Vectors a,b,c.
2. Input cell - this is the cell that was used for the findsym input.
3. Orthogonalized cell - this is an orthogonal coordinate system (not a unit cell in general) that has always z in direction of c, y in the plane bc and x orthogonal to z and y. This is an orthogonal coordinate system that is closest to the conventional unit cell. Not properly tested.