### Program for determining symmetry of linear transport coefficients ###

## How to install ##

Download the findsym code from the [isotropy suite](http://stokes.byu.edu/iso/isolinux.php). This is used for determining the symmetry. See the findsym.txt file for info on the input format for findsym.

Install python and sympy. Sympy is a symbolic library for python. Sympy is in the repositories of most distributions or you can download latest version from [github](https://github.com/sympy/sympy).

Download the latest version of the code using git:

```
#!bash

git clone git@bitbucket.org:zeleznyj/linear-response-symmetry.git
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
#!python

main.py v v < findsym.out

```


