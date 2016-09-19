### Program for determining symmetry of linear transport coefficients ###

## How to install ##
The program uses the findsym code from the [isotropy suite](http://stokes.byu.edu/iso/isolinux.php) for determining the symmetry. Findsym is now included in the repository so it does not have to be installed separately. See the findsym.txt file for info on the input format for findsym.

Install python, numpy and sympy. Sympy is a symbolic library for python. Sympy is in the repositories of most distributions or you can download the latest version from [github](https://github.com/sympy/sympy). If installing sympy from github, you will probably also need to install [mpmath](https://github.com/fredrik-johansson/mpmath#1-download--installation)

Download the latest version of the code using git:

```
#!bash

git clone https://zeleznyj@bitbucket.org/zeleznyj/linear-response-symmetry.git
```

This will create a folder called linear-response-symmetry.

If you want to later update the code, run:

```
#!bash

git pull origin master
```
from the folder where the code is located.

## How to use ##



Main executable is called main.py. It takes as input the type of operators, use v v for conductivity tensor, s v for non-equilibrium spin-polarization, specify a findsym input file with switch -f:

```
#!bash

main.py v v -f findsym.in
main.py s v -f findsym.in

```

You can also use a group name or number as an input:
```
#!bash

main.py v v -g Pmmm

```

Escape primes with a backslash:
```
#!bash

main.py v v -g Pmmm

```

Non-equilibrium spin-polarization can also be calculated projected on a particular atom, specify atom number with a switch -p or --projection. For example for atom with number 1 run:

```
#!bash

main.py s v -p 1 -f findsym.in

```
The atoms are numbered by order in the findsym input file, starting from 1. Note that it may happen that some of the atoms given in the input file are equivalent, then only the first can be used. List of all non-equivalent atoms is in the findsym output, where the positions are printed in the new (abc) basis. The first number is the atom number that the program uses.

The code can output the symmetrized linear response matrix in three different coordinate systems. Specify the basis using switch -b or --basis. Separate by commas (with no spaces!) for more than one basis. If no basis is specified, the input basis is used.

1. 'abc' - Conventional unit cell - this is the unit cell used in crystallographic tables. Vectors a,b,c.
2. 'i' - Input cell - this is the cell that was used for the findsym input.
3. 'cart' - A cartesian coordinate system - this is the cartesian coordinate system with respect to which the unit cell is defined in the findsym input.
4. 'abc_o' - Orthogonalized cell - this is an orthogonal coordinate system (not a unit cell in general) that has always z in direction of c, y in the plane bc and x orthogonal to z and y. This is an orthogonal coordinate system that is closest to the conventional unit cell. Not properly tested.

 For example:

```
#!bash

main.py s v -b cart -f findsym.out

```

### Equivalent magnetic configurations ###

The code can output the form of the tensor for all equivalent magnetic configurations (ie configurations related by some symmetry operation). To do that use the switch -e or --equiv

```
#!bash

main.py s v -e -f findsym.in

```

### Relation between two sublaticess ###

The code can try to find a relation between the tensors projected on two sublaticess. This will work if there is a symmetry operation relating the two. Use switch -p2 or --projection2:

```
#!bash

main.py s v -p 1 -p2 2 -f findsym.in

```

###Latex output###
The output can be also given in latex format, use --latex switch.

### Tests ###

There are tests in a folder tests/, you can see how they were created in a script run_tests.sh