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

If you want to later update the code, run:

```
#!bash

git pull origin master
```
from the folder where the code is located.

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


Non-equilibrium spin-polarization can also be calculated projected on a particular atom, specify atom number with a switch -p or --projection. Atom numbers start from 0. For example for atom with number 0 run:

```
#!bash

main.py s v -p 0 < findsym.out

```
Right now, the atoms are numbered according to the the findsym output - first atom in the findsym output will have number 0 etc. This will probably be changed in the future. 

The code can output the symmetrized linear response matrix in three different coordinate systems. Specify the basis using switch -b or --basis. Separate by commas (with no spaces!) for more than one basis. If no basis is specified, the input basis is used.

1. 'abc' - Conventional unit cell - this is the unit cell used in crystallographic tables. Vectors a,b,c.
2. 'i' - Input cell - this is the cell that was used for the findsym input.
3. 'abc_o' - Orthogonalized cell - this is an orthogonal coordinate system (not a unit cell in general) that has always z in direction of c, y in the plane bc and x orthogonal to z and y. This is an orthogonal coordinate system that is closest to the conventional unit cell. Not properly tested.

 For example:

```
#!bash

main.py s v -b i,abc < findsym.out

```

### Equivalent magnetic configurations ###

The code can output the form of the tensor for all equivalent magnetic configurations (ie configurations related by some symmetry operation). To do that, the symmetry operations of the nonmagnetic crystal are needed. Run findsym with the same input file, but change all the moments to zero. Do not use the nonmagnetic mode of findsym! Then specify the name of the findsym output using switch -e or --equiv:

```
#!bash

main.py s v -e nonmag_out < findsym.out

```

The output is always in the input basis!

### Relation between two sublaticess ###

The code can try to find a relation between the tensors projected on two sublaticess. This will work if there is a symmetry operation relating the two. Use switch -p2 or --projection2:

```
#!bash

main.py s v -p 0 -p 1 < findsym.out

```
The basis transformations work for this.

### Tests ###

There are tests in a folder tests/, you can see how they were created in a script run_tests.sh