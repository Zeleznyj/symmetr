Symmetr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Symmetr is a program for determining various symmetry properties of magnetic and nonmagnetic crystals. Currently the program can determine symmetry of various linear response properties and symmetry of classical magnetic interactions.

How to install
--------------

The easiest way how to install the program is using pip, simply run:

::
    
    pip install symmetr
    
    
The program currently only works on linux. If you use a different OS, you can run the program using docker, see `here <https://bitbucket.org/zeleznyj/linear-response-symmetry/wiki/docker>`__ . Although the program itself could be run on any operating system since it's written in python, most of the features rely on the findsym code from the `isotropy suite <http://stokes.byu.edu/iso/isolinux.php>`__, which only works on linux. Findsym is now included in the repository so it does not have
to be installed separately. See the symmetr/findsym/findsym.txt file for info on the input format for findsym or check the various input files present in tests/.

The program currently works both in Python2 and Python3

If you don't want to use pip, you can directly download the repository as a zip file from `this link <https://bitbucket.org/zeleznyj/linear-response-symmetry/get/7d569d9b5dbe.zip>`__ or use git:

::

    git clone https://zeleznyj@bitbucket.org/zeleznyj/linear-response-symmetry.git 

This will create a folder called linear-response-symmetry. To install you can then run:

::

    python setup.py

How to use
----------

The program is used from the command line. To run the code use the executable symmetr located in the exec directory. Alternatively you can use

::

    python $install_dir/symmetr

where $install_dir is the path to the directory where the code is located.

There are two modes, one for evaluating the linear response symmetry and the other for the symmetry of magnetic hamiltonians, which are used by:

::

    symmetr res ...
    symmetr mham ...

The syntax of the program is different for the two modes, but some arguments are shared. The main ones are the following:

Crystal structure input
^^^^^^^^^^^^^^^^^^^^^^^^

Crystal structure is specified using a findsym input file or directly by specifying the magnetic space group. Note, however, that many features only work when findsym input is used because the magnetic space group does not contain information about the sublatticess and magnetic structure. 

The findsym input file is specified with:

::

    -f findsym.in

where *findsym.in* is the name of the findsym input file.

See the `findsym/findsym.txt <https://bitbucket.org/zeleznyj/linear-response-symmetry/src/master/findsym/findsym.txt>`__ file for description of the findsym input format. Examples can be founds in the tests/ directory.

Note that even if you have nonmagnetic system, you have to define the system as magnetic with zero magnetic moments, otherwise the program will not work!

The magnetic space group is specified by the -g switch, for example:

::

    -g "P4mm"

This has to be one of the magnetic space groups specified either by its name or number as given in the `list <https://bitbucket.org/zeleznyj/linear-response-symmetry/src/master/findsym/mag_groups.txt>`__ of magnetic space groups distributed with findsym. The group name input can only be used in the response mode and only in some cases. 

Either the -f switch or -g switch have to always be present!

Coordinate system
^^^^^^^^^^^^^^^^^^

The program can return the result in different bases (coordinate systems). If findsym input is used then the following bases can be used.

1. **i** - Input cell - this is the conventional basis that was used for the findsym
   input, i.e. this is the basis defined on line 4 in the findsym input file in which the atomic positions are specified. 
#. **cart** - A cartesian coordinate system - this is the cartesian
   coordinate system with respect to which the conventional basis is defined in
   the findsym input. That is, if the conventional basis used in findsym input is given by three lattice vectors, then the cartesian basis vectors are simply: (1,0,0), (0,1,0), (0,0,1). If the conventional basis is given using the form 2 as a,b,c,alpha,beta,gamma, then the conventional basis vectors are defined only up to a rotation. We choose a convention such that the lattice vector c is along the (001) direction and the lattice vector b lies in the [100] plane. 
#. **abc** - Conventional crystallographic unit cell - this is the unit cell used in
   crystallographic tables. This is the basis in which findsym specifies the symmetry operations and is defined by vectors a,b,c in the findsym output.
#. **custom** - You can define a custom coordinate system in the findsym input file. Add a keyword *axes:* at the end of the findsym input file and define the three basis vectors each on separate line. These are given in the **cart** coordinate system.
#. **abc\_o** - Orthogonalized cell - this is an orthogonal coordinate
   system (not a unit cell in general) that has always z in direction of
   c, y in the plane bc and x orthogonal to z and y. This is an
   orthogonal coordinate system that is closest to the conventional crystallographic unit
   cell. Not properly tested and may be removed in future releases.

Specify basis by the -b switch, for example:

::

    -b abc

By default the 'cart' basis is used.

Response mode
---------------

In this mode the program finds the symmetry of response tensors. So far only linear response is implemented, but higher order response is planned.

The program returns tensor :math:`\chi` which describes linear response of some observable :math:`A` to field :math:`F`:

:math:`\delta A_i = \chi_{ij}F_j`.

The observable A can also be a tensor instead of vector and then the tensor :math:`\chi` has three components.

To run the code you have to specify the observable and the field, for example:

::

    symmetr res j E -f findsym.in

for response of current to electric field, i.e. the conductivity tensor.

So far the implemented **observables** are:

- *j* or *v* for current
- *s* for spin
- *t* for torque: Note that for the torque it is not taken into account that it must only exist for magnetic systems and must be perpendicular to the magnetic moments.
- *x* for position: I'm not sure whether this is actually useful for anything
- s.v for the spin current. You can try other combinations of operators, but I'm not sure if you get anything meaningful.

The implemented **fields**:

- *E* or *v* for electric field

The linear response has two fundamentally distinct components, which differ in their transformation under the time-reversal symmetry, i.e. one is even and the other is odd. The program returns the response tensors for both the even and the odd components. 

Note that the order in which the two components is returned is not fixed. This is because the program always returns first the component which corresponds to the Kubo formula component that includes a real part of the matrix elements and second the component which contains the imaginary part of the Kubo formula. Depending of the type of the response function these can be either even or odd, but always one is even and the other odd.

Additional features of the response mode:

Projection on atomic site
^^^^^^^^^^^^^^^^^^^^^^^^^^

Using a switch -p, atomic number can be specified. Then instead of the space group of the crystal, only the site point group is used, i.e. the group of symmetry operations of the whole crystal which leave the atomic site invariant. This is useful for evaluating response of local quantities such as the spin-polarization because then the response corresonds to response of spin on the chosen atomic site. I'm not sure if this would give anything meaningful for the conductivity tensor for example.

For example for projecting on site 2:

::

    symmetr res s E -f findsym.in -p 2 

The atomic sites are numbered by the order in which they appear in the findsym input file. Note that it can happen that the conventional unit cell you define in the findsym input file is larger than the crystallographic conventional cell. Then findsym uses a smaller number of atoms than in the input file. The numbering of the atoms in this case is still the same, but some of the atoms cannot be projected on. You can print the atomic sites used by findsym by

::

    --print-pos

This prints for each atomic site the index of the atom and its position in the crystallographic conventional basis used by findsym.

The projections are only possible with the findsym input as when the group input is used the information about lattice sites is not available.

If a projection on atomic site is done it is also possible to specify projection on a second site by -p2, for example:

::

    symmetr res s E -f findsym.in -p 2 -p2 3

This will try to find a symmetr operation which connects sites 2 and 3 and it if finds one, it will print out the relation between response tensors on sites 2 and 3. Note that this feature is not implemented for the three component response function (when A = s.j)

Equivalent magnetic configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a magnetic system it is often useful to know how are the response tensors related for all magnetic configuration which are connected by some symmetry operation (and are therefore equivalent). This can be done using a switch -e:


::

    symmetr s v -f findsym.in -e

This is only possible when the findsym input is used.

Expansions in the direction of magnetic moments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a magnetic system it is often useful to know how does the response tensor change when the direction of magnetic moments is rotated. For a collinear magnetic system:

:math:`\chi_{ij}(\hat{\mathbf{n}}) = \chi_{ij}^{(0)} + \chi_{ij,k}^{(1)} \hat{n}_k + \chi_{ij,kl}^{(2)} \hat{n}_k \hat{n}_l + \dots`

where :math:`\hat{n}` is the direction of the magnetic moments.

Symmetry of the expansion coefficients can be obtained using

::

    --exp n

where *n* is the order of the expansion.

This features is at the moment only implemented for ferromagnets. When you use findsym input file, the magnetic order is ignored and the output is always for a ferromagnetic system. In a two-sublattice collinear antiferromagnet, if a projection on one of the magnetic atoms is specified, then this expansion also applies. In general, however, for a collinear antiferromagnet the expansion is different from a ferromagnet and this is not implemented!

You can also specify a group as an input, however, you need to specify the correct nonmagnetic point group as this is not done automatically and is not checked by the program. That is, you need to find the point group of the nonmagnetic crystal since this is the group the determines the symmetry of the expansion.

Expansions are not implemented for the three component response.


Magnetic Hamiltonians mode
---------------------------

In this mode the program determines the symmetry of magnetic Hamiltonians, that is Hamiltonians of the form:

:math:`H = \sum_{ij} H_{ij}^{ab} M_i^a M_j^b`

and analogously for interaction of more magnetic momentes. Here :math:`a,b` denote different magnetic sublattices. 

To determine the symmetry of tensors :math:`H_{ij\dots}^{ab\dots}` use the *mham* mode. You need to use the findsym input with this mode. The magnetic order you specify in the input file is irrelevant. 

As an input you need to specify the sites for which you want to find  the interaction term using the -s switch:

::

    symmetr mham -s 1,2 -f findsym.in

Give the sites as a comma separated list with no space. Note that it doesn't matter whether you choose a site which is actually magnetic.

Equivalent sites
^^^^^^^^^^^^^^^^

You can also use the -e switch to determine the magnetic Hamiltonian for all sets of sites which are related by a symmetry operation. For example if sites 1 and 2 are related by a symmetry operation then the Hamiltonians for sites 1,1 and 2,2 are related.

Some other features
-------------------

Latex output
^^^^^^^^^^^^

The output can be also given in latex format, use --latex switch.

Selecting symmetries
^^^^^^^^^^^^^^^^^^^^

You can print the symmetries using --print-syms and you can choose to use only some symmetries using the --syms switch. Input the symmetries as a comma separated list and you can also use ranges:

::

    --syms 1-3,7,9-12

Symmetry without spin-orbit coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Symmetry without spin-orbit coupling is higher since the spin is then not coupled to the lattice directly. In a nonmagnetic system this means that any pure spin rotation is a symmetry. In a magnetic system a pure spin rotation, which in combination with some other symmetry operation leave the magnetic order invariant is a symmetry without spin-orbit coupling.

This is implemented for the response mode using the --noso switch.

However, this feature is very experimental!!! The alghorithm used is not guaranteed to obtain all of the symmetry operations and the feature is not tested very much so there may be bugs.

Tests
^^^^^^

There are tests in a folder tests/. You can run them using script run_tests.sh. This script demonstrates the main functionalities of the program. You can check that the results you obtained are correct using the script check_res.sh. If a difference is found, vimdiff is run, so that you can check whether there is an important difference.

License
-------

This code is distributed under the `MPL 2.0
license <https://www.mozilla.org/en-US/MPL/2.0/FAQ/>`__.

If the code is used in any research that leads to a publication, you should cite this website.

Contact
--------

If you have any issues or questions, you can contant me by email at jakub.zelezny@gmail.com or create issue here on bitbucket.