Installation
=============

The easiest way how to install the program is using pip, simply run:

::
    
    pip install symmetr
    
    
The program currently only works on linux. If you use a different OS, you can run the program using docker, see :doc:`docker`. Although the program itself could be run on any operating system since it's written in python, most of the features rely on the findsym code from the `isotropy suite <http://stokes.byu.edu/iso/isolinux.php>`__, which only works on linux. Findsym is now included in the repository so it does not have
to be installed separately. See the symmetr/findsym/findsym.txt file for info on the input format for findsym or check the various input files present in tests/.

Only Python3 is supported.

If you don't want to use pip, you can directly download the repository as a zip file from `this link <https://bitbucket.org/zeleznyj/linear-response-symmetry/get/7d569d9b5dbe.zip>`__ or use git:

::

    git clone https://zeleznyj@bitbucket.org/zeleznyj/linear-response-symmetry.git 

This will create a folder called linear-response-symmetry. To install you can then run:

::

    python setup.py


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   docker
