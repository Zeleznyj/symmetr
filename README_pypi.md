Python package for determining symmetry properties of crystals. See [bitbucket
repository](https://bitbucket.org/zeleznyj/linear-response-symmetry) for more details.

The code internally uses [findsym](https://stokes.byu.edu/iso/isolinux.php) to determine the symmetry of a given crystal. Since findsym only works on Linux, most features of Symmetr will also only work on Linux. However, it is possible to run the code on any platform through [Docker](https://bitbucket.org/zeleznyj/linear-response-symmetry/wiki/docker). The Windows Subsystem for Linux could likely be used for running the code, but this has not been tested. 

Currently both Python2 and Python3 are supported, however, the support for Python2 will likely be dropped in the future.

This code is still in development and may contain bugs! Use with caution!
