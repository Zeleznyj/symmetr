Docker
======

It is possible to run the code using Docker, which should work also on Windows and macOS.

To do so, install Docker for your operating system: https://www.docker.com/community-edition.

Then run the following command in terminal (powershell in Windows):

::

    docker run zeleznyj/symmetr


The program will be automatically download from the `Docker repository <https://hub.docker.com/r/zeleznyj/symmetr/>`_. The syntax is the same as for the normal version so, you can do for example:

::

    docker run zeleznyj/symmetr res s E -g P4mm

If you want to use findsym input, it is necessary to mount your working directory into docker, like this


::

    docker run -v $PWD:/workdir zeleznyj/symmetr res s E -f findsym.in

On linux $PWD contains the location of your current directory. On different OSs this may have to be modified or you can specify the path manually. The program will then look for the file findsym.in in your working directory.
