# ParticleModule
A python module for driving Lagrangian particles via vtk formatted data with collisions.

Still under heavy development.

Directory structure

* docs : Automatic documentation via sphinx.
* examples : Examples of code use
* fluidity_tools : Useful modules for coupling with Fluidity, http://fluidityproject.github.io
* particle_module : Main source directory
* particle_module/tests : Location for pytest tests.
* schemas : xml schemas for using the particle module with diamond.
* src : Source code for compiled modules
* CMakeLists.txt : cmake file
* setup.py.in : base template for disttools
* CONTRIBUTING.md : information on contributing
* README.md : This file.

Prerequisits
------------

* Python, currently only tested on Python 2.7
* VTK, version 6 or higher, including the python wrappings for Python 2.7.
* numpy
* scipy
* mpi4py

and optionally:

* libspud
* lxml

Building
--------

In the ParticleModule base directory run

```bash
cmake -DCMAKE_BUILD_TYPE=Release .
make
```

This assumes that cmake is installed on your system, and sanely configured.

Testing
-------

After building, the tests can be run using

```bash
py.test
```

in the base directory. This assumes pytest is installed.

Installation
------------

After building the code, run

```bash
sudo make install
```

This installs the package under your default path for python packages.
