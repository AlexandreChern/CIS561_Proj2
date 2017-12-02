# PETSc
This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

This package requires the MPI.jl package be installed.  Once it is installed you should be able to run both Julia and Petsc in parallel using MPI for all communication.  The testing verifies that PETSc can be used both serially and in parallel.


To run PETSC in parallel, do:

```
mpirun -np 4 julia ./name_of_file
```

for a 4 processor run.
Note that this launches 4 independent Julia processes.  They are not aware of each other using Julia's built-in parallelism, and MPI is used for all communications.  

To run in serial, do:
```
julia ./name_of_file
```

Even when running serially, the MPI.jl package must be installed.


## Notes on wrapping functions
  Wrappers generated by Clang.jl are in the src/auto directory.  Although not quite usable, the functions can be made useable with a few simple modifications:
  * Pass `comm.val` instead of `comm` itself for MPI communicators, and change the type to `comm_type`, which is typealiased to the the type used by the MPI.jl package
  * Pass obj.pobj instead of obj for Petsc objects as Ptr{Void}
  * For each Petsc object, you must create a type that holds a void pointer called pobj and use that in place of the (incorrect) type generated by Clang.jl
  * For every function you add, create a test


## Directory Structure
  `/src` : source files.  PETSc.jl is the main file containing initialization, with the functions for each type of Petsc object in its own file.  All constants are declared in `petsc_constants.jl`.

  `/src/auto`: auto generated wrappers from Clang.jl.  Not directly useful, but easy to modify to make useful

  `/test` : contains `runtest.jl`, which does some setup and runs all tests on the current PETSc installation.  Tests for each type of Petsc object (mirroring the files in `/src`) are contained in separate files.  The file `runtests.sh` builds PETSc and runs the tests on combinations of integer size, floating point precision, and type of scalar (real or complex).

  `/deps` : builds Petsc if needed.  See description below
  `/docs`: documentation (using Documenter.jl)

## Building Petsc (or not)
Upon installation of the package, the build script (`build.jl`) checks for the environmental variables `PETSC_DIR` and `PETSC_ARCH`.  If both are present, it does nothing, otherwise is downloads and builds Petsc using the script `install_petsc.sh`.  Note that the script assumes BLAS, LAPACK and MPI (currently only MPICH is supported) are already installed.  See `.travis.yml` to see the Ubuntu packages used for testing.  When complete, it generates two files, `use_petsc.sh` and `petsc_evars`, which contains the values of `PETSC_DIR` and `PETSC_ARCH` for the Petsc installation.

  At runtime, the module checks for the presence of the environmental variables and uses them if found, otherwise it uses the values in `petsc_evars`.  This enables you to use different Petsc installations if you choose.  When the module is imported in the user code, it auto-detects the size of the Petsc integers, the precision of floats, and whether scalars are real or complex.


## Installing MPI.jl
This package requires MPI.jl, although it is not listed in the REQUIRE file
because that would download the release version of MPI.jl, which does not work.
Instead, you must use an older verion.  
After you have an MPI implementation installed, Pkg.build("PETSc2") will
install MPI.jl and then Petsc, according to the description above.
If you wish to install MPI.jl manually, see `deps/build.jl`.

[![Build Status](https://travis-ci.org/OptimalDesignLab/PETSc2.svg?branch=master)](https://travis-ci.org/OptimalDesignLab/PETSc2.jl)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OptimalDesignLab.github.io/PETSc2.jl/latest)
