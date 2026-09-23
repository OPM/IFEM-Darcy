## Introduction

This repository contains applications for solving Darcy flow problems,
built using the [IFEM](https://github.com/OPM/IFEM) library.

### Getting all dependencies

1. Install IFEM from https://github.com/OPM/IFEM

### Getting the code

This is done by first navigating to a folder `<App root>` in which you want
the IFEM applications and typing

    git clone https://github.com/OPM/IFEM-Darcy

### Compiling the code

To compile, first navigate to the root catalogue `<App root>`.

    cd IFEM-Darcy
    mkdir Debug
    cd Debug
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    make

This will compile the library and the Darcy applications.
The executables can be found in the `bin` subfolder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

### Testing the code

IFEM is using cmake test system.
To compile and run all regression- and unit-tests, navigate to your build folder
(i.e., `<App root>/IFEM-Darcy/Debug`) and type

    make check

## Class overview

The integrands of the Darcy problem
are organized in the following class hierarchy:
```mermaid
graph TD;
    IFEM1([IntegrandBase]) --> IFEM2([HasGravityBase])
    IFEM2 --> DarcyBase
    DarcyBase --> Darcy
    DarcyBase --> DarcyAdvection
    Darcy --> DarcyTransport
    Darcy --> CompatibleDarcy
    IFEM1 --> DarcyTransportCorr
```

The simulation drivers are connected to the integrands as follows:
```mermaid
graph TD;
    IFEM1([SIM2D,SIM3D]) --> SIMDarcy
    IFEM2([SIMSolution]) --> SIMDarcy
    SIMDarcy -. has a .-> Darcy
    IFEM1 --> SIMDarcyAdvection
    IFEM2 --> SIMDarcyAdvection
    SIMDarcyAdvection -. has a .-> DarcyAdvection
    IFEM1 --> SIMDarcyTransportCorr
    SIMDarcyTransportCorr -. has a .-> DarcyTransportCorr
    IFEM3([SIMCoupled]) --> SIMDarcySchedule
    IFEM4([SIMadmin]) --> SIMDarcySchedule
    IFEM3 -. has a .-> SIMDarcy
    IFEM3 -. has a .-> SIMDarcyAdvection
```
