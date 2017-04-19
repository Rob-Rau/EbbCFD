# EbbCFD
EbbCFD is an unstructured, second order spatial, finite volume solver. It can solve both the Euler equations and the compressible Navier-Stokes equations.

Supported time integration:
* 1st order forward euler
* 2nd order total variation diminishing (TVD) Runge-Kutta
* 2nd order Runge-Kutta
* 4th order Runge-Kutta

## Features
* 2D second order spatial solver
* parallelization using OpenMPI
* clunky configuration files
* non-uniform interface that leaves your head scratching a bit (look at the plotting part, what am I doing)
* dependancy on cblas that isn't really used.

## Building
> dub build --config=solver

## Running
> ebb-solve -c configFile.json

## Plotting
> ebb-reconstruct configFile.json saveFile.esln
> PlotEbbSolution mach recon_meshFile.mmsh recon_saveFile.esln
