# Fluid Simulation using Navier-Stokes (C++)

![Fluid Simulator](FluidSim.gif)

This project implements a real-time 2D fluid simulation based on the incompressible Navier-Stokes equations.  
It uses a grid-based Stable Fluids solver and Raylib for visualization. The simulation models density and velocity fields and allows interactive injection of dye and momentum.

## Features

- Real-time fluid simulation
- Semi-Lagrangian advection
- Diffusion and viscosity modeling
- Pressure projection for incompressibility
- Interactive input (mouse-driven forces and dye)
- Boundary condition handling
- Iterative linear solver (Gauss-Seidel)

## Mathematical Model

Steps:

1. Diffuse velocity  
2. Project velocity  
3. Advect velocity  
4. Project again 
5. Diffuse density  
6. Advect density  

## Project Structure

main.cpp — Application loop and rendering  
FluidSimulation.h — Simulation interface and data structures  
FluidSimulation.cpp — Numerical solver implementation  

## Controls

Left mouse button: add density and velocity  
C key: clear simulation  

## Build

Requires C++17 and Raylib.

Linux / macOS:

g++ main.cpp FluidSimulation.cpp -o fluid -lraylib -std=c++17

Windows (MinGW):

g++ main.cpp FluidSimulation.cpp -o fluid.exe -lraylib -std=c++17

## References

- Jos Stam — Stable Fluids
- Real-Time Fluid Dynamics for Games
- Navier-Stokes equations