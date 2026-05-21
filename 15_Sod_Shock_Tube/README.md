Sod Shock Tube
==============

Classical 1D Riemann problem (Sod, 1978) solved with the compressible
module of code_saturne (constant gamma, Euler equations).

A diaphragm at x = 0 separates a high-pressure region (left) from a
low-pressure region (right). Removing the diaphragm produces three
characteristic wave structures: a left-moving rarefaction fan, a contact
discontinuity, and a right-moving shock wave.

Setup
-----

- Compressible model: constant gamma (gamma = 1.4)
- Turbulence: off (inviscid)
- Mesh: 2000 x 1 x 1 Cartesian (generated at runtime), x in [-1, 1]
- Boundary conditions: symmetry on all faces
- Time stepping: fixed dt, CFL ~ 0.54, 436 iterations (t_final ~ 0.2 s)

Initial conditions:

| Region | Density | Pressure | Velocity |
|--------|---------|----------|----------|
| Left   | 1.0     | 1.0      | 0.0      |
| Right  | 0.125   | 0.1      | 0.0      |

The user-defined file `CASE/SRC/cs_user_physical_properties.cpp` forces
molecular viscosity, volume viscosity, and thermal conductivity to zero
for a true inviscid (Euler) simulation.

Post-processing
---------------

Scripts in `POST/` require pvpython (ParaView environment):

- `exactSodTube.py`: generates the exact Riemann solution in VTU/PVD format.
- `postprocessSodTube.py`: extracts 1D profiles from code_saturne results
  and compares them against the exact solution.

Reference
---------

G. A. Sod, "A survey of several finite difference methods for systems of
nonlinear hyperbolic conservation laws", J. Comput. Phys., 27(1):1-31, 1978.

E. F. Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics",
3rd Edition, Springer, 2009.
