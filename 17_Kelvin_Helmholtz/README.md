Kelvin-Helmholtz Instability
============================

2D Kelvin-Helmholtz instability simulated with the VOF (Volume of Fluid)
method in code_saturne.

The instability develops at the interface between two fluid layers moving
in opposite directions. A sinusoidal perturbation of the interface grows
into the characteristic rolling vortex structures.

Setup
-----

- Multiphase model: VOF (Volume of Fluid)
- Turbulence: off (laminar)
- Mesh: 128 x 64 x 1 Cartesian (generated at runtime)
- Domain: [0, 1] x [0.25, 0.75] x [0, 1]
- Periodicity: x-direction
- Boundary conditions:
  - x0, x1: periodic
  - y0, y1: symmetry
  - z0, z1: symmetry (2D)
- Time stepping: fixed dt = 0.004 s, t_final = 10 s

Fluid properties:

| Property   | Phase 1 | Phase 2 |
|------------|---------|---------|
| Density    | 1.0     | 2.0     |
| Viscosity  | 1e-4    | 1e-4    |

Initial conditions (set in `cs_user_initialization.cpp`):

| Parameter                    | Value             |
|------------------------------|-------------------|
| Shear velocity U0            | 0.5 m/s           |
| Shear layer thickness delta  | 0.01 m            |
| Interface perturbation eps   | 0.03              |
| Perturbation wavelength      | Lx/2 (4 pi mode)  |

The velocity field is initialized with a hyperbolic tangent shear profile
and the VOF fraction with a smoothed sinusoidal interface.

Reference
---------

Lord Kelvin (W. Thomson), "Hydrokinetic solutions and observations",
Phil. Mag., 42:362-377, 1871.

H. von Helmholtz, "On discontinuous movements of fluids",
Phil. Mag., 36:337-346, 1868.
