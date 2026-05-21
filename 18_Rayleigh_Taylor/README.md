Rayleigh-Taylor Instability
===========================

2D Rayleigh-Taylor instability simulated with the VOF (Volume of Fluid)
method in code_saturne.

A heavy fluid rests on top of a light fluid under gravity. Any
perturbation of the interface grows exponentially, producing
characteristic mushroom-shaped structures as the heavy fluid sinks
and the light fluid rises.

Setup
-----

- Multiphase model: VOF (Volume of Fluid)
- Turbulence: off (laminar)
- Mesh: 64 x 192 x 1 Cartesian (generated at runtime)
- Domain: [0, 1] x [-1, 2] x [0, 1]
- Periodicity: x-direction
- Boundary conditions:
  - x0, x1: periodic
  - y0, y1, z0, z1: symmetry
- Gravity: g = -1 m/s^2 (y-axis)
- Time stepping: adaptive (CFL = 1), dt_ref = 0.0005 s, t_final = 4.4 s

Fluid properties:

| Property   | Light fluid | Heavy fluid |
|------------|-------------|-------------|
| Density    | 1.0         | 5.0         |
| Viscosity  | 0.001       | 0.001       |

Atwood number: At = (rho_h - rho_l) / (rho_h + rho_l) = 0.667

Initial conditions (set in `cs_user_initialization.cpp`):

| Parameter                     | Value           |
|-------------------------------|-----------------|
| Interface height y0           | 1.0 m           |
| Perturbation amplitude eps    | 0.03 m          |
| Perturbation wavelength       | Lx (2 pi mode)  |
| Interface smoothing thickness | 0.005 m         |

The VOF fraction is initialized with a smoothed sinusoidal interface.
Velocity is zero everywhere; the instability is driven purely by gravity.

Reference
---------

Lord Rayleigh, "Investigation of the character of the equilibrium of an
incompressible heavy fluid of variable density", Proc. London Math. Soc.,
14:170-177, 1883.

G. I. Taylor, "The instability of liquid surfaces when accelerated in a
direction perpendicular to their planes", Proc. R. Soc. Lond. A,
201(1065):192-196, 1950.
