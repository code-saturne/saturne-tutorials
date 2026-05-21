Sedov Blast Wave
================

2D Sedov blast wave problem solved with the compressible module of
code_saturne (constant gamma, Euler equations).

A point-like energy release at the origin creates a cylindrical shock
wave that propagates radially. The shock radius follows the self-similar
Sedov-Taylor solution: R(t) proportional to t^(2/5) in 2D (t^(1/2) for
the cylindrical symmetry used here with a quarter-domain).

Setup
-----

- Compressible model: constant gamma (gamma = 1.4)
- Turbulence: off (inviscid)
- Mesh: 200 x 200 x 1 Cartesian (generated at runtime), domain [0, 0.2]^2
- Boundary conditions:
  - x0, y0: symmetry (quarter-domain exploiting axial symmetry)
  - x1, y1: supersonic outlet
  - z0, z1: symmetry (2D)
- Time stepping: adaptive (CFL = 1), dt_ref = 5e-5 s, t_final = 0.5 s

Initial conditions (set in `cs_user_initialization.cpp`):

| Region    | Density | Pressure | Velocity |
|-----------|---------|----------|----------|
| r < 0.005 | 1.0     | 100.0    | 0.0      |
| r >= 0.005| 1.0     | 0.01     | 0.0      |

Reference
---------

L. I. Sedov, "Similarity and Dimensional Methods in Mechanics",
Academic Press, 1959.

G. I. Taylor, "The formation of a blast wave by a very intense explosion",
Proc. R. Soc. Lond. A, 201(1065):159-186, 1950.
