Mixer Vessel 2D (Frozen Rotor)
==============================

Quasi-2D stirred tank simulation using the frozen rotor turbomachinery
model of code_saturne.  A four-blade impeller rotates at 1000 RPM
inside a cylindrical vessel fitted with four radial baffles.

Setup
-----

- Turbomachinery model: frozen rotor
- Turbulence: k-epsilon with production limiter (scalable wall functions)
- Mesh: structured hexahedral (generated with Gmsh), 3072 cells
  - 32 blocks: 8 azimuthal sectors x 4 radial zones
  - Hub radius: 0.02 m, impeller tip: 0.05 m, MRF boundary: 0.06 m
  - Baffle tip: 0.07 m, tank radius: 0.1 m
  - Thickness: 0.01 m (1 cell in Z, quasi-2D)
- Thin walls (zero-thickness internal boundaries):
  - 4 blades at 0/90/180/270 deg, from hub to impeller tip (r < 0.05 m)
  - 4 baffles at 45/135/225/315 deg, from baffle tip to tank wall (0.07 < r < 0.1 m)
- Boundary conditions:
  - walls (hub, tank wall): no-slip wall
  - z0, z1: symmetry (quasi-2D)
- Velocity-pressure coupling: SIMPLEC
- Time stepping: fixed dt = 0.01 s, 100 iterations

Rotor configuration:

| Parameter       | Value                   |
|-----------------|-------------------------|
| Rotation speed  | 104.72 rad/s (1000 RPM) |
| Rotation axis   | (0, 0, 1)              |
| Rotation center | (0, 0, 0)              |
| Rotor cells     | cylinder r < 0.06 m    |

Physical properties:

| Property               | Value        |
|------------------------|--------------|
| Density                | 1 kg/m^3     |
| Kinematic viscosity    | 1e-5 m^2/s   |

Mesh
---------------

The mesh is a plain annular domain without blades or baffles.
Blades and baffles are defined as thin walls (zero-thickness internal
boundaries) in the `setup.xml`, using geometric selectors (plane
intersected with cylinder criteria).

