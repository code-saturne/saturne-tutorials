Mixer Vessel 2D (Transient Rotor / Sliding Mesh)
=================================================

Quasi-2D stirred tank simulation using the transient rotor (sliding
mesh) turbomachinery model of code_saturne.  A four-blade impeller
rotates at 60 RPM inside a cylindrical vessel fitted with four
radial baffles.

This case complements tutorial 19 (frozen rotor) by using a true
transient approach with sliding mesh at the rotor/stator interface.

Setup
-----

- Turbomachinery model: transient rotor (sliding mesh)
- Turbulence: laminar (no model)
- Mesh: two structured hexahedral meshes (Gmsh)
  - Rotor: hub (r=0.02 m) to interface (r=0.045 m), 768 cells
  - Stator: interface (r=0.045 m) to tank wall (R=0.1 m), 1920 cells
  - Thickness: 0.01 m (1 cell in Z, quasi-2D)
- Rotor/stator interface: face joining at r = 0.045 m
- Thin walls (zero-thickness internal boundaries):
  - 4 blades at 0/90/180/270 deg in rotor mesh (r < 0.045 m)
  - 4 baffles at 45/135/225/315 deg in stator mesh (0.07 < r < 0.1 m)
- Boundary conditions:
  - hub, tank_wall: no-slip wall
  - z0, z1: symmetry (quasi-2D)
- Velocity-pressure coupling: SIMPLEC
- Time stepping: fixed dt = 0.001 s, t_final = 5 s (5 revolutions)
- Output: every 100 time steps (transient connectivity)

Rotor configuration:

| Parameter       | Value                  |
|-----------------|------------------------|
| Rotation speed  | 6.2832 rad/s (60 RPM)  |
| Rotation axis   | (0, 0, 1)             |
| Rotation center | (0, 0, 0)             |
| Rotor cells     | cylinder r < 0.045 m  |

Physical properties:

| Property               | Value        |
|------------------------|--------------|
| Density                | 1 kg/m^3     |
| Kinematic viscosity    | 1e-5 m^2/s   |

Mesh generation
---------------

Two separate meshes are prefesent for stator and rotor.

The rotor mesh rotates as a rigid body at each time step.  The
rotor/stator interface is handled by a face joining that recomputes
connectivity at each iteration.
