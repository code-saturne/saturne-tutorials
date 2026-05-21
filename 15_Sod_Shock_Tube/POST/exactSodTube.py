#!/usr/bin/env python3
"""
Exact solution of the 1D Riemann problem for the Euler equations
of a perfect gas.

This script computes the analytical solution of the Riemann problem
(shock and/or rarefaction) for given left/right initial conditions.
The solution is sampled on a uniform 1D mesh and exported in VTU/PVD
format for visualization (ParaView).

Features:
- Perfect gas (constant gamma)
- Exact resolution of pressure and velocity in the star region
- Temporal sampling at constant time step
- ParaView-compatible export

The exact Riemann solver formulation for the 1D Euler equations of
a perfect gas follows:

[1] E. F. Toro,
    "Riemann Solvers and Numerical Methods for Fluid Dynamics",
    3rd Edition, Springer, 2009.
    Chapters 4 and 5 (exact solver, shocks and rarefactions).

Usage:
    pvpython exactSodTube.py

Author: fhermet
"""
import os
import math
import numpy as np
import xml.etree.ElementTree as ET

# ------------------------------------------------------------
# Physical and numerical parameters
# ------------------------------------------------------------
gamma = 1.4

xmin, xmax = -1.0, 1.0
nx = 32000
x0 = 0.0  # discontinuity position

# time (must match setup.xml: 436 iterations, dt = 4.5662e-4 s)
dt = 0.00045662100456621
n_iter = 436
t_final = dt * n_iter
output_stride = n_iter                # write a VTU every 'output_stride' steps (set to 1 for all)
output_dir = "."

n_steps = int(round(t_final / dt))
t_final_effective = n_steps * dt
print(f"n_steps = {n_steps}, t_final_effective = {t_final_effective}")

# ------------------------------------------------------------
# Exact Riemann functions (perfect gas)
# ------------------------------------------------------------
def sound_speed(rho, p, gamma=gamma):
    """
    Speed of sound for a perfect gas.
    """
    return math.sqrt(gamma * p / rho)

def f_shock(p, rhoK, pK, gamma=gamma):
    """
    Function f(p) for a shock wave in the exact Riemann solver.

    Parameters
    ----------
    p : float
        Tested pressure.
    rhoK : float
        Initial state density (left or right).
    pK : float
        Initial state pressure.
    gamma : float
        Ratio of specific heats.

    Returns
    -------
    float
        Shock contribution to the compatibility function.
    """
    AK = 2.0 / ((gamma + 1.0) * rhoK)
    BK = (gamma - 1.0) / (gamma + 1.0) * pK
    return (p - pK) * math.sqrt(AK / (p + BK))

def f_rare(p, rhoK, pK, gamma=gamma):
    """
    Function f(p) for a rarefaction wave in the exact Riemann solver.
    """
    aK = sound_speed(rhoK, pK, gamma)
    return (2.0 * aK / (gamma - 1.0)) * ((p / pK) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)

def f_total(p, rhol, ul, pl, rhor, ur, pr, gamma=gamma):
    """
    Total compatibility function for the Riemann problem.

    Corresponds to:
        f_L(p) + f_R(p) + (u_R - u_L)

    Used to determine the star pressure p*.
    """
    if p > pl:
        fL = f_shock(p, rhol, pl, gamma)
    else:
        fL = f_rare(p, rhol, pl, gamma)
    if p > pr:
        fR = f_shock(p, rhor, pr, gamma)
    else:
        fR = f_rare(p, rhor, pr, gamma)
    return fL + fR + (ur - ul)

def df_dp_side(p, rhoK, pK, gamma=gamma):
    """
    Derivative of f(p) with respect to p for one side (left or right).
    """
    if p > pK:
        AK = 2.0 / ((gamma + 1.0) * rhoK)
        BK = (gamma - 1.0) / (gamma + 1.0) * pK
        sqrt_term = math.sqrt(AK / (p + BK))
        return sqrt_term * (1.0 - 0.5 * (p - pK) / (p + BK))
    else:
        aK = sound_speed(rhoK, pK, gamma)
        return (1.0 / (rhoK * aK)) * (p / pK) ** (-(gamma + 1.0) / (2.0 * gamma))

def df_total(p, rhol, ul, pl, rhor, ur, pr, gamma=gamma):
    """
    Derivative of the total compatibility function with respect to pressure.
    """
    return df_dp_side(p, rhol, pl, gamma) + df_dp_side(p, rhor, pr, gamma)

def pressure_guess_pvrs(rhol, ul, pl, rhor, ur, pr, gamma=gamma):
    """
    Initial estimate of the star pressure p* using the PVRS approximation
    (Pressure-Velocity Riemann Solver), cf. Toro (2009), Section 4.6.

    Returns
    -------
    float
        Initial pressure for Newton iteration.
    """
    aL = sound_speed(rhol, pl, gamma)
    aR = sound_speed(rhor, pr, gamma)
    pPV = 0.5 * (pl + pr) - 0.5 * (ur - ul) * (rhol * aL + rhor * aR) / (rhol + rhor)
    return max(1e-6, pPV)

def solve_p_star(rhol, ul, pl, rhor, ur, pr, gamma=gamma,
                 tol=1e-10, max_iter=50):
    """
    Solve for the star pressure p* using Newton's method, cf. Toro (2009), Section 4.4.

    Parameters
    ----------
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    float
        Pressure in the star region.
    """
    p = pressure_guess_pvrs(rhol, ul, pl, rhor, ur, pr, gamma)
    for _ in range(max_iter):
        f = f_total(p, rhol, ul, pl, rhor, ur, pr, gamma)
        df = df_total(p, rhol, ul, pl, rhor, ur, pr, gamma)
        if abs(df) < 1e-20:
            break
        p_new = p - f / df
        if p_new < 1e-6:
            p_new = 1e-6
        if abs(p_new - p) < tol:
            return p_new
        p = p_new
    return p

def compute_u_star(pstar, rhol, ul, pl, rhor, ur, pr, gamma=gamma):
    """
    Compute the velocity in the star region u* from p*.
    """
    if pstar > pl:
        fL = f_shock(pstar, rhol, pl, gamma)
    else:
        fL = f_rare(pstar, rhol, pl, gamma)
    if pstar > pr:
        fR = f_shock(pstar, rhor, pr, gamma)
    else:
        fR = f_rare(pstar, rhor, pr, gamma)
    return 0.5 * (ul + ur) + 0.5 * (fR - fL)

def star_region_densities(pstar, rhol, pl, rhor, pr, gamma=gamma):
    """
    Compute left and right densities in the star region.

    Returns
    -------
    rhoL_star, rhoR_star : tuple(float, float)
        Left and right star densities.
    """
    if pstar > pl:
        p_ratio = pstar / pl
        rhoL_star = rhol * (p_ratio + (gamma - 1.0) / (gamma + 1.0)) / \
                    ((gamma - 1.0) / (gamma + 1.0) * p_ratio + 1.0)
    else:
        rhoL_star = rhol * (pstar / pl) ** (1.0 / gamma)
    if pstar > pr:
        p_ratio = pstar / pr
        rhoR_star = rhor * (p_ratio + (gamma - 1.0) / (gamma + 1.0)) / \
                    ((gamma - 1.0) / (gamma + 1.0) * p_ratio + 1.0)
    else:
        rhoR_star = rhor * (pstar / pr) ** (1.0 / gamma)
    return rhoL_star, rhoR_star

def sample_riemann(x, t, rhol, ul, pl, rhor, ur, pr,
                   gamma=gamma, x0=0.0):
    """
    Sample the exact Riemann problem solution at time t, cf. Toro (2009), Section 4.5.

    Parameters
    ----------
    x : array_like
        Spatial positions.
    t : float
        Time.
    x0 : float
        Initial discontinuity position.

    Returns
    -------
    rho, u, p : ndarray
        Density, velocity and pressure fields.
    """
    x = np.asarray(x)

    if t == 0.0:
        rho = np.where(x < x0, rhol, rhor)
        u   = np.where(x < x0, ul,   ur)
        p   = np.where(x < x0, pl,   pr)
        return rho, u, p

    # solve once for this problem (independent of t)
    pstar = solve_p_star(rhol, ul, pl, rhor, ur, pr, gamma)
    ustar = compute_u_star(pstar, rhol, ul, pl, rhor, ur, pr, gamma)
    rhoL_star, rhoR_star = star_region_densities(pstar, rhol, pl, rhor, pr, gamma)

    aL = sound_speed(rhol, pl, gamma)
    aR = sound_speed(rhor, pr, gamma)

    # wave speeds
    if pstar > pl:
        sL = ul - aL * math.sqrt((gamma + 1.0)/(2.0*gamma) * (pstar/pl - 1.0) + 1.0)
        headL = tailL = None  # not used
        aL_star = None
    else:
        headL = ul - aL
        aL_star = aL * (pstar / pl) ** ((gamma - 1.0)/(2.0*gamma))
        tailL = ustar - aL_star
        sL = None

    if pstar > pr:
        sR = ur + aR * math.sqrt((gamma + 1.0)/(2.0*gamma) * (pstar/pr - 1.0) + 1.0)
        headR = tailR = None
        aR_star = None
    else:
        headR = ur + aR
        aR_star = aR * (pstar / pr) ** ((gamma - 1.0)/(2.0*gamma))
        tailR = ustar + aR_star
        sR = None

    xi = (x - x0) / t

    rho = np.zeros_like(xi)
    u   = np.zeros_like(xi)
    p   = np.zeros_like(xi)

    for i, s in enumerate(xi):
        if s < ustar:
            # left side
            if pstar > pl:
                # left shock
                if s < sL:
                    rho[i] = rhol
                    u[i]   = ul
                    p[i]   = pl
                else:
                    rho[i] = rhoL_star
                    u[i]   = ustar
                    p[i]   = pstar
            else:
                # left rarefaction
                if s < headL:
                    rho[i] = rhol
                    u[i]   = ul
                    p[i]   = pl
                elif s > tailL:
                    rho[i] = rhoL_star
                    u[i]   = ustar
                    p[i]   = pstar
                else:
                    u_i = (2.0/(gamma+1.0)) * (aL + 0.5*(gamma-1.0)*ul + s)
                    a_i = aL - 0.5*(gamma-1.0)*(u_i - ul)
                    rho[i] = rhol * (a_i / aL) ** (2.0/(gamma-1.0))
                    u[i]   = u_i
                    p[i]   = pl * (a_i / aL) ** (2.0*gamma/(gamma-1.0))
        else:
            # right side
            if pstar > pr:
                # right shock
                if s > sR:
                    rho[i] = rhor
                    u[i]   = ur
                    p[i]   = pr
                else:
                    rho[i] = rhoR_star
                    u[i]   = ustar
                    p[i]   = pstar
            else:
                # right rarefaction
                if s > headR:
                    rho[i] = rhor
                    u[i]   = ur
                    p[i]   = pr
                elif s < tailR:
                    rho[i] = rhoR_star
                    u[i]   = ustar
                    p[i]   = pstar
                else:
                    u_i = (2.0/(gamma+1.0)) * (-aR + 0.5*(gamma-1.0)*ur + s)
                    a_i = aR + 0.5*(gamma-1.0)*(u_i - ur)
                    rho[i] = rhor * (a_i / aR) ** (2.0/(gamma-1.0))
                    u[i]   = u_i
                    p[i]   = pr * (a_i / aR) ** (2.0*gamma/(gamma-1.0))

    return rho, u, p

# ------------------------------------------------------------
# VTU / PVD writing
# ------------------------------------------------------------
def write_vtu(filename, x, rho, u, p, E):
    """
    Write a 1D VTU file containing the physical fields.

    Written fields:
    - rho: density
    - velocity: velocity
    - pressure: pressure
    - energy: total energy
    - temperature: temperature (perfect gas)

    Parameters
    ----------
    filename : str
        Path to the VTU file.
    """
    npts = len(x)
    ncells = npts - 1

    vtk = ET.Element("VTKFile", type="UnstructuredGrid",
                     version="0.1", byte_order="LittleEndian")
    grid = ET.SubElement(vtk, "UnstructuredGrid")
    piece = ET.SubElement(grid, "Piece",
                          NumberOfPoints=str(npts),
                          NumberOfCells=str(ncells))

    point_data = ET.SubElement(piece, "PointData", Scalars="rho")

    def add_array(name, data):
        arr = ET.SubElement(point_data, "DataArray",
                            type="Float64", Name=name, format="ascii")
        arr.text = " ".join(map(str, data))

    add_array("rho", rho)
    add_array("velocity", u)
    add_array("pressure", p)
    add_array("energy", E)
    add_array("temperature", p/(rho*287.1))


    points = ET.SubElement(piece, "Points")
    coords = ET.SubElement(points, "DataArray",
                           type="Float64", NumberOfComponents="3",
                           format="ascii")
    coords.text = " ".join(f"{xi} 0 0" for xi in x)

    cells = ET.SubElement(piece, "Cells")

    conn = ET.SubElement(cells, "DataArray",
                         type="Int32", Name="connectivity", format="ascii")
    conn.text = " ".join(f"{i} {i+1}" for i in range(ncells))

    offs = ET.SubElement(cells, "DataArray",
                         type="Int32", Name="offsets", format="ascii")
    offs.text = " ".join(str(2 * (i+1)) for i in range(ncells))

    types = ET.SubElement(cells, "DataArray",
                          type="UInt8", Name="types", format="ascii")
    types.text = " ".join("3" for _ in range(ncells))  # 3 = VTK_LINE

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    tree = ET.ElementTree(vtk)
    tree.write(filename)
    print(f"[VTU] Written: {filename}")

def write_pvd(filename, files, times):
    """
    Write a PVD file referencing a time series of VTU files.

    Parameters
    ----------
    files : list of str
        VTU file names.
    times : list of float
        Associated times.
    """
    vtk = ET.Element("VTKFile", type="Collection",
                     version="0.1", byte_order="LittleEndian")
    coll = ET.SubElement(vtk, "Collection")

    for f, t in zip(files, times):
        ET.SubElement(coll, "DataSet",
                      timestep=str(t), group="", part="0", file=f)

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    tree = ET.ElementTree(vtk)
    tree.write(filename)
    print(f"[PVD] Written: {filename}")

# ------------------------------------------------------------
# Main program
# ------------------------------------------------------------
os.makedirs(output_dir, exist_ok=True)

dx = (xmax - xmin) / nx
x = np.linspace(xmin + 0.5*dx, xmax - 0.5*dx, nx)

# Initial conditions
rhol, ul, pl = 1.0, 0.0, 1
rhor, ur, pr = 0.125, 0.0, 0.1

vtu_files = []
time_values = []

# time loop at constant dt (sampling)
for step in range(1, n_steps + 1):
    t = step * dt

    rho, u, p = sample_riemann(x, t, rhol, ul, pl, rhor, ur, pr, gamma, x0=x0)
    E = p/(gamma - 1.0) + 0.5*rho*u**2

    if (step % output_stride == 0) or (step == n_steps):
        base = f"exact_{step:06d}.vtu"
        fullpath = os.path.join(output_dir, base)
        write_vtu(fullpath, x, rho, u, p, E)
        vtu_files.append(base)
        time_values.append(t)

# PVD file
pvd_path = os.path.join(output_dir, "exact.pvd")
write_pvd(pvd_path, vtu_files, time_values)
