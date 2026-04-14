#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-
"""
ParaView post-processing (pvpython): 1D extraction and multi-resolution comparison.

This script opens several results (code_saturne .case files and/or reference
solutions in VTK .pvd/.vtu format), selects the available time step closest to
a user-specified time t_user, and samples the physical fields along a segment
P1->P2 using the ParaView PlotOverLine filter.

The extracted 1D profiles (density rho, pressure p, velocity u_x) are then
plotted with Matplotlib to compare different resolutions and/or an exact
solution. Zoom insets and arrow annotations highlight the shock tube structures
(rarefaction fan, contact surface, shock wave).

Inputs:
- code_saturne cases: CASE/RESU/*/postprocessing/RESULTS_FLUID_DOMAIN.case
- Exact reference: .pvd/.vtu file provided in CASES via the "path" key

Outputs:
- High-resolution PNG figure

Notes:
- Must be run with pvpython (ParaView environment) for paraview.simple.
- Explicit cleanup of ParaView objects (Delete) to limit memory usage.

Usage:
    pvpython postprocessSodTube.py

Author: fhermet
"""


import os, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from paraview.simple import OpenDataFile, PlotOverLine, GetTimeKeeper, Delete
from paraview.simple import servermanager
from vtk.util.numpy_support import vtk_to_numpy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# ---------------------------
# Matplotlib
# ---------------------------
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})

# ---------------------------
# Parameters
# ---------------------------
t_user = 0.2

# sampling segment
P1 = (-1.0, 0.0, 0.0)
P2 = ( 1.0, 0.0, 0.0)

ZOOMS = {
    "rarefaction": (-0.1, 0.1),
    "contact":     (0.10, 0.25),
    "contact_v":   (0.10, 0.30),
}

ARROWS = {
    "rho": dict(xy=(0.2, 0.32),  xytext=(-0.3, 0.15), text=r"Contact surface", fontsize=12),
    "p":   dict(xy=(0.38, 0.2),  xytext=(0.5, 0.3),   text=r"Shock wave",      fontsize=12),
    "u":   dict(xy=(-0.2, 0.20), xytext=(-0.75, 0.35),text=r"Rarefaction fan",   fontsize=12),
}

CASES = [
    dict(case_dir="../CASE", resolution=2000,   label=r"$N=2000$",   color="navy",  ls="-",  lw=1.4),
    dict(path="exact.pvd",        resolution=32000,   label=r"Exact", color="black", ls="--", lw=1.4),
]

FIELD_NAMES = {
    "P":   ["Pressure", "pressure"],
    "rho": ["Density", "rho"],
    "V":   ["Velocity", "velocity"],
}

# ---------------------------
# Utilities
# ---------------------------
def nearest_time(timesteps, t):
    """
    Return the available time step closest to the requested time.

    Parameters
    ----------
    timesteps : iterable of float
        Available time steps from the ParaView reader.
    t : float
        User-requested time.

    Returns
    -------
    float
        Actual time used (nearest neighbor).
    """
    ts = np.asarray(list(timesteps), dtype=float)
    return float(ts[np.argmin(np.abs(ts - float(t)))])

def reader_timesteps(reader):
    """
    Retrieve the list of time steps from a ParaView reader.

    Works with both:
    - code_saturne files (.case)
    - VTK files (.pvd, .vtu)

    Parameters
    ----------
    reader : paraview.simple.Proxy
        Already opened ParaView reader.

    Returns
    -------
    list of float
        Available time steps.
    """
    if getattr(reader, "TimestepValues", None):
        return list(reader.TimestepValues)
    return list(GetTimeKeeper().TimestepValues)

def resolve_saturne_case(case_dir, case_filename="RESULTS_FLUID_DOMAIN.case"):
    """
    Automatically resolve the path to a code_saturne .case file.

    Assumption:
        case_dir/RESU/*/postprocessing/RESULTS_FLUID_DOMAIN.case

    Parameters
    ----------
    case_dir : str
        Root directory of the code_saturne computation.
    case_filename : str
        Name of the .case file (default: standard code_saturne name).

    Returns
    -------
    str
        Full path to the unique .case file found.

    Raises
    ------
    FileNotFoundError
        If no .case file is found.
    RuntimeError
        If multiple .case files are found.
    """
    pattern = os.path.join(case_dir, "RESU", "*", "postprocessing", case_filename)
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No .case file found via:\n  {pattern}")
    if len(matches) > 1:
        raise RuntimeError("Multiple .case files found:\n  " + "\n  ".join(matches))
    return matches[0]

def get_np_any(dataset, names):
    """
    Extract a VTK field as a NumPy array.

    Searches successively in:
    - PointData
    - CellData

    Parameters
    ----------
    dataset : vtkDataSet
        Data retrieved via servermanager.Fetch().
    names : list of str
        Possible field names (e.g. ["Pressure", "pressure"]).

    Returns
    -------
    numpy.ndarray
        NumPy array of the field found.

    Raises
    ------
    RuntimeError
        If no matching field is found.
    """
    pd, cd = dataset.GetPointData(), dataset.GetCellData()
    for name in names:
        arr = pd.GetArray(name) or cd.GetArray(name)
        if arr is not None:
            return vtk_to_numpy(arr)
    raise RuntimeError(f"Array not found among: {names}")

def sample_line(reader, line, t):
    """
    Sample physical fields along a line at a given time.

    Extracted fields:
    - abscissa (arc_length or x coordinate)
    - velocity Ux
    - pressure P
    - density rho

    Compatible with:
    - code_saturne (.case)
    - exact solutions (.pvd / .vtu)

    Parameters
    ----------
    reader : paraview.simple.Proxy
        ParaView reader for the case.
    line : paraview.simple.Proxy
        Configured PlotOverLine filter.
    t : float
        Time at which to sample.

    Returns
    -------
    x : numpy.ndarray
        Abscissa along the line.
    Ux : numpy.ndarray
        X-component of velocity.
    P : numpy.ndarray
        Pressure.
    rho : numpy.ndarray
        Density.
    """
    reader.UpdatePipeline(time=t)
    line.UpdatePipeline(time=t)
    data = servermanager.Fetch(line)

    arc = data.GetPointData().GetArray("Point_X")
    x = vtk_to_numpy(arc) if arc is not None else vtk_to_numpy(data.GetPoints().GetData())[:, 0]

    P   = get_np_any(data, FIELD_NAMES["P"])
    rho = get_np_any(data, FIELD_NAMES["rho"])
    V   = get_np_any(data, FIELD_NAMES["V"])
    Ux  = V[:, 0] if (hasattr(V, "ndim") and V.ndim == 2) else V

    return x, Ux, P, rho

def add_zoom(ax, x_ref, y_ref, xlim, ylim=None, loc="upper right",
             width="38%", height="38%", borderpad=1.0):
    """
    Add a zoom inset on an existing matplotlib axis.

    All curves already plotted on the main axis are automatically
    copied into the inset with the same style.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Main axis.
    x_ref : numpy.ndarray
        Reference abscissa (used for auto-ylim).
    y_ref : numpy.ndarray
        Reference field for computing zoom bounds.
    xlim : tuple(float, float)
        X-bounds of the zoom.
    ylim : tuple(float, float) or None
        Y-bounds of the zoom (auto if None).
    loc : str
        Inset position.
    width, height : str
        Relative size of the inset.
    borderpad : float
        Padding between inset and main axis.

    Returns
    -------
    matplotlib.axes.Axes
        Created inset axis.
    """
    axins = inset_axes(ax, width=width, height=height, loc=loc, borderpad=borderpad)

    # copy all curves already plotted (same styles)
    for line in ax.get_lines():
        axins.plot(line.get_xdata(), line.get_ydata(),
                   color=line.get_color(), ls=line.get_linestyle(), lw=line.get_linewidth())

    axins.set_xlim(*xlim)

    if ylim is None:
        m = (x_ref >= xlim[0]) & (x_ref <= xlim[1])
        yy = y_ref[m]
        if yy.size:
            axins.set_ylim(yy.min() - 0.01, yy.max() + 0.01)
    else:
        axins.set_ylim(*ylim)

    axins.grid(True, alpha=0.25)
    axins.tick_params(labelsize=8)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.3", lw=0.8)
    return axins

def add_arrow(ax, xy, xytext, text, fontsize=12):
    """
    Add an arrow annotation on a plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    xy : tuple(float, float)
        Point the arrow points to.
    xytext : tuple(float, float)
        Text position.
    text : str
        Displayed text.
    fontsize : int
        Font size.
    """
    ax.annotate(
        text, xy=xy, xytext=xytext, textcoords="data", xycoords="data",
        ha="left", va="center", fontsize=fontsize,
        arrowprops=dict(arrowstyle="->", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8),
    )

# ---------------------------
# Extraction (case loop)
# ---------------------------
results = []

for k, c in enumerate(CASES):
    # ----------------------------------------------------------------
    # 1) Resolve the file path to open
    #    - either an exact file (.pvd/.vtu) provided directly
    #    - or a code_saturne computation (.case) resolved automatically
    # ----------------------------------------------------------------
    filepath = c.get("path") or resolve_saturne_case(c["case_dir"])
    print(f"[{k}] Opening: {filepath}")

    # ----------------------------------------------------------------
    # 2) Create ParaView reader and query metadata
    # ----------------------------------------------------------------
    reader = OpenDataFile(filepath)
    reader.UpdatePipelineInformation()

    # ----------------------------------------------------------------
    # 3) Retrieve available time steps
    # ----------------------------------------------------------------
    ts = reader_timesteps(reader)
    if not ts:
        raise RuntimeError(f"No time steps detected in {filepath}")

    # ----------------------------------------------------------------
    # 4) Select the time step closest to the requested one
    #    (important if t_user does not match an output exactly)
    # ----------------------------------------------------------------
    t = nearest_time(ts, t_user)
    print(f"Requested time={t_user} -> used={t}")

    # ----------------------------------------------------------------
    # 5) Create PlotOverLine filter
    #    - line defined by P1 -> P2
    #    - resolution controlled by the case (mesh / reference)
    # ----------------------------------------------------------------
    line = PlotOverLine(Input=reader,
                        Point1=list(P1),
                        Point2=list(P2),
                        Resolution=int(c["resolution"]))

    # ----------------------------------------------------------------
    # 6) Sample physical fields along the line
    #    - abscissa x (arc length or coordinate)
    #    - velocity Ux
    #    - pressure P
    #    - density rho
    # ----------------------------------------------------------------
    x, Ux, P, rho = sample_line(reader, line, t)

    # ----------------------------------------------------------------
    # 7) Store results for later plotting
    #    (plotting is intentionally separated from extraction)
    # ----------------------------------------------------------------
    results.append(dict(case=c, t=t, x=x, Ux=Ux, P=P, rho=rho))

    # ----------------------------------------------------------------
    # 8) Explicit cleanup of ParaView objects
    #    (important to avoid memory leaks with pvpython)
    # ----------------------------------------------------------------
    Delete(line)
    Delete(reader)

# Reference for auto-ylim of zooms
x_ref   = results[0]["x"]
Ux_ref  = results[0]["Ux"]
P_ref   = results[0]["P"]
rho_ref = results[0]["rho"]

# ---------------------------
# Plot
# ---------------------------
fig, (axR, axP, axU) = plt.subplots(3, 1, sharex=True, figsize=(10, 14))

for r in results:
    c = r["case"]
    style = dict(color=c.get("color"), ls=c.get("ls", "-"), lw=c.get("lw", 1.4), label=c.get("label", "case"))
    axU.plot(r["x"], r["Ux"], **style)
    axP.plot(r["x"], r["P"],  **style)
    axR.plot(r["x"], r["rho"],**style)

# Zooms
add_zoom(axR, x_ref, rho_ref, ZOOMS["contact"],     loc="upper right", borderpad=1.5)
add_zoom(axP, x_ref, P_ref,   ZOOMS["rarefaction"], loc="lower left", width="22%", height="22%", borderpad=2.5)
add_zoom(axU, x_ref, Ux_ref,  ZOOMS["contact_v"],   loc="upper right", width="22%", height="22%", borderpad=2.5)

# Arrows
add_arrow(axR, **ARROWS["rho"])
add_arrow(axP, **ARROWS["p"])
add_arrow(axU, **ARROWS["u"])

# Labels
axU.set_ylabel(r"$u_x$", fontsize=16)
axP.set_ylabel(r"$p$",   fontsize=16)
axR.set_ylabel(r"$\rho$",fontsize=16)

for ax in (axR, axP, axU):
    ax.set_xlabel(r"$x$", fontsize=16)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.3)

# Single legend without duplicates (on axP)
handles, labels = axP.get_legend_handles_labels()
uniq = {}
for h, lab in zip(handles, labels):
    uniq.setdefault(lab, h)
axP.legend(uniq.values(), uniq.keys(), fancybox=True, frameon=True, framealpha=0.9, fontsize=12)

fig.suptitle(rf"Comparison at $t = {results[0]['t']:.2f}\,\mathrm{{s}}$", fontsize=20, y=0.96)
fig.tight_layout(rect=[0, 0, 1, 0.96])

out = f"results_t{t_user:.3f}.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"Figure saved: {out}")
