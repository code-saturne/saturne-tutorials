/*============================================================================
 * User initialization for the Kelvin-Helmholtz instability (VOF).
 *
 * A shear layer between two fluid phases with a sinusoidal perturbation
 * of the interface triggers the Kelvin-Helmholtz rolling vortex instability.
 *============================================================================*/

/* code_saturne version 9.1 */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.cpp
 *
 * \brief Kelvin-Helmholtz instability initialization (VOF).
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define initial conditions for the Kelvin-Helmholtz problem.
 *
 * Sets a hyperbolic tangent shear velocity profile and a smoothed
 * sinusoidal VOF interface between two fluid phases.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t  *domain)
{
  const cs_lnum_t n_cells = domain->mesh->n_cells;

  /* Do not reinitialize on restart */
  if (domain->time_step->nt_prev > 0)
    return;

  /* ---- Parameters (adapt if needed) ---- */

  const cs_real_t Lx = 1.0;         /* domain length in x */
  const cs_real_t y0 = 0.5;         /* mean interface position */
  const cs_real_t eps = 0.03;       /* interface perturbation amplitude */
  const cs_real_t U0 = 0.5;         /* shear velocity magnitude */
  const cs_real_t delta = 0.01;     /* shear layer thickness */
  const cs_real_t thickness = 0.01; /* interface smoothing thickness */

  /* ---- Get mesh coordinates (v9.x: array [n_cells][3]) ---- */

  const cs_real_t (*cell_cen)[3] =
    domain->mesh_quantities->cell_cen;

  /* ---- Get fields ---- */

  /* Volume fraction field (VOF) */
  cs_field_t *vf = cs_field_by_name_try("void_fraction");
  if (vf == nullptr)
    vf = cs_field_by_name_try("volume_fraction");
  if (vf == nullptr)
    vf = cs_field_by_name_try("alpha");

  if (vf == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "VOF fraction field not found (tried: void_fraction, volume_fraction, alpha).");

  /* Velocity field */
  cs_real_t *vel = CS_F_(vel)->val;

  /* ---- Initialization loop ---- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t x = cell_cen[c_id][0];
    const cs_real_t y = cell_cen[c_id][1];

    /* 1) Wavy interface */
    const cs_real_t y_int =
      y0 + eps * sin(4.0 * M_PI * x / Lx);

    /* 2) Smoothed VOF using tanh */
    const cs_real_t s = (y - y_int) / thickness;

    cs_real_t alpha =
      0.5 * (1.0 + tanh(s));

    /* Clamp safety */
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    vf->val[c_id] = alpha;

    /* 3) Smooth shear velocity profile */
    const cs_real_t r = (y - y0) / delta;

    const cs_real_t u =
      U0 * tanh(r);

    vel[3*c_id + 0] = u;     /* u */
    vel[3*c_id + 1] = 0.0;   /* v */
    vel[3*c_id + 2] = 0.0;   /* w */
  }
}

END_C_DECLS