/*============================================================================
 * User initialization for the Rayleigh-Taylor instability (VOF).
 *
 * A heavy fluid (rho=5) rests on top of a light fluid (rho=1) under
 * gravity. A sinusoidal perturbation of the interface triggers the
 * instability, producing characteristic mushroom-shaped structures.
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
 * \brief Rayleigh-Taylor instability initialization (VOF).
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define initial conditions for the Rayleigh-Taylor problem.
 *
 * Sets zero velocity and a smoothed sinusoidal VOF interface between
 * a heavy fluid on top and a light fluid below, under gravity.
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

  /* ---- Parameters ---- */

  const cs_real_t Lx = 1.0;        /* domain length in x */
  const cs_real_t y0 = 1.0;        /* interface height (middle of Ly=2) */
  const cs_real_t eps = 0.03;      /* perturbation amplitude */
  const cs_real_t thickness = 0.005; /* interface smoothing thickness */

  /* ---- Mesh coordinates (v9.x format) ---- */

  const cs_real_t (*cell_cen)[3] =
    domain->mesh_quantities->cell_cen;

  /* ---- Get VOF field ---- */

  cs_field_t *vf = cs_field_by_name_try("void_fraction");
  if (vf == nullptr)
    vf = cs_field_by_name_try("volume_fraction");
  if (vf == nullptr)
    vf = cs_field_by_name_try("alpha");

  if (vf == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "VOF fraction field not found.");

  /* ---- Velocity field ---- */

  cs_real_t *vel = CS_F_(vel)->val;

  /* ---- Initialization loop ---- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t x = cell_cen[c_id][0];
    const cs_real_t y = cell_cen[c_id][1];

    /* 1) Perturbed interface */
    const cs_real_t y_int =
      y0 + eps * sin(2.0 * M_PI * x / Lx);

    /* 2) Smoothed VOF (heavy fluid above) */
    const cs_real_t s = (y - y_int) / thickness;

    cs_real_t alpha =
      0.5 * (1.0 + tanh(s));

    /* Clamp */
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    vf->val[c_id] = alpha;

    /* 3) Initial velocity = 0 everywhere */
    vel[3*c_id + 0] = 0.0;
    vel[3*c_id + 1] = 0.0;
    vel[3*c_id + 2] = 0.0;
  }
}

END_C_DECLS