/*============================================================================
 * User initialization for the Sedov blast wave problem.
 *
 * A high-pressure region (r < r0) surrounded by a low-pressure ambient
 * creates a cylindrical shock wave propagating radially outward.
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
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define initial conditions for variables.
 *
 * This function is not yet called by code_saturne, use
 * `cs_user_initialization` for the time being
 * This function is called before the beginning of the computation
 * allowing an overload of the GUI defined initialization (called just
 * after cs_gui_initial_conditions).
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t *domain)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_3_t *cell_cen =
    (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

  /* Access field pointers */
  cs_real_t *rho = CS_F_(rho)->val;
  cs_real_t *p   = CS_F_(p)->val;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

  /* Sedov blast wave parameters:
   *   rho0   - uniform ambient density
   *   pinfty - ambient pressure (very low)
   *   phot   - high pressure inside the blast region
   *   r0     - initial blast radius */
  const cs_real_t rho0   = 1.0;
  const cs_real_t pinfty = 0.01;
  const cs_real_t phot   = 100.0;
  const cs_real_t r0     = 0.005;

  /* Loop over all cells: set uniform density and zero velocity,
   * then apply a pressure discontinuity at radius r0 */
  for (cs_lnum_t c = 0; c < n_cells; c++) {
    double x = cell_cen[c][0];
    double y = cell_cen[c][1];
    double r = sqrt(x*x + y*y);

    rho[c]    = rho0;
    vel[c][0] = 0.0;
    vel[c][1] = 0.0;
    vel[c][2] = 0.0;

    if (r < r0)
      p[c] = phot;     /* hot region: energy source */
    else
      p[c] = pinfty;   /* ambient: near-vacuum pressure */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
