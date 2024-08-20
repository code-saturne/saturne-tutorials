/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particle location.
 *
 * This function aims at modifying injection coordinates, particle properties
 * and cell_id depending on the position are updated based on the modified
 * position after this function and before cs_user_lagr_in.
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modified by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               zone injection set data
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  NULL otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_in_force_coords(cs_lagr_particle_set_t         *particles,
                             const cs_lagr_injection_set_t  *zis,
                             const cs_lnum_t                 particle_range[2],
                             const cs_lnum_t                 particle_face_id[],
                             const cs_real_t                 visc_length[])
{
  /* Loop on new particles
     --------------------- */

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    /* particles are injected along pipe axis (x=0,y=0) */

    cs_real_t *part_coords = (cs_real_t *)cs_lagr_particles_attr(particles,
                                                                 p_id,
                                                                 CS_LAGR_COORDS);

    part_coords[0] = 0.; // x = 0
    part_coords[1] = 0.; // y = 0
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
