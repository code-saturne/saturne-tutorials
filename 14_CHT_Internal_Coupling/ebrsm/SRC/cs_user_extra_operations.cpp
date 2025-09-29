/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
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

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  CS_UNUSED(domain);

  if (cs_log_default_is_active()) {

    cs_real_t volume_heat_flux
      = cs_notebook_parameter_value_by_name("v_heat_flux");

    cs_lnum_t *b_face_list;
    CS_MALLOC(b_face_list, domain->mesh->n_b_faces, cs_lnum_t);

    const cs_real_t *b_face_surf = domain->mesh_quantities->b_face_surf;

    const char *names[] = {"auto:internal_coupling_0_solid",
                           "auto:internal_coupling_0_fluid"};

    cs_function_t *th_flux_func = cs_function_by_name("boundary_thermal_flux");

    cs_real_t wall_flux[2] = {0, 0};

    for (int i = 0; i < 2; i++) {
      cs_lnum_t n_sel_faces = 0;
      cs_selector_get_b_face_list(names[i], &n_sel_faces, b_face_list);

      /* Get thermal flux */

      cs_real_t *th_flux_density;
      CS_MALLOC(th_flux_density, n_sel_faces, cs_real_t);

      cs_function_evaluate(th_flux_func,
                           domain->time_step,
                           CS_MESH_LOCATION_BOUNDARY_FACES,
                           n_sel_faces,
                           b_face_list,
                           th_flux_density);

      for (cs_lnum_t idx = 0; idx < n_sel_faces; idx++) {
        cs_lnum_t f_id = b_face_list[idx];
        wall_flux[i] += th_flux_density[idx] * b_face_surf[f_id];
      }

      CS_FREE(th_flux_density);
    }

    CS_FREE(b_face_list);

    cs_parall_sum(2, CS_REAL_TYPE, wall_flux);

    bft_printf("\n"
               " Thermal fluxes\n"
               " --------------\n"
               "  volume injection:  %g\n"
               "  solid to fluid:    %g\n"
               "  fluid to solid:    %g\n\n",
               volume_heat_flux, wall_flux[0], wall_flux[1]);

  }

  /* Check velocity is zero in solid zone
     ------------------------------------ */

  const cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  const cs_zone_t *z = cs_volume_zone_by_name("solid");

  cs_real_t max_vel = 0;
  for (cs_lnum_t idx = 0; idx < z->n_elts; idx++) {
    cs_lnum_t c_id = z->elt_ids[idx];
    max_vel = cs_math_fmax(max_vel, cs_math_3_square_norm(vel[c_id]));
  }
  cs_parall_max(1, CS_REAL_TYPE, &max_vel);
  max_vel = sqrt(max_vel);

  if (max_vel > 1e-12)
    bft_error(__FILE__, __LINE__, 0,
              "%s: nonzero velocity in solid zone (max. %g).",
              __func__, max_vel);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
