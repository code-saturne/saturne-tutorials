/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <stdio.h>

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
cs_user_extra_operations(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  const cs_lnum_t n_cells = cs_glob_mesh-> n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)cs_glob_mesh->i_face_cells;

  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_real_3_t *cell_cen
    = (cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;
  const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_t *i_face_surf = cs_glob_mesh_quantities->i_face_surf;

  const cs_real_t *cvar_temp = CS_F_(t)->val;
  const cs_real_t *cvar_rho = CS_F_(rho)->val;

  const cs_real_t *f_b_temp
    = (const cs_real_t *)cs_field_by_name("boundary_temperature")->val;

  int iflmas = cs_field_get_key_int(CS_F_(vel), cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  const cs_real_t *coefap = CS_F_(t)->bc_coeffs->a;

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {

    /* define geometrical data */
    const cs_real_t dz = 1.e-2;

    /* get physical values */
    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
    const cs_real_t cp0 = cs_glob_fluid_properties->cp0;
    const cs_real_t lambda0 = cs_glob_fluid_properties->lambda0;

    /* Collect information about the case, geometrical and physical */
    cs_real_t min_T = 1.e15, max_T = 0.;
    cs_real_t min_X = 1.e15, max_X = 0.;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      min_T = fmin(f_b_temp[face_id],min_T);
      max_T = fmax(f_b_temp[face_id],max_T);

      min_X = fmin(b_face_cog[face_id][0],min_X);
      max_X = fmax(b_face_cog[face_id][0],max_X);
    }

    /* MPI communications */
    if (cs_glob_rank_id >= 0) {
      cs_parall_min(1, CS_REAL_TYPE, &min_T);
      cs_parall_max(1, CS_REAL_TYPE, &max_T);
      cs_parall_min(1, CS_REAL_TYPE, &min_X);
      cs_parall_max(1, CS_REAL_TYPE, &max_X);
    }

    /* Position of mid planes and length of the domain */
    const cs_real_t xy_midplane = 0.5*(max_X - min_X);
    const cs_real_t L = max_X - min_X;

    /* Compute Nusselt nb on the left vertical wall */
    cs_real_t Nu0 = 0.;
    const cs_zone_t *zone = cs_boundary_zone_by_name("hot_wall");
    for (cs_lnum_t elt_id = 0; elt_id < zone->n_elts; elt_id++) {
      cs_lnum_t face_id = zone->elt_ids[elt_id];
      cs_lnum_t c_id = b_face_cells[face_id];

      Nu0 += b_face_surf[face_id] * (cvar_temp[c_id] - coefap[face_id])
           / (cell_cen[c_id][0] - b_face_cog[face_id][0]);
    }

    /* MPI communications */
    if (cs_glob_rank_id >= 0)
      cs_parall_sum(1, CS_REAL_TYPE, &Nu0);

    Nu0 = Nu0 / ((min_T - max_T) * dz);

    /* Compute Nusselt nb and the max vel along x in the vertical mid-plane */
    cs_real_t Nu12 = 0., max_vel_x = -1., y_vel_x = -1.;

    cs_lnum_t  nelts = 0;
    cs_lnum_t *lstelt = NULL;

    BFT_MALLOC(lstelt, n_i_faces, cs_lnum_t);    
    cs_selector_get_i_face_list("(x > 0.499) and (x < 0.501)", &nelts, lstelt);

    for (cs_lnum_t elt_id = 0; elt_id < nelts; elt_id++) {
      cs_lnum_t face_id = lstelt[elt_id];
      cs_lnum_t c_id1 = i_face_cells[face_id][0];
      cs_lnum_t c_id2 = i_face_cells[face_id][1];

      /* factor is used to correct contributions from internal faces  */
      /* that belongs to two MPI parts (i.e. c_id in halo) */
      cs_real_t factor = 1.;
      if ((cs_glob_rank_id >= 0) && ((c_id1 > n_cells) || (c_id2 > n_cells)))
        factor = 0.5;

      /* Compute values of density, temperature and velocity at faces */
      cs_real_t rho_face = 0.5*(cvar_rho[c_id1] + cvar_rho[c_id2]);
      cs_real_t temp_face = 0.5*(cvar_temp[c_id1] + cvar_temp[c_id2]);
      /* flux has to be divided by the surface */
      cs_real_t vel_x_face = i_mass_flux[face_id]
                           / (rho_face * i_face_surf[face_id]);

      if (fabs(vel_x_face) > max_vel_x) {
        max_vel_x = fabs(vel_x_face);
        y_vel_x = cell_cen[c_id1][1];
      }

      Nu12 += factor * i_face_surf[face_id] *
            (ro0 * cp0 * vel_x_face * (temp_face - min_T) / lambda0
            - (cvar_temp[c_id2] - cvar_temp[c_id1]) / 
              (cell_cen[c_id2][0]- cell_cen[c_id1][0]));
    }

    /* MPI communications */
    if (cs_glob_rank_id >= 0) {
      cs_lnum_t rank_id_max;
      cs_lnum_t node = 0; 
      cs_parall_sum(1, CS_REAL_TYPE, &Nu12);
      /* rank id of the max is required to cast the coordinate */
      cs_parall_min_id_rank_r(&node, &rank_id_max, -max_vel_x);
      cs_parall_max(1, CS_REAL_TYPE, &max_vel_x);
      cs_parall_bcast(rank_id_max, 1, CS_REAL_TYPE, &y_vel_x);
    }

    Nu12 = Nu12 / ((max_T - min_T) * dz);
    /* non dimensional velocity */
    max_vel_x = max_vel_x * ro0 * cp0 * L / lambda0;
    y_vel_x = y_vel_x / L;

    /* Compute the max vel along y in the horizontal mid-plane */
    cs_real_t max_vel_y = -1., x_vel_y = -1.;

    cs_selector_get_i_face_list("(y > 0.499) and (y < 0.501)", &nelts, lstelt);

    for (cs_lnum_t elt_id = 0; elt_id < nelts; elt_id++) {
      cs_lnum_t face_id = lstelt[elt_id];
      cs_lnum_t c_id1 = i_face_cells[face_id][0];
      cs_lnum_t c_id2 = i_face_cells[face_id][1];

      if (((cell_cen[c_id1][1] < xy_midplane) && (cell_cen[c_id2][1] > xy_midplane)) ||
          ((cell_cen[c_id1][1] > xy_midplane) && (cell_cen[c_id2][1] < xy_midplane))) {
        /* Compute values of density, temperature and velocity at faces */
        cs_real_t rho_face = 0.5*(cvar_rho[c_id1] + cvar_rho[c_id2]);
        /* flux has to be divided by the surface */
        cs_real_t vel_y_face = i_mass_flux[face_id]
                             / (rho_face * i_face_surf[face_id]);

        if (fabs(vel_y_face) > max_vel_y) {
          max_vel_y = fabs(vel_y_face);
          x_vel_y = cell_cen[c_id1][0];
        }
      }
    }

    /* MPI communications */
    if (cs_glob_rank_id >= 0) {
      cs_lnum_t rank_id_max;
      cs_lnum_t node = 0;
      /* rank id of the max is required to cast the coordinate */
      /* face_id is not usefull in this case so it is set to 0 */
      cs_parall_min_id_rank_r(&node, &rank_id_max, -max_vel_y);
      cs_parall_max(1, CS_REAL_TYPE, &max_vel_y);
      cs_parall_bcast(rank_id_max, 1, CS_REAL_TYPE, &x_vel_y);
    }

    /* non dimensional velocity */
    max_vel_y = max_vel_y * ro0 * cp0 * L / lambda0;
    x_vel_y = x_vel_y / L;

    bft_printf(" Nu0 %f Nu1/2 %f Ux_max %f y_Ux_max %f Uy_max %f x_Uy_max %f\n",
               Nu0, Nu12, max_vel_x, y_vel_x, max_vel_y, x_vel_y);

    BFT_FREE(lstelt);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
