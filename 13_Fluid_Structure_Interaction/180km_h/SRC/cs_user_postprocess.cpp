/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, NULL otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts)
{
  CS_NO_WARN_IF_UNUSED(probes);
  CS_NO_WARN_IF_UNUSED(n_vertices);
  CS_NO_WARN_IF_UNUSED(vertex_list);

  /* Output of the vorticity
     ----------------------- */

  if (cat_id == CS_POST_MESH_VOLUME) { /* filter: only for volume
                                          postprocessing mesh */

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

    cs_real_33_t *gradv;
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    cs_field_gradient_vector(CS_F_(vel),
                             false,   /* use_previous_t, */
                             1,       /* inc */
                             gradv);


    cs_real_t *vorticity;
    BFT_MALLOC(vorticity, n_cells, cs_real_t);
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      cs_lnum_t c_id = cell_list[i];
      cs_real_t v[3] = {gradv[c_id][2][1] - gradv[c_id][1][2],
                        gradv[c_id][0][2] - gradv[c_id][2][0],
                        gradv[c_id][1][0] - gradv[c_id][0][1]};
      vorticity[i] = cs_math_3_square_norm(v);
    }

    BFT_FREE(gradv);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "vorticity",                    /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      true,                           /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      vorticity,                      /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    BFT_FREE(vorticity);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
