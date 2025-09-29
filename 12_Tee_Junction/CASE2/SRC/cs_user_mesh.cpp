/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between code_saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
 *============================================================================*/

/* VERS */

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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.cpp
 *
 * \brief Definition and modification of the calculation mesh.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  {
    /* Define extrusion parameters for each face */

    int n_zones = 1;
    const char *sel_criteria[] = {"walls"};
    const int zone_layers[] = {4};
    const double zone_thickness[] = {-1};
    const float zone_expansion[] = {0.8};

    cs_mesh_extrude_face_info_t *efi = cs_mesh_extrude_face_info_create(mesh);

    cs_lnum_t n_faces;
    cs_lnum_t *face_list;

    CS_MALLOC(face_list, mesh->n_b_faces, cs_lnum_t);

    for (int z_id = 0; z_id < n_zones; z_id++) {

      cs_selector_get_b_face_list(sel_criteria[z_id], &n_faces, face_list);

      cs_mesh_extrude_set_info_by_zone(efi,
                                       zone_layers[z_id],
                                       zone_thickness[z_id],
                                       zone_expansion[z_id],
                                       n_faces,
                                       face_list);
    }

    CS_FREE(face_list);

    /* Determine vertex values for extrusion */

    cs_mesh_extrude_vectors_t *e = cs_mesh_extrude_vectors_create(efi);

    /* Insert boundary layer */

    cs_mesh_extrude_face_info_destroy(&efi);
    cs_mesh_boundary_layer_insert(mesh, e, 0.2, false, 0, NULL);
    cs_mesh_extrude_vectors_destroy(&e);

  }

  mesh->modified |= CS_MESH_MODIFIED;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
