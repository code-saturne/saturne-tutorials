/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* code_saturne version 9.0 */

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

#include "base/cs_defs.h"

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
 * \file cs_user_extra_operations.cpp
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

  cs_probe_set_t *pset = cs_probe_set_get("probes");
  const cs_lnum_t *elt_ids = cs_probe_set_get_elt_ids(pset, CS_MESH_LOCATION_CELLS);

  if (elt_ids != NULL) {
    /* probe at output */
    cs_lnum_t iel0 = elt_ids[0];
    cs_real_t val_temp_out = CS_F_(t)->val[iel0];
    cs_notebook_parameter_set_value("OUTPUT", val_temp_out);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
