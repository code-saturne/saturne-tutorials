/*============================================================================
 * User functions for input of calculation parameters.
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
#include <string.h>

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
 * \file cs_user_parameters.cpp
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup([[maybe_unused]] cs_domain_t   *domain)
{
  cs_property_t *lambda = cs_property_by_name(CS_THERMAL_LAMBDA_NAME);

  cs_real_t lambda_disk1[3][3] = {{25, 0, 0}, {0, 25, 0}, {0, 0, 25}};
  cs_property_def_aniso_by_value(lambda, "disk1", lambda_disk1);
  cs_property_def_aniso_by_value(lambda, "sheath", lambda_disk1);

  cs_real_t lambda_disk2[3][3] = {{25, 0, 0}, {0, 5, 0}, {0, 0, 1}};
  cs_property_def_aniso_by_value(lambda, "disk2", lambda_disk2);

  // Rotation matrix (45 degrees aroound z axis)
  cs_real_t cos_t = cos(cs_math_pi/4.0);
  cs_real_t sin_t = sin(cs_math_pi/4.0);
  cs_real_t r[3][3] = {{cos_t, -sin_t, 0}, {sin_t, cos_t, 0}, {0, 0, 1}};
  cs_real_t lambda_disk3[3][3];

  // Tensor rotation: lambda_disk3 == lambda_disk3, rotated by R
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
       lambda_disk3[i][j] = 0.;
       for (int m = 0; m < 3; m++) {
         for (int n = 0; n < 3; n++) {
           lambda_disk3[i][j] += r[i][m]*r[j][n]*lambda_disk2[m][n];
         }
       }
    }
   }
  cs_property_def_aniso_by_value(lambda, "disk3", lambda_disk3);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
