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
 * \brief Modification of the calculation of the particle relaxation time
 *  with respect to the chosen formulation for the drag coefficient
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 *      \f$\tau_c = \frac{m_p{C_p}_p}{PId_p^2h_e}\f$
 *
 *      \f$\tau_c\f$  : Thermal relaxation time (value to be computed)
 *
 *      \f$m_p\f$     : Particle mass
 *
 *      \f${C_p}_p\f$ : Particle specific heat
 *
 *      \f$d_p\f$     : Particle diameter
 *
 *      \f$h_e\f$     : Coefficient of thermal exchange
 *
 *  The coefficient of thermal exchange is calculated from a Nusselt number,
 *  itself evaluated by a correlation (Ranz-Marshall by default)
 *
 *      \f$\nu = \frac{h_ed_p}{\lambda} = 2 + 0.55{\Re_e}_p^{0.5}P_{rt}^{0.33}\f$
 *
 *      \f$\lambda\f$ : Thermal conductivity of the carrier field
 *
 *      \f${\Re_e}_p\f$     : Particle Reynolds number
 *
 *     \f$ P_{rt}\f$    : Prandtl number
 *
 * \param[in]   id_p   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[out]  taup   thermal relaxation time
 * \param[in]   dt     time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_rt(cs_lnum_t        id_p,
                cs_real_t        re_p,
                cs_real_t        uvwr,
                cs_real_t        rho_f,
                cs_real_t        rho_p,
                cs_real_t        nu_f,
                cs_real_t        taup[],
                const cs_real_t  dt[])
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  unsigned char *particle = p_set->p_buffer + p_am->extents * id_p;
  cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_DIAMETER);

  /*==========================================================================
   * Relaxation time with Haider and Levenspiel drag coefficient
   * (Drag coefficient and terminal velocity of spherical and
   * nonspherical particles,
   * Powder technology 58 (1989))
   *==========================================================================*/

  /* Sphericity parameter */

  cs_real_t phi = 0.9;

  /* Intermediate calculation */

  cs_real_t cte1 = exp(2.3288-6.4581*phi+2.4486*(phi*phi));
  cs_real_t cte2 = 0.0964+0.5565*phi;
  cs_real_t cte3 = exp(4.9050-13.8944*phi+18.4222*(phi*phi)
                       -10.2599*(phi*phi*phi));
  cs_real_t cte4 = exp(1.4681+12.2584*phi-20.7322*(phi*phi)
                       +15.8855*(phi*phi*phi));

  /* Drag coefficient (Haider and Levenspiel) */
  cs_real_t cd = 24.0/re_p*(1+cte1*pow(re_p,cte2))+cte3*re_p/(cte4+re_p);

  /* Intermediate calculation 2 */
  cs_real_t fdr = 3.0*cd*uvwr/p_diam;

  taup[id_p] = rho_p / rho_f / fdr;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
