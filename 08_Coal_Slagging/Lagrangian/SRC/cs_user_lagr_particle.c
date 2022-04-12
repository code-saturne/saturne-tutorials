/*============================================================================
 * Methods for lagrangian module
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_geom.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of an external force field acting on the particles.
 *
 * It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
 * By default gravity and drag force are the only forces acting on the particles
 * (the gravity components gx gy gz are assigned in the GUI or in usipsu)
 *
 * \param[in]     dt_p     time step (for the cell)
 * \param[in]     taup     particle relaxation time
 * \param[in]     tlag     relaxation time for the flow
 * \param[in]     piil     term in the integration of the sde
 * \param[in]     bx       characteristics of the turbulence
 * \param[in]     tsfext   infos for the return coupling
 * \param[in]     vagaus   Gaussian random variables
 * \param[in]     gradpr   pressure gradient
 * \param[in]     gradvf   gradient of the flow velocity
 * \param[in,out] rho_p     particle density
 * \param[out]    fextla   user external force field (m/s^2)$
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_ef(cs_real_t            dt_p,
                const cs_real_t      taup[],
                const cs_real_3_t    tlag[],
                const cs_real_3_t    piil[],
                const cs_real_33_t   bx[],
                const cs_real_t      tsfext[],
                const cs_real_33_t   vagaus[],
                const cs_real_3_t    gradpr[],
                const cs_real_33_t   gradvf[],
                cs_real_t            rho_p[],
                cs_real_3_t          fextla[])
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function (non-mandatory intervention)
 *
 * User-defined modifications on the variables at the end of the
 * Lagrangian time step and calculation of user-defined
 * additional statistics on the particles.
 *
 * About the user-defined additional statistics, we recall that:
 *
 *   isttio = 0 : unsteady Lagrangian calculation
 *          = 1 : steady Lagrangian calculation
 *
 *   isuist : Restart of statistics calculation if >= 1, else no stats
 *
 *   idstnt : Number of the time step for the start of the statistics calculation
 *
 *   nstist : Number of the time step of the start of the steady computation
 *
 * \param[in]  dt      time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_extra_operations(const cs_real_t  dt[])
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Impose the motion of a particle flagged CS_LAGR_PART_IMPOSED_MOTION.
 *
 * User-defined modifications on the particle position and its
 * velocity.
 *
 * \param[in]   coords    old particle coordinates
 * \param[in]   dt        time step (per particle)
 * \param[out]  disp      particle dispacement
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_imposed_motion(const cs_real_t  coords[3],
                            cs_real_t        dt,
                            cs_real_t        disp[3])
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particles.
 *
 * This function is called after the initialization of the new particles in
 * order to modify them according to new particle profiles (injection
 * profiles, position of the injection point, statistical weights,
 * correction of the diameter if the standard-deviation option is activated).
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modidied by this function.
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
cs_user_lagr_in(cs_lagr_particle_set_t         *particles,
                const cs_lagr_injection_set_t  *zis,
                const cs_lnum_t                 particle_range[2],
                const cs_lnum_t                 particle_face_id[],
                const cs_real_t                 visc_length[])
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the calculation of the particle relaxation time
 *  with respect to the chosen formulation for the drag coefficient
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 *                m   Cp
 *                 p    p
 *       Tau = ---------------
 *          c          2
 *                PI d    h
 *                    p    e
 *
 *      Tau  : Thermal relaxation time (value to be computed)
 *         c
 *
 *      m    : Particle mass
 *       p
 *
 *      Cp   : Particle specific heat
 *        p
 *
 *      d    : Particle diameter
 *       p
 *
 *      h    : Coefficient of thermal exchange
 *       e
 *
 *  he coefficient of thermal exchange is calculated from a Nusselt number,
 *  itself evaluated by a correlation (Ranz-Marshall by default)
 *
 *             h  d
 *              e  p
 *      Nu = --------  = 2 + 0.55 Re **(0.5) Prt**(0.33)
 *            Lambda                p
 *
 *      Lambda : Thermal conductivity of the carrier field
 *
 *      Re     : Particle Reynolds number
 *        p
 *
 *      Prt    : Prandtl number
 *
 * \param[in]   id_p   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[in]   cp_f   specific heat of the fluid at particle position
 * \param[in]   k_f    diffusion coefficient of the fluid at particle position
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
/*!
 * \brief Modification of the computation of the thermal relaxation time
 *        of the particles with respect to the chosen formulation of
 *        the Nusselt number.
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 * \param[in]   id_p   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[in]   cp_f   specific heat of the fluid at particle position
 * \param[in]   k_f    diffusion coefficient of the fluid at particle position
 * \param[out]  tauc   thermal relaxation time
 * \param[in]   dt     time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_rt_t(cs_lnum_t        id_p,
                  cs_real_t        re_p,
                  cs_real_t        uvwr,
                  cs_real_t        rho_f,
                  cs_real_t        rho_p,
                  cs_real_t        nu_f,
                  cs_real_t        cp_f,
                  cs_real_t        k_f,
                  cs_real_t        tauc[],
                  const cs_real_t  dt[])
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User integration of the SDE for the user-defined variables.
 *
 * The variables are constant by default. The SDE must be of the form:
 *
 * \f[
 *    \frac{dT}{dt}=\frac{T - PIP}{Tca}
 * \f]
 *
 * T:   particle attribute representing the variable
 * Tca: characteristic time for the sde
 *      to be prescribed in the array auxl1
 * PIP: coefficient of the SDE (pseudo RHS)
 *      to be prescribed in the array auxl2.
 *      If the chosen scheme is first order (nordre=1) then, at the first
 *      and only call pip is expressed as a function of the quantities of
 *      the previous time step (contained in the particle data).
 *      If the chosen scheme is second order (nordre=2)
 *      then, at the first call (nor=1) pip is expressed as a function of
 *      the quantities of the previous time step, and at the second passage
 *      (nor=2) pip is expressed as a function of the quantities of the
 *      current time step.
 *
 * \param[in]  dt      time step (per cell)
 * \param[in]  taup    particle relaxation time
 * \param[in]  tlag    relaxation time for the flow
 * \param[in]  tempct  characteristic thermal time and implicit source
 *                     term of return coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_sde(const cs_real_t  dt[],
                 cs_real_t        taup[],
                 cs_real_3_t      tlag[],
                 cs_real_t        tempct[])
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
