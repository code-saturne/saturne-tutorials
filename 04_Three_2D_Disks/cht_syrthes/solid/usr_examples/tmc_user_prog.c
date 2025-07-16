/*-----------------------------------------------------------------------

                         SYRTHES version 6.0
                         -------------------

     This file is part of the SYRTHES Kernel, element of the
     thermal code SYRTHES.

     Copyright (C) 2019 EDF S.A., France

     contact: syrthes-support@edf.fr


     The SYRTHES Kernel is free software; you can redistribute it
     and/or modify it under the terms of the GNU General Public License
     as published by the Free Software Foundation; either version 3 of
     the License, or (at your option) any later version.

     The SYRTHES Kernel is distributed in the hope that it will be
     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.


     You should have received a copy of the GNU General Public License
     along with the SYRTHES Kernel; if not, write to the
     Free Software Foundation, Inc.,
     51 Franklin St, Fifth Floor,
     Boston, MA  02110-1301  USA

-----------------------------------------------------------------------*/

#include "tmc.h"
#include "tmc_user_prog.h"
#include <sdis.h>

#define UNKNOWN_TEMPERATURE -DBL_MAX

/* Dummy function definitions to avoid link errors */

res_T
tmc_user_prog_init(struct tmc* tmc)
{
  (void)tmc;
  /* TODO initialise global user prog context */
  return RES_OK;
}

void
tmc_user_prog_release(void)
{
  /* TODO release global user prog context */
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer le coefficient d'echange d'une condition d'échange en paroi  |
  |======================================================================| */
double
tmc_user_prog_interface_h
  (const struct sdis_interface_fragment* frag,
   struct tmc_user_prog_interface_data* data)
{
  double x, y, z, temps;
  int ref;
  double coef_ech = 0;
  ASSERT(frag && data);

  x = frag->P[0];
  y = frag->P[1];
  z = frag->P[2];
  temps = frag->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0 && y < 0 && z < 0) {
      coef_ech = 100;
    } else {
      coef_ech = 1000;
    }
  } else {
    coef_ech = 10;
  } */

  return coef_ech;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer le flux en paroi                                             |
  |======================================================================| */
double
tmc_user_prog_interface_flux
  (const struct sdis_interface_fragment* frag,
   struct tmc_user_prog_interface_data* data)
{
  double x, y, z, temps;
  int ref;
  double flux = SDIS_FLUX_NONE;
  ASSERT(frag && data);

  x = frag->P[0];
  y = frag->P[1];
  z = frag->P[2];
  temps = frag->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0 && y < 0 && z < 0) {
      flux = 100;
    } else {
      flux = 1000;
    }
  } else {
    flux = 10;
  } */

  return flux;
}
/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la temperature d'une condition de Dirichlet                  |
  |======================================================================| */
double
tmc_user_prog_interface_temperature
  (const struct sdis_interface_fragment* frag,
   struct tmc_user_prog_interface_data* data)
{
  double x, y, z, temps;
  int ref;
  double temperature = UNKNOWN_TEMPERATURE;
  ASSERT(frag && data);

  x = frag->P[0];
  y = frag->P[1];
  z = frag->P[2];
  temps = frag->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    temperature = 100;
  } else if(x < 0) {
    temperature = 10;
  } else if(y < 0) {
    temperature = 20;
  } else if(z < 0) {
    temperature = 30;
  } else {
    temperature = UNKNOWN_TEMPERATURE;
  } */

  if(temperature != UNKNOWN_TEMPERATURE)
    temperature = tmc_celsius_to_kelvin(temperature);

  return temperature;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Définir la temperature de référence                                  |
  |======================================================================| */
double
tmc_user_prog_interface_ray_tref
  (const struct sdis_interface_fragment* frag,
   struct tmc_user_prog_interface_data* data)
{
  double x, y, z, temps;
  int ref;
  double temperature = UNKNOWN_TEMPERATURE;
  ASSERT(frag && data);

  x = frag->P[0];
  y = frag->P[1];
  z = frag->P[2];
  temps = frag->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    temperature = 100;
  } else if(x < 0) {
    temperature = 10;
  } else if(y < 0) {
    temperature = 20;
  } else if(z < 0) {
    temperature = 30;
  } else {
    temperature = UNKNOWN_TEMPERATURE;
  } */

  if(temperature != UNKNOWN_TEMPERATURE)
    temperature = tmc_celsius_to_kelvin(temperature);

  return temperature;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la temperature d'une condition d'échange en paroi            |
  |======================================================================| */
double
tmc_user_prog_fluid_temperature
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double temperature = UNKNOWN_TEMPERATURE;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    temperature = 100;
  } else if(x < 0) {
    temperature = 10;
  } else if(y < 0) {
    temperature = 20;
  } else if(z < 0) {
    temperature = 30;
  } else {
    temperature = UNKNOWN_TEMPERATURE;
  } */

  if(temperature != UNKNOWN_TEMPERATURE)
    temperature = tmc_celsius_to_kelvin(temperature);

  return temperature;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer des flux volumiques (W/m3)                                   |
  |======================================================================| */
double
tmc_user_prog_volumic_power
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double volpow = SDIS_VOLUMIC_POWER_NONE;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      volpow = 10000;
    } else {
      volpow = 500;
    }
  } */

  return volpow;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la masse volumique (kg/m3)                                   |
  |======================================================================| */
double
tmc_user_prog_solid_rho
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double rho = 0;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      rho = 1000;
    } else {
      rho = 7000;
    }
  } */

  return rho;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la chaleur spécifique (J/kg/K)                               |
  |======================================================================| */
double
tmc_user_prog_solid_cp
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double cp = 0;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      cp = 460;
    } else {
      cp = 500;
    }
  } */

  return cp;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la conductivité thermique (W/m/K)                            |
  |======================================================================| */
double
tmc_user_prog_solid_lambda
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double lambda = 0;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      lambda = 5;
    } else {
      lambda = 25;
    }
  } */

  return lambda;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la conductivité thermique (W/m/K)                            |
  |======================================================================| */
double
tmc_user_prog_solid_init_temperature
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double temperature = 0;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  if(temps > 0) return UNKNOWN_TEMPERATURE;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      temperature = 5;
    } else {
      temperature = 25;
    }
  } */

  if(temperature != UNKNOWN_TEMPERATURE)
    temperature = tmc_celsius_to_kelvin(temperature);

  return temperature;
}

/*|======================================================================|
  | SYRTHES 6.0                                       COPYRIGHT EDF 2022 |
  |======================================================================|
  | AUTEURS  : MESO-STAR                                                 |
  |======================================================================|
  | Imposer la conductivité thermique (W/m/K)                            |
  |======================================================================| */
double
tmc_user_prog_delta_solid
  (const struct sdis_rwalk_vertex* vert,
   struct tmc_user_prog_medium_data* data)
{
  double x, y, z, temps;
  int ref;
  double delta = 0;
  ASSERT(vert && data);

  x = vert->P[0];
  y = vert->P[1];
  z = vert->P[2];
  temps = vert->time;
  ref = data->ref;

  /* Avoid unused variable warnings */
  (void)x, (void)y, (void)z, (void)temps, (void)ref;

  /* TODO Write the function body as for instance the following code
  if(ref == 6) {
    if(x < 0) {
      delta = 1.e-4;
    } else {
      delta = 1.e-5;
    }
  } */

  return delta;
}

