/*-----------------------------------------------------------------------

                         SYRTHES version 5.0
                         -------------------

     This file is part of the SYRTHES Kernel, element of the
     thermal code SYRTHES.

     Copyright (C) 2009 EDF S.A., France

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

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

#include "syr_usertype.h"
#include "syr_tree.h"
#include "syr_bd.h"
#include "syr_meteo.h"
#include "syr_parall.h"
#include "syr_abs.h"
#include "syr_option.h"
#include "syr_const.h"
#include "syr_fluid0d.h"

#include "syr_user_ray.h"
#include "syr_init_radiation.h"
#include "syr_utilitaire_parall.h"

#ifdef _SYRTHES_MPI_
#include "mpi.h"
static MPI_Status status;
#endif

extern int user_stop;     /* pour generer un arret prmature mais propre du code */
extern char nommeteo[CHLONG];
extern int nbVar;


/* Dans les fonctions qui suivent, on pourra acceder aux grandeurs suivantes :
/    - maillnodes.coords[j][i] : coordonnee j du noeud i
/    - maillnodes.nodes[j][i]  : noeud j de l'element i
/    - tmps[i]                 : temperature au noeud i
/    - maillnodes.nrefe[i]     : reference de l'element i
/
/    - maillnodes.ndim  : dimension du probleme
/    - maillnodes.nelem : nombre d'elements du maillage
/    - maillnodes.npoin : nombre de noeud du maillage
/
/ fonction data_element_moy permet de recuperer les donnees moyennes sur l'element
/ -------------------------
/   --> data_element_moy(i,maillnodes,var,&nrefe,&x,&y,&z,&v);
/   avec en entree :
/            i =numero de l'element
/          var =champ sur les noeuds (en general temperature)
/   et en retour :
/        nrefe = numero de materiau de l'element,
/        x,y,z = coordonnees barycentriques de l'element
/            v = valeur moyenne de la variable sur l'element
/
/ fonction interpol_table1D permet d'interpoler dans une table de dimension 1
/ -------------------------
/ (TABLE (tabX) = tabFX)
/   --> y = interpol_table1D(double x,int nb,double *tabX,double *tabFX)
/   avec en entree :
/            nb = la longueur de la table
/          tabX = la liste des entrees de la table (tabX[nb])
/         tabFX = la liste des valeurs de la variable (tabFX[nb])
/   et en retour :
/             y = valeur interpolee
/
*/

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Proprietes physiques et conditions aux limites rayonnement    |
  |======================================================================| */
void user_ray(struct Maillage maillnodray,
              double *tmpray,struct ProphyRay phyray,
              struct PropInfini *propinf,struct Clim fimpray,struct Vitre vitre,
              struct PasDeTemps pasdetemps)
{
  int i,n,nr,ngfac;
  double ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  int mesproprietes=0;

  if (mesproprietes)
    {
      /* current time (s) */
      ct=pasdetemps.tempss;

      /* Definitions des proprietes physiques
       * Emissivite, reflectivite, absorptivite, transmittivite
       * ------------------------------------------------------  */
      for (i=0;i<maillnodray.nelem;i++)
        {
          nr=maillnodray.nrefe[i];  /* reference de la facette en cours */
          for (n=0;n<phyray.bandespec.nb;n++)
            {
              phyray.emissi[n][i]=1.;
              phyray.absorp[n][i]=1.;
              phyray.reflec[n][i]=0.;
              phyray.transm[n][i]=0.;
            }
        }

      if (vitre.actif && !vitre.prop_glob)
        proprvitre(vitre,maillnodray,phyray);


      /* Definitions de la temperature imposee */
      for (i=0;i<maillnodray.nelem;i++)
        tmpray[i]=20;


      /* Definitions des flux imposes par bande */
      for (i=0;i<fimpray.nelem;i++)
        for (n=0;n<phyray.bandespec.nb;n++)
          {
            ngfac=fimpray.numf[i];
            fimpray.val1[n][i]=0.;
            fimpray.val2[n][i]=0.;
          }


    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |  Calcul personnel des flux solaires directs et diffus                |
  |======================================================================| */
void user_solaire(double *fdirect,double *fdiffus, int nbande,
                  struct PasDeTemps pasdetemps)
{
  int i;
  double ct;

  /* si on souhaite programmer ses propres formules,
   * mettre la variable ci dessous a 1, puis fournir les formules.
   * Il faut ici fournir les flux solaires directs et diffus recus
   * par une surface horizontale                                    */

  int mesproprietes=0;

  if (mesproprietes)
    {
      /* current time (s) */
      ct=pasdetemps.tempss;

      for (i=0;i<nbande;i++)
        {
          fdirect[i]=0;
          fdiffus[i]=0;
        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |  Calcul personnel des proprietes physiques radiatives en fonction    |
  |  de l'angle d'incidence                                              |
  |======================================================================| */
void user_propincidence(int n,int nbande,double teta,
                        double **e,double **r,double **a,double **t,
                        double **emissi,double **reflec,
                        double **absorp,double **transm)
{
  int i;
  double tetadeg;
  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  int mesproprietes=0;

  if (mesproprietes)
    {
      tetadeg=teta*180/Pi;

      for (i=0;i<nbande;i++)
        {
          e[i][n]=emissi[i][n];
          r[i][n]=reflec[i][n];
          a[i][n]=absorp[i][n];
          t[i][n]=transm[i][n];
        }
    }

}
