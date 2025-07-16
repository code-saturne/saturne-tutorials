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

#include "syr_fluid0d.h"
#include "syr_user_fluid0d.h"
#include "syr_utilitaire_parall.h"

#ifdef _SYRTHES_MPI_
#include "mpi.h"
static MPI_Status status;
#endif

extern int user_stop;     /* pour generer un arret prmature mais propre du code */
extern char nommeteo[CHLONG];
extern int nbVar;
extern FILE *syr_fadd;

/* Dans les fonctions qui suivent, on pourra acceder aux grandeurs suivantes :
/
/ fluid0d.cavite[i].Tf     temperature courante de la cavite i
/ fluid0d.cavite[i].rho    masse volumique du fluide de la cavite i
/ fluid0d.cavite[i].cp     chaleur specifique du fluide de la cavite i
/ fluid0d.cavite[i].volume volume de fluide dans la cavite i
/
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
/
/ cas des calculs paralleles
/ --------------------------
/ les fonctions suivantes peuvent etre appelees pour realiser des operations
/ sur l'ensemble des processeurs
/ s=somme_int_parall(n)  ou s=somme_double_parall(x) : sommer une variable
/ mx=max_int_parall(n)   ou mx=max_double_parall(x)   : calculer le max d'une variable
/ mn=min_int_parall(n)   ou mn=min_double_parall(x0   : calculer le min d'une variable
/
/ Rq : ces fontions fonctionnent egalement si le calcul est lance en sequentiel
/
/ Pour n'avoir les impressions qu'a 1 seul exemplaire dans le fichier listing,
/ utiliser le test suivant :
/ if (syrglob_nparts==1 ||syrglob_rang==0)
/    printf("mon impression\n");
*/



/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2015 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |  Temperature initiale    modele fluide 0D                            |
  |======================================================================| */
void user_fluid0d_cini(struct Fluid0d fluid0d,
                       struct Meteo meteo,struct Myfile myfile)
{
  int n;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mescondinit a 1, puis fournir les formules */

  int mescondinit=0;

  if (mescondinit)
    {
      /* on parcourt les cavites fluides */
      for (n=0;n<fluid0d.nb_cavite;n++)
        {
          fluid0d.cavite[n].Tini=20.;
        }
    }
}


/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2015 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Proprietes physiques      modele fluide 0D                    |
  |======================================================================| */
void user_fluid0d_cphy(struct Fluid0d fluid0d,
                       struct PasDeTemps *pasdetemps,
                       struct Meteo meteo,struct Myfile myfile)
{
  int n;
  double ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mesproprietes a 1,
  /   puis fournir les formules */
  int mesproprietes=0;

  if (mesproprietes)
    {
      /* current time (s) */
      ct=pasdetemps->tempss;
      /* on parcourt les cavites fluides */
      for (n=0;n<fluid0d.nb_cavite;n++)
        {
          fluid0d.cavite[n].rho=1.;
          fluid0d.cavite[n].cp=1017;
        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Conditions aux limites                                        |
  |======================================================================| */
void user_fluid0d_clim(struct Maillage maillnodes,struct MaillageBord maillnodebord,
                       struct MaillageCL maillnodeus,
                        double *t, struct Clim0D ech0D,
                        struct PasDeTemps *pasdetemps,
                        struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nr,np;
  double x,y,z,xt,ct;

  /* remarque : on fait ici des boucles sur les noeuds des faces de bord    */
  /* A partir du numero global du noeud (np) on accede directement          */
  /* a ses coordonnees (x,y,z) et a sa temperature (xt)                     */


  /* si on souhaite programmer ses propres formules,
  /   mettre la(les) variables ci dessous a 1,
  /   puis fournir les formules */
  int mescoeffech=0;


  /* current time (s) */
  ct=pasdetemps->tempss;


  /* Conditions d'echange */
  /* -------------------- */
  if (mescoeffech)
    {
      /* pour chaque noeud de chaque element de bord de type echange */
      for (j=0;j<ech0D.ndmat;j++)
        for (i=0;i<ech0D.nelem;i++)
          {
            nr=maillnodeus.nrefe[ech0D.numf[i]];            /* reference de l'element de bord */
            np=maillnodeus.node[j][ech0D.numf[i]];          /* numero global du noeud courant */

            x=maillnodes.coord[0][np];                       /* coordonnees du noeud courant   */
            y=maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z=maillnodes.coord[2][np];

            xt=t[np];                                       /* temperature du noeud courant   */

            ech0D.val2[j][i]=0;  /* coeff d'echange  */
          }
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2015 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Flux volumiques sur la temperature fluide 0D                  |
  |======================================================================| */
void user_fluid0d_fluv(struct Fluid0d fluid0d,
                       struct PasDeTemps *pasdetemps,
                       struct Meteo meteo,struct Myfile myfile)
{
  int n;
  double tf,ct;

  /* Attention : on fait ici des boucles sur les ELEMENTS               */
  /* Il faut utiliser la fonction 'data_element_moy' pour avoir acces   */
  /* aux coordonnees du barycentre de chaque element (x,y,z)            */
  /* et a la tyemperature moyenne de l'element (tmoy)                   */


/* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  int mesfluxvol=0;

  /* Conditions de flux volumiques */
  /* ----------------------------- */
  if (mesfluxvol==0)
    return;
  else
    {
      /* current time (s) */
      ct=pasdetemps->tempss;

      for (n=0;n<fluid0d.nb_cavite;n++)
        {
          tf=fluid0d.cavite[n].Tf; /* temperature fluide de la cavite */

          fluid0d.cavite[n].phiva=0.;
          fluid0d.cavite[n].phivb=0.;
        }

    }

}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2019 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Ventilation : transferts entre cavites                               |
  |======================================================================| */
void user_fluid0d_vent(struct Fluid0d fluid0d,
                       struct PasDeTemps *pasdetemps,
                       struct Meteo meteo,struct Myfile myfile)
{
  int i,j,numcav,numcav1,numcav2;
  double tf,ct;

  /* A noter : si besoin, on peut acceder aux donnees de la cavite    */
  /* fluid0d.Tf[i];  temperature courante de la cavite i              */
  /* fluid0d.rho[i]; masse volumique du fluide de la cavite i         */
  /* fluid0d.cp[i];  chaleur specifique du fluide de la cavite i      */


  /* si on souhaite programmer ses propres formules, mettre la variable ci dessous a 1 */
  /* puis fournir les formules                                                         */
  int vent=0;

  /* Conditions de flux volumiques */
  /* ----------------------------- */
  if (vent==0)
    return;
  else
    {
      /* !!! attention, dans les tableaux de debit et de temperature, les cavites sont numerotees de 0 a n-1 !!! */
      /* !!! flowrate value (in or out)  is always > 0                                                           */


      /* example : flowrate going into cavity 3, then going into cavity 5, then going into cavity 1, and go out (from cavity 1) */

      /* current time (s) */
      /* ct=pasdetemps->tempss; */

      /* ex : incoming flowrate (QmIn, kg/s) into cavity 3, inlet temperature (TmIn, deg C) */
      /* numcav=3; */
      /* fluid0d.cavite[numcav-1].QmIn = 0.7;   */
      /* fluid0d.cavite[numcav-1].TmIn = 25.; */

      /* ex : flowrate (kg/s) from cavity 3 to cavity 5 and then to cavity 1 */
      /* numcav1=3; numcav2=5;  fluid0d.Qm[numcav1-1][numcav2-1] = 0.7; */
      /* numcav1=5; numcav2=1;  fluid0d.Qm[numcav1-1][numcav2-1] = 0.7; */



      /* flowrate (kg/s) going out of cavity 1 */
      /* numcav=1;   fluid0d.cavite[numcav-1].QmOut = 0.7; */
    }

}
