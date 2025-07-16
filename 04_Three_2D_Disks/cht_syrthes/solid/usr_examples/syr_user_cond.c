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

#include "syr_user_cond.h"
#include "syr_utilitaire.h"
#include "syr_utilitaire_parall.h"
#include "syr_contact.h"

#ifdef _SYRTHES_MPI_
#include "mpi.h"
static MPI_Status status;
#endif

extern int user_stop;     /* pour generer un arret prmature mais propre du code */
extern char nommeteo[CHLONG];
extern int nbVar;
extern FILE *syr_fadd;

/* Dans les fonctions qui suivent, on pourra acceder aux grandeurs suivantes :
/    - maillnodes.coords[j][i] : coordonnee j du noeud i
/    - maillnodes.nodes[j][i]  : noeud j de l'element i
/    - t[i]                    : temperature au noeud i
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
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |  Temperature initiale                                                |
  |======================================================================| */
void user_cini(struct Maillage maillnodes,double *t,struct PasDeTemps *pasdetemps,
               struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nr;
  double x,y,z,T;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mesproprietes a 1, puis fournir les formules */

  int mescondinit=0;

  if (mescondinit)
    {

      /* exemple 1 : on parcourt l'ensemble des noeuds du maillage */
      for (i=0;i<maillnodes.npoin;i++)
        {
          /*      x=maillnodes.coord[0][i]; */
          /*      y=maillnodes.coord[1][i]; */
          /*      if (maillnodes.ndim==3) z=maillnodes.coord[2][i]; */
          t[i]=20.;    /* temperature initiale */
        }

      /* exemple 2 : on parcourt tous les elements et suivant leur reference,  */
      /*             on impose une valeur initiale sur les noeuds de l'element */
      /* for (i=0;i<maillnodes.nelem;i++)
        {
          data_element_moy(i,maillnodes,t,&nr,&x,&y,&z,&T);
          if (nr== 6 || nr== 7 || nr== 8) {
            for (j=0;j<maillnodes.ndmat;j++)
              t[maillnodes.node[j][i]] = 22;
          }
        }
       */
      /* fin de l'exemple 2 */

    }
}


/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Proprietes physiques                                          |
  |======================================================================| */
void user_cphyso(struct Maillage maillnodes,
                 double *t,struct Prophy physol,struct PasDeTemps *pasdetemps,
                 struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nrefe,ng,nr;
  double x,y,z,tmoy,ct;
  static int prem;

  /* Attention : on fait ici des boucles sur les ELEMENTS               */
  /* Il faut utiliser la fonction 'data_element_moy' pour avoir acces   */
  /* aux coordonnees du barycentre de chaque element (x,y,z)            */
  /* et a la tyemperature moyenne de l'element (tmoy)                   */


  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mesproprietes a 1,
  /   puis fournir les formules */

  int mesproprietes=0;

  if (mesproprietes)
    {
      /* current time (s) */
      ct=pasdetemps->tempss;

      /* masse volumique et chaleur specifique */
      /*--------------------------------------*/
      for (i=0;i<physol.nelem;i++) /* boucle sur tous les elements */
        {
          /* si besoin, valeurs sur l'element (i est le numero globalde l'element) */
          /* data_element_moy(i,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);*/

          physol.rho[i]=7700;
          physol.cp[i]=460;
        }

      /* Conductivite */
      /*--------------*/
      /* Si conductivite constante, remplir le mot-cle CONDUCTIVITE VARIABLE= NON */
      /* conductivite des elements isotropes */
      for (i=0;i<physol.kiso.nelem;i++)
        {
          /* ng=physol.kiso.ele[i]    */                                    /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng] */                                     /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);*/  /* valeurs moy sur l'elt */
          physol.kiso.k[i]=25.;
        }

      /* conductivite des elements orthotropes */
      for (i=0;i<physol.kortho.nelem;i++)
        {
          /* ng=physol.kortho.ele[i]    */                                   /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng]   */                                    /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);*/   /* valeurs moy sur l'elt */
          physol.kortho.k11[i]=25.;
          physol.kortho.k22[i]=25.;
          if (maillnodes.ndim==3) physol.kortho.k33[i]=25.;
        }

      /* conductivite des elements anisotropes */
      for (i=0;i<physol.kaniso.nelem;i++)
        {
          /* ng=physol.kaniso.ele[i]    */                                   /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng]    */                                   /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);*/   /* valeurs moy sur l'elt */

          physol.kaniso.k11[i]=25.;
          physol.kaniso.k22[i]=25.;
          physol.kaniso.k12[i]=0.;
          if (maillnodes.ndim==3)
            {
              physol.kaniso.k33[i]=25.;
              physol.kaniso.k13[i]=0.;
              physol.kaniso.k23[i]=0.;
            }
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
void user_limfso(struct Maillage maillnodes,struct MaillageBord maillnodebord,
                 struct MaillageCL maillnodeus,
                 double *t,struct Clim diric,struct Clim flux,
                 struct Clim echang,struct Clim rayinf,struct PasDeTemps *pasdetemps,
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
  int mesflux=0;
  int mesdirichlet=0;
  int mescoeffech=0;
  int mesrayinf=0;


  /* current time (s) */
  ct=pasdetemps->tempss;


  /* Conditions de flux */
  /* ------------------ */
  if (mesflux)
    {
      /* pour chaque noeud de chaque element de bord de type flux */
      for (j=0;j<flux.ndmat;j++)
        for (i=0;i<flux.nelem;i++)
          {
            nr=maillnodeus.nrefe[flux.numf[i]];             /* reference de l'element de bord */
            np=maillnodeus.node[j][flux.numf[i]];           /* numero global du noeud courant */

            x=maillnodes.coord[0][np];                      /* coordonnees du noeud courant   */
            y=maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z=maillnodes.coord[2][np];

            xt=t[np];                                      /* temperature du noeud courant   */

            flux.val1[j][i]=0;                             /* valeur du flux a imposer */
          }
    }

  /* Conditions de dirichlet */
  /* ----------------------- */
  if (mesdirichlet)
    {
      /* pour chaque noeud de chaque element de bord de type dirichlet */
      for (j=0;j<diric.ndmat;j++)
        for (i=0;i<diric.nelem;i++)
          {
            nr=maillnodebord.nrefe[diric.numf[i]];             /* reference de l'element de bord  */
            np=maillnodebord.node[j][diric.numf[i]];            /* numero global du noeud courant */

            x=maillnodes.coord[0][np];                       /* coordonnees du noeud courant    */
            y=maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z=maillnodes.coord[2][np];

            xt=t[np];                                       /* temperature du noeud courant     */

            diric.val1[j][i]=0;                             /* valeur du dirichlet a imposer    */
          }
    }


  /* Conditions d'echange */
  /* -------------------- */
  if (mescoeffech)
    {
      /* pour chaque noeud de chaque element de bord de type echange */
      for (j=0;j<echang.ndmat;j++)
        for (i=0;i<echang.nelem;i++)
          {
            nr=maillnodeus.nrefe[echang.numf[i]];            /* reference de l'element de bord */
            np=maillnodeus.node[j][echang.numf[i]];          /* numero global du noeud courant */

            x=maillnodes.coord[0][np];                       /* coordonnees du noeud courant   */
            y=maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z=maillnodes.coord[2][np];

            xt=t[np];                                       /* temperature du noeud courant   */

            echang.val1[j][i]=0;  /* Temperature en C */
            echang.val2[j][i]=0;  /* coeff d'echange  */
          }
    }


  /* Conditions de rayonnement infini */
  /* -------------------------------- */
  if (mesrayinf)
    {
      for (j=0;j<rayinf.ndmat;j++)
        for (i=0;i<rayinf.nelem;i++)
          {
            nr=maillnodeus.nrefe[rayinf.numf[i]];            /* reference de l'element de bord */
            np=maillnodeus.node[j][rayinf.numf[i]];          /* numero global du noeud courant */

            x=maillnodes.coord[0][np];                       /* coordonnees du noeud courant   */
            y=maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z=maillnodes.coord[2][np];

            xt=t[np];                                        /* temperature du noeud courant   */

            rayinf.val1[j][i]=0;                             /* Temperature de l'infini en C   */
            rayinf.val2[j][i]=0;                             /* emissivite                     */
          }
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Flux volumiques sur la temperature                            |
  |======================================================================| */
void user_cfluvs(struct Maillage maillnodes,
                 double *t,struct Cvol fluxvol,struct PasDeTemps *pasdetemps,
                 struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nrefe,ng,nr;
  double x,y,z,tmoy,ct;

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


      /*---------------------------------------------------------------------------- */
      /* exemple 1 : imposer une puissance de 300 W sur les elements de reference 13 */
      /*---------------------------------------------------------------------------- */
      if (1==2){
        /* etape 1 : calculer le volume discretise correspondant a la zone */
        double vol, vol_tot;
        for (vol=0.,i=0;i<fluxvol.nelem;i++)
          {
            ng=fluxvol.nume[i];                                            /* numero global de l'element       */
            nr=maillnodes.nrefe[ng];                                       /* numero de reference de l'element */
            if (nr==13) vol += maillnodes.volume[ng];
          }
        vol_tot=somme_double_parall(vol);  /* volume total = somme sur les processeurs en cas de parallelisme */

        /* etape 2 : imposer le flux volumique correspondant a la puissance */
        for (i=0;i<fluxvol.nelem;i++)
          {
            ng=fluxvol.nume[i];                                            /* numero global de l'element       */
            nr=maillnodes.nrefe[ng];                                       /* numero de reference de l'element */
            if (nr==13) {
              fluxvol.val1[i] = 0;            /* W/m3/C  */
              fluxvol.val2[i] = 300/vol_tot;  /* W/m3    */
            }
          }
      }
      /*---------------------------------------------------------------------------- */
      /* exemple 1 : FIN                                                             */
      /*---------------------------------------------------------------------------- */




      /*---------------------------------------------------------------------------- */
      /* exemple 2 : imposer un flux volumique partout                               */
      /*---------------------------------------------------------------------------- */
      if (1==2){
        for (i=0;i<fluxvol.nelem;i++)
          {
            ng=fluxvol.nume[i];                                            /* numero global de l'element       */
            nr=maillnodes.nrefe[ng];                                       /* numero de reference de l'element */
            data_element_moy(ng,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);       /* valeurs moy sur l'elt            */

            /* flux vol du type = val1 * T + val2 */
            fluxvol.val1[i]=0;  /* W/m3/C  */
            fluxvol.val2[i]=0;  /* W/m3    */
          }
      }
      /*---------------------------------------------------------------------------- */
      /* exemple 2 : FIN                                                             */
      /*---------------------------------------------------------------------------- */
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Conditions aux limites - resistances de contact               |
  |======================================================================| */
void user_rescon(struct Maillage maillnodes,struct MaillageCL maillnodeus,
                 double *t,double *tcor,struct Contact rescon,struct PasDeTemps *pasdetemps,
                 struct SDparall sdparall)
{
  int i,j,nr,ne,num;
  double t1,t2,ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1, puis fournir les formules */
  int mesrescon=0;


  if (mesrescon==0)
    return;
  else
    {
      /* current time (s) */
      ct=pasdetemps->tempss;

      /* recuperation de la temperature pour les couples de RC - ne pas toucher */
      prepare_paires_rc(maillnodes,rescon,t,tcor,sdparall);


      for (i=0;i<rescon.nelem;i++)     /* pour chaque facette avec RC */
        {
          ne=rescon.numf[i];           /* numero de la facette    */
          nr=maillnodeus.nrefe[ne];    /* reference de la facette */

          for (j=0;j<rescon.ndmat;j++) /* pour chacun des noeuds de la facette */
            {
              num=maillnodeus.node[j][i];   /* numero du noeud      */
              t1=t[num];                    /* temperature du noeud */
              t2=tcor[num];                 /* temperature du noeud correspondant */

              rescon.g[ADR_T][j][i]=0;      /* resistance de contact a remplir  */
            }
        }

    }

}
