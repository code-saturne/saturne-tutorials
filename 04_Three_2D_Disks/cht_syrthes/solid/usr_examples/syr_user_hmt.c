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
#include "syr_hmt_bd.h"
#include "syr_abs.h"
#include "syr_option.h"
#include "syr_const.h"
#include "syr_hmt_libmat.h"
#include "syr_parall.h"
#include "syr_utilitaire_parall.h"

#include "syr_user_hmt.h"
#include "syr_contact.h"

/* Dans les fonctions qui suivent, on pourra acceder aux grandeurs suivantes :
/    - coords[j][i] : coordonnee j du noeud i
/    - nodes[j][i]  : noeud j de l'element i
/    - tmps[i]      : temperature au noeud i
/    - pv[i]        : pression de vapeur au noeud i
/    - pt[i]        : pression totale au noeud i
/    - nrefe[i]     : reference de l'element i
/
/    - maillnodes.ndim  : dimension du probleme
/    - maillnodes.nelem : nombre d'elements du maillage
/    - maillnodes.npoin : nombre de noeud du maillage
/
/
/  Utilisation d'un fichier meteo
/       1ere ligne = nbre de variables (colonnes)
/         meteo.var[nbvar][nelem]
*/

extern int user_stop;     /* pour generer un arret premature mais propre du code */
extern int nbVar;

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP,  C. PENIGUEL                                     |
  |======================================================================|
  | Affectation des materiaux dans le cas du bati                        |
  |======================================================================| */
void user_hmt_affectmat(struct Maillage maillnodes,struct Humid humid)
{
  int i,j;

  /*  Voir la liste des materiaux disponibles dans hmt_libmat.h */

  /* initialisation par defaut au beton */
  /* Exemple 1 : cas ou il n'y a que du beton */
/*   for (i=0;i<maillnodes.nelem;i++)  */
/*     { */
/*       humid.mat[i]=MAT_BETON; */
/*     } */




  /* Exemple 2 : cas avec 3 materiaux differents */
  /*  for (i=0;i<maillnodes.nelem;i++)
      {
      if (maillnodes.nrefe[i]==1)
        humid.mat[i]=MAT_BETON;
      else if (maillnodes.nrefe[i]==2)
        humid.mat[i]=MAT_PLATRE;
      else if (maillnodes.nrefe[i]==3)
        humid.mat[i]=MAT_BOIS_PIN;
      }
  */


}


/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Programmation des conditions initiales par l'utilisateur dans le cas |
  | des transferts couples                                               |
  |======================================================================| */
void user_hmt_cini(struct Maillage maillnodes,
                   double *tmps,double *pv,double *pt,struct Humid humid,
                   struct PasDeTemps *pasdetemps,
                   struct Meteo meteo,struct Myfile myfile)
{
  int i,j;

  /* si on souhaite programmer ses propres valeurs,
  /   mettre la variable mescondinit 1, puis fournir les formules */
  int mescondinit=0;

  if (mescondinit)
    {
      /* initialisation de la temperature */
      for (i=0;i<maillnodes.npoin;i++)
        {
          tmps[i]=20.;
        }

      /* initialisation de PV */
      for (i=0;i<maillnodes.npoin;i++)
        {
          pv[i]=2316.;
        }

      /* initialisation de PT */
      for (i=0;i<maillnodes.npoin;i++)
        {
          pt[i]=101300.;
        }


    } /* fin de mescondinit */

}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Conditions aux limites pour les transferts couples            |
  |======================================================================| */
void user_hmt_limfso(struct Maillage maillnodes,struct MaillageCL maillnodeus,
                     double *t,double *pv,double *pt,
                     struct HmtClimhhh hmtclimhhh,
                     struct Clim rayinf,struct PasDeTemps *pasdetemps,
                     struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nr,np;
  double x,y,z, xt,xpv,xpt,ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la(les) variables ci dessous a 1,
  /   puis fournir les formules */
  int meshhh=0;

  /* le temps physique (en s) est donne par 'tempss' */

  /* Conditions d'echange sur T, PV et PT */
  /* ------------------------------------ */
  if (meshhh)
    {
      /* current time (s) */
      ct=pasdetemps->tempss;


      /* pour chaque noeud de chaque element de bord portant une condition hhh */
      for (j=0;j<hmtclimhhh.ndmat;j++)
        for (i=0;i<hmtclimhhh.nelem;i++)
          {
            nr = maillnodeus.nrefe[hmtclimhhh.numf[i]];             /* reference de l'element de bord */
            np = maillnodeus.node[j][hmtclimhhh.numf[i]];           /* numero global du noeud courant */

            x = maillnodes.coord[0][np];                            /* coordonnees du noeud courant   */
            y = maillnodes.coord[1][np];
            if (maillnodes.ndim==3) z = maillnodes.coord[2][np];

            xt  = t[np];              /* temperature du noeud courant   */
            xpv = pv[np];            /* pression de vapeur             */
            xpt = pt[np];            /* pression totale                */

            hmtclimhhh.t_ext[j][i]=0;   /* temperature exterieure */
            hmtclimhhh.t_h  [j][i]=0;   /* coefficient d'echangee */

            hmtclimhhh.pv_ext[j][i]=0;  /* pv exterieure          */
            hmtclimhhh.pv_h  [j][i]=0;  /* coefficient d'echangee */

            if (nbVar==3) /* pour le modele a 3 equations */
              {
                hmtclimhhh.pt_ext[j][i]=0;  /* pt exterieure          */
                hmtclimhhh.pt_h  [j][i]=0;  /* coefficient d'echangee */
              }
          }
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Flux volumiques sur T,PV,PT                                   |
  |======================================================================| */
void user_hmt_cfluvs(struct Maillage maillnodes,
                     double *t,  struct Cvol fluxvol_t,
                     double *pv, struct Cvol fluxvol_pv,
                     double *pt, struct Cvol fluxvol_pt,
                     struct PasDeTemps *pasdetemps,struct Humid humid,
                     struct Meteo meteo,struct Myfile myfile)
{
  int i,j,nrefe,ng,nr;
  double x,y,z,tmoy,ct;

  /* si on souhaite programmer ses propres formules, mettre la(les) */
  /* variable(s) ci dessous a 1, puis fournir les formules          */
  int mesfluxvol_t=0;
  int mesfluxvol_pv=0;
  int mesfluxvol_pt=0;

  /* current time (s) */
  ct=pasdetemps->tempss;


  /* Flux volumiques sur la temperature flux= val1 *T + val2 */
  /* ------------------------------------------------------ */
  if (mesfluxvol_t)
    {
      for (i=0;i<fluxvol_t.nelem;i++)
        {
          /* ng=fluxvol_t.nume[i]        */                               /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng]   */                                 /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,t,&nrefe,&x,&y,&z,&tmoy);*/   /* valeurs moy sur l'elt            */

          fluxvol_t.val1[i]=0;  /* W/m3/C */
          fluxvol_t.val2[i]=0;  /* W/m3   */
        }
    }


  /* Flux volumiques sur PV  */
  /* ----------------------- */
  if (mesfluxvol_pv)
    {
      for (i=0;i<fluxvol_pv.nelem;i++)
        {
          /* ng=fluxvol_pv.nume[i]        */                                /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng]   */                                   /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,pv,&nrefe,&x,&y,&z,&tmoy);*/    /* valeurs moy sur l'elt            */

          fluxvol_pv.val1[i]=0;
          fluxvol_pv.val2[i]=0;
        }
    }

  /* Flux volumiques sur PT  */
  /* ----------------------- */
  if (nbVar==3 && mesfluxvol_pt)
    {
      for (i=0;i<fluxvol_pt.nelem;i++)
        {
          /* ng=fluxvol_pt.nume[i]        */                                /* numero global de l'element       */
          /* nr=maillnodes.nrefe[ng]   */                                   /* numero de reference de l'element */
          /* data_element_moy(ng,maillnodes,pv,&nrefe,&x,&y,&z,&tmoy);*/    /* valeurs moy sur l'elt            */

          fluxvol_pt.val1[i]=0;
          fluxvol_pt.val2[i]=0;
        }
    }

}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |        Conditions aux limites - resistances de contact               |
  |======================================================================| */
void user_hmt_rescon(struct Maillage maillnodes,struct MaillageCL maillnodeus,
                     double *t,double *pv,double *pt,
                     double *tcor,double *pvcor,double *ptcor,
                     struct Contact rescon,struct PasDeTemps *pasdetemps,
                     struct SDparall sdparall)
{
  int i,j,nr,ne,num;
  double t1,t2,pv1,pv2,pt1,pt2,ct;

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
      prepare_paires_rc(maillnodes,rescon,pv,pvcor,sdparall);
      if (nbVar==3)
        prepare_paires_rc(maillnodes,rescon,pt,ptcor,sdparall);


      /* remplir les valeurs de RC ci-dessous */
      /* ************************************ */
      for (i=0;i<rescon.nelem;i++)     /* pour chaque facette avec RC */
        {
          ne=rescon.numf[i];           /* numero de la facette    */
          nr=maillnodeus.nrefe[ne];    /* reference de la facette */

          for (j=0;j<rescon.ndmat;j++) /* pour chacun des noeuds de la facette */
            {
              num=maillnodeus.node[j][i];   /* numero du noeud      */

              /* valeur de la RC sur T */
              /* --------------------- */
              t1=t[num];                    /* temperature du noeud */
              t2=tcor[num];                 /* temperature du noeud correspondant */

              rescon.g[ADR_T][j][i]=0;      /* resistance de contact a remplir  */


              /* valeur de la RC sur PV */
              /* --------------------- */
              pv1=pv[num];                    /* pression de vapeur du noeud */
              pv2=pvcor[num];                 /* pression de vapeur du noeud correspondant */

              rescon.g[ADR_PV][j][i]=0;      /* resistance de contact a remplir  */


              /* si modele a 3 equations */
              if (nbVar==3){
                /* valeur de la RC sur PT */
                /* --------------------- */
                pt1=pt[num];                    /* pression totale du noeud */
                pt2=ptcor[num];                 /* pression totale du noeud correspondant */

                rescon.g[ADR_PT][j][i]=0;      /* resistance de contact a remplir  */
              }
            }
        }

    }

}
