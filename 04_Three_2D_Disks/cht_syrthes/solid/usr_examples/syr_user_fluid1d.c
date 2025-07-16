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
#include "syr_parall.h"
#include "syr_abs.h"
#include "syr_option.h"
#include "syr_const.h"
#include "syr_fluid1d.h"

#include "syr_user_fluid1d.h"
#include "syr_utilitaire_parall.h"

#ifdef _SYRTHES_MPI_
#include "mpi.h"
static MPI_Status status;
#endif

extern int user_stop;     /* pour generer un arret prmature mais propre du code */
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
  | AUTEURS  : I. RUPP, C. PENIGUEL, D. HERMOUET (BT)                    |
  |======================================================================|
  |        Proprietes physiques                                          |
  |======================================================================| */
void user_fluid1d_cini(struct Maillage maillnodef1d,
                       double *tf,double *vf,double *pf,
                       double *dtf,double *tempsf,
                       struct Myfile myfile)
{
  int i,nr,np0,np1;
  double x,y,z;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mesproprietes a 1,
  /   puis fournir les formules */

  /* dtf = pas de temps fluide */
  /* tempsf = temps physique fluide */
  /* rq : pf[] non utilise pour le moment */

  int mescondinit=0;

  if (mescondinit)
    {
      /* temperature initiale */
      for (i=0;i<maillnodef1d.nelem;i++)
        {
          /*     nr=maillnodef1d.nrefe[i];      */   // reference de l'element fluide 1d
          /*     np0=maillnodef1d.node[0][i];   */   // numero des 2 noeuds de l'element
          /*     np1=maillnodef1d.node[1][i];   */

          /*      x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); */   // coordonnees du milieu de l'element fluide
          /*      y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]); */
          /*      z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]); */

          tf[i]=20.;            // temperature initiale
          vf[i]=0;              // vitesse  initiale
        }
    }
}


/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL , D. HERMOUET (BT)                   |
  |======================================================================|
  |        Proprietes physiques                                          |
  |======================================================================| */
void user_fluid1d_cphy(struct Maillage maillnodef1d, struct Prophyf1d phyf1d,
                       double *tf,double *vf,double *pfv,
                       double *dtf,double *tempsf,
                       struct Myfile myfile)
{
  int i,nr,np0,np1;
  double x,y,z,xt,ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable mesproprietes a 1,
  /   puis fournir les formules */
  /* rq : pfv[] non utilise pour le moment */

  int mesproprietes=0;

  if (mesproprietes)
    {

      /* masse volumique et chaleur specifique */
      /*--------------------------------------*/
      for (i=0;i<maillnodef1d.nelem;i++) /* boucle sur tous les elements */
        {
          nr=maillnodef1d.nrefe[i];          // reference de l'element
          np0=maillnodef1d.node[0][i];       // numero des 2 noeuds
          np1=maillnodef1d.node[1][i];
          x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); // coordonnees du barycentre de l'element
          y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]);
          z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]);

          //  note :  temperature=tf[i]   vitesse=vf[i]

          // phyf1d.rho[i]=851.773+1.28658*(tf[i]+273.15)-0.002*(tf[i]+273.15)*(tf[i]+273.15); /* Masse volumique [kg/m3] Eau */
          // phyf1d.rho[i]=1.293*273/(273+tf[i]);//Air
          // phyf1d.rho[i]=1000;/*Masse volumique [kg/m3]*/

          // phyf1d.cp[i]=4185; /*Chaleur specifique massique [J/kg.K]*/
          // phyf1d.lambda[i]=1;/*Conductivite thermique [W/m]*/
          // phyf1d.mu[i]=0.001;/*Viscocite dynamique [kg/m.s]*/
          // phyf1d.dh[i]=0.001;/*Diametre hydraulique[m]*/
          // phyf1d.section[i]=0.001;/*Section de passage [m^2]*/
          // phyf1d.source1[i]=1;/*Terme source1 [kg/s]*/
          // phyf1d.source2[i]=2;/*Terme source2 [kg/s]*/

        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL , D. HERMOUET (BT)                   |
  |======================================================================|
  |        Flux volumiques sur la temperature                            |
  |======================================================================| */
void user_fluid1d_fluv(struct Maillage maillnodef1d,struct Climvolf1d  cvolf1d,
                       double *tf,double *vf,double *pf,
                       double *dtf,double *tempf,struct Myfile myfile)
{
  int i,j,nr,np0,np1,nume;
  double x,y,z,T,P,V;

  // dtf = pas de temps fluide
  // tempsf = temps physique fluide

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  /* rq : pf[] non utilise pour le moment */

  int mesfluxvol=0;

  /* Conditions de flux volumiques */
  /* ----------------------------- */
  if (mesfluxvol==0)
    return;
  else
    {
      /* current time (s) */
      for (i=0;i<cvolf1d.nelem;i++)
        {
          nume=cvolf1d.nume[i];        // numero de l'element
          nr=maillnodef1d.nrefe[nume];          // reference de l'element
          np0=maillnodef1d.node[0][nume];       // numero des 2 noeuds
          np1=maillnodef1d.node[1][nume];
          x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); // coordonnees du barycentre de l'element
          y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]);
          z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]);
          // temperature=tf[i]   vitesse=vf[i]

          /*Flux volumique en AT+B [W/m^3]*/
          cvolf1d.val1[i]=0;
          cvolf1d.val2[i]=0;
        }
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL , D. HERMOUET (BT)                   |
  |======================================================================|
  |       Entrees Conditions aux limites                                 |
  |======================================================================| */
void user_fluid1d_entrees(struct Maillage maillnodef1d, struct Entreef1d entreef1d,
                          double *dtf,double *tempsf,struct Myfile myfile)
{
  int i;
  double x,y,z;

  /* si on souhaite programmer ses propres formules, mettre la(les) variables ci dessous a 1,
     puis fournir les formules.
     *dtf = pas de temps fluide
     *tempsf = temps physique fluide
 */

  int mesentrees=0;

  if (mesentrees)
    {
      for (i=0;i<entreef1d.nelem;i++)   /* pour chaque entree */
        {
          /* coordonnees du point proche de l'entree */
          x=entreef1d.x[i];
          y=entreef1d.y[i];
          z=entreef1d.z[i];

          if (entreef1d.type[i]==0)        /* cette fibre 1D est une boucle */
            {
              entreef1d.vx[i]=0;          /* vecteur sens du mouvement [x] */
              entreef1d.vy[i]=0;          /* vecteur sens du mouvement [y] */
              entreef1d.vz[i]=0;          /* vecteur sens du mouvement [z] */
              entreef1d.q[i]=0;           /* Debit en entree [kg/s]     */
            }
          else if (entreef1d.type[i]==1)   /* pour cette entree on donne le debit et la temperature */
            {
              entreef1d.q[i]=0;          /* Debit en entree [kg/s]     */
              entreef1d.t[i]=0;          /*  Temperature d'entree [°C] */
            }
          else if (entreef1d.type[i]==2)  /* pour cette entree on donne le deltaP et la temperature */
            {
              entreef1d.Pin[i]=0;        /* Pression a l'entree       */
              entreef1d.Pout[i]=0;       /* Pression a la sortie      */
              entreef1d.t[i]=0;          /* Temperature d'entree [°C] */
            }

        }
    }


}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL , D. HERMOUET (BT)                   |
  |======================================================================|
  |   Conditions aux limites                                             |
  |======================================================================| */
void user_fluid1d_clim(struct Maillage maillnodef1d, struct Climvolf1d  echf1d,struct Climvolf1d fimpf1d,
                       double *tf,double *vf,double *pf,double *dtf,double *tempsf,struct Myfile myfile)
{
  int i,j,nr,np,nume,np0,np1;
  double x,y,z,xt,ct;

  /* tf[i] = temperature fluide */
  /* vf[i]    = vitesse fluide */
  /* pf[i]    = pression fluide - non utilise pour le moment*/
  /* dtf   = pas de temps fluide */
  /* tempsf   =  temps fluide */

  /* si on souhaite programmer ses propres formules,
  /   mettre la(les) variables ci dessous a 1,
  /   puis fournir les formules */
  int mesech=0;
  int mesflux=0;


  /* Conditions d'echanghe */
  /* ------------------ */
  if (mesech)
    {
      /*pour chaque noeud de chaque element de bord de type flux*/
      for (i=0;i<echf1d.nelem;i++)
        {
          nume=echf1d.nume[i];   /* numero de l'elt fluide */
          nr=maillnodef1d.nrefe[nume];          // reference de l'element
          np0=maillnodef1d.node[0][nume];       // numero des 2 noeuds
          np1=maillnodef1d.node[1][nume];
          x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); // coordonnees du barycentre de l'element
          y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]);
          z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]);

          echf1d.val1[i]=50;   /*Temperature d'echange [°C]*/
          echf1d.val2[i]=10;   /*Coeffecient d'echange [W/m^2.K]*/
        }
    }
  /* Conditions de flux */
  /* ------------------ */
 if (mesflux)
    {
      /*pour chaque noeud de chaque element de bord de type flux*/
      for (i=0;i<fimpf1d.nelem;i++)
        {
          nume=fimpf1d.nume[i];   /* numero de l'elt fluide */
          nr=maillnodef1d.nrefe[nume];          // reference de l'element
          np0=maillnodef1d.node[0][nume];       // numero des 2 noeuds
          np1=maillnodef1d.node[1][nume];
          x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); // coordonnees du barycentre de l'element
          y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]);
          z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]);

          fimpf1d.val1[i]=50;/*Temperature d'echange [W]*/
        }
    }

}

/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL , D. HERMOUET (BT)                   |
  |======================================================================|
  |        Coefficient Echange entre Fluide et Solide                    |
  |======================================================================| */
void user_fluid1d_ech(struct Maillage maillnodef1d, struct Couplef1ds f1dcoups,
                      double *tf,double *vf,double *pf,
                      double *dtf,double *tempsf,struct Myfile myfile)
{
  int i,j,ne,nr,np0,np1;
  double x,y,z,T,V,P,ct;

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  /* rq : pf[] non utilise pour le moment */

  int mescoefech=0;

  /* Coeffcient d'echange */
  /* ----------------------------- */
  if (mescoefech==0)
    return;
  else
    {
      /* current time (s) */
      for (i=0;i<f1dcoups.nelem;i++)
        {
          if (f1dcoups.type[i]==1) /* si c'est pas colburn */
            {
              ne=f1dcoups.numglob[i];
              nr=maillnodef1d.nrefe[ne];
              np0=maillnodef1d.node[0][ne];       // numero des 2 noeuds
              np1=maillnodef1d.node[1][ne];
              x=0.5*(maillnodef1d.coord[0][np0]+maillnodef1d.coord[0][np1]); // coordonnees du barycentre de l'element
              y=0.5*(maillnodef1d.coord[1][np0]+maillnodef1d.coord[1][np1]);
              z=0.5*(maillnodef1d.coord[2][np0]+maillnodef1d.coord[2][np1]);
              T=tf[ne];
              V=vf[ne];

              f1dcoups.h[i]=200;/*Coefficent d'echange [W/m²]*/
            }
        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2015 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Coefficient de pertes de charge singulieres                          |
  |======================================================================| */
void user_fluid1d_pdcsing(int debut,int fin,struct Maillage maillnodef1d,struct Prophyf1d phyf1d,
                          double*vf,double q,double *dtf,double *tempsf,struct Myfile myfile)
{
  int i,j,nrefe,ng,nr;
  double x,y,z,tmoy,ct,a,b,re;

  /* tempf[i] = temperature fluide */
  /* vf[i]    = vitesse fluide */
  /* q        = debit fluide  */
  /* dtf      = pas de temps fluide */
  /* tempsf   = temps fluide */

  /* si on souhaite programmer ses propres formules,
  /   mettre la variable ci dessous a 1,
  /   puis fournir les formules */
  int mespdc=0;

  /* Coeffcient d'echange */
  /* ----------------------------- */
  if (mespdc==0)
    return;
  else
    {
      a=0.; b=1.;

      for (i=debut;i<fin;i++)
        {
          re=phyf1d.reynolds[i];
          phyf1d.pdcsing[i]=a*pow(re,-b);    /* valeur de la perte de charge */
        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2015 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Coefficient de pertes de charge lineaires                            |
  |======================================================================| */
void user_fluid1d_pdclin(int debut,int fin,struct Maillage maillnodef1d,struct Prophyf1d phyf1d,
                          double*vf,double q,double *dtf,double *tempsf,struct Myfile myfile)
{
  int i,j,nrefe,ng,nr;
  double x,y,z,tmoy,ct,a,b,re;

  /* tempf[i] = temperature fluide */
  /* vf[i]    = vitesse fluide */
  /* q        = debit fluide  */
  /* dtf      = pas de temps fluide */
  /* tempsf   = temps fluide */
  /* phyf1d.reynolds[i]   = reynolds */


  /* Coeffcient de perte de charge lineaire */
  /* par defaut c'est une loi de Colbrook en turbulent */

  /* si on souhaite programmer une autre loi, mettre  mespdc=1 */
  /* et programmer les formules ci-dessous                     */

  int mespdc=0;

  /* -------------------------------------- */
  if (mespdc==0)
    return;
  else
    {
      for (i=debut ; i <fin ; i++)
        {
          if (phyf1d.reynolds[i] < 2000)
            {
              /* On est en laminaire */
              phyf1d.pdclin[i] = 0;
            }
          else
            {
              /* en turbulent */
              phyf1d.pdclin[i]=0.;
            }
        }
    }
}
