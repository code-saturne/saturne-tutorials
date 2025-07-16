/*-----------------------------------------------------------------------

                         SYRTHES version 5.0
                         -------------------

     This file is part of the SYRTHES Kernel, element of the
     thermal code SYRTHES.

     Copyright (C) 2009 EDF S.A., France

     contact: syrthes-support@edf.fr


     The SYRTHES Kernel is free software; you can redistribute it
     and/or modify it under the terms of the GNU General Public License
     as published by the Free Software Foundation; either version 2 of
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
# include <string.h>

#include "syr_usertype.h"
#include "syr_tree.h"
#include "syr_bd.h"
#include "syr_meteo.h"
#include "syr_hmt_bd.h"
#include "syr_parall.h"
#include "syr_abs.h"
#include "syr_option.h"
#include "syr_const.h"
#include "syr_fluid1d.h"
#include "syr_fluid0d.h"
#include "syr_cfd.h"

#include "syr_user.h"
#include "syr_generic.h"
#include "syr_utilitaire.h"
#include "syr_utilitaire_parall.h"
#include "syr_ecrire_fichier.h"

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
/ somme_int_parall(n) ou somme_double_parall(x) : sommer une variable
/ max_int_parall(n)   ou max_double_parall(x)   : calculer le max d'une variable
/ min_int_parall(n)   ou min_double_parall(x0   : calculer le min d'une variable
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
  | Lecture specifique d'un fichier utilisateur                          |
  |======================================================================| */
int user_read_myfile(struct Myfile *myfile)
{
  int i,j,n;
  double *p,t,fg,fd,te,tab[8];
  FILE *fmyfile;
  static int prem=1;
  int my_file;
  char *token;

  /* --------------------------------------------------- */
  /* -1- modify to activate reading of the user file      */
  /* --------------------------------------------------- */

  /* pour programmer la lecture d'un fichier personnel, mettre my_file=1;*/
  my_file=0;

  /* longueur maximale d'une ligne du fichier (nbre de caracteres) */
#define MYCHLONG 1000
  char ch[MYCHLONG];


  /* name of the file */
  char filename[CHLONG]="nom_de_fichier";

  /* values of the file will fill the array :           */
  /*     myfile->var[myfile->nbvar][myfile->nelem]      */
  /*     myfile->var[columns      ][rows         ]      */


  /* number of variables to read (number of columns) */
  myfile->nbvar=2;

  /* number of values per variable (number of rows) */
  myfile->nelem=10;

  /* ------------------------------------------------------- */
  /* -2- don't modify the code below                         */
  /* ------------------------------------------------------- */

  myfile->actif = my_file;

  if (prem){
    prem=0;
    return my_file;
  }

  if (my_file)
    {
      if ((fmyfile=fopen(filename,"r")) == NULL)
        {
          if (SYRTHES_LANG == FR)
            printf("user_lire_myfile : Impossible d'ouvrir le fichier personnel :%s\n",filename);
          else if (SYRTHES_LANG == EN)
            printf("user_lire_myfile : Impossible to open the personal file :%s\n",filename);
          syrthes_exit(1) ;
        }

      myfile->var=(double**)malloc(myfile->nbvar*sizeof(double*));
      for (n=0;n<myfile->nbvar;n++)
        myfile->var[n]=(double*)malloc(myfile->nelem*sizeof(double));

      verif_alloue_double2d(myfile->nbvar,"user_read_myfile",myfile->var);


  /* ------------------------------------------------------- */
  /* -3- modify to read the file                             */
  /* ------------------------------------------------------- */

      /* reading of the file - method 1 : automatic, nothing to do */
      /* reading of all the columns of the file                    */
      /* --------------------------------------------------------- */

      for (i=0;i<myfile->nelem;i++) // rows
        {
          fgets(ch,MYCHLONG,fmyfile);
          token = strtok(ch," \t");
          j=0; // columns
          while(token != NULL)
            {
              sscanf(token,"%lf",&(myfile->var[j][i]));
              token = strtok (NULL, " \t");
              j++;
            }
        }


      /* reading of the file - method 2 : to choose columns to be read */
      /* ------------------------------------------------------------- */
      /* to program ... */

      /* example : reading of a file with nelem rows, and 4 columns but only 2 required */
      /*           --> to get values of columns 1 and 4 only                                 */

      /* double NoNeed;      */
      /* for (i=0;i<myfile->nelem;i++){ */
      /*  fscanf(fmyfile,"%lf%lf%lf%lf",&(myfile->var[0][i]),&NoNeed,&NoNeed,&(myfile->var[1][i])); */
      /*  fgets(ch,CHLONG,fmyfile); */
      /*   }   */

  /* -------------------------------------------------------- */
  /* --- end                                                  */
  /* -------------------------------------------------------- */


      fclose(fmyfile);
    } /* fin if my_file */

  return my_file;

}


/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Ecriture de variables sur le fichier additionnel                     |
  |======================================================================| */
void user_add_var_in_file(struct Maillage maillnodes,
                          struct Cvol *fluxvol,struct Variable variable,
                          struct Prophy physol,struct PasDeTemps pasdetemps)
{
  int i;
  double *trav;

  /* l'ecriture sur fichier se fait en appelant add_var_in_file                    */
  /* Remarque : l'appel a add_var_in_file peut etre fait a                          */
  /*            n'importe quel endroit de la boucle en temps                        */
  /*            donc en particulier dans les autres sous-programmes                 */
  /*            utilisateurs                                                        */
  /* Dans la structure Variable on trouve :                                         */
  /*                variable.var[variable.adr_t] = temperature[maillnodes.npoins]       */
  /*                   => variable.var[variable.adr_t][i]= temperature du noeud i       */
  /*                variable.var[variable.adr_pv] = pv[maillnodes.npoins]               */
  /*                   => variable.var[variable.adr_pv][i]= pression vapeur du noeud i  */
  /*                variable.var[variable.adr_pt] = pt[maillnodes.npoins]               */
  /*                   => variable.var[variable.adr_pt][i]= pression totale du noeud i  */
  /* les parametres de la fonction d'ecriture sont les suivants :                   */
  /* add_var_in_file(nb,var,nomvar,idiscr)                                          */
  /*    - nb    = nbre valeurs du tableau (int)                                     */
  /*    - var   = variable (double*)                                                */
  /*    - nomvar= variable (char* de 12 carateres max)                              */
  /*    - idiscr= variable sur : 1->les elts de bord                                */
  /*                             2->les elements                                    */
  /*                             3-> les noeuds                                     */


  /* exemple de sorties de variables supplementaires        */
  /* on demande ici une ecriture sur fichier des proprietes */
  /* des materiaux                                          */

  int addvar=0;

  if (addvar)
    {
      /* exemple 1 : ecriture de rho,cp et k */
      add_var_in_file(physol.nelem,physol.rho,"RHO",2);
      add_var_in_file(physol.nelem,physol.cp,"CP",2);
      add_var_in_file(physol.kiso.nelem,physol.kiso.k,"K",2);

      /* exemple 2 : adimensionnalisation de la temperature */
      double tmin=1e6,tmax=-273.;
      for (i=0;i<maillnodes.npoin;i++) {
        tmin=min(tmin,variable.var[variable.adr_t][i]);       /* variable.var[variable.adr_t] --> variable "temperature" */
        tmax=max(tmax,variable.var[variable.adr_t][i]);
      }
      trav=(double*)malloc(maillnodes.npoin*sizeof(double));  /* alloue un tableau de longueur "npoin" */
      for (i=0;i<maillnodes.npoin;i++) trav[i]=(variable.var[variable.adr_t][i]-tmin)/(tmax-tmin);
      add_var_in_file(maillnodes.npoin,trav,"T_ADIM",3);
      free(trav);
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2009 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  |  Donnees de la transformation geometrique pour la periodicite        |
  |  dans le cas ou ce n'est ni une translation ni une rotation          |
  |======================================================================| */
void user_transfo_perio(int ndim,int num,
                        double x,double y,double z,
                        double *xt, double *yt, double *zt)
{
  int i;

  int matransfo=0;
  /* si on programme une transformation, mettre "matransfo=1"              */

  /*    a partir des coordonnees du point (x,y,z)                          */
  /*    et de la reference auquel il appartient (une des nb1 references    */
  /*    de la liste list_ref1) fournir, les coordonnees du point           */
  /*    apres transformation periodique (xt,yt,zt)                         */


  if (matransfo)
    {
      if (num==1)
        {
          *xt = x+4;
          *yt = y+4;
          if (ndim==3) *zt = z+4;
        }
      else if (num==2)
        {
          *xt = x+2;
          *yt = y+2;
          if (ndim==3) *zt = z+2;
        }
    }
}
/*|======================================================================|
  | SYRTHES 5.0                                       COPYRIGHT EDF 2019 |
  |======================================================================|
  | AUTEURS  : I. RUPP, C. PENIGUEL                                      |
  |======================================================================|
  | Modification des maillages volumiques et de bord                     |
  | fonction pour modifier manuellement le maillage : maillnodes         |
  |  et/ou le maillage de bord : maillnodebord                           |
  |======================================================================| */
void user_modif_mesh(struct Maillage *maillnodes,
                     struct MaillageBord *maillnodebord)
{
  int i,n;
  int nr,np,ne;
  double x,y,z;

  /* mettre mamodif=1 si on programme cette fonction */
  /* =============================================== */
  int mamodif=0;
  if (!mamodif) return;


  if (syrglob_nparts==1 ||syrglob_rang==0)
    {
      printf("\n -------------------------------------- \n");
      printf(" WARNING : mesh modified by user function \n");
      printf(" ---------------------------------------- \n");
    }


  /* exemple 1 : passage des coordonnees de mm em metre */
  if (1==0){
    for (i=0;i<maillnodes->ndim;i++)
      for (n=0;n<maillnodes->npoin;n++)
        maillnodes->coord[i][n] *= 0.001;
  }

  /* exemple 2 : changement de references de volume en fonction de la coordonnee z */
  if (1==0){
    for (n=0;n<maillnodes->nelem;n++)
      {
        if (maillnodes->nrefe[n] == 2)                  /* si l'element porte la reference 2 */
          {
            data_element_bary(n,*maillnodes,&x,&y,&z);  /* calcul des coord barycentriques de l'element */
            if (z>0.) maillnodes->nrefe[n] = 22;       /* si le z du barycentre de l'element est positif, on change la reference en 22 */
          }
      }
  }

  /* exemple 3 : changement d'une reference de bord en fonction de la coordonnee z */
  if (1==0){
    for (n=0;n<maillnodebord->nelem;n++)
      {
        if (maillnodebord->nrefe[n] == 5)                  /* si l'element porte la reference 5 */
          {
            data_elementBord_bary(n,*maillnodes,*maillnodebord,&x,&y,&z);  /* calcul des coord barycentriques de l'element */
            if (z>0.) maillnodebord->nrefe[n] = 7;        /* si le z du barycentre de l'element est positif, on change la reference en 7 */
          }
      }
  }

  /* pour ecrire le maillage modifie dans un fichier : décommenter ci-dessous en mettant le nom du fichier que vous souhaitez */
  char nom[300];
  if (syrglob_nparts>1)
     sprintf(nom,"PART/maillmodif_%05dpart%05d.syr",syrglob_nparts,syrglob_rang);
  else
     sprintf(nom,"maillmodif.syr");

  /* a decommenter pour ecrire le maillage modifieee         */
  /* ecrire_geom(nom,maillnodes,maillnodebord); */

}
