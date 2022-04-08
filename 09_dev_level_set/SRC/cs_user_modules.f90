!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.2-alpha
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> See \subpage cs_user_modules for examples.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module user_module

  !=============================================================================

  implicit none

  !=============================================================================

  !> -- Définition du domaine géométrique et temporel --

  ! Nombre de volume de contrôle par dimension
  integer   nligne, ncolonne

  parameter(nligne = 128, ncolonne = 64)!FIXME remove

  double precision Lx, Ly, dx, dy, tf

  parameter(Lx = 1d0, Ly = 2d0)
  parameter(dx = Lx/dble(ncolonne), dy = Ly/dble(nligne))
  parameter(tf = 3d0)

  !=============================================================================

  !> -- Définition des paramètres du transport de Level Set --

  ! Nombre de courant et maximums des vitesses (stabilité numérique)
  double precision xcfl, maxu, maxv

  parameter (xcfl = 0.1d0, maxu = 0.5d0, maxv = 0.5d0)

  !=============================================================================

  !> -- Définition des paramètres physiques des deux phases --

  ! Masses volumiques et viscosités
  double precision ro_neg, ro_pos, mu_neg, mu_pos

  parameter(ro_neg = 100.d0, ro_pos = 1000.d0)
  parameter(mu_neg = 1.d0, mu_pos = 10.d0)

  !=============================================================================

  !> -- Définition des paramètres de l'interface --

  ! Coefficient de lissage
  double precision alpha

  parameter(alpha=1.5d0)

  !=============================================================================

  !> -- Définition des fonctions particulières

  ! Fonction Heaviside
  double precision heaviside (nligne*ncolonne)!FIXME allocate a field

  ! Fonctioin Dirac
  double precision dirac (nligne*ncolonne)!FIXME allocate a field


  !=============================================================================

end module user_module
