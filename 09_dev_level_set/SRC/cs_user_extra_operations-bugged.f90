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

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          cpt, iel
double precision centre_masse, vit_ascension, xrtp

double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:,:), pointer :: cvar_vel

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

cpt = 0
centre_masse = 0d0
vit_ascension = 0d0

open(11, file="Surface_indicator", form="formatted")
open(12, file="Centre_of_Mass", form="formatted")
open(13, file="Vertical_velocity", form="formatted")

!===============================================================================
! User operations
!===============================================================================

call field_get_val_s_by_name("level_set", cvar_scal)
call field_get_val_v(ivarfl(iu), cvar_vel)

do iel = 1, ncel

  xrtp = cvar_scal(iel)

  ! Comptage de noeuds où la Level Set est négative
  if (xrtp.lt.0d0) then
    cpt = cpt + 1
    centre_masse = centre_masse + xyzcen(2,iel)
    vit_ascension = vit_ascension + cvar_vel(2,iel)
  endif
enddo

centre_masse = centre_masse / cpt
vit_ascension = vit_ascension / cpt

write(11,*) cpt
write(12,*) centre_masse
write(13,*) vit_ascension

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_f_user_extra_operations
