!-------------------------------------------------------------------------------

!                      VERS
!                      --------------------------
! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_user_initialization.f90
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

subroutine cs_user_f_initialization &
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
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use vof

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel, iutile
integer, allocatable, dimension(:) :: lstelt
double precision R0, X0, Y0, R, width, xrtp
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:), pointer :: heavyside, dirac

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(ncel)) ! temporary array for cells selection

!===============================================================================
! Variables initialization:
!
!   cvar_scal(iel) is the value of the variable in cell number iel.
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then
  call field_get_val_s_by_name("level_set", cvar_scal) !FIXME to bug "lev_set"
  R0 = 0.25d0
  X0 = 0.5d0
  Y0 = 0.5d0

  do iel = 1, ncel

    ! Initialisation de la variable scalaire level set

    R = sqrt((xyzcen(1,iel)-X0)**2+(xyzcen(2,iel)-Y0)**2)
    cvar_scal(iel) = R-R0

  enddo

  call field_get_val_s(icrom, crom)

  width = 1.5d0/64.d0

  call field_get_val_s_by_name("heavyside", heavyside)
  call field_get_val_s_by_name("dirac", dirac)

  do iel = 1, ncel

    xrtp = cvar_scal(iel)

    if (xrtp.lt.(-width)) then
      heavyside(iel) = 0.d0
      dirac(iel) = 0.d0
    else if (xrtp.gt.width) then
      heavyside(iel) = 1.d0
      dirac(iel) = 0.d0
    else
      heavyside(iel) = 0.5d0*(1.d0+xrtp/width+sin(pi*xrtp/width)/pi)
      dirac(iel) = 0.5d0/width*(1.d0+cos(pi*xrtp/width))
    endif

    ! rho1 is density of the liquid
    ! rho2 is density of the vapor
    crom(iel) = rho2 + (rho1 - rho2)*heavyside(iel)
  enddo
endif

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_f_initialization
