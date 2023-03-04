!================================================================================!
! This file is part of gfn0.
!
! Copyright (C) 2022-2023 Philipp Pracht
!
! gfn0 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfn0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfn0. If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module cn_module
   use iso_fortran_env, only: wp=>real64
   implicit none
   private

   public :: cnType, getCoordinationNumber, cutCoordinationNumber


   interface getCoordinationNumber
   !   module procedure :: getCoordinationNumberWrap
      module procedure :: getCoordinationNumberLP
   end interface getCoordinationNumber


   !> Possible counting functions for calculating coordination numbers
   type :: TCNTypeEnum

      !> Counting function not specified
      integer :: invalid = 0

      !> Original DFT-D3 coordination number
      integer :: exp = 1

      !> Faster decaying error function CN, better for dense systems
      integer :: erf = 2

      !> Error function CN with covalency correction
      integer :: cov = 3

      !> Particular long-ranged version of the DFT-D3 coordination number
      integer :: gfn = 4

   end type TCNTypeEnum

   !> Enumerator for different coordination number types
   type(TCNTypeEnum), parameter :: cnType = TCNTypeEnum()


   abstract interface
   !> Abstract interface for the counting function (and its derivative)
   pure function countingFunction(k, r, r0)
      import :: wp

      !> Constant for counting function
      real(wp), intent(in) :: k

      !> Actual distance
      real(wp), intent(in) :: r

      !> Critical distance
      real(wp), intent(in) :: r0

      !> Value of the counting function in the range of [0,1]
      real(wp) :: countingFunction

   end function countingFunction
   end interface

!========================================================================================!
!> PARAMETER block

   !> convert bohr (a.u.) to Ångström and back
   real(wp),private,parameter :: autoaa = 0.52917726_wp
   real(wp),private,parameter :: aatoau = 1.0_wp/autoaa

   !> Pi
   real(wp),private,parameter :: pi=3.1415926535897932385_wp

   !> Parameter for electronegativity scaling
   real(wp),parameter :: k4=4.10451_wp

   !> Parameter for electronegativity scaling
   real(wp),parameter :: k5=19.08857_wp

   !> Parameter for electronegativity scaling
   real(wp),parameter :: k6=2*11.28174_wp**2


   !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
   !> 188-197), values for metals decreased by 10%.
   real(wp), parameter :: covalentRadD3(1:118) = [ &
      & 0.32_wp,0.46_wp, & ! H,He
      & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
      & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
      & 1.76_wp,1.54_wp, & ! K,Ca
      &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
      &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
      &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
      & 1.89_wp,1.67_wp, & ! Rb,Sr
      &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
      &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
      &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
      & 2.09_wp,1.76_wp, & ! Cs,Ba
      &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
      &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
      &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
      &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
      &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
      & 2.01_wp,1.81_wp, & ! Fr,Ra
      &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
      &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp] & ! Nh-Og
      & * aatoau * 4.0_wp / 3.0_wp

   !> Pauling electronegativities, used for the covalent coordination number.
   real(wp), parameter :: paulingEN(1:118) = [ &
      & 2.20_wp,3.00_wp, & ! H,He
      & 0.98_wp,1.57_wp,2.04_wp,2.55_wp,3.04_wp,3.44_wp,3.98_wp,4.50_wp, & ! Li-Ne
      & 0.93_wp,1.31_wp,1.61_wp,1.90_wp,2.19_wp,2.58_wp,3.16_wp,3.50_wp, & ! Na-Ar
      & 0.82_wp,1.00_wp, & ! K,Ca
      &                 1.36_wp,1.54_wp,1.63_wp,1.66_wp,1.55_wp, & ! Sc-
      &                 1.83_wp,1.88_wp,1.91_wp,1.90_wp,1.65_wp, & ! -Zn
      &                 1.81_wp,2.01_wp,2.18_wp,2.55_wp,2.96_wp,3.00_wp, & ! Ga-Kr
      & 0.82_wp,0.95_wp, & ! Rb,Sr
      &                 1.22_wp,1.33_wp,1.60_wp,2.16_wp,1.90_wp, & ! Y-
      &                 2.20_wp,2.28_wp,2.20_wp,1.93_wp,1.69_wp, & ! -Cd
      &                 1.78_wp,1.96_wp,2.05_wp,2.10_wp,2.66_wp,2.60_wp, & ! In-Xe
      & 0.79_wp,0.89_wp, & ! Cs,Ba
      &         1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp, & ! La-Eu
      &         1.20_wp,1.21_wp,1.22_wp,1.23_wp,1.24_wp,1.25_wp,1.26_wp, & ! Gd-Yb
      &                 1.27_wp,1.30_wp,1.50_wp,2.36_wp,1.90_wp, & ! Lu-
      &                 2.20_wp,2.20_wp,2.28_wp,2.54_wp,2.00_wp, & ! -Hg
      &                 1.62_wp,2.33_wp,2.02_wp,2.00_wp,2.20_wp,2.20_wp, & ! Tl-Rn
      ! only dummies below
      & 1.50_wp,1.50_wp, & ! Fr,Ra
      &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Ac-Am
      &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Cm-No
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-Cn
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp ] ! Nh-Og


!========================================================================================!
contains


!> Geometric fractional coordination number, supports both error function
!  and exponential counting functions.
subroutine getCoordinationNumberWrap(nat, at, xyz, cf, cn, dcndr, cutoff)

   !> Source for error creation
   character(len=*), parameter :: source = &
      & 'getCoordinationNumberWrap'

   !> Molecular structure information
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   !> Coordination number type (by counting function)
   integer, intent(in) :: cf

   !> Error function coordination number
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Real space cutoff for the coordination number
   real(wp), intent(in), optional :: cutoff

   logical :: exitRun
   real(wp) :: rCutoff

   if (present(cutoff)) then
      rCutoff = cutoff
   else
      rCutoff = 40.0_wp
   end if

   !> Actual call to the lattice point version of the CN evaluation
   call getCoordinationNumberLP(nat, at, xyz, rCutoff, cf, cn, dcndr)

end subroutine getCoordinationNumberWrap



!> Geometric fractional coordination number, supports both error function
!  and exponential counting functions.
subroutine getCoordinationNumberLP(nat, at, xyz, cutoff, cf, cn, dcndr)

   !> Molecular structure information
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Coordination number type (by counting function)
   integer, intent(in) :: cf

   !> Error function coordination number
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(out) :: dcndr(:, :, :)

   real(wp), parameter :: kcn_exp = 16.0_wp
   real(wp), parameter :: kcn_erf = 7.5_wp
   real(wp), parameter :: kcn_gfn = 10.0_wp

   select case(cf)
   case(cnType%exp)
      call ncoord(nat, at, xyz, cutoff, kcn_exp, expCount, dexpCount, &
         & .false., covalentRadD3, paulingEN, cn, dcndr)
   case(cnType%erf)
      call ncoord(nat, at, xyz, cutoff, kcn_erf, erfCount, derfCount, &
         & .false., covalentRadD3, paulingEN, cn, dcndr)
   case(cnType%cov)
      call ncoord(nat, at, xyz, cutoff, kcn_erf, erfCount, derfCount, &
         & .true., covalentRadD3, paulingEN, cn, dcndr)
   case(cnType%gfn)
      call ncoord(nat, at, xyz, cutoff, kcn_gfn, gfnCount, dgfnCount, &
         & .false., covalentRadD3, paulingEN, cn, dcndr )
   end select

end subroutine getCoordinationNumberLP


!> Actual implementation of the coordination number, takes a generic counting
!  function to return the respective CN.
subroutine ncoord(nat, at, xyz, cutoff, kcn, cfunc, dfunc, enscale, &
      & rcov, en, cn, dcndr)

   !> Molecular structure information
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Function implementing the counting function
   procedure(countingFunction) :: cfunc

   !> Function implementing the derivative of counting function w.r.t. distance
   procedure(countingFunction) :: dfunc

   !> Use a covalency criterium by Pauling EN's
   logical, intent(in) :: enscale

   !> Steepness of counting function
   real(wp), intent(in) :: kcn

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   integer :: iat, jat, ati, atj, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), stress(3, 3), den, cutoff2

   cn = 0.0_wp
   dcndr = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) private(den) shared(enscale, rcov, en)&
   !$omp reduction(+:cn, dcndr) shared(nat, at, xyz, kcn, cutoff2) &
   !$omp private(jat, itr, ati, atj, r2, rij, r1, rc, countf, countd)
   do iat = 1, nat
      ati = at(iat)
      do jat = 1, iat
         atj = at(jat)

         if (enscale) then
            den = k4*exp(-(abs(en(ati)-en(atj)) + k5)**2/k6)
         else
            den = 1.0_wp
         end if

         rij = xyz(:, iat) - xyz(:, jat)
         r2 = sum(rij**2)
         if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
         r1 = sqrt(r2)

         rc = rcov(ati) + rcov(atj)

         countf = den * cfunc(kcn, r1, rc)
         countd = den * dfunc(kcn, r1, rc) * rij/r1

         cn(iat) = cn(iat) + countf
         if (iat /= jat) then
            cn(jat) = cn(jat) + countf
         end if

         dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
         dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
         dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
         dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

      end do
   end do
   !$omp end parallel do

end subroutine ncoord


!> Error function counting function for coordination number contributions.
pure function erfCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

end function erfCount


!> Derivative of the counting function w.r.t. the distance.
pure function derfCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp), parameter :: sqrtpi = sqrt(pi)

   real(wp) :: count

   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

end function derfCount


!> Exponential counting function for coordination number contributions.
pure function expCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))

end function expCount


!> Derivative of the counting function w.r.t. the distance.
pure function dexpCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count
   real(wp) :: expterm

   expterm = exp(-k*(r0/r-1._wp))

   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))

end function dexpCount


!> Exponential counting function for coordination number contributions.
pure function gfnCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = expCount(k, r, r0) * expCount(2*k, r, r0+2)

end function gfnCount


!> Derivative of the counting function w.r.t. the distance.
pure function dgfnCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = dexpCount(k, r, r0) * expCount(2*k, r, r0+2) &
      &  + expCount(k, r, r0) * dexpCount(2*k, r, r0+2)

end function dgfnCount


!> Cutoff function for large coordination numbers
pure subroutine cutCoordinationNumber(nAtom, cn, dcndr, maxCN)

   !> number of atoms
   integer, intent(in) :: nAtom

   !> on input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(:)

   !> on input derivative of CN w.r.t. cartesian coordinates,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(:, :, :)

   !> maximum CN (not strictly obeyed)
   real(wp), intent(in), optional :: maxCN

   real(wp) :: cnmax
   integer :: iAt

   if (present(maxCN)) then
      cnmax = maxCN
   else
      cnmax = 4.5_wp
   end if

   if (cnmax <= 0.0_wp) return

   if (present(dcndr)) then
      do iAt = 1, nAtom
         dcndr(:, :, iAt) = dcndr(:, :, iAt) * dCutCN(cn(iAt), cnmax)
      end do
   end if

   do iAt = 1, nAtom
      cn(iAt) = cutCN(cn(iAt), cnmax)
   end do

end subroutine cutCoordinationNumber


!> Cutting function for the coordination number.
elemental function cutCN(cn, cut) result(cnp)

   !> Current coordination number.
   real(wp), intent(in) :: cn

   !> Cutoff for the CN, this is not the maximum value.
   real(wp), intent(in) :: cut

   !> Cuting function vlaue
   real(wp) :: cnp

   cnp = log(1.0_wp + exp(cut)) - log(1.0_wp + exp(cut - cn))

end function cutCN


!> Derivative of the cutting function w.r.t. coordination number
elemental function dCutCN(cn, cut) result(dcnpdcn)

   !> Current coordination number.
   real(wp), intent(in) :: cn

   !> Cutoff for the CN, this is not the maximum value.
   real(wp), intent(in) :: cut

   !> Derivative of the cutting function
   real(wp) :: dcnpdcn

   dcnpdcn = exp(cut)/(exp(cut) + exp(cn))

end function dCutCn


end module cn_module
