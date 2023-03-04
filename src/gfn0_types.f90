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
! along with gfn0.  If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfn0_types
!> Types for storing xTB method parametrization data
   use iso_fortran_env, only : wp=>real64
   implicit none
   private

   public :: TxTBData_mod,TxTBParameter
   public :: TRepulsionData, TCoulombData, THamiltonianData
   public :: TDispersionData,TShortRangeData
   public :: newData, getData
   public :: gfn0_partials

   interface newData
      module procedure :: newAtomicData
      module procedure :: newShellData
   end interface newData


   interface getData
      module procedure :: getAtomicData
      module procedure :: getShellData
   end interface getData

!========================================================================================!

   !> Data for the repulsion contribution
   type :: TRepulsionData

      !> Repulsion exponent for heavy elements
      real(wp) :: kExp

      !> Repulsion exponent for light elements
      real(wp) :: kExpLight

      !> Repulsion exponent
      real(wp) :: rExp

      !> Electronegativity scaling of repulsion
      real(wp) :: enScale

      !> Exponents of repulsion term
      real(wp), allocatable :: alpha(:)

      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:)

      !> Electronegativitity for scaling of repulsion
      real(wp), allocatable :: electronegativity(:)

      !> FIXME: real space cutoff should not be part of data
      real(wp) :: cutoff

   end type TRepulsionData

!========================================================================================!

   !> Data for the evaluation of the xTB core Hamiltonian
   type :: THamiltonianData

      !> Scaling factors for different interacting shells
      real(wp) :: kScale(0:3, 0:3)

      !> Scaling factor for diffuse or polarisation function
      real(wp) :: kDiff

      !> Shell dependence of the EN polynom
      real(wp) :: enScale(0:3, 0:3)

      !> Quartic contribution to EN polynom
      real(wp) :: enscale4

      !> Exponent for shell exponent weighting
      real(wp) :: wExp

      !> Principal quantum number of each shell
      integer, allocatable :: principalQuantumNumber(:, :)

      !> Angular momentum of each shell
      integer, allocatable :: angShell(:, :)

      !> Valence character of each shell
      integer, allocatable :: valenceShell(:, :)

      !> Number of primitives for expansion of Slater functions
      integer, allocatable :: numberOfPrimitives(:, :)

      !> Exponent of the Slater function
      real(wp), allocatable :: slaterExponent(:, :)

      !> Atomic level information
      real(wp), allocatable :: selfEnergy(:, :)

      !> Reference occupation of the atom
      real(wp), allocatable :: referenceOcc(:, :)

      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kCN(:, :)

      !> Electronegativity used in the shell polynomials
      real(wp), allocatable :: electronegativity(:)

      !> Atomic radii used in the shell polynomials
      real(wp), allocatable :: atomicRad(:)

      !> Shell polynomials to scale Hamiltonian elements
      real(wp), allocatable :: shellPoly(:, :)

      !> Pair parameters to scale Hamiltonian elements
      real(wp), allocatable :: pairParam(:, :)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQShell(:, :)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQAtom(:)

   end type THamiltonianData

!========================================================================================!

   !> Data for the evalutation of the Coulomb interactions
   type :: TCoulombData

      !> Use electronegativity equlibration for reference density
      logical :: enEquilibration

      !> Include second order electrostatics
      logical :: secondOrder

      !> Include third order electrostatics
      logical :: thirdOrder

      !> Third order electrostatics is shell resolved
      logical :: shellResolved

      !> Exponent of the generalized gamma function
      real(wp) :: gExp

      !> Atomic hardnesses used in second order electrostatics
      real(wp), allocatable :: chemicalHardness(:)

      !> Scaling factors for shell electrostatics
      real(wp), allocatable :: shellHardness(:, :)

      !> Third order Hubbard derivatives
      real(wp), allocatable :: thirdOrderAtom(:)

      !> Shell resolved third order Hubbard derivatives
      real(wp), allocatable :: thirdOrderShell(:, :)

      !> Charge widths for EEQ model
      real(wp), allocatable :: chargeWidth(:)

      !> Electronegativity for EEQ model
      real(wp), allocatable :: electronegativity(:)

      !> Coordination number dependence of the EN
      real(wp), allocatable :: kCN(:)

   end type TCoulombData

!========================================================================================!
   !> Data for the dispersion contribution
   type :: TDispersionData

      !> Damping parameters
      real(wp) :: s6  = 1.0_wp
      real(wp) :: s8  = 0.0_wp
      real(wp) :: s10 = 0.0_wp
      real(wp) :: a1  = 0.0_wp
      real(wp) :: a2  = 0.0_wp
      real(wp) :: s9  = 0.0_wp
      integer  :: alp = 16

      !> Weighting factor for Gaussian interpolation
      real(wp) :: wf

      !> Charge steepness
      real(wp) :: g_a

      !> Charge height
      real(wp) :: g_c

      !> Reference data for the dispersion
      integer, allocatable :: atoms(:)
      integer, allocatable :: nref(:)
      integer, allocatable :: ncount(:, :)
      real(wp), allocatable :: cn(:, :)
      real(wp), allocatable :: q(:, :)
      real(wp), allocatable :: alpha(:, :, :)
      real(wp), allocatable :: c6(:, :, :, :)

   end type TDispersionData

!========================================================================================!
   !> Short range basis correction
   type TShortRangeData

      !> Additional offset for the reference bond lengths
      real(wp) :: shift

      !> Scaling factor for the energy contribution
      real(wp) :: prefactor

      !> Steepness of the EN dependence
      real(wp) :: steepness

      !> Scaling factor for electronegativity differences
      real(wp) :: enScale

   end type TShortRangeData


!========================================================================================!

   !> Parametrisation data for the xTB method
   type :: TxTBData_mod

      !> Name of the parametrisation
      character(len=:), allocatable :: name

      !> Reference to the publication
      character(len=:), allocatable :: doi

      !> Internal version number
      integer :: level

      !> accuaracy selection
      real(wp) :: acc = 1.0_wp

      !> electronic temperature
      real(wp) :: etemp = 300.0_wp

      !> Number of shells
      integer, allocatable :: nShell(:)

      !> Parametrisation data for repulsive interactions
      type(TRepulsionData) :: repulsion

      !> Parametrisation data for core Hamiltonian
      type(THamiltonianData) :: hamiltonian

      !> Parametrisation data for Coulombic interactions
      type(TCoulombData) :: coulomb

      !> Parametrization for Dispersion interactions
      type(TDispersionData) :: dispersion

      !> Parametrisation data for the short range basis correction (optional)
      type(TShortRangeData), allocatable :: srb

   contains

      !> Write informative printout for the parametrisation data
      procedure :: writeInfo

   end type TxTBData_mod

!========================================================================================!

   !> global default parameter data
   type :: TxTBParameter
      real(wp) :: kshell(0:3) = 0.0_wp
      real(wp) :: ksp = 0.0_wp
      real(wp) :: ksd = 0.0_wp
      real(wp) :: kpd = 0.0_wp
      real(wp) :: kdiff = 0.0_wp
      real(wp) :: kdiffa = 0.0_wp
      real(wp) :: kdiffb = 0.0_wp
      real(wp) :: enshell(0:3) = 0.0_wp
      real(wp) :: enscale4 = 0.0_wp
      real(wp) :: cnshell(2, 0:3) = 0.0_wp
      real(wp) :: gam3shell(2, 0:3) = 0.0_wp
      real(wp) :: srbshift = 0.0_wp
      real(wp) :: srbpre = 0.0_wp
      real(wp) :: srbexp = 0.0_wp
      real(wp) :: srbken = 0.0_wp
      real(wp) :: wllscal = 0.0_wp
      real(wp) :: gscal = 0.0_wp
      real(wp) :: zcnf = 0.0_wp
      real(wp) :: tscal = 0.0_wp
      real(wp) :: kcn = 0.0_wp
      real(wp) :: fpol = 0.0_wp
      real(wp) :: ken = 0.0_wp
      real(wp) :: lshift = 0.0_wp
      real(wp) :: lshifta = 0.0_wp
      real(wp) :: split = 0.0_wp
      real(wp) :: zqf = 0.0_wp
      real(wp) :: alphaj = 0.0_wp
      real(wp) :: kexpo = 0.0_wp
      real(wp) :: dispa = 0.0_wp
      real(wp) :: dispb = 0.0_wp
      real(wp) :: dispc = 0.0_wp
      real(wp) :: dispatm = 0.0_wp
      real(wp) :: xbdamp = 0.0_wp
      real(wp) :: xbrad = 0.0_wp
      real(wp) :: aesshift = 0.0_wp
      real(wp) :: aesexp = 0.0_wp
      real(wp) :: aesrmax = 0.0_wp
      real(wp) :: aesdmp3 = 0.0_wp
      real(wp) :: aesdmp5 = 0.0_wp
      real(wp) :: ipeashift = 0.0_wp
      real(wp) :: renscale = 0.0_wp
   end type TxTBParameter

!========================================================================================!
   type :: gfn0_partials
    real(wp),allocatable :: cn(:)
    real(wp),allocatable :: dcndr(:,:,:)
    real(wp),allocatable :: qat(:)
    real(wp),allocatable :: dqdr(:,:,:)
    real(wp),allocatable :: selfEnergy(:,:)
    real(wp),allocatable :: dSEdcn(:,:)
    real(wp),allocatable :: dSEdq(:,:)
   end type gfn0_partials 

!========================================================================================!

contains
!========================================================================================!
!========================================================================================!

subroutine writeInfo(self, unit, num)

   !> Instance of the parametrisation data
   class(TxTBData_mod), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Atomic numbers
   integer, intent(in), optional :: num(:)

   character(len=*), parameter :: rnum = '(f12.6)'
   character(len=*), parameter :: offset = '(8x)'
   character(len=*), parameter :: head = '(6x,"*",1x,a,":")'
   character(len=*), parameter :: rfmt = '('//offset//',a,t36,'//rnum//')'
   character(len=*), parameter :: afmt = '('//offset//',a,t36,4x,a)'
   character(len=:), allocatable :: name
   integer :: ii, jj

   write(unit, '(a)')
   if (allocated(self%name)) then
      allocate(character(len=2*len(self%name)-1) :: name)
      name = repeat(' ', len(name))
      do ii = 1, len(self%name)
         jj = 2*ii-1
         name(jj:jj) = self%name(ii:ii)
      end do
   else
      name = repeat(' ', 10)
      write(name, '(i0)') self%level
      name = 'xTB level '//trim(name)
   end if
   call generic_header(unit, name, 49, 10)

   write(unit, '(a)')
   if (allocated(self%doi)) then
      write(unit, afmt) "Reference", self%doi
   end if

   write(unit, head) "Hamiltonian"
   write(unit, rfmt, advance='no') "H0-scaling (s, p, d)"
   do ii = 0, 2
      write(unit, rnum, advance='no') self%hamiltonian%kScale(ii, ii)
   end do
   write(unit, '(a)')
   write(unit, rfmt) "zeta-weighting", self%hamiltonian%wExp

   write(unit, head) "Dispersion"
   write(unit, rfmt) "s8", self%dispersion%s8
   write(unit, rfmt) "a1", self%dispersion%a1
   write(unit, rfmt) "a2", self%dispersion%a2
   write(unit, rfmt) "s9", self%dispersion%s9

   write(unit, head) "Repulsion"
   write(unit, rfmt, advance='no') "kExp", self%repulsion%kExp
   if (self%repulsion%kExpLight /= self%repulsion%kExp) then
      write(unit, rnum, advance='no') self%repulsion%kExpLight
   end if
   write(unit, '(a)')
   write(unit, rfmt) "rExp", self%repulsion%rExp

   write(unit, head) "Coulomb"
   write(unit, rfmt) "alpha", self%coulomb%gExp
   if (allocated(self%coulomb%thirdOrderShell)) then
      write(unit, afmt) "third order", "shell-resolved"
   else if (allocated(self%coulomb%thirdOrderAtom)) then
      write(unit, afmt) "third order", "atomic"
   else
      write(unit, afmt) "third order", "false"
   end if

   if (allocated(self%srb)) then
      write(unit, head) "Polar bond correction"
      write(unit, rfmt) "rad-shift", self%srb%shift
      write(unit, rfmt) "strength", self%srb%prefactor
      write(unit, rfmt) "en-exp", self%srb%steepness
      write(unit, rfmt) "en-scale", self%srb%enScale
   end if

   write(unit, '(a)')

contains
   subroutine generic_header(iunit,string,width,offset)
   implicit none
   integer,intent(in) :: iunit
   integer,intent(in) :: offset
   integer,intent(in) :: width
   character(len=*),intent(in) :: string
   character(len=width) :: dum1,dum2
   character(len=2*width) :: outstring
   character(len=width) :: formatstr
   integer :: strlen,ifront,iback
   strlen = len(string)
   ifront = (width - strlen)/2
   iback  = width - ifront - strlen
   write(dum1,*) width
   write(dum2,*) offset
   write(formatstr,'(i0,"x,a,",i0,"x")') ifront,iback
   write(outstring,'("|",'//formatstr//',"|")') string
   write(iunit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
   write(iunit,'('//dum2//'x,a)') trim(outstring)
   write(iunit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
   end subroutine generic_header
end subroutine writeInfo


subroutine newAtomicData(vec, num, data)

   real(wp), allocatable, intent(out) :: vec(:)

   integer, intent(in) :: num(:)

   real(wp), intent(in) :: data(:)

   allocate(vec(size(num)))
   call getAtomicData(vec, num, data)

end subroutine newAtomicData


subroutine getAtomicData(vec, num, data)

   real(wp), intent(out) :: vec(:)

   integer, intent(in) :: num(:)

   real(wp), intent(in) :: data(:)

   integer :: ii, izp

   do ii = 1, size(vec, dim=1)
      izp = num(ii)
      vec(ii) = data(izp)
   end do

end subroutine getAtomicData


subroutine newShellData(vec, num, nshell, data)

   real(wp), allocatable, intent(out) :: vec(:, :)

   integer, intent(in) :: num(:)

   integer, intent(in) :: nshell(:)

   real(wp), intent(in) :: data(:, :)

   allocate(vec(maxval(nshell), size(num)))
   call getShellData(vec, num, nshell, data)

end subroutine newShellData


subroutine getShellData(vec, num, nshell, data)

   real(wp), intent(out) :: vec(:, :)

   integer, intent(in) :: num(:)

   integer, intent(in) :: nshell(:)

   real(wp), intent(in) :: data(:, :)

   integer :: ii, ish, izp

   vec(:, :) = 0.0_wp
   do ii = 1, size(vec, dim=2)
      izp = num(ii)
      do ish = 1, nshell(izp)
         vec(ish, ii) = data(ish, izp)
      end do
   end do

end subroutine getShellData


end module gfn0_types
