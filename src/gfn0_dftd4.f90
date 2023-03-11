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
!> Copyright (C) 2019-2020 Sebastian Ehlert, Stefan Grimme, Eike Caldeweyher
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfn0_dftd4
   use iso_fortran_env, only: wp=>real64,stdout=>output_unit
   use gfn0_types, only: TDispersionData
   use gfn0_math_wrapper, only: gemv,contract
   use dftd4param, gam=>chemical_hardness
   implicit none
   private

   public :: newD4Model
   public :: d4_gradient
  
   real(wp),private,parameter :: autoaa = 0.52917726_wp
   real(wp),private,parameter :: aatoau = 1.0_wp/autoaa
   real(wp),private,parameter :: pi=3.1415926535897932385_wp

contains

subroutine newD4Model(dispm,g_a,g_c,mode)
   type(TDispersionData), intent(inout) :: dispm
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode

   integer  :: i,ia,is,icn,j,ii,jj
   integer  :: cncount(0:18)
   real(wp) :: sec_al(23),iz,c6,alpha(23)
   real(wp) :: tmp_hq(7,118)

   intrinsic :: nint

   !call init(dispm) !< see "initGFN0Dispersion"

   secq = 0.0_wp
   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
!     print'(1x,''* using PBE0/def2-TZVP Hirshfeld charges'')'
      refq = dftq
      refh = dfth
      secq = dfts
!  case(2)
!     refq = pbcq
!     refh = pbch
!     secq = pbcs
   case(p_refq_gasteiger)
!     print'(1x,''* using classical Gasteiger charges'')'
      refq = gffq
      refh = gffh
      secq = gffs
   case(p_refq_goedecker)
      refq = clsq
      refh = clsh
      secq = clss
   case(p_refq_gfn2xtb_gbsa_h2o)
!     print'(1x,''* using GFN2-xTB//GBSA(H2O) charges'')'
      refq = solq
      refh = solh
      secq = sols
   end select

   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
      dispm%q = dftq
      tmp_hq = dfth
   case(p_refq_gasteiger)
      dispm%q = gffq
      tmp_hq = gffh
   case(p_refq_goedecker)
      dispm%q = clsq
      tmp_hq = clsh
   case(p_refq_gfn2xtb_gbsa_h2o)
      dispm%q = solq
      tmp_hq = solh
   case default
      dispm%q = refq
      tmp_hq = refh
   end select

   dispm%atoms = 0
   dispm%nref = 0
  
   do ia = 1, 118
      cncount = 0
      cncount(0) = 1
      dispm%nref(ia) = refn(ia)
      do j = 1, refn(ia)
         is = refsys(j,ia)
         iz = zeff(is)
         sec_al = sscale(is)*secaiw(:,is) &
            &  * zeta(g_a,gam(is)*g_c,secq(is)+iz,tmp_hq(j,ia)+iz)
         dispm%cn(j,ia) = refcovcn(j,ia)
         icn =nint(refcn(j,ia))
         cncount(icn) = cncount(icn) + 1
         dispm%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia)-hcount(j,ia)*sec_al),0.0_wp)
      enddo
      do j = 1, refn(ia)
         icn = cncount(nint(refcn(j,ia)))
         dispm%ncount(j,ia) = icn*(icn+1)/2
      enddo
   end do

   ! integrate C6 coefficients
   do i = 1, 118
      do j = 1, i
         do ii = 1, dispm%nref(i)
            do jj = 1, dispm%nref(j)
               alpha = dispm%alpha(:,ii,i)*dispm%alpha(:,jj,j)
               c6 = thopi * trapzd(alpha)
               dispm%c6(jj,ii,j,i) = c6
               dispm%c6(ii,jj,i,j) = c6
            enddo
         enddo
      enddo
   enddo

end subroutine newD4Model

subroutine d4dim(dispm,nat,at,ndim)
   type(TDispersionData), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   integer, intent(out) :: ndim

   integer :: i

   ndim = 0

   do i = 1, nat
      ndim = ndim + dispm%nref(at(i))
   enddo

end subroutine d4dim

subroutine prmolc6(molc6,molc8,molpol,nat,at,  &
      & cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   real(wp),intent(in)  :: molc6
   real(wp),intent(in)  :: molc8
   real(wp),intent(in)  :: molpol
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in),optional :: cn(nat)
   real(wp),intent(in),optional :: covcn(nat)
   real(wp),intent(in),optional :: q(nat)
   real(wp),intent(in),optional :: qlmom(3,nat)
   real(wp),intent(in),optional :: c6ab(nat,nat)
   real(wp),intent(in),optional :: alpha(nat)
   real(wp),intent(in),optional :: rvdw(nat)
   real(wp),intent(in),optional :: hvol(nat)
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(stdout,'(a)')
   write(stdout,'(''   #   Z   '')',advance='no')
   if(present(cn))   write(stdout,'(''        CN'')',advance='no')
   if(present(covcn))write(stdout,'(''     covCN'')',advance='no')
   if(present(q))    write(stdout,'(''         q'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(stdout,'(''      C6AA'')',advance='no')
   if(present(alpha))write(stdout,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(stdout,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(stdout,'(''    relVol'')',advance='no')
   write(*,'(a)')
   do i=1,nat
      write(*,'(i4,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),toSymbol(at(i))
      if(present(cn))   write(stdout,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(stdout,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(stdout,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(stdout,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(stdout,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(stdout,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(stdout,'(f10.3)',advance='no')hvol(i)
      write(*,'(a)')
   enddo
   endif
   write(stdout,'(/,1x,''Mol. C6AA /au·bohr⁶  :'',f18.6,'// &
   &         '/,1x,''Mol. C8AA /au·bohr⁸  :'',f18.6,'// &
   &         '/,1x,''Mol. α(0) /au        :'',f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

pure elemental function zeta(a,c,qref,qmod)
   !$acc routine seq
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: zeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1._wp - exp( c * ( 1._wp - qref/qmod ) ) ) )
   endif

end function zeta

pure elemental function dzeta(a,c,qref,qmod)
   !$acc routine seq
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: dzeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      dzeta = 0._wp
   else
      dzeta = - a * c * exp( c * ( 1._wp - qref/qmod ) ) &
      &           * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

pure function trapzd(pol)
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp)            :: tmp1, tmp2
   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)
!  just precalculate all weights and get the job done
   real(wp),parameter :: weights(23) = 0.5_wp * (/ &
&  ( freq (2) - freq (1) ),  &
&  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
&  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
&  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
&  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
&  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
&  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
&  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
&  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
&  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
&  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
&  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
&  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
&  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
&  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
&  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
&  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
&  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
&  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
&  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
&  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
&  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
&  ( freq(23) - freq(22) ) /)

!!  do average between trap(1)-trap(22) .and. trap(2)-trap(23)
!   tmp1 = 0.5_wp * ( &
!&  ( freq (2) - freq (1) ) * ( pol (2) + pol (1) )+ &
!&  ( freq (4) - freq (3) ) * ( pol (4) + pol (3) )+ &
!&  ( freq (6) - freq (5) ) * ( pol (6) + pol (5) )+ &
!&  ( freq (8) - freq (7) ) * ( pol (8) + pol (7) )+ &
!&  ( freq(10) - freq (9) ) * ( pol(10) + pol (9) )+ &
!&  ( freq(12) - freq(11) ) * ( pol(12) + pol(11) )+ &
!&  ( freq(14) - freq(13) ) * ( pol(14) + pol(13) )+ &
!&  ( freq(16) - freq(15) ) * ( pol(16) + pol(15) )+ &
!&  ( freq(18) - freq(17) ) * ( pol(18) + pol(17) )+ &
!&  ( freq(20) - freq(19) ) * ( pol(20) + pol(19) )+ &
!&  ( freq(22) - freq(21) ) * ( pol(22) + pol(21) ))
!   tmp2 = 0.5_wp * ( &
!&  ( freq (3) - freq (2) ) * ( pol (3) + pol (2) )+ &
!&  ( freq (5) - freq (4) ) * ( pol (5) + pol (4) )+ &
!&  ( freq (7) - freq (6) ) * ( pol (7) + pol (6) )+ &
!&  ( freq (9) - freq (8) ) * ( pol (9) + pol (8) )+ &
!&  ( freq(11) - freq(10) ) * ( pol(11) + pol(10) )+ &
!&  ( freq(13) - freq(12) ) * ( pol(13) + pol(12) )+ &
!&  ( freq(15) - freq(14) ) * ( pol(15) + pol(14) )+ &
!&  ( freq(17) - freq(16) ) * ( pol(17) + pol(16) )+ &
!&  ( freq(19) - freq(18) ) * ( pol(19) + pol(18) )+ &
!&  ( freq(21) - freq(20) ) * ( pol(21) + pol(20) )+ &
!&  ( freq(23) - freq(22) ) * ( pol(23) + pol(22) ))

   trapzd = sum(pol*weights)

end function trapzd

pure elemental function cngw(wf,cn,cnref)
   !$acc routine seq
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight

   intrinsic :: exp

   cngw = exp ( -wf * ( cn - cnref )**2 )

end function cngw

pure elemental function dcngw(wf,cn,cnref)
   real(wp),intent(in) :: wf,cn,cnref
   real(wp) :: dcngw

   dcngw = 2*wf*(cnref-cn)*cngw(wf,cn,cnref)

end function dcngw

!* BJ damping function ala DFT-D3(BJ)
!  f(n,rab) = sn*rab**n/(rab**n + R0**n)  w/ R0 = a1*sqrt(C6/C8)+a2
!  see: https://doi.org/10.1002/jcc.21759
pure elemental function fdmpr_bj(n,r,c) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: fdmp
   fdmp = 1.0_wp / ( r**n + c**n )
end function fdmpr_bj
pure elemental function fdmprdr_bj(n,r,c) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: dfdmp
   dfdmp = -n*r**(n-1) * fdmpr_bj(n,r,c)**2
end function fdmprdr_bj

!* original DFT-D3(0) damping
!  f(n,rab) = sn/(1+6*(4/3*R0/rab)**alp)  w/ R0 of unknown origin
pure elemental function fdmpr_zero(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (c/r)**(n+alp)))
end function fdmpr_zero
pure elemental function fdmprdr_zero(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(c/r)**(alp)) &
             - alp*r**n/c*six*(c/r)**(alp-1) ) &
           * fdmpr_zero(n,r,c,alp)**2
!  fdmp = 1.0_wp / (r**n*(1 + 6.0_wp * (c/r)**(n+alp)))
end function fdmprdr_zero

!* fermi damping function from TS and MBD methods
!  f(n,rab) = sn/(1+exp[-alp*(rab/R0-1)]) w/ R0 as experimenal vdW-Radii
pure elemental function fdmpr_fermi(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1.0_wp+exp(-alp*(r/c - 1.0))))
end function fdmpr_fermi
pure elemental function fdmprdr_fermi(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = -(-alp/c*r**n*exp(-alp*(r/c - 1.0)) &
             + n*r**(n-1)*(1.0_wp+exp(-alp*(r/c - 1.0)))) &
             * fdmpr_fermi(n,r,c,alp)**2
end function fdmprdr_fermi

!* optimized power zero damping (M. Head-Gordon)
!  f(n,rab) = sn*rab**(n+alp)/(rab**(n+alp) + R0**(n+alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.7b00176
pure elemental function fdmpr_op(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmpr_op
pure elemental function fdmprdr_op(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = (alp*r*(alp-1) - (n+alp)*r**alp*r**(n+alp-1)) &
           * fdmpr_op(n,r,c,alp)**2
!  fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmprdr_op

!* Sherrill's M-zero damping function
!  f(n,rab) = sn/(1+6*(4/3*R0/rab+a2*R0)**(-alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.6b00780
pure elemental function fdmpr_zerom(n,r,c,rsn,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (r/c+rsn*c)**(-alp)))
end function fdmpr_zerom
pure elemental function fdmprdr_zerom(n,r,c,rsn,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(r/c+rsn*c)**(-alp)) &
              - alp*r**n/c*six*(r/c+rsn*c)**(-alp-1) ) &
           * fdmpr_zerom(n,r,c,rsn,alp)**2
end function fdmprdr_zerom


!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(dispm, nat, atoms, g_a, g_c, wf, q, cn, zeff, gam, &
      &                      zetavec, zerovec, zetadcn, zerodcn, zetadq)
   type(TDispersionData), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Charge scaling height.
   real(wp), intent(in) :: g_a
   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)
   real(wp), intent(in) :: zeff(:)
   real(wp), intent(in) :: gam(:)
   !> weighting and scaling function for the atomic reference systems
   real(wp), intent(out) :: zetaVec(:, :)
   !> weighting and scaling function for the atomic reference systems for q=0
   real(wp), intent(out) :: zeroVec(:, :)
   !> derivative of the weight'n'scale function w.r.t. the partial charges
   real(wp), intent(out) :: zetadq(:, :)
   !> derivative of the weight'n'scale function w.r.t. the coordination number
   real(wp), intent(out) :: zetadcn(:, :)
   !> derivative of the weight'n'scale function w.r.t. the CN for q=0
   real(wp), intent(out) :: zerodcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, twf, gw, expw, expd, gwk, dgwk
   real(wp) :: gi, zi

   zetavec = 0.0_wp
   zerovec = 0.0_wp
   zetadcn = 0.0_wp
   zerodcn = 0.0_wp
   zetadq  = 0.0_wp

   !$omp parallel do shared(zetavec, zetadcn, zetadq, zerodcn) &
   !$omp shared(nat, atoms, dispm, cn, q, g_a, g_c, wf, zerovec) &
   !$omp private(iat, ati, zi, gi, norm, dnorm, iref, icount, twf, gw, expw, &
   !$omp& expd, gwk, dgwk)
   do iat = 1, nat
      ati = atoms(iat)

      zi = zeff(ati)
      gi = g_c * gam(ati)

      norm = 0.0_wp
      dnorm = 0.0_wp
      ! acc loop vector
      do iref = 1, dispm%nref(ati)
         ! acc loop seq private(twf, gw)
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            norm = norm + gw
            dnorm = dnorm + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo
      end do
      norm = 1.0_wp / norm
      ! acc loop vector private(expw, expd)
      do iref = 1, dispm%nref(ati)
         expw = 0.0_wp
         expd = 0.0_wp
         ! acc loop seq private(twf, gw)
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            expw = expw + gw
            expd = expd + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(dispm%cn(:dispm%nref(ati), ati)) &
               & == dispm%cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zetavec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerovec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         zetadcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * dgwk
         zetadq(iref, iat) = dzeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerodcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * dgwk

      end do
   end do

end subroutine weight_references

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(dispm, nat, atoms, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)
   type(TDispersionData), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting and scaling function for the atomic reference systems
   real(wp), intent(in) :: zetaVec(:, :)
   !> derivative of the weight'n'scale function w.r.t. the partial charges
   real(wp), intent(in) :: zetadq(:, :)
   !> derivative of the weight'n'scale function w.r.t. the coordination number
   real(wp), intent(in) :: zetadcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)
   !> derivative of the C6 w.r.t. the partial charge
   real(wp), intent(out) :: dc6dq(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp
   dc6dq = 0.0_wp

   !$omp parallel do default(none) shared(c6, dc6dcn, dc6dq) &
   !$omp shared(nat, atoms, dispm, zetavec, zetadcn, zetadq) &
   !$omp private(iat, ati, jat, atj, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj, &
   !$omp& iref, jref, refc6)
   do iat = 1, nat
      do jat = 1, nat
         if (jat > iat) cycle
         ati = atoms(iat)
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         dc6dqi = 0.0_wp
         dc6dqj = 0.0_wp
         do iref = 1, dispm%nref(ati)
            do jref = 1, dispm%nref(atj)
               refc6 = dispm%c6(iref, jref, ati, atj)
               dc6 = dc6 + zetavec(iref, iat) * zetavec(jref, jat) * refc6
               dc6dcni = dc6dcni + zetadcn(iref, iat) * zetavec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + zetavec(iref, iat) * zetadcn(jref, jat) * refc6
               dc6dqi = dc6dqi + zetadq(iref, iat) * zetavec(jref, jat) * refc6
               dc6dqj = dc6dqj + zetavec(iref, iat) * zetadq(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
         dc6dq(iat, jat) = dc6dqi
         dc6dq(jat, iat) = dc6dqj
      end do
   end do
end subroutine get_atomic_c6

!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_gradient(nat, at, xyz, dispm, cutoff, &
      &  cn, dcndr, q, dqdr, energy, gradient)

   !> Molecular Structure information.
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   !> D4 and damping parameters
   type(TDispersionData), intent(in) :: dispm

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)

   !> Derivative of partial charges w.r.t. atomic coordinates.
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   integer :: max_ref

   !> Charge scaling height.
   real(wp) :: g_a

   !> Charge scaling steepness.
   real(wp) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp) :: wf

   !> DFT-D parameters
   real(wp) :: a1,a2,s6,s8,s10

   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:), dEdq(:)
   integer :: i

   g_a = dispm%g_a
   g_c = dispm%g_c
   wf  = dispm%wf 
   a1 = dispm%a1
   a2 = dispm%a2
   s6 = dispm%s6
   s8 = dispm%s8
   s10= dispm%s10

   max_ref = maxval(dispm%nref(at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), dEdcn(nat), dEdq(nat), source=0.0_wp)

   call weight_references(dispm, nat, at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, at, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)

   call disp_gradient(nat, at, xyz, a1,a2,s6,s8,s10, cutoff, sqrtZr4r2, &
      & c6, dc6dcn, dc6dq, energies, gradient, dEdcn, dEdq)

   !call gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (present(dqdr)) then
      !call gemv(dqdr, dEdq, gradient, beta=1.0_wp)
       call contract(dqdr, dEdq, gradient, beta=1.0_wp)
   end if

   energy = energy + sum(energies)

   deallocate( dEdq, dEdcn, energies, &
   &           dc6dq, dc6dcn, c6, zerodq, zerodcn, &
   &           zerovec, zetadq, zetadcn, zetavec)
end subroutine d4_gradient


subroutine disp_gradient &
      & (nat, at, xyz, a1,a2,s6,s8,s10, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      &  energies, gradient, dEdcn, dEdq)

   !> Molecular Structure information.
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   !> D4 damping parameters
   real(wp),intent(in) :: a1,a2 
   real(wp),intent(in) :: s6,s8,s10

   real(wp), intent(in) :: r4r2(:)

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   real(wp), intent(in) :: c6(:, :)

   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(in) :: dc6dq(:, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energies(:)

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   real(wp), intent(inout) :: dEdcn(:)

   real(wp), intent(inout) :: dEdq(:)

   integer :: iat, jat, ati, atj, itr

   real(wp) :: cutoff2
   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   cutoff2 = cutoff**2
   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, dEdcn, dEdq) &
   !$omp shared(nat, at, xyz, cutoff2, a1,a2,s6,s8,s10, r4r2, c6, dc6dcn, dc6dq) &
   !$omp private(iat, jat, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG)
   do iat = 1, nat
      ati = at(iat)
      do jat = 1, iat
         atj = at(jat)

         r4r2ij = 3.0_wp*r4r2(ati)*r4r2(atj)
         r0 = a1*sqrt(r4r2ij) + a2
         rij = xyz(:, iat) - xyz(:, jat)
         r2 = sum(rij**2)
         if (r2 > cutoff2 .or. r2 < 1.0e-10_wp) cycle

         t6 = 1.0_wp/((r2**3)+(r0**6))
         t8 = 1.0_wp/((r2**4)+(r0**8))
         t10 = 1.0_wp/((r2**5)+(r0**10))

         d6 = (-6.0_wp)*(r2**2)*(t6**2)
         d8 = (-8.0_wp)*(r2**3)*(t8**2)
         d10 = (-10.0_wp)*(r2**4)*(t10**2)

         disp = s6*t6 + s8*r4r2ij*t8 &
            &  + s10*49.0_wp/40.0_wp*r4r2ij**2*t10
         ddisp= s6*d6 + s8*r4r2ij*d8 &
            & + s10*49.0_wp/40.0_wp*r4r2ij**2*d10

         dE = -c6(iat, jat)*disp * 0.5_wp
         dG = -c6(iat, jat)*ddisp*rij

         energies(iat) = energies(iat) + dE
         dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
         dEdq(iat) = dEdq(iat) - dc6dq(iat, jat) * disp
         if (iat /= jat) then
            energies(jat) = energies(jat) + dE
            dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
            dEdq(jat) = dEdq(jat) - dc6dq(jat, iat) * disp
            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG
         end if

      end do
   end do
   !$omp end parallel do

end subroutine disp_gradient

end module gfn0_dftd4
