!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


subroutine compute_interpolated_dva(displ,veloc,accel,NGLOB_AB, &
                        ispec,NSPEC_AB,ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

  implicit none
  include 'constants.h'

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB) :: displ,veloc,accel
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,iglob

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO
  vxd = ZERO
  vyd = ZERO
  vzd = ZERO
  axd = ZERO
  ayd = ZERO
  azd = ZERO

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! displacement
    dxd = dble(displ(1,iglob))
    dyd = dble(displ(2,iglob))
    dzd = dble(displ(3,iglob))
    ! velocity
    vxd = dble(veloc(1,iglob))
    vyd = dble(veloc(2,iglob))
    vzd = dble(veloc(3,iglob))
    ! acceleration
    axd = dble(accel(1,iglob))
    ayd = dble(accel(2,iglob))
    azd = dble(accel(3,iglob))

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

          ! displacement
          dxd = dxd + dble(displ(1,iglob))*hlagrange
          dyd = dyd + dble(displ(2,iglob))*hlagrange
          dzd = dzd + dble(displ(3,iglob))*hlagrange
          ! velocity
          vxd = vxd + dble(veloc(1,iglob))*hlagrange
          vyd = vyd + dble(veloc(2,iglob))*hlagrange
          vzd = vzd + dble(veloc(3,iglob))*hlagrange
          ! acceleration
          axd = axd + dble(accel(1,iglob))*hlagrange
          ayd = ayd + dble(accel(2,iglob))*hlagrange
          azd = azd + dble(accel(3,iglob))*hlagrange

        enddo
      enddo
    enddo

  endif

end subroutine compute_interpolated_dva

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_interpolated_dva_ac(displ_element,veloc_element,&
                        potential_dot_dot_acoustic,potential_dot_acoustic,&
                        potential_acoustic,NGLOB_AB, &
                        ispec,NSPEC_AB,ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)

! acoustic elements
! returns displacement/velocity/pressure (dxd,..,vxd,..,axd,.. ) at receiver location

  implicit none
  include 'constants.h'

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_acoustic

  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,iglob

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO
  vxd = ZERO
  vyd = ZERO
  vzd = ZERO
  axd = ZERO
  ayd = ZERO
  azd = ZERO

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    ! displacement
    dxd = displ_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    dyd = displ_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    dzd = displ_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))
    ! velocity
    vxd = veloc_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    vyd = veloc_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    vzd = veloc_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))

    ! global index
    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! x component -> acoustic potential
    axd = potential_acoustic(iglob)
    ! y component -> first time derivative of potential
    ayd = potential_dot_acoustic(iglob)
    ! z component -> pressure
    azd = - potential_dot_dot_acoustic(iglob)

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

          ! displacement
          dxd = dxd + hlagrange*displ_element(1,i,j,k)
          dyd = dyd + hlagrange*displ_element(2,i,j,k)
          dzd = dzd + hlagrange*displ_element(3,i,j,k)
          ! velocity
          vxd = vxd + hlagrange*veloc_element(1,i,j,k)
          vyd = vyd + hlagrange*veloc_element(2,i,j,k)
          vzd = vzd + hlagrange*veloc_element(3,i,j,k)

          ! x component -> acoustic potential
          axd = axd + hlagrange*potential_acoustic(iglob)
          ! y component -> first time derivative of potential
          ayd = ayd + hlagrange*potential_dot_acoustic(iglob)
          ! z component -> pressure
          azd = azd - hlagrange*potential_dot_dot_acoustic(iglob)

        enddo
      enddo
    enddo

  endif

end subroutine compute_interpolated_dva_ac
!
!Carene compute the divergence and curl 

subroutine compute_interpolated_dva_div(displ,veloc,accel,NGLOB_AB, &
                        ispec,NSPEC_AB,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
!                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
!                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        jacobian, &
                        ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,div,curlx,curly,curlz)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

  use specfem_par, only : zstore,xstore,ystore
  implicit none
  include 'constants.h'

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd
  double precision,intent(out) :: div,curlx,curly,curlz

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB) :: displ,veloc,accel
! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,l,iglob
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL)::  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz
!

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO
  vxd = ZERO
  vyd = ZERO
  vzd = ZERO
  axd = ZERO
  ayd = ZERO
  azd = ZERO

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! displacement
    dxd = dble(displ(1,iglob))
    dyd = dble(displ(2,iglob))
    dzd = dble(displ(3,iglob))
    ! velocity
    vxd = dble(veloc(1,iglob))
    vyd = dble(veloc(2,iglob))
    vzd = dble(veloc(3,iglob))
    ! acceleration
    axd = dble(accel(1,iglob))
    ayd = dble(accel(2,iglob))
    azd = dble(accel(3,iglob))

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

          ! displacement
          dxd = dxd + dble(displ(1,iglob))*hlagrange
          dyd = dyd + dble(displ(2,iglob))*hlagrange
          dzd = dzd + dble(displ(3,iglob))*hlagrange
          ! velocity
          vxd = vxd + dble(veloc(1,iglob))*hlagrange
          vyd = vyd + dble(veloc(2,iglob))*hlagrange
          vzd = vzd + dble(veloc(3,iglob))*hlagrange
          ! acceleration
          axd = axd + dble(accel(1,iglob))*hlagrange
          ayd = ayd + dble(accel(2,iglob))*hlagrange
          azd = azd + dble(accel(3,iglob))*hlagrange

        enddo
      enddo
    enddo

  endif

!Calculation of derivatives - non fast scheme a priori
!NB we want div (v)!!!! 

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           iglob = ibool(i,j,k,ispec)
!Carene TEST gradients
!            dummyx_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
!            dummyy_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
!            dummyz_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
!
            dummyx_loc(i,j,k) = veloc(1,iglob)
            dummyy_loc(i,j,k) = veloc(2,iglob)
            dummyz_loc(i,j,k) = veloc(3,iglob)
         enddo
      enddo
  enddo
  
  i = 5! nint(xi_r)
  j = 5! nint(eta_r)
  k = 5! nint(gamma_r)
!Debug
!write(*,*) 'I J K',i,j,k

  tempx1 = 0._CUSTOM_REAL
  tempx2 = 0._CUSTOM_REAL
  tempx3 = 0._CUSTOM_REAL

  tempy1 = 0._CUSTOM_REAL
  tempy2 = 0._CUSTOM_REAL
  tempy3 = 0._CUSTOM_REAL

  tempz1 = 0._CUSTOM_REAL
  tempz2 = 0._CUSTOM_REAL
  tempz3 = 0._CUSTOM_REAL

  do l=1,NGLLX
    hp1 = hprime_xx(i,l)
    tempx1 = tempx1 + dummyx_loc(l,j,k)*hp1
    tempy1 = tempy1 + dummyy_loc(l,j,k)*hp1
    tempz1 = tempz1 + dummyz_loc(l,j,k)*hp1

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ
    hp2 = hprime_yy(j,l)
    tempx2 = tempx2 + dummyx_loc(i,l,k)*hp2
    tempy2 = tempy2 + dummyy_loc(i,l,k)*hp2
    tempz2 = tempz2 + dummyz_loc(i,l,k)*hp2

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ
    hp3 = hprime_zz(k,l)
    tempx3 = tempx3 + dummyx_loc(i,j,l)*hp3
    tempy3 = tempy3 + dummyy_loc(i,j,l)*hp3
    tempz3 = tempz3 + dummyz_loc(i,j,l)*hp3
  enddo

  ! get derivatives of ux, uy and uz with respect to x, y and z
  xixl = xix(i,j,k,ispec)
  xiyl = xiy(i,j,k,ispec)
  xizl = xiz(i,j,k,ispec)
  etaxl = etax(i,j,k,ispec)
  etayl = etay(i,j,k,ispec)
  etazl = etaz(i,j,k,ispec)
  gammaxl = gammax(i,j,k,ispec)
  gammayl = gammay(i,j,k,ispec)
  gammazl = gammaz(i,j,k,ispec)
  jacobianl = jacobian(i,j,k,ispec)

!Debug
!        print*, xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  duxdxl = xixl*tempx1 + etaxl*tempx2 + gammaxl*tempx3
  duxdyl = xiyl*tempx1 + etayl*tempx2 + gammayl*tempx3
  duxdzl = xizl*tempx1 + etazl*tempx2 + gammazl*tempx3

  duydxl = xixl*tempy1 + etaxl*tempy2 + gammaxl*tempy3
  duydyl = xiyl*tempy1 + etayl*tempy2 + gammayl*tempy3
  duydzl = xizl*tempy1 + etazl*tempy2 + gammazl*tempy3

  duzdxl = xixl*tempz1 + etaxl*tempz2 + gammaxl*tempz3
  duzdyl = xiyl*tempz1 + etayl*tempz2 + gammayl*tempz3
  duzdzl = xizl*tempz1 + etazl*tempz2 + gammazl*tempz3

! outputs
  div=duxdxl+duydyl+duzdzl
!Need expression curl in 3D
  curlx = 0.d0
  curly = 0.0
  curlz = 0.0

!Debug
  iglob = ibool(5,5,5,ispec)
!  write(1023,*) ispec,' UX ',zstore(iglob)/1000.d0,duxdxl,duxdyl,duxdzl
!  write(1023,*) ispec,' UY ',zstore(iglob)/1000.d0,duydxl,duydyl,duydzl
!  write(1023,*) ispec,' UZ ',zstore(iglob)/1000.d0,duzdxl,duzdyl,duzdzl


end subroutine compute_interpolated_dva_div

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_interpolated_dva_ac_div(displ_element,veloc_element,&
                        potential_dot_dot_acoustic,potential_dot_acoustic,&
                        potential_acoustic,NGLOB_AB, &
                        ispec,NSPEC_AB, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        jacobian, ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,div)

! acoustic elements
! returns displacement/velocity/pressure (dxd,..,vxd,..,axd,.. ) at receiver location
!  use specfem_par, only: zstore
  use specfem_par, only: zstore,xstore,ystore
  implicit none
  include 'constants.h'
  

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,div

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_acoustic
!  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: zstore

! arrays with mesh parameters per slice
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,l, iglob
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL)::  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
!  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3,zzl

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO
  vxd = ZERO
  vyd = ZERO
  vzd = ZERO
  axd = ZERO
  ayd = ZERO
  azd = ZERO

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    ! displacement
    dxd = displ_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    dyd = displ_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    dzd = displ_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))
    ! velocity
    vxd = veloc_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    vyd = veloc_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    vzd = veloc_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))

    ! global index
    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! x component -> acoustic potential
    axd = potential_acoustic(iglob)
    ! y component -> first time derivative of potential
    ayd = potential_dot_acoustic(iglob)
    ! z component -> pressure
    azd = - potential_dot_dot_acoustic(iglob)

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

          ! displacement
          dxd = dxd + hlagrange*displ_element(1,i,j,k)
          dyd = dyd + hlagrange*displ_element(2,i,j,k)
          dzd = dzd + hlagrange*displ_element(3,i,j,k)
          ! velocity
          vxd = vxd + hlagrange*veloc_element(1,i,j,k)
          vyd = vyd + hlagrange*veloc_element(2,i,j,k)
          vzd = vzd + hlagrange*veloc_element(3,i,j,k)

          ! x component -> acoustic potential
          axd = axd + hlagrange*potential_acoustic(iglob)
          ! y component -> first time derivative of potential
          ayd = ayd + hlagrange*potential_dot_acoustic(iglob)
          ! z component -> pressure
          azd = azd - hlagrange*potential_dot_dot_acoustic(iglob)

        enddo
      enddo
    enddo

  endif

!Carene ??? pressure that should contains the update of internal force
!  div = P/rho = azd/rho

!Calculation of derivatives - non fast scheme a priori
!NB we want div (v)!!!! 
!LR: put electron density profile (a guaussian, for a start) + constant B (-55
!degre inclination angle, 0 degre inclination angle)

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
! Carene TEST gradient
!            dummyx_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
!            dummyy_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
!            dummyz_loc(i,j,k) = 8.00*xstore(iglob)+5.00*ystore(iglob)+10.0*zstore(iglob)
            dummyx_loc(i,j,k) = veloc_element(1,i,j,k)
            dummyy_loc(i,j,k) = veloc_element(2,i,j,k)
            dummyz_loc(i,j,k) = veloc_element(3,i,j,k)
            !dummyx_loc(i,j,k) = 0.d0 !-veloc_element(1,i,j,k)*1.d+09*dexp(-((dble(zstore(iglob))*1.d-3-250.)/50.)**2)
            !dummyy_loc(i,j,k) = -veloc_element(2,i,j,k)*cos(55*3.1416/180.)*1.d+09*dexp(-((dble(zstore(iglob))*1.d-3-250.)/50.)**2)
            !dummyz_loc(i,j,k) = veloc_element(3,i,j,k)*sin(55*3.1416/180.)*1.d+09*dexp(-((dble(zstore(iglob))*1.d-3-250.)/50.)**2)
         enddo
      enddo
  enddo

! Fast computation - place in the middle 
  i = 5
  j = 5
  k = 5

  tempx1 = 0._CUSTOM_REAL
  tempx2 = 0._CUSTOM_REAL
  tempx3 = 0._CUSTOM_REAL

  tempy1 = 0._CUSTOM_REAL
  tempy2 = 0._CUSTOM_REAL
  tempy3 = 0._CUSTOM_REAL

  tempz1 = 0._CUSTOM_REAL
  tempz2 = 0._CUSTOM_REAL
  tempz3 = 0._CUSTOM_REAL


  do l=1,NGLLX
    hp1 = hprime_xx(i,l)
    tempx1 = tempx1 + dummyx_loc(l,j,k)*hp1
    tempy1 = tempy1 + dummyy_loc(l,j,k)*hp1
    tempz1 = tempz1 + dummyz_loc(l,j,k)*hp1

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ
    hp2 = hprime_yy(j,l)
    tempx2 = tempx2 + dummyx_loc(i,l,k)*hp2
    tempy2 = tempy2 + dummyy_loc(i,l,k)*hp2
    tempz2 = tempz2 + dummyz_loc(i,l,k)*hp2

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ
    hp3 = hprime_zz(k,l)
    tempx3 = tempx3 + dummyx_loc(i,j,l)*hp3
    tempy3 = tempy3 + dummyy_loc(i,j,l)*hp3
    tempz3 = tempz3 + dummyz_loc(i,j,l)*hp3
  enddo

  ! get derivatives of ux, uy and uz with respect to x, y and z
  xixl = xix(i,j,k,ispec)
  xiyl = xiy(i,j,k,ispec)
  xizl = xiz(i,j,k,ispec)
  etaxl = etax(i,j,k,ispec)
  etayl = etay(i,j,k,ispec)
  etazl = etaz(i,j,k,ispec)
  gammaxl = gammax(i,j,k,ispec)
  gammayl = gammay(i,j,k,ispec)
  gammazl = gammaz(i,j,k,ispec)
  jacobianl = jacobian(i,j,k,ispec)

!Cross-terms = 0 
  duxdxl = xixl*tempx1 + etaxl*tempx2 + gammaxl*tempx3
!  duxdyl = xiyl*tempx1 + etayl*tempx2 + gammayl*tempx3
!  duxdzl = xizl*tempx1 + etazl*tempx2 + gammazl*tempx3

!  duydxl = xixl*tempy1 + etaxl*tempy2 + gammaxl*tempy3
  duydyl = xiyl*tempy1 + etayl*tempy2 + gammayl*tempy3
!  duydzl = xizl*tempy1 + etazl*tempy2 + gammazl*tempy3

!  duzdxl = xixl*tempz1 + etaxl*tempz2 + gammaxl*tempz3
!  duzdyl = xiyl*tempz1 + etayl*tempz2 + gammayl*tempz3
  duzdzl = xizl*tempz1 + etazl*tempz2 + gammazl*tempz3

! outputs
  div=duxdxl+duydyl+duzdzl

!Debug
  iglob = ibool(5,5,5,ispec)
!  write(1022,*) ispec,zstore(iglob)/1000.d0,duxdxl,duydyl,duzdzl
!  write(1022,*) ispec,' DEBUG ',veloc_element(1,5,5,5),xstore(iglob),ystore(iglob),zstore(iglob)

end subroutine compute_interpolated_dva_ac_div
