program get_tomo
!C. Larmat March 2015
use tomography_par
use generate_databases_par, only: TOMOGRAPHY_PATH,undef_mat_prop,nundefMat_ext_mesh
implicit none
!
integer:: ier,getcwd
integer ::nrecord_max=8120601
double precision :: x,y,z
real(kind=CUSTOM_REAL) :: rho,vp,vs,qkappa,qmu
double precision :: Z_BEGIN,Z_END
integer :: NZO,iz
!parameter ii = 1 ask point by point output; ii=2 profile
integer, parameter :: ii=2
!
nundefMat_ext_mesh=1
allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
undef_mat_prop(1,1)="-1"
undef_mat_prop(2,1)="tomography"
undef_mat_prop(3,1)="acoustic"
!undef_mat_prop(4,1)="atmo_model_SAM_800km_4km.xyz"
undef_mat_prop(5,1)="1"

!ASK USER FILE TOMO TO READ
write(*,*) 'Enter filename of the tomo file'
read(*,*) undef_mat_prop(4,1)
ier = getcwd(TOMOGRAPHY_PATH)
ier = len_trim(TOMOGRAPHY_PATH)
TOMOGRAPHY_PATH(ier+1:ier+1)="/"
!TOMOGRAPHY_PATH="/panfs/scratch3/vol7/carene/run_SPECFEM3D_CUBIT/SAM_4km/DATA/tomo_files/"
!
NFILES_TOMO=1
allocate(ORIG_X(NFILES_TOMO),ORIG_Y(NFILES_TOMO),ORIG_Z(NFILES_TOMO),stat=ier)
allocate(SPACING_X(NFILES_TOMO),SPACING_Y(NFILES_TOMO),SPACING_Z(NFILES_TOMO),stat=ier)
allocate(NX(NFILES_TOMO),NY(NFILES_TOMO),NZ(NFILES_TOMO),stat=ier)
allocate(VP_MIN(NFILES_TOMO),VS_MIN(NFILES_TOMO),RHO_MIN(NFILES_TOMO), &
       VP_MAX(NFILES_TOMO),VS_MAX(NFILES_TOMO),RHO_MAX(NFILES_TOMO),stat=ier)
allocate(nrecord(NFILES_TOMO),stat=ier)
nrecord(1)=nrecord_max
allocate(vp_tomography(NFILES_TOMO,nrecord_max), &
       vs_tomography(NFILES_TOMO,nrecord_max), &
       rho_tomography(NFILES_TOMO,nrecord_max), &
       z_tomography(NFILES_TOMO,nrecord_max),stat=ier)

!
write(*,*) shape(undef_mat_prop)
!
call read_model_tomography(0)
write(*,*) 'Reading model Done'
write(*,*) 'Stats ',minval(vp_tomography),maxval(vp_tomography)
write(*,*) 'Stats ',minval(vs_tomography),maxval(vs_tomography)
write(*,*) 'Stats ',minval(rho_tomography),maxval(rho_tomography)
write(*,*) 'Stats ',minval(z_tomography),maxval(z_tomography)
!
!ASK USER HORIZONTAL COORDINATES
write(*,*) 'Enter horizontal coordinates'
read(*,*) x
read(*,*) y
!x=(-28959.+771041.0)/2.0
!y=(-2069620.0-1269620.0)/2.0

if (ii == 1) then
write(*,*) 'altitude in m ,enter -1 to stop'
read(*,*) z
do while ( z>-1) 
 call model_tomography(x,y,z,rho,vp,vs,qkappa,qmu,-1)
 write(*,*) 'z, rho, vp ',z,rho,vp
 write(*,*) 'Enter altitude in m, enter -1 to stop'
 read(*,*) z
enddo
endif

if (ii== 2) then
Z_BEGIN = 0.d0
Z_END = 400000.d0
NZO=801
do iz = 0, NZO
 z=dble(iz)*(Z_END-Z_BEGIN)/dble(NZO)
 call model_tomography(x,y,z,rho,vp,vs,qkappa,qmu,-1)
 write(35,*) z/1000.0,vp,rho
enddo
endif
!
end program get_tomo
