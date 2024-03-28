MODULE INPUT_OUTPUT

! contains subroutines:
! readin (program input file)
! read_psf (CHARMM psf file)
! read_dcd_header (read NAMD DCD file header)
! get_bin_record_dcd (read NAMD DCD file record)
! read_definitions (read solute site definitions file)
! read_dlpoly2_header (get dfr and natms from Ward's dlpoly trj file)
! get_dlpoly2_record (read Ward's DLPOLY trj file record)

USE types
USE declarations,only : nmax

implicit none

contains

SUBROUTINE readin(nunit,dcdfile,psffile,sitedeffile,vector,solutesegname,watertype,htype,otype,nmaxsol,icentr,roughcut,dt,restart)

character(len=40),intent(out)     :: dcdfile,psffile,sitedeffile
character(len=4),intent(out)      :: vector   ! OH or vector orthogonal to water plane
character(len=3),intent(out)      :: solutesegname,watertype ! SPC or TIP
character(len=4),intent(out)      :: htype,otype
integer,intent(out)               :: nmaxsol
integer,intent(in)                :: nunit
integer,intent(out)       :: icentr
real(kind=4),intent(out) :: roughcut
real(kind=4),intent(out)          :: dt
integer,intent(out)          :: restart

read(nunit,*) !XPLOR style psf
read(nunit,*) psffile
read(nunit,*) !NAMD dcd trajectory
read(nunit,*) dcdfile
read(nunit,*) !file with H-bond criteria for water and each site
read(nunit,*) sitedeffile
read(nunit,*) !only OH for the moment
read(nunit,*) vector      
read(nunit,*) 
read(nunit,*) solutesegname
read(nunit,*) !TIP or SPC (used to define r_oh distance)
read(nunit,*) watertype
read(nunit,*) 
read(nunit,*) htype
read(nunit,*)
read(nunit,*) otype
read(nunit,*) !number of solvent vectors in system, used for memory allocation
read(nunit,*) nmaxsol
read(nunit,*) !set to zero to avoid using it
read(nunit,*) icentr
read(nunit,*)
read(nunit,*) roughcut
read(nunit,*) !MD timestep (in fs) 
read(nunit,*) dt
read(nunit,*) !restart ? 1=yes, 2=no
read(nunit,*) restart

END SUBROUTINE readin


SUBROUTINE read_dcd_header(ninput,nfr,natms,dfr)

    ! Reads the header of a dcd trajectory by NAMD

    integer,intent(in)  :: ninput 
    integer,intent(out) :: nfr ! number of frames
    integer,intent(out) :: natms 
    integer,intent(out) :: dfr ! number of timesteps between frames

    ! Local variables

    integer           :: ntitle,nset,istart,icntrl(6),jcntrl(9)
    character(len=80) :: title(10)
    character(len=4)  :: hdr
    real(kind=8)      :: delta
    integer           :: i

    read(ninput) hdr,nfr,istart,dfr,(icntrl(i),i=1,6),delta,&
         &(jcntrl(i),i=1,9)
    
    read(ninput) ntitle,(title(i),i=1,ntitle)
    read(ninput) natms

  ! Do not close input file !!!

END SUBROUTINE read_dcd_header

SUBROUTINE get_bin_record_dcd(ninput,natms,cell,x,y,z)

    integer,intent(in) :: ninput
    integer,intent(in) :: natms 
    real(kind=4),intent(out)       :: x(natms),y(natms),z(natms) 
    real(kind=4)       :: xx(natms),yy(natms),zz(natms) 
    real(kind=4),intent(out) :: cell(3)


    ! Local variables

 
    integer      :: i
    real(kind=8) :: delta,xtable(6)
    logical      :: lerr

    read(ninput,ERR=10,END=10) (xtable(i),i=1,6)
    read(ninput,ERR=10,END=10) (xx(i),i=1,natms)
    read(ninput,ERR=10,END=10) (yy(i),i=1,natms)
    read(ninput,ERR=10,END=10) (zz(i),i=1,natms)

!    cell(1)=dble(xtable(1))
!    cell(2)=dble(xtable(3))
!    cell(3)=dble(xtable(6))

!    x = dble(xx) - nint(dble(xx)/cell(1))*cell(1)
!    y = dble(yy) - nint(dble(yy)/cell(2))*cell(2)
!    z = dble(zz) - nint(dble(zz)/cell(3))*cell(3)
      
    cell(1)=(xtable(1))
    cell(2)=(xtable(3))
    cell(3)=(xtable(6))

    x = xx - nint(xx/cell(1))*cell(1)
    y = yy - nint(yy/cell(2))*cell(2)
    z = zz - nint(zz/cell(3))*cell(3)

    lerr=.false.

    return

10  lerr = .true.

END SUBROUTINE get_bin_record_dcd

SUBROUTINE read_xyz(nunit,natms,x,y,z,atomnames)

integer,intent(in)                        :: nunit
integer,intent(in)                        :: natms
character(len=4),dimension(natms),intent(out)         :: atomnames
real(kind=4),dimension(natms),intent(out)             :: x,y,z

integer :: i

read(nunit,*)
read(nunit,*)

do i=1,natms
  read(nunit,*) atomnames(i),x(i),y(i),z(i)
enddo

END SUBROUTINE read_xyz

SUBROUTINE read_prmtop(n,natms,resindex,resname,atomname,atomtype)

integer,intent(in)                                :: n,natms
character(len=4),dimension(natms),intent(out)     :: resname
character(len=4),dimension(natms),intent(out)     :: atomname
character(len=4),dimension(natms),intent(out)     :: atomtype
integer,dimension(natms),intent(out)              :: resindex

character(len=15) :: cflag
character(len=6) :: frmt
integer :: i,j,idum,nres,mod20
integer,dimension(:),allocatable :: respointer
character(len=4),dimension(:),allocatable :: resnametemp

do while (cflag.ne.'POINTERS       ')
  read(n,'(6x,a15)') cflag
enddo

read(n,*)
read(n,*)

read(n,'(2i8)') idum,nres

allocate(respointer(nres+1),resnametemp(nres))

cflag='XXXXXXXXX      '

! atomnames

do while (cflag.ne.'ATOM_NAME      ')
  read(n,'(6x,a15)') cflag
enddo

read(n,*)

mod20=mod(natms,20)

do i=1,natms-mod20,20
  read(n,'(20a4)') (atomname(j),j=i,i+19)
enddo

write(frmt,'("(",i2,a2,")")') mod20,'a4'
read(n,frmt) (atomname(j),j=natms-mod20+1,natms)

! resnames and resids

do while (cflag.ne.'RESIDUE_LABEL  ')
  read(n,'(6x,a15)') cflag
enddo

read(n,*)

mod20=mod(nres,20)

do i=1,nres-mod20,20
  read(n,'(20a4)') (resnametemp(j),j=i,i+19)
enddo

write(frmt,'("(",i2,a2,")")') mod20,'a4'
read(n,frmt) (resnametemp(j),j=nres-mod20+1,nres)

read(n,*)
read(n,*)

do i=1,nres-mod20,20
  read(n,'(10i8)') (respointer(j),j=i,i+19)
enddo

write(frmt,'("(",i1,a2,")")') mod20,'i8'
read(n,frmt) (respointer(j),j=nres-mod20+1,nres)

respointer(nres+1)=natms+1
do i=1,nres
  do j=respointer(i),respointer(i+1)-1
    resname(j)=resnametemp(i)
    resindex(j)=i
  enddo
enddo

deallocate(respointer,resnametemp)

! atomtypes

do while (cflag.ne.'AMBER_ATOM_TYPE')
  read(n,'(6x,a15)') cflag
enddo

read(n,*)

mod20=mod(natms,20)

do i=1,natms-mod20,20
  read(n,'(20a4)') (atomtype(j),j=i,i+19)
enddo

write(frmt,'("(",i2,a2,")")') mod20,'a4'  
read(n,frmt) (atomtype(j),j=natms-mod20+1,natms)

END SUBROUTINE read_prmtop

SUBROUTINE read_psf(n,natms,segname,resindex,resname,atomname,atomtype)
    
    integer,intent(in)                                :: n,natms
    character(len=4),dimension(natms),intent(out)     :: resname
    character(len=4),dimension(natms),intent(out)     :: atomname
    character(len=4),dimension(natms),intent(out)     :: atomtype
    character(len=3),dimension(natms),intent(out)     :: segname
    integer,dimension(natms),intent(out)              :: resindex            
    
    
    !--- local variables -------------------
    
    integer                           :: i,j,ntitle,check
    character                         :: str

    read(n,*)
    read(n,*)
    read(n,*) ntitle
    
    do  j=1,ntitle
       read(n,*)
    enddo

    read(n,*)
    read(n,*) check

    if (check.ne.natms) write(*,*) 'Warning: mismatch between no. atoms in dcd and in psf'
    
    do j=1,natms
       read(n,1000,err=20,end=20) i,segname(i),resindex(i),resname(i),atomname(i),atomtype(i)
    enddo

20 continue
1000 format(3x,i5,1x,a3,2x,i5,a4,x,a4,1x,a4) 

END SUBROUTINE read_psf

SUBROUTINE read_definitions(nunit,sitedata,nsitetypes,jump_criterion)

integer,intent(out) :: nsitetypes
integer,intent(in)  :: nunit
type(site),dimension(nmax),intent(out) :: sitedata
type(criteria),intent(out) :: jump_criterion

integer :: ierror,i

i=0
read(nunit,*)
read(nunit,*) jump_criterion%factor_ad,jump_criterion%factor_ah,jump_criterion%factor_hda

do
  i = i + 1
  read(nunit,*,iostat=ierror) sitedata(i)%siteatomtype,sitedata(i)%siteresname,&
                              sitedata(i)%iacceptor,sitedata(i)%d_ad_acc,sitedata(i)%d_ah_acc,&
                              sitedata(i)%th_hda_acc,sitedata(i)%idonor,sitedata(i)%d_ad_don,&
                              sitedata(i)%d_ah_don,sitedata(i)%th_hda_don,&
                              sitedata(i)%ihydrophobe,sitedata(i)%d_hydr
  if (ierror.lt.0) exit
  if (sitedata(i)%iacceptor == 1) then
    sitedata(i)%sitetype = 'acc'
  elseif (sitedata(i)%idonor == 1) then
    sitedata(i)%sitetype = 'don'
  elseif (sitedata(i)%ihydrophobe == 1) then
    sitedata(i)%sitetype = 'hpb'
!  else
!    write(*,*) 'Error in site definitions file, line ',i
  endif
  sitedata(i)%d_ad_acc = sitedata(i)%d_ad_acc**2.0
  sitedata(i)%d_ah_acc = sitedata(i)%d_ah_acc**2.0
  sitedata(i)%d_ad_don = sitedata(i)%d_ad_don**2.0
  sitedata(i)%d_ah_don = sitedata(i)%d_ah_don**2.0
  sitedata(i)%d_hydr   = sitedata(i)%d_hydr**2.0
! ### maybe turn angles into dotproducts here too
enddo

nsitetypes=i-1

END SUBROUTINE read_definitions

SUBROUTINE get_sites_by_index(nin,nsites,soluteatom,nmaxsol)

! for defining sites by index
! for use when only looking at a subsection of sites in the solute
! or when solute cannot be identified by segname, resname etc.

integer,intent(in) :: nin,nmaxsol
integer,intent(out) :: nsites
integer,intent(out),dimension(nmaxsol) :: soluteatom ! nmaxsol far too big, but usually allocated nmaxsol

integer s

read(nin,*) nsites
do s=1,nsites
  read(nin,*) soluteatom(s)
enddo

END SUBROUTINE get_sites_by_index


END MODULE INPUT_OUTPUT
