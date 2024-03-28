PROGRAM hydr_shell

USE declarations
USE sys_def
USE input_output

implicit none

! tcf - reorientation time correlation function (2nd order)
! jcf - jump time correlation function
! tcffr - frame reorientation tcf (between jumps)
! jcf(i1,12,i3)
! i1 = time
! i2 = index of the site
! i3 = 0 for jump towards bulk oxygen or 1 for jump towards an(other) acceptor site, 2 for total jcf

! ** read input files **

open(nin,file='input-hydrshell',readonly,status='old')
open(nrep,file='calc_details.out')

read(nin,*)
read(nin,*) inputformat

  call readin(nin,dcdfile,psffile,sitedeffile,vector,solutesegname,watertype,htype,otype,nmaxsol,icentr,roughcut,tmax,&
              tdecorr,nblock,tfitmin,tfitmax,tfitmin2,tfitmax2,dt,nbins,tdistmin,tdistmax)

roughcut = (roughcut+4.5)**2.0

close(nin)



  open(ndcd,file=dcdfile,form='unformatted',status='old',readonly)
  call read_dcd_header(ndcd,nstep,natms,dfr)
  write(nrep,*) 'Trajectory file name: ',dcdfile
  write(nrep,*) 'No. of frames in trj file = ',nstep,'no. of timesteps between frames = ',dfr
  write(nrep,*) 'No. of atoms = ',natms
  allocate (x(natms),y(natms),z(natms))
  dt = dt*dble(dfr) 
  write(nrep,*) 'time between frames in fs = ',dt

if (inputformat.eq.'NAMD  ') then
  allocate(resindex(natms),resname(natms),segname(natms),atomname(natms),atomtype(natms))
elseif (inputformat.eq.'DLPOLY') then
  allocate(resindex(natms),resname(natms),segname(natms),atomname(natms))
endif

allocate(beta_reor(natms))
beta_reor = 0.0

! psf for namd

  open(npsf,file=psffile,status='old')
  call read_psf(npsf,natms,segname,resindex,resname,atomname,atomtype)
  close(npsf)

      call get_bin_record_dcd(ndcd,natms,cell,x,y,z)
open(npdbc,file='before.pdb')

do i=1,natms
    write(npdbc,200) 'ATOM ',i,atomname(i),resname(i),resindex(i),x(i),y(i),z(i),1.00,beta_reor(i),segname(i)
enddo

close(npdbc)

call recentre_box_cut(natms,segname,resindex,x,y,z,cell,icentr)

open(npdbc,file='after.pdb')

do i=1,natms
    write(npdbc,200) 'ATOM ',i,atomname(i),resname(i),resindex(i),x(i),y(i),z(i),1.00,beta_reor(i),segname(i)
enddo

close(npdbc)

200 FORMAT (A5,I6,1X,A4,1X,A4,1X,I4,4X,3F8.3,2F6.2,7X,A3)
300 FORMAT (I4,1X,I6,1X,A4,1X,A4,1X,A3,6F10.4)
400 FORMAT (I4,1X,I6,1X,A4,1X,A4,1X,A3,8F10.4)
500 FORMAT (I4,1X,I6,1X,A4,1X,A4,1X,A3,4F10.4)

END PROGRAM hydr_shell
