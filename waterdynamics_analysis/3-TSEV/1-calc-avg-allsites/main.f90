PROGRAM tsev

! for use with NAMD dcd and psf file

USE declarations
USE sys_def
USE input_output

implicit none


! ** read input files for system setup **
open(nin,file='input',readonly,status='old')
open(nrep,file='calc_details.out')

read(nin,*)
read(nin,*) inputformat

if (inputformat.eq.'NAMD  ') then
  call readin(nin,dcdfile,psffile,sitedeffile,vector,solutesegname,watertype,htype,otype,nmaxsol,icentr,roughcut,dt,restart)
else
  write(*,*) 'Input files should be in NAMD format.'
  write(*,*) 'Use the appropriate main program file for the input format.'
endif

roughcut = (roughcut+4.5)**2.0

close(nin)

! ** read input files for jump tcf calculation **
open(nin,file='input_tsev',readonly,status='old')
call readin_tsev(nin,exclvolfile,r_ts,theta_ts,nbins,tdecorr)
close(nin)

write(nrep,*) 'For tsev calculation: r_ts=',r_ts,'theta_ts=',theta_ts,'nbins=',nbins


!conversion in fs
tdecorr=tdecorr*1000

!radius calculation for tsev
radius = sin(theta_ts*pi/180.0) * r_ts
write(nrep,*) 'radius=',radius

if (watertype.eq.'TIP') then ! TIP3P,TIP4P,TIP4P/Ice,TIP4P/2005,TIP5P and others (not TIP4P/2005f or TIP3P/Fw)
  r_oh = 0.9572
elseif (watertype.eq.'SPC') then ! SPC,SPC/E and others (not SPC/Fw)
  r_oh = 1.000
else
  write(*,*) 'Warning: water type in input file not recognised!'
endif

invr_oh = 1.0/r_oh

natmsfix=0

open(ndcd,file=dcdfile,form='unformatted',status='old',readonly)
call read_dcd_header(ndcd,nstep,natms,dfr)



write(nrep,*) 'Trajectory file name: ',dcdfile
write(nrep,*) 'No. of frames in trj file = ',nstep,'no. of timesteps between frames = ',dfr
write(*,*) "nstep=",nstep
write(nrep,*) 'No. of atoms = ',natms
write(*,*) "natms=",natms
allocate (x(natms),y(natms),z(natms))
dt = dt*dble(dfr) 
write(nrep,*) 'time between frames in fs = ',dt

allocate(resindex(natms),resname(natms),segname(natms),atomname(natms),atomtype(natms))

ntdecorr= int(tdecorr/dt)
print *,"ntdecorr=",ntdecorr

! psf for namd
open(npsf,file=psffile,status='old')
call read_psf(npsf,natms,segname,resindex,resname,atomname,atomtype)
close(npsf)

! !!!!! nmax hard coded in parameter file. Max=800 by default.
allocate(sitedata(nmax))

! site definitions file

open(ndef,file=sitedeffile,status='old')
call read_definitions(ndef,sitedata,nsitetypes,jump_criterion)
close(ndef)

! ** find solute and water **
allocate(sitetypeindex(nmaxsol),soluteatom(nmaxsol),hwatatom(nmaxsol),owatatom(nmaxsol),donorlist(nmaxsol),pdonorlist(nmaxsol))

call define_system(natms,nsitetypes,solutesegname,segname,resname,atomtype,&
                         sitedata,sitetypeindex,soluteatom,nsites,htype,hwatatom,nhyd,&
                         otype,owatatom,nox,donorlist,pdonorlist,nmaxsol)

! setup excluded volume calculation
allocate(exclvoldata(800))
open(ndef,file=exclvolfile,status='old',readonly)
call read_exclvol(ndef,exclvoldata,nsolutetypes)
close(ndef)

allocate(solutetypeindex(nmaxsol),soluteatomall(nmaxsol))
call exclvol_setup(natms,nsolutetypes,solutesegname,segname,resname,atomtype,&
                       exclvoldata,solutetypeindex,soluteatomall,nsolatoms,nmaxsol)

!initialization
allocate (r_ooinit(3),r_osol(3),bin(nbins,3),centrep(3))
allocate (neighbour(nsolatoms),rho_tsev(nsites),irho_tsev(nsites),rho_tsevsq(nsites),rho_tsevstdev(nsites),poprho_tsev(nsites))
poprho_tsev = 0
rho_tsevstdev = 0.0
rho_tsevsq = 0.0
irho_tsev = 0.0
rho_tsev = 0.0
neighbour = .false.


!find heavy atom for each hydrogen on solute donor sites
call get_bin_record_dcd(ndcd,natms,cell,x,y,z)
print *, "cell=",cell

call pdonorlist_fill(natms,nmaxsol,nsites,pdonorlist,x,y,z,cell,sitedata,sitetypeindex,soluteatom)

write(nrep,*) 'Solute: ',solutesegname
write(nrep,*) 'Solute: ',nsites,' sites'
write(nrep,*) 'Water: ',nhyd,' ',vector,' vectors'
write(nrep,*) 'Watertype: ',watertype

open(nout,file='site_details.out')
write(nout,*) '# site index, atom index, atomtype, resname, sitetype'
do i=1,nsites
  j = sitetypeindex(i)
  write(nout,*) i,soluteatom(i),sitedata(j)%siteatomtype,sitedata(j)%siteresname,&
             sitedata(j)%sitetype
enddo
close(nout)

allocate (site_wat(nhyd),state(nhyd),hbacc(nhyd))




write(*,*) "allocations done"

!open files with site data
open(nsit,file="site_t.out",FORM="UNFORMATTED",status="old",readonly)

write(*,*) "files opened"

! ** loop over dcd **
do istep=1,nstep 
   if (istep.ne. 1) then
        call get_bin_record_dcd(ndcd,natms,cell,x,y,z)
   end if

   read(nsit) (site_wat(k),k=1,nhyd)

   call assign_hbacc(hbacc,site_wat,nhyd,nox,sitetypeindex,sitedata,donorlist,&
                        hwatatom,owatatom,x,y,z,natms,cell,r_oh,soluteatom,nmaxsol)


   if (mod(istep,ntdecorr).ne.0) cycle


   !debug
   write(*,*) "istep=", istep

  ! loop over OH vectors
  do h=1,nhyd
    if ((site_wat(h).eq.0).or.(hbacc(h).eq.0).or.(site_wat(h).le.0)) cycle ! OH vectors not in the hydration layer, or in the hydration layer but not donating a H bond

    neighbour = .false.
    ! O* O_init vector
    io = donorlist(h) ! O*
    ia = hbacc(h) ! O_init
    call vec_distance(io,ia,x,y,z,natms,cell,r_ooinit(1),r_ooinit(2),r_ooinit(3),d_ooinit2)
    d_ooinit = sqrt(d_ooinit2)
    r_ooinit = r_ooinit/d_ooinit
    ! point at centre of circle
    temp = sqrt(r_ts**2 - radius**2)
    centrep(1) = x(io) + r_ooinit(1)*temp
    centrep(2) = y(io) + r_ooinit(2)*temp
    centrep(3) = z(io) + r_ooinit(3)*temp
    ! construct circle of bins
    ! construct vector perp to r_ooinit of length radius
    ! Gramm-Schmidt using vector (1,0,0)
    bin(1,1) = 1.0 - r_ooinit(1)**2
    bin(1,2) = -r_ooinit(1)*r_ooinit(2)
    bin(1,3) = -r_ooinit(1)*r_ooinit(3)
    length = sqrt((bin(1,1))**2 + (bin(1,2))**2 + (bin(1,3))**2)
    bin(1,:) = bin(1,:)/length*radius
    bin(1,1) = bin(1,1) + centrep(1)     
    bin(1,2) = bin(1,2) + centrep(2)     
    bin(1,3) = bin(1,3) + centrep(3)
    theta = 2*pi/float(nbins) 
    do i=2,nbins
      call rotate(bin(i-1,1),bin(i-1,2),bin(i-1,3),centrep(1),centrep(2),centrep(3),&
                r_ooinit(1),r_ooinit(2),r_ooinit(3),theta,&
                bin(i,1),bin(i,2),bin(i,3)) 
    enddo
    ! find solute atoms which could lie on ring
    do i=1,nsolatoms
      m = solutetypeindex(i)             
      is=soluteatomall(i)                   
      call vec_distance(is,io,x,y,z,natms,cell,r_osol(1),r_osol(2),r_osol(3),r_osol2)
      call vec_distance(is,ia,x,y,z,natms,cell,dx,dy,dz,r_oinitsol2)
      if ((r_osol2.le.((r_ts + exclvoldata(m)%r_excl)**2)).and.(r_oinitsol2.le.((radius+exclvoldata(m)%r_excl)**2))) then 
        neighbour(i) = .true.                  
      end if
    enddo

    f = 0
    outer: do i=1,nbins
      inner: do j=1,nsolatoms
       if (neighbour(j)) then 
        m = solutetypeindex(j)
        is=soluteatomall(j)
        call vec_distance_points(x(is),y(is),z(is),bin(i,1),bin(i,2),bin(i,3),cell,dx,dy,dz,r2)
        if (r2.le.exclvoldata(m)%r2_excl) then
          f = f + 1
          cycle outer
        endif
       end if !neighbor
      enddo inner
    enddo outer

    s = site_wat(h)
    irho_tsev(s) = float(f)/float(nbins)
    if (irho_tsev(s).eq.1.0) then ! ### treat this properly
      irho_tsev(s) = 0.0
    else
      irho_tsev(s) = 1.0/(1.0-irho_tsev(s))
      poprho_tsev(s) = poprho_tsev(s) + 1
    endif
    rho_tsev(s) = rho_tsev(s) + irho_tsev(s)
    rho_tsevsq(s) = rho_tsevsq(s) + (irho_tsev(s))**2
  enddo ! loop over OH vectors
enddo ! loop over steps

where (poprho_tsev.gt.0)
  rho_tsev = rho_tsev/float(poprho_tsev)
  rho_tsevsq = rho_tsevsq/float(poprho_tsev)
  rho_tsevstdev = rho_tsevsq - rho_tsev**2
elsewhere
  rho_tsevstdev = 0.0
endwhere
where (rho_tsevstdev.ge.0)
  rho_tsevstdev = sqrt(rho_tsevstdev)
elsewhere
  rho_tsevstdev = -1.0
endwhere


open(nout,file ='tsev_factor.out')
do s=1,nsites
  j = sitetypeindex(s)
  write(nout,100) s,soluteatom(s),sitedata(j)%siteatomtype,sitedata(j)%siteresname,rho_tsev(s),rho_tsevstdev(s),poprho_tsev(s)
enddo

100 FORMAT (2I5,X,A4,X,A4,X,F14.6,X,F14.6,X,I9)


END PROGRAM tsev 
