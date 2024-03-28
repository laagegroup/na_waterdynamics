PROGRAM reor_sites

!write proper description here !!!!!!!!
! for use with NAMD dcd and psf files
! calculates the reorientation tcf at each site, reading site data printed by save_acc code.

!$ use OMP_LIB

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

! ** read input files for reorientation tcf calculation **
open(nin,file='input_reor',readonly,status='old')
call readin_reor(nin,tmax,tdecorr,nblock)
close(nin)


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

!test
!nstep=100 


write(nrep,*) 'Trajectory file name: ',dcdfile
write(nrep,*) 'No. of frames in trj file = ',nstep,'no. of timesteps between frames = ',dfr
write(nrep,*) 'No. of atoms = ',natms
allocate (x(natms),y(natms),z(natms))
dt = dt*dble(dfr) 
write(nrep,*) 'time between frames in fs = ',dt

allocate(resindex(natms),resname(natms),segname(natms),atomname(natms),atomtype(natms))

write(nrep,*) 'Times length of tcf  = ',tmax,'and decorrelation time = ',tdecorr
!conversion in fs and into number of frames
ntmax = int(dble(tmax)*1000/dt)
ntdecorr= int(dble(tdecorr)*1000/dt)

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

!find heavy atom for each hydrogen on solute donor sites
call get_bin_record_dcd(ndcd,natms,cell,x,y,z)

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


! Memory allocation for tcf
nstepblock = int(dble(nstep)/dble(nblock))
write(nrep,*) 'no. of frames per block = ',nstepblock
popcut=int(dble(nstepblock)/(ntdecorr))
if (popcut.lt.250) then
  write(*,*) 'Warning: poor statistics expected. Use longer trajectory.'
endif
popcut=min(popcut,250)
write(nrep,*) 'population cutoff = ',popcut    

nt0max = int(dble(nstepblock-1)/dble(ntdecorr))+1

allocate(nt0(nhyd),site0(nt0max,nhyd),time0(nt0max,nhyd),dx0(nt0max,nhyd),dy0(nt0max,nhyd),dz0(nt0max,nhyd),&
         acc0(nt0max,nhyd),state0(nt0max,nhyd),abs(nt0max,nhyd)) 
allocate(tcf(ntmax,nsites+1),norm(ntmax,nsites+2),tcfblk(ntmax,nsites+1),normblk(ntmax,nsites+1),&
         poptcf(nsites+1),popblktcf(nsites+1))
!allocate(jumped(nt0max,nhyd))

tcf = 0.0
norm = 0.0
poptcf = 0.0

!open files with site data
open(nstat,file="state_t.out",FORM="UNFORMATTED",status="old")
open(nhb,file="hbacc_t.out",FORM="UNFORMATTED",status="old")
open(nsit,file="site_t.out",FORM="UNFORMATTED",status="old")


! ** loop over dcd **
do iblock=1,nblock
  charblock = intchar2(iblock)
  ! initialisations
  popblktcf = 0 
  tcfblk = 0.0
  normblk = 0.0
  nt0 = 0

  do istep=1,nstepblock
    write(*,*) istep

    !read trajectory to get coordinates for reorientation
    call get_bin_record_dcd(ndcd,natms,cell,x,y,z)

    !read site data
    read(nstat) (state(h),h=1,nhyd)
    read(nsit) (site_wat(h),h=1,nhyd)
    read(nhb) (hbacc(h),h=1,nhyd)
   
    ! loop over vectors
    do h=1,nhyd
      t0 = nt0(h)

      ! new time origin for tcf ?
      new = .false.
      if (t0.eq.0) then ! first time origin for this H
        if (site_wat(h).gt.0) new = .true.
      elseif (istep-time0(t0,h).gt.ntdecorr) then 
        if (site_wat(h).gt.0) new = .true.
      endif
      if ((t0.gt.0).or.new) then
        ih = hwatatom(h)
        io = donorlist(h)
        call vec_distance_nodr2(io,ih,x,y,z,natms,cell,dx,dy,dz) !save times by not calculating dr2
        if (watertype.ne.'SPC') then
          dx = dx*invr_oh
          dy = dy*invr_oh
          dz = dz*invr_oh
        endif
      endif

      ! take new origin
      if (new) then
        t0 = t0+1
        nt0(h) = t0
        time0(t0,h) = istep
        site0(t0,h) = site_wat(h) ! i.e. each time origin is assigned to a site
        popblktcf(site_wat(h)) = popblktcf(site_wat(h)) + 1
        ! for reorientation tcf
        dx0(t0,h) = dx
        dy0(t0,h) = dy
        dz0(t0,h) = dz
      endif

      ! increment tcfs for all time origins, for all sites
      do t=1, t0
        delt = istep-time0(t,h)+1
        if (delt.le.ntmax) then
          s = site0(t,h)
          ! for tcf
          tcfblk(delt,s) = tcfblk(delt,s)+1.5*(dx*dx0(t,h)+dy*dy0(t,h)+dz*dz0(t,h))**2-0.5
          normblk(delt,s) = normblk(delt,s)+1.0
          tcfblk(delt,nsites+1) = tcfblk(delt,nsites+1)+1.5*(dx*dx0(t,h)+dy*dy0(t,h)+dz*dz0(t,h))**2-0.5
          normblk(delt,nsites+1) = normblk(delt,nsites+1)+1.0
        endif
      enddo ! loop over time origins

    enddo ! loop over vectors
  enddo ! loop over steps

  do s=1,nsites+1
    do t=1, ntmax
      ! normalise block tcf   
      if (normblk(t,s).gt.0.0) tcfblk(t,s) = tcfblk(t,s)/normblk(t,s)
      ! update global tcf 
      if (tcfblk(t,s).gt.0) then
        tcf(t,s) = tcf(t,s) + tcfblk(t,s)
      endif
    enddo
  enddo


  ! output of block tcf   
  do s=1,nsites+1

    if ((popblktcf(s).lt.popcut).and.(s.ne.(nsites+1))) cycle

    charsite=intchar4(s)
    open(nout,file='tcf_block'//charblock//'site'//charsite//'.out')
    write(nout,"('# Block between ',f15.3,' and ',f15.3)")dt*dble((iblock-1)*nstepblock)/1000,dt*dble(iblock*nstepblock)/1000
    do t=1, ntmax
      write(nout,'(2f20.10)') dble(t-1)*dt/1000,tcfblk(t,s)
    enddo
    close(nout)
  enddo

  poptcf = poptcf + popblktcf

enddo ! loop over blocks

close(ndcd)

open(nout,file='pop.dat')
write(nout,*) 'Total population: '
write(nout,*) 'Site index, reor tcf population'
do s=1,nsites
  write(nout,*) s,poptcf(s)
enddo
close(nout)

! normalize global tcf and determine std deviation
tcf = tcf/float(nblock)

! output of total tcf
do s=1,nsites+1
    if ((poptcf(s).lt.popcut).and.(s.ne.(nsites+1))) cycle
    charsite=intchar4(s)
   open(nout,file='tcf.out'//charsite)
   do i=1, ntmax
       write(nout,'(2f20.10)') dble(i-1)*dt/1000,tcf(i,s)
   enddo
   close(nout)
enddo


END PROGRAM reor_sites
