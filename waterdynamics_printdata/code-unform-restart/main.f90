PROGRAM hydr_shell_hbprint

!write proper description here !!!!!!!!

!$ use OMP_LIB

USE declarations
USE sys_def
USE input_output

implicit none


! ** read input files **

open(nin,file='input',readonly,status='old')
open(nrep,file='calc_details.out')

read(nin,*)
read(nin,*) inputformat

if (inputformat.eq.'NAMD  ') then
  call readin(nin,dcdfile,psffile,sitedeffile,vector,solutesegname,watertype,htype,otype,nmaxsol,icentr,roughcut,dt,restart)
else
  write(*,*) 'Input files should be in NAMD/AMBER format.'
  write(*,*) 'Use the appropriate main program file for the input format.'
endif

roughcut = (roughcut+4.5)**2.0

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


write(nrep,*) 'Trajectory file name: ',dcdfile
write(nrep,*) 'No. of frames in trj file = ',nstep,'no. of timesteps between frames = ',dfr
write(nrep,*) 'No. of atoms = ',natms
allocate (x(natms),y(natms),z(natms))
dt = dt*dble(dfr) 
write(nrep,*) 'time between frames in fs = ',dt

allocate(resindex(natms),resname(natms),segname(natms),atomname(natms),atomtype(natms))

ntmax = int(tmax/dt)
ntdecorr= int(tdecorr/dt)

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

allocate(nt0(nhyd))

if (restart.eq.0) then  ! if first run
        print *,"not restart"
        open(nsit,file="site_t.out",form="unformatted")
        open(nstat,file="state_t.out",form="unformatted")
        open(nhb,file="hbacc_t.out",form="unformatted")
        write(nrep,*) "NOT RESTARTING from previous simulation"
        nskip=0
        nt0=0
else if (restart.eq.1) then  ! restarting previous calculation
        open(nsit,file="site_t.out",form="unformatted",status="old",access="append")
        open(nstat,file="state_t.out",form="unformatted",status="old",access="append")
        open(nhb,file="hbacc_t.out",form="unformatted",status="old",access="append")
        open(nrest,file="restartnt0",form="formatted",status="old")
        read(nrest,*) (nskip,nt0(h),h=1,nhyd)
        write(nrep,*) "RESTART --- restarting from step", nskip, "next step", nskip+1
        close(nrest)
end if 


! ** loop over dcd **


  do istep=1,nstep
    write(*,*) istep
    if (istep.ne.1) call get_bin_record_dcd(ndcd,natms,cell,x,y,z) !bcs already read the first frame earlier

    if (istep.gt.nskip) then !skip first steps already done
        call assign_sites_h(natms,x,y,z,nhyd,cell,icentr,roughcut,r_oh,hwatatom,soluteatom,sitetypeindex,sitedata,&
                        site_wat,nsites,pdonorlist,donorlist,nmaxsol)
        call assign_states(state,hbacc,site_wat,nhyd,nox,sitetypeindex,sitedata,donorlist,hwatatom,owatatom,&
                x,y,z,natms,cell,r_oh,jump_criterion,nt0,soluteatom,nmaxsol)


        do h=1,nhyd
            if (site_wat(h).gt.0) then
                nt0(h)=nt0(h)+1
            end if
        end do

        !!! PRINT RESTART FILES
        !print restart file for nt0, needed for assign_states
        open(nout,file="restartnt0")
        write(nout,'(17512I8)') (istep,nt0(h),h=1,nhyd)
        close(nout)

        ! unformatted writing to save space
        write(nsit) (site_wat(h),h=1,nhyd)
        write(nstat) (state(h),h=1,nhyd)
        write(nhb) (hbacc(h),h=1,nhyd)


    end if !istep>nskip
   end do !steps

close(ndcd)
close(nsit)
close(nstat)
close(nhb)


END PROGRAM hydr_shell_hbprint
