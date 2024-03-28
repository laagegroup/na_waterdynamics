MODULE sys_def

! contains subroutines:
! define_system - identifies and lists sites in solute, and solvent atoms
! assign_sites_h - finds which OH vectors are in hydration layer and assigns site to each; for following OH vector
! assign_states - finds Hbond acceptors for OH vectors which are or were in the hydration layer; uses stricter jcf criteria

USE types
USE calc
USE declarations,only : maxpossh,nmax

implicit none

contains

SUBROUTINE define_system(natms,nsitetypes,solutesegname,segname,resname,atomtype,&
                         sitedata,sitetypeindex,soluteatom,nsites,htype,hwatatom,nhyd,&
                         otype,owatatom,nox,donorlist,pdonorlist,nmaxsol)

!!!!!! in reality for the site_definition file with DNA, Aoife based everything on atom name and not atom type !!!
!!!!!! we cheat by switching the two columns for the solute in the PSF file

!dummy
integer,intent(in)    :: natms,nsitetypes,nmaxsol
integer,intent(out)   :: nsites,nhyd,nox
integer,dimension(nmaxsol),intent(out)         :: sitetypeindex,soluteatom,hwatatom,owatatom,donorlist,pdonorlist
character(len=4),dimension(natms),intent(in) :: resname,atomtype
character(len=3),dimension(natms),intent(in) :: segname
type(site),dimension(nmax),intent(in) :: sitedata 
character(len=3),intent(in) :: solutesegname
character(len=4),intent(in) :: htype,otype

!local
integer :: i,j,k,io
logical :: found

sitetypeindex = 0
nsites = 0
nhyd = 0
nox = 0

write(*,*) "natms in routine",natms

do i=1,natms
  !write(*,*) i
  found = .false.
  if (segname(i).eq.solutesegname) then
    do j=2,nsitetypes ! because sitedata(1) is water
      if ((atomtype(i).eq.sitedata(j)%siteatomtype).and.(resname(i).eq.sitedata(j)%siteresname)) then
        found = .true.
        if (sitedata(j)%iacceptor==1.or.sitedata(j)%ihydrophobe==1) then
          nsites = nsites + 1
          sitetypeindex(nsites) = j
          soluteatom(nsites) = i
        endif
        if (sitedata(j)%idonor==1) then
          nsites = nsites + 1
          sitetypeindex(nsites) = j
          soluteatom(nsites) = i
          ! assumes that for each H the donor atom is the last non-H immediately before 
          k = i - 1
          do
            if (atomtype(k)(1:1).ne.'H') exit
            k = k - 1
          enddo
          pdonorlist(nsites) = k 
        endif
        cycle
      endif
    enddo
    if ((.not.found).and.(atomtype(i)(1:1).ne.'H')) then
      write(*,*) 'Warning: Atom ',i,' of solute unaccounted for in define_system'
    endif
  elseif (atomtype(i).eq.otype) then
   ! write(*,*) "O"
    nox = nox + 1
    owatatom(nox) = i
  elseif (atomtype(i).eq.htype) then
    ! assumes that for each H the donor O is the last O immediately before
    !write(*,*) "H"
    nhyd = nhyd + 1
    hwatatom(nhyd) = i
    donorlist(nhyd) = owatatom(nox)
  !else
     !   write(*,*) "other"
  endif
enddo

write(*,*) "end define system"


END SUBROUTINE define_system

SUBROUTINE pdonorlist_fill(natms,nmaxsol,nsites,pdonorlist,x,y,z,cell,sitedata,sitetypeindex,soluteatom)

!finds nearest heavy atom for the H of each donor site

integer,intent(in)                             :: natms,nmaxsol,nsites
integer,dimension(nmaxsol),intent(in)          :: sitetypeindex,soluteatom
type(site),dimension(nmax),intent(in)           :: sitedata
real(kind=4),intent(in)                        :: cell(3)
real(kind=4),dimension(natms),intent(in)       :: x,y,z
integer,dimension(nmaxsol),intent(out)         :: pdonorlist


integer i,s,m,j,mini
real(kind=4) :: dx,dy,dz,r2,mindist

do s=1,nsites
  m=sitetypeindex(s)
  if (sitedata(m)%sitetype.ne.'don') cycle
  j = soluteatom(s)
  mini=1
  mindist=100.0 
  do i=1,natms
    if (i.eq.j) cycle
!    if ( (( x(i)-x(j) )**2).gt.4.0 ) cycle ###add cell
    call vec_distance(i,j,x,y,z,natms,cell,dx,dy,dz,r2)
    if (r2.lt.mindist) then
      mindist=r2
      mini=i
    endif
  enddo
  pdonorlist(s)=mini
enddo

END SUBROUTINE pdonorlist_fill

SUBROUTINE assign_sites_h(natms,x,y,z,nhyd,cell,icentr,cutsq,&
       r_oh,hwatatom,soluteatom,sitetypeindex,sitedata,site_wat,&
       nsites,pdonorlist,donorlist,nmaxsol)

! For each OH water bond, finds out which sites it belongs to (with criteria for don, acc and hpb sites), with priority to acc, don and then hpb. If site_wat=0, then OH not in the hydration shell.

integer,intent(in)      :: natms,nsites,nhyd,icentr,nmaxsol
real(kind=4),intent(in) :: cell(3),cutsq,r_oh
integer,dimension(nmaxsol),intent(in) :: hwatatom
real(kind=4),dimension(natms),intent(in)       :: x,y,z
integer,dimension(nhyd),intent(out) :: site_wat
type(site),dimension(nmax),intent(in) :: sitedata
integer,dimension(nmaxsol),intent(in) :: sitetypeindex
integer,dimension(nmaxsol),intent(in) :: soluteatom
integer,dimension(nmaxsol),intent(in) :: pdonorlist
integer,dimension(nmaxsol),intent(in) :: donorlist

integer :: d,h,i,j,k,l,m,n,nposs
real(kind=4) :: dx,dy,dz,rtemp,rtemp2,theta,rtemp_hd,rtemp_do,rtemp_do2,rtemp_hd2,r2,rij
integer,dimension(maxpossh)  :: poss_site,rsave

poss_site = 0

!$OMP PARALLEL DO &
!$OMP PRIVATE(k,i,nposs,dx,dy,dz,r2,l,j,m,n,rtemp2,rij,rtemp,theta,poss_site,rsave,d,rtemp_do2,rtemp_hd2,rtemp_hd,rtemp_do) &
!$OMP SHARED(nhyd,hwatatom,icentr,nmaxsol,x,y,z,natms,cell,cutsq,site_wat,nsites,soluteatom,sitetypeindex),& 
!$OMP& SHARED(sitedata,donorlist,r_oh,pdonorlist) 

!loop over H atoms
do k=1,nhyd
  i = hwatatom(k)
  nposs = 0 ! for counting possible sites
  if (icentr.gt.0) then
    ! get rid of molecules not within large sphere around surface of solute ### save even more time here by only looking at one H of the water molecule
    ! for system not containing globular solute, set icentr to 0
    call vec_distance(i,icentr,x,y,z,natms,cell,dx,dy,dz,r2)
    if (r2.gt.cutsq) then
      !$OMP CRITICAL      
      site_wat(k) = 0
      !$OMP END CRITICAL      
      cycle
    endif
  endif
  ! then find which sites are H-bonded to OH vector or within hydrophobic cutoff
  do l=1,nsites
    j = soluteatom(l)
    m = sitetypeindex(l)
    select case (sitedata(m)%sitetype)    
      case ('acc')
        call vec_distance(i,j,x,y,z,natms,cell,dx,dy,dz,r2)
        ! * criterion: A...H distance *
        if (r2.gt.(sitedata(m)%d_ah_acc)) cycle
        ! * criterion: A...O distance *
        n = donorlist(k) ! O
        call vec_distance(n,j,x,y,z,natms,cell,dx,dy,dz,rtemp2)
        if (rtemp2.gt.(sitedata(m)%d_ad_acc)) cycle
        ! * criterion: HOA angle
        rij = sqrt(r2) 
        rtemp = sqrt(rtemp2)
        call calc_angle(r_oh,rtemp,rij,theta)  ! i.e. r_oh,r_oa,r_ha
        if (theta.gt.(sitedata(m)%th_hda_acc)) cycle
        ! all criteria satisfied 
        nposs = nposs + 1
        poss_site(nposs) = l
        rsave(nposs) = r2
      case ('don') ! ### save time here by not checking this a second time for the second H of a water
        ! * criterion: H(of site)...O distance *
        n = donorlist(k) ! O
        call vec_distance(n,j,x,y,z,natms,cell,dx,dy,dz,rtemp2)
        if (rtemp2.gt.(sitedata(m)%d_ah_don)) cycle
        ! * criterion: D...O distance *
        d = pdonorlist(l) ! D of D-H site
        call vec_distance(d,n,x,y,z,natms,cell,dx,dy,dz,rtemp_do2)
        if (rtemp_do2.gt.(sitedata(m)%d_ad_don)) cycle
        ! * criterion: HDO angle *
        call vec_distance(j,d,x,y,z,natms,cell,dx,dy,dz,rtemp_hd2)  ! H-D distance
        rtemp_hd = sqrt(rtemp_hd2)
        rtemp_do = sqrt(rtemp_do2)
        rtemp = sqrt(rtemp2)
        call calc_angle(rtemp_hd,rtemp_do,rtemp,theta)  ! i.e. r_hd,r_do,r_oh (h of site!)
        if (theta.gt.(sitedata(m)%th_hda_don)) cycle
        ! all criteria satisfied
        nposs = nposs + 1
        poss_site(nposs) = l
        rsave(nposs) = rtemp_do2
      case ('hpb')
        ! * criterion: C...O distance
        n = donorlist(k) ! O
        call vec_distance(n,j,x,y,z,natms,cell,dx,dy,dz,r2)
        if (r2.gt.(sitedata(m)%d_hydr)) then
          cycle
        else 
          nposs = nposs + 1
          poss_site(nposs) = l
          rsave(nposs) = r2
        endif
      case default
        write(*,*) 'Problem with sitetype ',m
    end select
  enddo
  ! choose nearest/most important site
  !$OMP CRITICAL      
  if (nposs.gt.0) then
    site_wat(k) = poss_site(1)
    if (nposs.gt.1) then
      do l = 2,nposs
        m = sitetypeindex(site_wat(k))
        n = sitetypeindex(poss_site(l))
        if (sitedata(m)%sitetype.lt.sitedata(n)%sitetype) cycle ! works because a < d < h 
        if (sitedata(m)%sitetype.gt.sitedata(n)%sitetype) then 
          site_wat(k) = poss_site(l)
          cycle
        endif
        if (sitedata(m)%sitetype.eq.sitedata(n)%sitetype) then
          if (rsave(l).gt.rsave(l-1)) then
            cycle
          else
            site_wat(k) = poss_site(l)
          endif
        endif
      enddo
    endif
  else
    site_wat(k) = 0 ! i.e. OH vector k is not in hydration shell
  endif
  !$OMP END CRITICAL      
enddo

!$OMP END PARALLEL DO

END SUBROUTINE assign_sites_h

SUBROUTINE assign_states(state,hbacc,site_wat,nhyd,nox,sitetypeindex,&
                   sitedata,donorlist,hwatatom,owatatom,x,y,z,natms,cell,&
                   r_oh,jump_criterion,nt0,soluteatom,nmaxsol)

integer,intent(in) :: nhyd,nox,natms,nmaxsol
integer,dimension(nhyd),intent(inout) :: hbacc
integer,dimension(nhyd),intent(in) :: site_wat,nt0
integer,dimension(nhyd),intent(out) :: state
type(site),dimension(nmax),intent(in) :: sitedata
integer,dimension(nmaxsol),intent(in) :: sitetypeindex
integer,dimension(nmaxsol),intent(in) :: donorlist,soluteatom,hwatatom,owatatom
real(kind=4),dimension(natms),intent(in) :: x,y,z
real(kind=4),dimension(3),intent(in) :: cell
real(kind=4),intent(in) :: r_oh
type(criteria),intent(in) :: jump_criterion


integer :: ia,ih,io,io2,j,k,m
real(kind=4) :: dx,dy,dz,rtemp,rtemp2,theta,rtemp_oo,rtemp_ho,rtemp_oo2,rtemp_ho2
character(len=3) :: ctype
real(kind=4) :: r2,r

! state = 1 if OH is H-bonded to an acceptor site of the solute
! state = 0 if OH is donating a H bond to a water oxygen
! state = -1 if there is no HB acc with strict criteria
! if there is no HB acc, hbacc(k) retains its previous value
! state is used to differentiate jumps paths (towards another site or towards a water oxygen)
! state is used to decide whether a hbond with a tcf time origin also meets the stricter criteria to be a jcf time origin
! hbacc is used to decide whether a jump has occurred

!$OMP PARALLEL DO &
!$OMP PRIVATE(k,m,ctype,ih,ia,dx,dy,dz,r2,r,io,rtemp2,rtemp,theta,j,io2,rtemp_ho2,rtemp_oo2,rtemp_ho,rtemp_oo) &
!$OMP SHARED(nhyd,nt0,nmaxsol,site_wat,sitedata,sitetypeindex,hwatatom,soluteatom,x,y,z,natms,cell,jump_criterion,state), & 
!$OMP& SHARED(donorlist,r_oh,nox,hbacc,owatatom) 

do k=1,nhyd
  if ((nt0(k).eq.0).and.(site_wat(k).eq.0)) cycle ! want hbacc for all H with time origins, and also for H which will have their first time origin this step
  if (site_wat(k).gt.0) then
     m = sitetypeindex(site_wat(k))
     ctype = sitedata(m)%sitetype
  else
     ctype = 'xxx'
  endif
  if (ctype.eq.'acc') then ! must apply stronger criteria to existing hbond
    ih = hwatatom(k)
    ia = soluteatom(site_wat(k)) 
    ! * criterion: A...H distance *
    call vec_distance(ih,ia,x,y,z,natms,cell,dx,dy,dz,r2) ! r2 has already been calculated but large array r2(natms,natms) contributes to memory problems and time gain is minimal
    if (r2.gt.(sitedata(m)%d_ah_acc*jump_criterion%factor_ah)) then ! doesn't meet stricter criteria, hbacc remains as before or undefined
      !$OMP CRITICAL      
      state(k) = -1
      !$OMP END CRITICAL      
      cycle
    endif
    ! * criterion: A...O distance *
    io = donorlist(k) ! O of water molecule
    call vec_distance(ia,io,x,y,z,natms,cell,dx,dy,dz,rtemp2) ! rtemp2 has already been calculated but... see above
    if (rtemp2.gt.(sitedata(m)%d_ad_acc*jump_criterion%factor_ad)) then
      !$OMP CRITICAL      
      state(k) = -1
      !$OMP END CRITICAL      
      cycle
    endif
    ! * criterion: HOA angle *
    r = sqrt(r2) 
    rtemp = sqrt(rtemp2)
    call calc_angle(r_oh,rtemp,r,theta)  ! i.e. r_oh,r_oa,r_ha
    if (theta.gt.(sitedata(m)%th_hda_acc*jump_criterion%factor_hda)) then
      !$OMP CRITICAL      
      state(k) = -1
      !$OMP END CRITICAL      
      cycle
    endif
    ! * all criteria met *
    !$OMP CRITICAL      
    state(k) = 1 ! hbacc is solute site
    hbacc(k) = site_wat(k)+nox ! indices 1 to nox refer to water oxygen acceptor
    !$OMP END CRITICAL      
  else ! hb acceptor not yet known
    ih = hwatatom(k)
    io = donorlist(k)
    !$OMP CRITICAL      
    state(k) = -1 
    !$OMP END CRITICAL      
    do j = 1,nox ! hb acceptor has to be a water bcs otherwise OH would already be assigned to that site
      io2 = owatatom(j)
      if (io2.eq.io) cycle
      ! * criterion: H...O2 distance *
      call vec_distance(ih,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_ho2)
      if (rtemp_ho2.gt.(sitedata(1)%d_ah_acc*jump_criterion%factor_ah)) cycle
      ! * criterion: O...O2 distance *
      call vec_distance(io,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_oo2)
      if (rtemp_oo2.gt.(sitedata(1)%d_ad_acc*jump_criterion%factor_ad)) cycle
      ! * criterion: H O O2 angle (HDA angle) *
      rtemp_ho = sqrt(rtemp_ho2)
      rtemp_oo = sqrt(rtemp_oo2)
      call calc_angle(r_oh,rtemp_oo,rtemp_ho,theta) ! i.e. r_oh,r_o-o2,r_o2-h
      if (theta.gt.(sitedata(1)%th_hda_acc*jump_criterion%factor_hda)) cycle
      !$OMP CRITICAL      
      hbacc(k) = j ! ox index, not atom index
      state(k) = 0 ! hbacc is water oxygen
      !$OMP END CRITICAL      
      exit ! ### only take first acceptor found, fine since criteria are strict anyway
      ! and what if no acceptor is found? in that case hbacc(k) retains its previous value and state = -1
    enddo
  endif
enddo

!$OMP END PARALLEL DO

END SUBROUTINE assign_states

SUBROUTINE assign_hbacc(hbacc,site_wat,nhyd,nox,sitetypeindex,sitedata,donorlist,&
                  hwatatom,owatatom,x,y,z,natms,cell,r_oh,soluteatom,nmaxsol)

integer,intent(in) :: nhyd,nox,natms,nmaxsol
integer,dimension(nhyd),intent(out) :: hbacc
integer,dimension(nhyd),intent(in) :: site_wat
type(site),dimension(nmax),intent(in) :: sitedata
integer,dimension(nmaxsol),intent(in) :: sitetypeindex
integer,dimension(nmaxsol),intent(in) :: donorlist,soluteatom,hwatatom,owatatom
real(kind=4),dimension(natms),intent(in) :: x,y,z
real(kind=4),dimension(3),intent(in) :: cell
real(kind=4),intent(in) :: r_oh

integer :: ia,ih,io,io2,j,k,m
real(kind=4) :: dx,dy,dz,rtemp,rtemp2,theta,rtemp_oo,rtemp_ho,rtemp_oo2,rtemp_ho2,r2,r
character(len=3) :: ctype

hbacc = 0

do k=1,nhyd
  if (site_wat(k).eq.0) cycle
  m = sitetypeindex(site_wat(k))
  ctype = sitedata(m)%sitetype
  if (ctype.eq.'acc') then
    hbacc(k) = soluteatom(site_wat(k)) ! not the same indexing as used in assign_states
  else ! hb acceptor not yet known
    ih = hwatatom(k)
    io = donorlist(k)
    do j = 1,nox ! hb acceptor has to be a water bcs otherwise OH would already be assigned to that site
      io2 = owatatom(j)
      if (io2.eq.io) cycle
      ! * criterion: H...O2 distance *
      call vec_distance(ih,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_ho2)
      if (rtemp_ho2.gt.(sitedata(1)%d_ah_acc)) cycle
      ! * criterion: O...O2 distance *
      call vec_distance(io,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_oo2)
      if (rtemp_oo2.gt.(sitedata(1)%d_ad_acc)) cycle
      ! * criterion: H O O2 angle (HDA angle) *
      rtemp_ho = sqrt(rtemp_ho2)
      rtemp_oo = sqrt(rtemp_oo2)
      call calc_angle(r_oh,rtemp_oo,rtemp_ho,theta) ! i.e. r_oh,r_o-o2,r_o2-h
      if (theta.gt.(sitedata(1)%th_hda_acc)) cycle
      hbacc(k) = io2 ! atom index, not ox index, i.e. not the same indexing as used in assign_states
      exit ! only take first acceptor found
      ! and what if no acceptor is found? in that case hbacc(k)=0
    enddo
  endif
enddo

END SUBROUTINE assign_hbacc

SUBROUTINE find_hbacc(state,hbacc,nhyd,nox,donorlist,hwatatom,owatatom,x,y,z,natms,cell,r_oh,sitedata,jump_criterion,nmaxsol)

! find hb acceptor for an individual OH vector
! uses stricter stable-state criteria

! state = 1 if OH is H-bonded to an acceptor site of the solute
! state = 0 if OH is donating a H bond to a water oxygen
! state = -1 if there is no HB acc
! if there is no HB acc, hbacc(k) retains its previous value
! state is used to differentiate jumps paths (towards another site or towards an oxygen)
! state is used to decide whether a hbond with a tcf time origin also meets the stricter criteria to be a jcf time origin
! hbacc is used to decide whether a jump has occurred

! ### for the moment only used for bulk calculation, so only loop over nox, must add loop over acceptor sites for expanded use

integer,intent(in) :: nhyd,nox,natms,nmaxsol
integer,dimension(nhyd),intent(out) :: hbacc
integer,dimension(nhyd),intent(out) :: state
integer,dimension(nmaxsol),intent(in) :: donorlist,hwatatom,owatatom
type(site),dimension(nmax),intent(in) :: sitedata
real(kind=4),dimension(natms),intent(in) :: x,y,z
real(kind=4),dimension(3),intent(in) :: cell
real(kind=4),intent(in) :: r_oh
type(criteria),intent(in) :: jump_criterion

integer :: ia,ih,io,io2,j,k,m
real(kind=4) :: dx,dy,dz,rtemp,rtemp2,theta,rtemp_oo,rtemp_ho,rtemp_oo2,rtemp_ho2
real(kind=4) :: r2,r

do k=1,nhyd
  ih = hwatatom(k)
  io = donorlist(k)
  state(k) = -1
  do j = 1,nox
      io2 = owatatom(j)
      if (io2.eq.io) cycle
      ! * criterion: H...O2 distance *
      call vec_distance(ih,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_ho2)
      if (rtemp_ho2.gt.(sitedata(1)%d_ah_acc*jump_criterion%factor_ah)) cycle
      ! * criterion: O...O2 distance *
      call vec_distance(io,io2,x,y,z,natms,cell,dx,dy,dz,rtemp_oo2)
      if (rtemp_oo2.gt.(sitedata(1)%d_ad_acc*jump_criterion%factor_ad)) cycle
      ! * criterion: H O O2 angle (HDA angle) *
      rtemp_ho = sqrt(rtemp_ho2)
      rtemp_oo = sqrt(rtemp_oo2)
      call calc_angle(r_oh,rtemp_oo,rtemp_ho,theta) ! i.e. r_oh,r_o-o2,r_o2-h
      if (theta.gt.(sitedata(1)%th_hda_acc*jump_criterion%factor_hda)) cycle
      hbacc(k) = j ! ox index, not atom index
      state(k) = 0 ! hbacc is water oxygen
      exit ! ### only take first acceptor found, f
      ! and what if no acceptor is found? in that case hbacc(k) retains its previous value and state = -1
  enddo
enddo

END SUBROUTINE find_hbacc


SUBROUTINE count_wat_neighbours(nwatneighbours,nox,owatatom,io,sitedata,nmaxsol,natms,x,y,z,cell)

integer,intent(in) :: io,nox,nmaxsol,natms
integer,dimension(nmaxsol),intent(in)    :: owatatom
type(site),dimension(nmax),intent(in)    :: sitedata
real(kind=4),intent(in)                  :: cell(3)
real(kind=4),dimension(natms),intent(in) :: x,y,z
real(kind=4),intent(out) :: nwatneighbours

real(kind=4) :: r_oo2
integer :: o,io2

nwatneighbours = 0.0

do o=1,nox
  io2 = owatatom(o)
  if (io.eq.io2) cycle
  call distance(io,io2,x,y,z,natms,cell,r_oo2)
  if (r_oo2.le.(sitedata(1)%d_ad_acc)) nwatneighbours = nwatneighbours + 1.0
enddo

END SUBROUTINE count_wat_neighbours

SUBROUTINE exclvol_setup(natms,nsolutetypes,solutesegname,segname,resname,atomtype,&
                       exclvoldata,solutetypeindex,soluteatomall,nsolatoms,nmaxsol)

integer,intent(in)    :: natms,nsolutetypes,nmaxsol
integer,intent(out)   :: nsolatoms
integer,dimension(nmaxsol),intent(out)       :: solutetypeindex,soluteatomall
character(len=4),dimension(natms),intent(in) :: resname,atomtype
character(len=3),dimension(natms),intent(in) :: segname
type(exclvol),dimension(nmax),intent(in)      :: exclvoldata
character(len=3),intent(in) :: solutesegname

!local
integer :: i,j,k
logical :: found

solutetypeindex = 0
nsolatoms = 0

do i=1,natms
  found = .false.
  if (segname(i).eq.solutesegname) then
    do j=1,nsolutetypes
      if ((atomtype(i).eq.exclvoldata(j)%soluteatomtype).and.(resname(i).eq.exclvoldata(j)%soluteresname)) then
        found = .true.
        nsolatoms = nsolatoms + 1
        solutetypeindex(nsolatoms) = j
        soluteatomall(nsolatoms) = i
        cycle
      endif
    enddo
    if (.not.found) then
        write(*,*) 'Warning: no data in exclvolfile for solute atom',i
        !print *,"segname == solutesegname"
        print *,"atomtype =",atomtype(i)
        !print *,"resname =",resname(i)
    end if
  endif
enddo

END SUBROUTINE exclvol_setup


END MODULE sys_def


