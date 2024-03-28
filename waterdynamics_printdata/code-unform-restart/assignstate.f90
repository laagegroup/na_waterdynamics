
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



