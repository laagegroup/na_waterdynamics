
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

