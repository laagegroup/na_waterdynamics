

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



END MODULE sys_def


