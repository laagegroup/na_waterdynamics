subroutine o4setup(natoms,no4,nmax,no4index)

!find O4' atoms

integer,intent(in) :: natoms,nmax?
integer,intent(out) :: no4
integer,intent(out),dimension(nmax?) :: no4index

integer :: i

no4=0

do i=1,natoms
  if (atomname(i).eq.'O4\'') then
    no4=no4+1
    no4index(no4)=i
  endif
enddo

end subroutine o4setup

subroutine get_o4avedist(nsites,no4,siteindex

integer,intent(in) :: no4,nsites

! define average O4'-O4' distance for all sites 
! for use when the only sites defined are those in the narrow groove

real(kind=4) :: mindist1,mindist2

do s=1,nsites
  m = sitetypeindex(s)
  i = soluteatom(s)
  mindist1=100
  mindist2=100
  do k=1,no4
    j=no4index(k)
    dist ij
    mindist1=
  enddo
enddo

end subroutine get_o4avedist
