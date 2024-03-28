MODULE MISC

implicit none

!global 

integer
real(kind=8)

contains

SUBROUTINE update_nbrlist(nbrlist,nbrlist_pointers,natms,nbrlist_cutoff,x,y,z)

! updates neighbour list for each atom
! includes cell index method

integer,intent(inout),allocatable :: nbrlist,nbrlist_pointers
integer,intent(in)                :: natms
real(kind=8),intent(in)           :: cell(3)
real(kind=8),intent(in)             :: nbrlist_cutoff
real(kind=8),intent(in),allocatable :: x,y,z

integer  :: i,j,nbrlistsize
real(kind=8) :: dr,dx,dy,dz,r2

nbrlistsize=10000*natms
r2=nbrlist_cutoff**2.0

allocate(x(natms),y(natms),z(natms),nbrlist(nbrlistsize),nbrlist_pointers(natms))

do i=1,natms
  do j = 1,natms
    call distance(i,j,x,y,z,cell,dr2,dx,dx,dz)
    if (dr2.lt.r2) then

    endif
  enddo
enddo


END SUBROUTINE update_nbr_list

SUBROUTINE

END SUBROUTINE

END MODULE
