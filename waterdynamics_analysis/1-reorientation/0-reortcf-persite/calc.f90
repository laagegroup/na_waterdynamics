MODULE calc

! contains subroutines:
! distribution - construct distribution of times
! vec_distance - vector between atoms i and j

implicit none

contains

SUBROUTINE vec_distance(i,j,x,y,z,natms,cell,dx,dy,dz,dr2)

! vector from i to j
! outputs ij distance-squared !!
! vector not normalised !!

integer,intent(in)            :: i,j,natms
real(kind=4),intent(in)       :: x(natms),y(natms),z(natms),cell(3)
real(kind=4),intent(out)      :: dx,dy,dz,dr2

dx = x(j) - x(i)
dy = y(j) - y(i)
dz = z(j) - z(i)
dx = dx - anint(dx/cell(1))*cell(1)
dy = dy - anint(dy/cell(2))*cell(2)
dz = dz - anint(dz/cell(3))*cell(3)
dr2 = dx**2 + dy**2 + dz**2

END SUBROUTINE vec_distance

SUBROUTINE vec_distance_nodr2(i,j,x,y,z,natms,cell,dx,dy,dz)

! vector from i to j
! vector not normalised !!
! saves time by not calculating dr2 when not needed

integer,intent(in)            :: i,j,natms
real(kind=4),intent(in)       :: x(natms),y(natms),z(natms),cell(3)
real(kind=4),intent(out)      :: dx,dy,dz

dx = x(j) - x(i)
dy = y(j) - y(i)
dz = z(j) - z(i)
dx = dx - anint(dx/cell(1))*cell(1)
dy = dy - anint(dy/cell(2))*cell(2)
dz = dz - anint(dz/cell(3))*cell(3)

END SUBROUTINE vec_distance_nodr2

SUBROUTINE distance(i,j,x,y,z,natms,cell,dr2)

! outputs ij distance-squared !!

integer,intent(in)            :: i,j,natms
real(kind=4),intent(in)       :: x(natms),y(natms),z(natms),cell(3)
real(kind=4),intent(out)      :: dr2

real(kind=4)                  :: dx,dy,dz

dx = x(j) - x(i)
dy = y(j) - y(i)
dz = z(j) - z(i)
dx = dx - anint(dx/cell(1))*cell(1)
dy = dy - anint(dy/cell(2))*cell(2)
dz = dz - anint(dz/cell(3))*cell(3)
dr2 = dx**2 + dy**2 + dz**2

END SUBROUTINE distance

SUBROUTINE vec_distance_points(x1,y1,z1,x2,y2,z2,cell,dx,dy,dz,dr2)

! outputs ij distance-squared !!
! vector not normalised !!

real(kind=4),intent(in)       :: x1,y1,z1,x2,y2,z2,cell(3)
real(kind=4),intent(out)      :: dx,dy,dz,dr2

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
dx = dx - anint(dx/cell(1))*cell(1)
dy = dy - anint(dy/cell(2))*cell(2)
dz = dz - anint(dz/cell(3))*cell(3)
dr2 = dx**2 + dy**2 + dz**2

END SUBROUTINE vec_distance_points

SUBROUTINE vec_distance_points_nopbc(x1,y1,z1,x2,y2,z2,dx,dy,dz,dr2)

! outputs ij distance-squared !!
! vector not normalised !!

real(kind=4),intent(in)       :: x1,y1,z1,x2,y2,z2
real(kind=4),intent(out)      :: dx,dy,dz,dr2

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
dr2 = dx**2 + dy**2 + dz**2

END SUBROUTINE vec_distance_points_nopbc

SUBROUTINE calc_angle(rij,rjk,rik,theta)

! angle /_ ijk between atoms i, j and k
! returns angle in DEGREES

real(kind=4),intent(in)  :: rij,rjk,rik
real(kind=4),intent(out) :: theta

real,parameter :: pi = 3.14159265


theta = ((rij**2 + rjk**2 - rik**2)/(2*rij*rjk))
if (theta.gt.1.0) theta = 1.0 ! ### quick fix for numerical accuracy problem for linear H-bonds, better work in dotproduct
theta = acos(theta)
theta = theta*180/pi

END SUBROUTINE calc_angle

function intchar4(i)
  
  implicit none
  integer i
  character*4 intchar4
  character*1 char1
  character*2 char2
  character*3 char3
  
  if (i.lt.10) then
     write(char1,'(i1)') i
     intchar4='000'//char1
  elseif (i.lt.100) then
     write(char2,'(i2)') i
     intchar4='00'//char2
  elseif (i.lt.1000) then
     write(char3,'(i3)') i
     intchar4='0'//char3
  elseif (i.lt.10000) then
     write(intchar4,'(i4)') i
  else
     write(*,*) 'error in intchar4: i>9999'
     stop
  endif

  return
end function intchar4

function intchar3(i)
  
  implicit none
  integer i
  character*3 intchar3
  character*1 char1
  character*2 char2
  
  if (i.lt.10) then
     write(char1,'(i1)') i
     intchar3='00'//char1
  elseif (i.lt.100) then
     write(char2,'(i2)') i
     intchar3='0'//char2
  elseif (i.lt.1000) then
     write(intchar3,'(i3)') i
  else
     write(*,*) 'error in intchar3: i>999'
     stop
  endif

  return
end function intchar3

function intchar2(i)
  
  implicit none
  integer i
  character*2 intchar2
  character*1 char1
  
  if (i.lt.10) then
     write(char1,'(i1)') i
     intchar2='0'//char1
  elseif (i.lt.100) then
     write(intchar2,'(i2)') i
  else
     write(*,*) 'error in intchar2: i>100'
     stop
  endif

  return
end function intchar2

SUBROUTINE constr_hist(tau,ntimes,nbins,tdistmin,tdistmax,hist,histcoord,weight)

! constructs histogram of times, weighted by population (not yet)

integer,intent(in) :: ntimes,nbins
real(kind=4),intent(in)       :: tdistmin,tdistmax 
real(kind=4),dimension(ntimes),intent(in) :: tau
real(kind=4),dimension(nbins),intent(out) :: hist
real(kind=4),dimension(nbins),intent(out) :: histcoord
real(kind=4),dimension(ntimes),intent(in) :: weight

integer :: i,j,bin
real(kind=4) :: binwidth,halfbinwidth,sumhist

binwidth = (tdistmax-tdistmin)/float(nbins)
halfbinwidth = binwidth*0.5

hist=0.0
histcoord=0.0
sumhist=0.0

do i=1,nbins
  histcoord(i) = (i-1)*binwidth + halfbinwidth + tdistmin
enddo

do i=1,ntimes
  if (tau(i).gt.0.0) then
    bin = ceiling( (tau(i)-tdistmin)/binwidth )
    if ((bin.ge.1).and.(bin.le.nbins)) then
      hist(bin) = hist(bin) + weight(i)
      sumhist = sumhist + weight(i)
    else
      write(*,*) 'Warning: for site ',i,' tau = ',tau(i)
      write(*,*) ' which does not fall within histogram range ',tdistmin,' to ',tdistmax,' ps'
    endif
  endif
enddo

hist = hist/sumhist

END SUBROUTINE constr_hist

SUBROUTINE rotate(x,y,z,a,b,c,u,v,w,theta,xr,yr,zr)

! rotates point (x,y,z) about the line through (a,b,c) parallel to
! <u,v,w> by the angle theta

real x,y,z,a,b,c,u,v,w,theta,xr,yr,zr
real uvwsqr, sintheta, costheta,vwsqr,uwsqr,uvsqr

sintheta = sin(theta)
costheta = cos(theta)
uvwsqr = u**2.0 + v**2.0 + w**2.0
vwsqr = v**2.0 + w**2.0
uwsqr = u**2.0 + w**2.0
uvsqr = u**2.0 + v**2.0

xr = ( a*vwsqr + u*(-b*v-c*w+u*x+v*y+w*z) + ( (x-a)*vwsqr + u*(b*v+c*w-v*y-w*z) )*costheta + &
        sqrt(uvwsqr)*(b*w-c*v-w*y+v*z)*sintheta )/uvwsqr
yr = ( b*uwsqr + v*(-a*u-c*w+u*x+v*y+w*z) + ( (y-b)*uwsqr + v*(a*u+c*w-u*x-w*z) )*costheta + &
        sqrt(uvwsqr)*(-a*w+c*u+w*x-u*z)*sintheta )/uvwsqr
zr = ( c*uvsqr + w*(-a*u-b*v+u*x+v*y+w*z) + ( (z-c)*uvsqr + w*(a*u+b*v-u*x-v*y) )*costheta + &
        sqrt(uvwsqr)*(a*v-b*u-v*x+u*y)*sintheta )/uvwsqr

END SUBROUTINE rotate

SUBROUTINE recentre_box(natms,segname,resindex,x,y,z,cell,icentr)

integer,intent(in) :: natms,icentr
integer,intent(in) :: resindex(natms)
real,intent(in)    :: cell(3)
real,intent(inout) :: x(natms),y(natms),z(natms)
character(len=3),intent(in) :: segname(natms)

integer i
real xifirst,yifirst,zifirst

i=1
do while (i.le.natms) 
  xifirst = x(i)
  yifirst = y(i)
  zifirst = z(i)
  x(i) = x(i) - cell(1)*anint( (xifirst-x(icentr))/cell(1) )
  y(i) = y(i) - cell(2)*anint( (yifirst-y(icentr))/cell(2) )
  z(i) = z(i) - cell(3)*anint( (zifirst-z(icentr))/cell(3) )
  i=i+1
  if (i.gt.natms) exit
  if ( (resindex(i).ne.resindex(i-1)) .or. (segname(i).ne.segname(i-1)) ) then
    cycle
  else
    do while ( (resindex(i).eq.resindex(i-1)) .and. (segname(i).eq.segname(i-1)) )
      x(i) = x(i) - cell(1)*anint( (xifirst-x(icentr))/cell(1) )
      y(i) = y(i) - cell(2)*anint( (yifirst-y(icentr))/cell(2) )
      z(i) = z(i) - cell(3)*anint( (zifirst-z(icentr))/cell(3) )
      i=i+1
      if (i.gt.natms) exit
    enddo
  endif
enddo 

END SUBROUTINE recentre_box

SUBROUTINE recentre_box_cut(natms,x,y,z,cell,icentr)

integer :: natms,icentr
integer :: resindex(natms)
real,intent(in)    :: cell(3)
real,intent(inout) :: x(natms),y(natms),z(natms)
character(len=3) :: segname(natms)

integer i

i=1
do while (i.le.natms) 
  x(i) = x(i) - cell(1)*anint( (x(i)-x(icentr))/cell(1) )
  y(i) = y(i) - cell(2)*anint( (y(i)-y(icentr))/cell(2) )
  z(i) = z(i) - cell(3)*anint( (z(i)-z(icentr))/cell(3) )
  i=i+1
enddo 

END SUBROUTINE recentre_box_cut

SUBROUTINE bulk_test(code,x,y,z,natms,cell,atom,soluteatom,nsites,r0max)

! checks whether atom is more than r0max from surface of solute

implicit none
integer,intent(in)                   :: atom,natms,nsites
real(kind=4),dimension(natms),intent(in)     :: x,y,z
real(kind=4),intent(in)                      :: r0max,cell(3)
logical,intent(out)                  :: code
integer,dimension(nsites),intent(in) :: soluteatom !original allocated size: nmaxsol

real(kind=4) dx,dy,dz,rd2,rdnearest,rmax2
integer i,is

rmax2 = r0max**2.0
code = .false.

rdnearest = 1000.0

do i=1,nsites
  is=soluteatom(i)
  call vec_distance(is,atom,x,y,z,natms,cell,dx,dy,dz,rd2)
  rdnearest=min(rdnearest,rd2)
enddo

if (rdnearest.gt.rmax2) code=.true.

END SUBROUTINE bulk_test


END MODULE calc
