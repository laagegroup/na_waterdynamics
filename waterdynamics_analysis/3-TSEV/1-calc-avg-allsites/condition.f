      subroutine condition
     x     (x,y,sig,nfit,acfavg,acfstdev,ndata,itfitmin,
     x     ifitmax,dt)

      implicit none

      integer nfit,ndata,itfitmin,ifitmax
      real dt
      real, dimension(nfit):: x,y,sig
      real, dimension(ndata):: acfavg,acfstdev

      integer idat,ifit

      do ifit=1, ifitmax

         idat=itfitmin+ifit
         x(ifit) = dt*dble(idat-1)
         y(ifit) = acfavg(idat)
         sig(ifit) = acfstdev(idat)
         
      enddo

      return
      end
