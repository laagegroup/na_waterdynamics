      SUBROUTINE fit(x,y,ndim,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER mwt,ndim
      REAL*8 a,b,chi2,q,siga,sigb,sig(ndim),x(ndim),y(ndim)
C     USES gammq
      INTEGER i
      REAL*8 sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0d0
      sy=0d0
      st2=0d0
      b=0d0
      if(mwt.ne.0) then
        ss=0d0
        do 11 i=1,ndim
          wt=1d0/(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndim
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=dble(ndim)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndim
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndim
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1d0+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1d0/st2)
      chi2=0d0
      if(mwt.eq.0) then
        do 15 i=1,ndim
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=1d0
        sigdat=sqrt(chi2/(ndim-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndim
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        q=gammq(0.5d0*(ndim-2),0.5d0*chi2)
        if (q.ge.huge(1d0)) goto 100
      endif
 100  continue
      return
      END
