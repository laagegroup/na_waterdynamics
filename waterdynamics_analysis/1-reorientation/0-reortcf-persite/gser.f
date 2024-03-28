      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*4 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL*4 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)then
           write(*,*) 'x < 0 in gser'
           stop
        endif
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
 11   continue
      write(*,*) 'a too large, ITMAX too small in gser'
      stop
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
