      FUNCTION gammq(a,x)
      REAL*8 a,gammq,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln

      if(x.lt.0d0.or.a.le.0d0)then
         write(*,*) 'bad arguments in gammq',x,a
         gammq = huge(1d0)
         goto 100
      endif
      if(x.lt.a+1d0)then
        call gser(gamser,a,x,gln)
        gammq=1d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
 100  return
      END
