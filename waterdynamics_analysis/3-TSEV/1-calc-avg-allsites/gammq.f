      FUNCTION gammq(a,x)
      REAL*4 a,gammq,x
CU    USES gcf,gser
      REAL*4 gammcf,gamser,gln

      if(x.lt.0.0.or.a.le.0.0)then
         write(*,*) 'bad arguments in gammq'
         stop
      endif
      if(x.lt.a+1.0)then
        call gser(gamser,a,x,gln)
        gammq=1.0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
