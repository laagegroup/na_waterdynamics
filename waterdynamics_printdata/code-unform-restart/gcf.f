      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*4 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=500,EPS=1e-10,FPMIN=1e-30)
CU    USES gammln
      INTEGER i
      REAL*4 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.0-a
      c=1.0/FPMIN
      d=1.0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.0
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.0/d
        del=d*c
        h=h*del
        if(abs(del-1.0).lt.EPS)goto 1
11    continue
      write(*,*) 'a too large, ITMAX too small in gcf'
      write(*,*) a,ITMAX
      stop	
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
