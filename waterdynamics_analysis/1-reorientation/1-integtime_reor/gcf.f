      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=500,EPS=1d-10,FPMIN=1d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1d0-a
      c=1d0/FPMIN
      d=1d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2d0
        d=an*d+b
        if(dabs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(dabs(c).lt.FPMIN)c=FPMIN
        d=1d0/d
        del=d*c
        h=h*del
        if(dabs(del-1d0).lt.EPS)goto 1
11    continue
      write(*,*) 'a too large, ITMAX too small in gcf'
      write(*,*) a,ITMAX
      gammcf=huge(1d0)
      goto 2
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
 2    return
      END
