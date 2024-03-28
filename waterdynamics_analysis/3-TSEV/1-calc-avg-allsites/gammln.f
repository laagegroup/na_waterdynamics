      FUNCTION gammln(xx)
      REAL*4 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146,-86.50532032941677,
     *24.01409824083091,-1.231739572450155,.1208650973866179e-2,
     *-.5395239384953e-5,2.5066282746310005/
      x=xx
      y=x
      tmp=x+5.5
      tmp=(x+0.5)*log(tmp)-tmp
      ser=1.000000000190015
      do 11 j=1,6
        y=y+1.0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
