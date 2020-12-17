      subroutine bsec(f,xmin,xmax,accur,xzero)
      implicit double precision (a-h,o-z)
      icount=0
c     write(2,*)
c     write(2,*) 'Entering BSEC' 
c     write(2,*) xmin,f(xmin),xmax,f(xmax)
c     stop
      if(f(xmin)*f(xmax).ge.0.d0)go to 1
    2 icount=icount+1
      x2=(xmax+xmin)/2.d0
      if(f(xmin)*f(x2).lt.0.d0)xmax=x2
      if(f(xmax)*f(x2).lt.0.d0)xmin=x2
c     if(icount.le.10)write(2,*) icount,xmin,xmax,f(xmin),f(xmax)
c     if(icount.gt.50)stop
      if((xmax-xmin).lt.accur)go to 3
      go to 2
    1 stop '>>> CHANGE XMIN, XMAX FOR BSEC'
    3 xzero=x2
      return
      end
