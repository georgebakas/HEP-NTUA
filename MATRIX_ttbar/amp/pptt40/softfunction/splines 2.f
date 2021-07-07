CC    Routines for splines interpolation
CC    Taken from Numerical Recipes

      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      implicit none
      INTEGER m,n,NN
      double precision x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100) ! Maximum expected value of n and m.
C     USES spline
C     Given an m by n tabulated function ya(1:m,1:n), and tabulated independent variables x2a(1:n), this routine constructs one-dimensional natural cubic splines of the rows of ya and returns the second-derivatives in the array y2a(1:m,1:n). (The array x1a is included in the argument list merely for consistency with routine splin2.)
      INTEGER j,k
      double precision y2tmp(NN),ytmp(NN)
      do j=1,m
      do k=1,n
      ytmp(k)=ya(j,k)
      enddo
      call spline_new(x2a,ytmp,n,1.d30,1.d30,y2tmp)
      do k=1,n
      y2a(j,k)=y2tmp(k)
      enddo
      enddo
      return
      END


      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      implicit none
      INTEGER m,n,NN
      double precision x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100) ! Maximum expected value of n and m.
C     USES spline,splint
C      Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by that routine; and given a desired interpolating point x1,x2; this routine returns an interpolated function value y by bicubic spline interpolation.

      INTEGER j,k
      double precision y2tmp(NN),ytmp(NN),yytmp(NN)
      do j=1,m
      do k=1,n
      ytmp(k)=ya(j,k)
      y2tmp(k)=y2a(j,k)
      enddo
      call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
      enddo
      call spline_new(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      end

      SUBROUTINE spline_new(x,y,n,yp1,ypn,y2)
      implicit none
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
       y2(1)=0d0
       u(1)=0d0
      else
      y2(1)=-0.5d0
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99d30) then
       qn=0d0
       un=0d0
      else
       qn=0.5d0
       un=(3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END


       SUBROUTINE splint(xa,ya,y2a,n,x,y)
       implicit none
       INTEGER n
       double precision x,y,xa(n),y2a(n),ya(n)
C     Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xai’s in order), and given the array y2a(1:n), which is the output from spline above, and given a value of x, this routine returns a cubic-spline interpolated value y.
       INTEGER k,khi,klo
       double precision a,b,h
       klo=1
       khi=n
 1     if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).gt.x)then
       khi=k
       else
       klo=k
       endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.) then
          write(*,*)"Bad xa input in splint" ! The xa’s must be distinct.
          stop
       endif
       a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated.
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+
     &   ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0

       return
       END
