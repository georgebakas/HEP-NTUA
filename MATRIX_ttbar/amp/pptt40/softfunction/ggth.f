C     Threshold expansion of <M0|M2>+cc
C     gg channel

      double precision function ggth(beta,ctheta,nl)
      implicit none
      integer nh
      double precision beta,ctheta,pi,nl
      double precision CF,CA,lb,log2,ggth1,ggth8s,ggth8a
      double precision norm1,norm8s,norm8a
      double precision delta0g1,delta2g1,delta0g8s,delta2g8s
      double precision c0g1(1:10),c0g8s(1:10),c2g1(1:4),c2g8s(1:4)
      data c0g1/0d0,-136.86306d0,25.21974837d0,
     &     -21.52182675d0,0d0,0d0,0.1795129111d0,
     &     6.632216653d0,-1.331755821d0,-6.965100683d0/ 

      data c0g8s/12.64743334d0,-8.317638322d0,-179.8492935d0,
     &     43.04365351d0,-5.762337505d0,1.250977414d0,
     &     19.45833374d0,-11.07432229d0,14.36406511d0,
     &     24.28165052d0/

      data c2g1/0d0,-16.71069286d0,33.42138573d0,-16.71069286d0/

      data c2g8s/0d0,0d0,-8.355346432d0,33.42138573d0/
      
      pi=3.14159265359d0
      CF=4d0/3
      CA=3d0
      nh=1
      
      delta0g1=8*(27*c0g1(1)+3*c0g1(2)+c0g1(3)/3+c0g1(4)/27
     &            +9*(nl*c0g1(5)+nh*c0g1(6))+
     &     nl*c0g1(7)+nh*c0g1(8)+1/9d0*(nl*c0g1(9)+nh*c0g1(10)))

      
      delta2g1=8*(27*c2g1(1)+3*c2g1(2)+c2g1(3)/3+c2g1(4)/27)
      
      delta0g8s=8*(27*c0g8s(1)+3*c0g8s(2)+c0g8s(3)/3+c0g8s(4)/27
     &            +9*(nl*c0g8s(5)+nh*c0g8s(6))+
     &     nl*c0g8s(7)+nh*c0g8s(8)+1/9d0*(nl*c0g8s(9)+nh*c0g8s(10)))
      
      delta2g8s=8*(27*c2g8s(1)+3*c2g8s(2)+c2g8s(3)/3+c2g8s(4)/27)

      norm1=32d0/3
      norm8s=80d0/3
      norm8a=48*beta**2*ctheta**2
      
      log2=0.693147180560d0
      lb=dlog(beta)

C     Eq. (3.10)

      ggth1=norm1*CF*Pi**2*
     &    (-1d0/beta**2*CF*(lb**2+2*log2*lb+log2**2-pi**2/12)
     &   +1/beta*((CA*(-11d0/3-4*log2)+2*nl/3d0)*lb
     &   +CA*(49d0/18-11d0/3*log2-6*log2**2+pi**2/3)+CF*(-5+pi**2/4)
     &   +nl*(-5d0/9+2d0/3*log2))-4*CF*lb**2
     &   -(4*CA+CF*(4+8*log2)+CF*ctheta**2)*lb)
     &   +delta0g1+delta2g1*ctheta**2

C     Eq. (3.11)

      ggth8s=norm8s*(CF-CA/2)*pi**2*
     &      (-1d0/beta**2*(CF-CA/2)*(lb**2+2*log2*lb+log2**2-pi**2/12)
     &     +1d0/beta*((CA*(-11d0/3-2*log2)+2*nl/3d0)*lb
     &     +CA*(67d0/18-8d0/3*log2-4*log2**2+5d0*pi**2/24)
     &     +CF*(-5+pi**2/4)+nl*(-5d0/9+2d0/3*log2))
     &     -4*(CF-CA/2)*lb**2
     &     -(CA*(2d0-4*log2)+CF*(4d0+8*log2)+(CF-CA/2)*ctheta**2)*lb)
     &     +delta0g8s+delta2g8s*ctheta**2

C    Eq. (3.12)

      ggth8a=norm8a*(CF-CA/2)**2*pi**2
     &      *1d0/beta**2*(-lb**2+2*(1-log2)*lb+log2*(2-log2)+pi**2/12)

      ggth=ggth1+ggth8s+ggth8a

      return
      end
