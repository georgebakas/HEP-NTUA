      double precision function qqth(beta,ctheta,nl)
      implicit none 
      integer nh
      double precision beta,ctheta,pi,nl
      double precision CF,CA,lb,log2      
      double precision delta0q,delta1q,delta2q
      double precision c0q(1:10),c1q(1:3),c2q
      data c0q/21.39026702d0,41.88839539d0,-53.87050002d0,
     &       -1.736456791d0,-2.451777776d0,-0.5233480317d0,
     &       -13.46305281d0,-2.175776838d0,0.3322931029d0,1.580246914d0/

      data c1q/0d0,-19.23445615d0,76.93782458d0/

      data c2q/0.8658454253d0/

      pi=3.14159265359d0
      CF=4d0/3
      CA=3d0
      nh=1

      delta0q=16d0*(9*c0q(1)+c0q(2)+c0q(3)/9+3*(nl*c0q(4)+nh*c0q(5))
     &       +(nl*c0q(6)+nh*c0q(7))/3d0
     &       +nl**2*c0q(8)+nl*nh*c0q(9)+nh**2*c0q(10))

      delta1q=16d0*(9*c1q(1)+c1q(2)+c1q(3)/9)

      delta2q=16d0/9*c2q

      log2=0.693147180560d0
      lb=dlog(beta)

      qqth=16*pi**2*(CF-CA/2)*(
     &    -1d0/beta**2*(CF-CA/2)*(lb**2+2*log2*lb+log2**2-pi**2/12)
     &    +1d0/beta*((CA*(-31d0/6+2*log2)+CF*(3-4*log2)+4*nl/3d0)*lb
     &    +CA*(131d0/18-43d0/6*log2+2*log2**2-pi**2/4)
     &    +CF*(-8+6*log2-6*log2**2+7*pi**2/12d0)-8d0*nh/9
     &    +nl*(-10d0/9+2*log2))
     &    -0.5d0*(CF-CA/2)*(3+ctheta**2)*lb**2+(CA*(-41d0/12+3d0/2*log2)
     &    +CF*(-7d0/6-3*log2)+3*nh+(-2*CA+16d0/3*CF)*ctheta
     &    +0.5d0*(CF-CA/2)*(1-2*log2)*ctheta**2)*lb)
     &    +delta0q+delta1q*ctheta+delta2q*ctheta**2

      return
      end
