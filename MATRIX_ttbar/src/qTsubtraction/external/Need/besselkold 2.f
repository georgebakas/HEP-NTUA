c.....Function BK(n,z)
c.....BK(n,z) is the n-derivative of BesselK[nu,z]
c.....with respect to nu in nu=1

c.....Itilde defined as in the paper

      function Itilde(m)
      implicit none
      external BK
      real *8 BK,Itilde,argum,dloqt,a_param,b0p
      integer m
C      include 'scales.h'
C      include 'const.h'
C      common/a_param/a_param,b0p

      real *8 Eulergamma,b0,pi,z3,logx
      real *8 xmio 
      common/xmio/xmio

      Eulergamma=0.577215664902d0
      pi=3.14159265359d0
      z3=1.20205690316d0
      b0=2*dexp(-Eulergamma)

C

      argum=b0*xmio
      logx=dlog(xmio)
      
      if (m.eq.1) then
         Itilde=-b0/xmio*BK(0,argum)
      elseif (m.eq.2) then
         Itilde=2*b0/xmio*(BK(0,argum)*logx-BK(1,argum))
      endif


      return
      end



      function BK(n,z)
      implicit none
      external fb
      real *8 bk,fb,errest,z,zz,max,adpint
      real *8 loz,zhalf,egamma
      integer n,nn,ifail
      common/nuorder/nn
      common/zz/zz
      nn=n
      zz=z
      max=10d0
      egamma=0.577215664902d0
      zhalf=z/2d0
      loz=dlog(zhalf)
C     Use approximated form for n=0,1 and z<1.5d0
      if((n.eq.0).and.(z.lt.1.5d0)) then

      bk=1d0/z+zhalf*(loz-0.5d0*(1-2*egamma))
     & +0.5d0*zhalf**3*(loz-0.5d0*(2.5d0-2*egamma))
     & +zhalf**5/12d0*(loz-0.5d0*(10d0/3-2*egamma))
     & +zhalf**7/144d0*(loz-0.5d0*(47d0/12-2*egamma))
     & +zhalf**9/2880d0*(loz-0.5d0*(131d0/30-2*egamma))

      elseif((n.eq.1).and.(z.lt.1.5d0)) then

      bk=-(loz+egamma)-zhalf**2*(loz-1+egamma)
     & -0.25d0*zhalf**4*(loz-1.5d0+egamma)
     & -zhalf**6/36d0*(loz-11d0/6+egamma)
     & -zhalf**8/576d0*(loz-25d0/12+egamma)
      bk=bk/z
      else
      bk=adpint(fb,0d0,max,1d-10,1d-5,errest,ifail)
      endif
      return
      end
      
      
      function fb(t)
      implicit none
      integer nn,nu
      real *8 fb,t,zz
      common/nuorder/nn
      common/zz/zz
      nu=1
      if(nn.eq.0) then
         fb=dexp(-zz*dcosh(t))*dcosh(nu*t)
      elseif(nn.eq.1) then
         fb=dexp(-zz*dcosh(t))*t*dsinh(nu*t)
      elseif(nn.eq.2) then
         fb=dexp(-zz*dcosh(t))*t*t*dcosh(nu*t)
      elseif(nn.eq.3) then
         fb=dexp(-zz*dcosh(t))*t*t*t*dsinh(nu*t)
      endif      
      return
      end
      
