      subroutine xtest()
      implicit none
      double precision H1

      H1 = 100.d0

      return
      end


      subroutine callmyli3(x, result)
      implicit none
      double precision x
      double precision myLI3
      double precision result

      result = myLI3(x)

      return
      end


      
       function myLI3(x)
       implicit none
       double precision myLI3,xlog,x,PI,Z3
       PI=3.14159265358979312D0
       Z3=1.20205690315959429D0
  
       if (x.lt.0.5d0) then
       xlog=dlog(1d0-x)
       myLI3=  -xlog -(3*xlog**2)/8.-(17*xlog**3)/216.-(5*xlog**4)/576
     #   - (7*xlog**5)/54000. + (7*xlog**6)/86400. + 19*xlog**7/5556600
     #   - xlog**8/752640-11*xlog**9/127008000+11*xlog**10/435456000
       elseif (x.lt.1d0) then
       xlog=dlog(x)
       myLI3=Z3+(Pi**2*xlog)/6+(3d0/4-dlog(-xlog)/2)*xlog**2 
     # -xlog**3/12-xlog**4/288+xlog**6/86400-xlog**8/10160640  
       elseif (x.eq.1d0) then
         myLI3=1.20205690315959429D0
       else
        write(6,*)'wrong argument of Li3!!' 
       endif
       return
       end   


      subroutine initialization_momenta(pin, mu_Q, m_t)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      double precision pin(20)
      double precision pout(mxpart,4)
      double precision mu_Q
      double precision m_t

      integer i,j

c      do i = 1, 20
c         print *,i,pin(i)
c      enddo

      do i = 1, 2
         pout(i,4) = -pin((i - 1) * 5 + 1)
         do j = 1, 3
            pout(i,j) = -pin((i - 1) * 5 + j + 1)
         enddo
      enddo
      do i = 3, 4
         pout(i,4) = pin((i - 1) * 5 + 1)
         do j = 1, 3
            pout(i,j) = pin((i - 1) * 5 + j + 1)
         enddo
      enddo

      scale = mu_Q
      mt = m_t

c      do i = 1, 4
c         do j = 1, 4
c            print *,'pout',i,j,pout(i,j)
c         enddo
c      enddo

      call setup
      call col_operators(pout)

      return
      end



      subroutine fourcorrelators(B4, i1, i2, i3, i4, channel)
      implicit none
      include 'born_col_correl.f'
      integer i1, i2, i3, i4
      double precision B4
      integer channel;

      if (channel.eq.1) then
         B4 = Tgg4(i1, i2, i3, i4)
      elseif (channel.eq.2) then
         B4 = Tqq4(i1, i2, i3, i4)
      endif

      return
      end


      subroutine deltaext(pin, mu_Q, m_t, H1, H1_T34, H1_T13, H1_T23, 
     &     channel)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
c      double precision pin(mxpart*4)
      double precision pin(20)
c      double precision pout(5,4)
      double precision pout(mxpart,4)
      double precision H1qqdelta, T23H1qq, T34H1qq, T13H1qq
      double precision H1ggdelta, T23H1gg, T34H1gg, T13H1gg
      double precision H1, H1_T34, H1_T13, H1_T23
      double precision mu_Q
c mt could be set only once !!!
      double precision m_t

      integer i,j
      integer channel;

c      do i = 1, 20
c         print *,i,pin(i)
c      enddo

      do i = 1, 2
         pout(i,4) = -pin((i - 1) * 5 + 1)
         do j = 1, 3
            pout(i,j) = -pin((i - 1) * 5 + j + 1)
         enddo
      enddo
      do i = 3, 4
         pout(i,4) = pin((i - 1) * 5 + 1)
         do j = 1, 3
            pout(i,j) = pin((i - 1) * 5 + j + 1)
         enddo
      enddo

      scale = mu_Q
      mt = m_t

c      do i = 1, 4
c         do j = 1, 4
c            print *,'pout',i,j,pout(i,j)
c         enddo
c      enddo

c      H1 = 100.d0

c      H1qqdelta(p)

      call setup
      call col_operators(pout)
      if (channel.eq.1) then
         H1 = H1ggdelta(pout)
         H1_T34 = T34H1gg(pout)
         H1_T13 = T13H1gg(pout)
c         H1_T23 = T23H1gg(pout)
      elseif (channel.eq.2) then
         H1 = H1qqdelta(pout)
         H1_T34 = T34H1qq(pout)
         H1_T13 = T13H1qq(pout)
         H1_T23 = T23H1qq(pout)
      endif

      return
      end



ch contributions from the azimuthal average

ch  G*G interference for MUNICH
c      subroutine callmsqGGav(pin, msqGGav)
      subroutine callmsqggav(pin, m_t, msqGGav)
c     double precision function msqGGav(pin)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'rescoeff.f'
      integer i,j,k,naem,i1,j1
      double precision pin(20)
      double precision msqGGav
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,u1,factor
      double precision m_t
c      do i = 1, 20
c         print *,i,pin(i)
c      enddo

      call setup

      mt = m_t
      
      do i = 1, 2
         p(i,4) = -pin((i - 1) * 5 + 1)
         do j = 1, 3
            p(i,j) = -pin((i - 1) * 5 + j + 1)
         enddo
      enddo
      do i = 3, 4
         p(i,4) = pin((i - 1) * 5 + 1)
         do j = 1, 3
            p(i,j) = pin((i - 1) * 5 + j + 1)
         enddo
      enddo


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
      v34=dsqrt(1d0-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)

c      print *, 'gsq = ', gsq
c      print *, 'xn = ', xn
c      print *, 'mt = ', mt
c      print *, 'm_t = ', m_t
      
      msqGGav=(-4d0*gsq**2*(-1d0 + xn**2)*
     .    ((-1d0 + xn**2)*s12**2 + 2d0*xn**2*s12*u1 + 2d0*xn**2*u1**2)*
     .    (2d0*mt**4*s12**2 + 
     .     2d0*mt**2*s12*(2d0*pt**2*s12 + u1*(s12 + u1)) + 
     .    (pt**2*s12 + u1*(s12 + u1))*(3d0*pt**2*s12 + u1*(s12 + u1))))/
     .    (xn*s12**2*u1**2*(s12 + u1)**2)

c      print *, 'msqGGav = ',msqGGav
c      print *, 'avegg = ',avegg

      msqGGav=msqGGav*avegg
c      print *, 'msqGGav = ',msqGGav
 
      return
      end




ch  G*G interference for MUNICH
      subroutine callmsqdgav(pin, m_t, msqDGav)
ch  G*D interference
c      double precision function msqDGav(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'rescoeff.f'
            integer i,j,k,naem,j1
      double precision pin(20)
      double precision msqDGav
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,u1,I1,I2,I3,I3c,I31,I32,a
      double precision bbb

      double precision m_t
c      do i = 1, 20
c         print *,i,pin(i)
c      enddo

      call setup

      mt = m_t
      
      do i = 1, 2
         p(i,4) = -pin((i - 1) * 5 + 1)
         do j = 1, 3
            p(i,j) = -pin((i - 1) * 5 + j + 1)
         enddo
      enddo
      do i = 3, 4
         p(i,4) = pin((i - 1) * 5 + 1)
         do j = 1, 3
            p(i,j) = pin((i - 1) * 5 + j + 1)
         enddo
      enddo

      
      s12=2d0*dot(p,1,2)
      s34=2d0*dot(p,3,4)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)
      a=pt/mt
      if(1d0-4d0*mtrans**2/q2.le.0d0)then
        I3c=0d0
      else
        I3c=dsqrt(1d0-4d0*mtrans**2/q2)
      endif
      

      I1=(1d0 - 1d0/a**2*dlog(1d0 + a**2))/4d0

      I2=(-1d0 + (a**2 + 1d0)/a**2*dlog(1d0 + a**2))/4d0

      I31=dlog((1d0 + v34)/(1d0 - v34))/v34/s34

      I32=((1d0 - dsqrt(1d0 - v34**2))/2d0/v34*
     .      dlog((1d0 + v34)/(1d0 - v34))-
     .     dlog(mtrans**2/mt**2)-
     .     I3c*dlog((1d0 + I3c)/(1d0 - I3c)))/pt**2

      I3=s34*(I31 - I32)/4d0

      msqDGav=(-8d0*gsq**2*mt**2*(-1d0 + xn**2)*pt**2*
     .    ((I3 + 2d0*I2*xn**2*(-2d0 + xn**2) - 
     .           2d0*I1*(-1d0 + xn**2)**2)*s12**2 -
     .      2d0*xn**2*(-2d0*I1 + I3 + 2d0*(I1 - I2)*xn**2)*s12*u1 - 
     .      2d0*xn**2*(-2d0*I1 + I3 + 2d0*(I1 - I2)*xn**2)*u1**2))/
     .  (xn**2*u1**2*(s12 + u1)**2)

      msqDGav=-msqDGav*avegg

      return
      end
      



      


c      double precision function L34(q2,beta34,mtrans)
      double precision function L34(p)
      implicit none
      include 'masses.f'
      include 'npart.f'
      include 'constants.f'
      integer i
      double precision q2,beta34,mtrans,beta34t,c,ct
      double precision p(mxpart,4),pt,L34new,yrap
      double precision myli2,dot,yraptwo,y34,y,y3
      double precision r34(1:2)
      external yraptwo,yrap

      beta34t=dsqrt(1d0-(2*mtrans**2/(q2-2*mtrans**2))**2)
      y=yraptwo(3,4,p)
      y3=yrap(3,p)
      c=dsqrt((1d0-beta34)/(1d0+beta34))
      ct=dsqrt((1d0-beta34t)/(1d0+beta34t))
      do i=1,2
       r34(i)=dot(p,i,3)/dot(p,i,4)
      enddo
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)
      beta34=dsqrt(1d0-mt**4/(dot(p,3,4)**2))
      y34=(y3-y)*2d0
    

      L34=0.5d0*dlog((1d0+beta34)/(1d0-beta34))*dlog(mtrans**4/mt**4)
     .    -2d0*myli2(2d0*beta34/(1d0+beta34))
     .    -0.25d0*dlog((1d0+beta34)/(1d0-beta34))**2
      do i=1,2
       L34=L34+myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*r34(i))
     .        +myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))/r34(i))
     .        +0.5d0*dlog(r34(i))**2
      enddo

      L34new=0.5d0*dlog((1d0+beta34)/(1d0-beta34))*dlog(mtrans**4/mt**4)
     .    -2d0*myli2(2d0*beta34/(1d0+beta34))
     .    -0.25d0*dlog((1d0+beta34)/(1d0-beta34))**2
       L34new=L34new+
     .         2d0*myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*dexp(y34))
     .       +2d0*myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*dexp(-y34))
     .        +y34**2

ch      L34=L34new


      return
      end




      double precision function H1qqdelta(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision saveFacDif



      double precision sn,tn,FinVirtqq
      integer i

      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))

      t2=1d0-t1
      
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)

      qqQQv_1=-64d0*gsq**2/4d0
      
      PoleEp2=-qqQQv_0/2d0*(Tqq(1,1)+Tqq(2,2))
      PoleEp1=-qqQQv_0/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))


c      write (*,*) 'gsq',gsq
c      do i = 1, 4
c         do j = 1, 4
c            print *,'p',i,j,p(i,j)
c         enddo
c      enddo
c      write (*,*) 'V',V
c      write (*,*) 's12',s12
c      write (*,*) 'qqQQv_0',qqQQv_0

c      write (*,*) 'PoleEp1',PoleEp1

ch    this is the additional 1/ep pole
     .  -qqQQv_1/2d0*(Tqq(1,1)+Tqq(2,2))

ch    Begin The IR finite terms

      Facdif=
     .       -qqQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tqq(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tqq(3,3)+Tqq(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +qqQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tqq(i1,j1)
      enddo
      enddo

ch    End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif

ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-qqQQv_0/2d0*(
     .  (-2d0*B1q
     .  +Tqq(3,3)+Tqq(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tqq(1,1)+Tqq(2,2))
     .  *0.5d0*(dlog(scale**2/q2)**2))

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -qqQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tqq(i1,j1)
      enddo
      enddo

ch    pi**2/12 term
      Facdif=Facdif+qqQQv_0/2d0
     .        *pisq/12d0*(Tqq(1,1)+Tqq(2,2))
     

      Facdif=Facdif
     .     -qqQQv_1/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))
     .  -qqQQv_1/2d0*(Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      H1qqdelta=FinVirtqq(sn,tn)

      H1qqdelta=gsq**2/4d0*H1qqdelta
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*PoleEp1
      
      H1qqdelta=(H1qqdelta-Facdif)*aveqq
csk      H1qqdelta=(0.d0-saveFacdif)*aveqq
     
      return
      end

      double precision function T34H1qq(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision saveFacDif



      double precision sn,tn,FinVirtT34qq


      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
      qqQQv_1=-64d0*gsq**2/4d0
      
      PoleEp2=-qqQQv_0/2d0*(Tqq(1,1)+Tqq(2,2))
      PoleEp1=-qqQQv_0/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))

ch    this is the additional 1/ep pole
     .  -qqQQv_1/2d0*(Tqq(1,1)+Tqq(2,2))



ch    Begin The IR finite terms

      Facdif=
     .       -qqQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tqq(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tqq(3,3)+Tqq(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +qqQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tqq(i1,j1)
      enddo
      enddo

ch    End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif


ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-qqQQv_0/2d0*(
     .  (-2d0*B1q
     .  +Tqq(3,3)+Tqq(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tqq(1,1)+Tqq(2,2))
     .  *0.5d0*(dlog(scale**2/q2)**2))

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -qqQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tqq(i1,j1)
      enddo
      enddo

ch    pi**2/12 term

      Facdif=Facdif+qqQQv_0/2d0
     .        *pisq/12d0*(Tqq(1,1)+Tqq(2,2))


      Facdif=Facdif
     .     -qqQQv_1/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))
     .  -qqQQv_1/2d0*(Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T34H1qq=FinVirtT34qq(sn,tn)


      T34H1qq=gsq**2/4d0*T34H1qq
     .        +Tqq(3,4)*(
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1)

      T34H1qq=(T34H1qq-Tqq(3,4)*Facdif)*aveqq
csk      T34H1qq=(0.d0-Tqq(3,4)*saveFacdif)*aveqq

      return
      end

      double precision function T13H1qq(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision saveFacDif



      double precision sn,tn,FinVirtT13qq,FinVirtqq



      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
      qqQQv_1=-64d0*gsq**2/4d0
      
      PoleEp2=-qqQQv_0/2d0*Tqq(1,3)*(Tqq(1,1)+Tqq(2,2))
      PoleEp1=-qqQQv_0/2d0*(
     .   +Tqq(1,3)*(Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(1,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(1,4,1,3)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(2,3,1,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(2,4,1,3))

ch    this is the additional 1/ep pole
     .  -qqQQv_1/2d0*Tqq(1,3)*(Tqq(1,1)+Tqq(2,2))


ch    Begin The IR finite terms

      Facdif=
     .       -qqQQv_0/2d0*Tqq(1,3)*(
     .      1d0/v34*L34(p)*Tqq(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tqq(3,3)+Tqq(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +qqQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tqq4(1,3,i1,j1)
      enddo
      enddo

ch    End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif


ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-qqQQv_0/2d0*Tqq(1,3)*(
     .  (-2d0*B1q
     .  +Tqq(3,3)+Tqq(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tqq(1,1)+Tqq(2,2))
     .  *0.5d0*(dlog(scale**2/q2)**2))

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -qqQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tqq4(1,3,i1,j1)
      enddo
      enddo

ch    pi**2/12 term

      Facdif=Facdif+qqQQv_0/2d0*Tqq(1,3)
     .        *pisq/12d0*(Tqq(1,1)+Tqq(2,2))


      Facdif=Facdif
     .     -qqQQv_1/2d0*(
     .   +Tqq(1,3)*(Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(1,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(1,3,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(1,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(1,3,2,4))
     .  -qqQQv_1/2d0*Tqq(1,3)*
     .     (Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T13H1qq=FinVirtT13qq(sn,tn)


      T13H1qq=gsq**2/4d0*T13H1qq
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1

      T13H1qq=(T13H1qq-Facdif)*aveqq
csk      T13H1qq=(0.d0-saveFacdif)*aveqq

      return
      end

      double precision function T23H1qq(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision saveFacDif



      double precision sn,tn,FinVirtT23qq,FinVirtqq,
     &                 FinVirtT34qq,FinvirtT13qq



      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
      qqQQv_1=-64d0*gsq**2/4d0
      
      PoleEp2=-qqQQv_0/2d0*Tqq(2,3)*(Tqq(1,1)+Tqq(2,2))
      PoleEp1=-qqQQv_0/2d0*(
     .   +Tqq(2,3)*(Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(1,3,2,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(1,4,2,3)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(2,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(2,4,2,3))

ch    this is the additional 1/ep pole
     .  -qqQQv_1/2d0*Tqq(2,3)*(Tqq(1,1)+Tqq(2,2))




ch    Begin The IR finite terms

      Facdif=
     .       -qqQQv_0/2d0*Tqq(2,3)*(
     .      1d0/v34*L34(p)*Tqq(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tqq(3,3)+Tqq(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +qqQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tqq4(2,3,i1,j1)
      enddo
      enddo


ch    End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif


ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-qqQQv_0/2d0*Tqq(2,3)*(
     .  (-2d0*B1q
     .  +Tqq(3,3)+Tqq(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tqq(1,1)+Tqq(2,2))
     .  *0.5d0*(dlog(scale**2/q2)**2))

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -qqQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tqq4(2,3,i1,j1)
      enddo
      enddo


ch    pi**2/12 term

      Facdif=Facdif+qqQQv_0/2d0*Tqq(2,3)
     .        *pisq/12d0*(Tqq(1,1)+Tqq(2,2))


      Facdif=Facdif
     .     -qqQQv_1/2d0*(
     .   +Tqq(2,3)*(Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(2,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(2,3,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(2,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(2,3,2,4))
     .  -qqQQv_1/2d0*Tqq(2,3)*
     .     (Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T23H1qq=FinVirtT23qq(sn,tn)

      T23H1qq=gsq**2/4d0*T23H1qq
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1

      T23H1qq=(T23H1qq-Facdif)*aveqq
csk      T23H1qq=(0.d0-saveFacdif)*aveqq

      return
      end

      double precision function H1ggdelta(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'rescoeff.f'
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision ggQQv_0,ggQQv_1,ggQQv_2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision M0T34eps,M0T13eps,M0T23eps,M0sqeps
      double precision sn,tn,FinVirtgg
      double precision saveFacDif


      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps), _2 for O(eps^2)
      ggQQv_0=gsq**2/16d0/pi**2*M0sq

      ggQQv_1=-16d0/xn*V*(t1**2+t1*t2+t2**2)
     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
     . *gsq**2/4d0

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0


ch    colour-correlated O(eps) piece 
      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)*gsq**2
      M0T34eps=Tggeps(3,4)*gsq**2
      M0T13eps=Tggeps(1,3)*gsq**2
      M0T23eps=Tggeps(2,3)*gsq**2

      PoleEp2=-ggQQv_0/2d0*(Tgg(1,1)+Tgg(2,2))
      PoleEp1=
     .   -ggQQv_0/2d0*(
     .   +Tgg(3,3)+Tgg(4,4)
     .  -2d0*B1g
     .  -(Tgg(1,1)+Tgg(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tgg(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tgg(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tgg(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tgg(2,4))
      
ch    this is the additional 1/ep pole
     .  -ggQQv_1/2d0*(xn+xn)




ch    Begin the IR finite terms

      Facdif=
     .       -ggQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tgg(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tgg(3,3)+Tgg(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +ggQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tgg(i1,j1)
      enddo
      enddo


ch   End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif


ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-ggQQv_0/2d0*(
     .  (-2d0*B1g
     .  +Tgg(3,3)+Tgg(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tgg(1,1)+Tgg(2,2))
     .  *0.5d0*(dlog(scale**2/q2))**2)

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -ggQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tgg(i1,j1)
      enddo
      enddo


ch    pi**2/12 term
      Facdif=Facdif+ggQQv_0/2d0
     .        *pisq/12d0*(Tgg(1,1)+Tgg(2,2))


      Facdif=Facdif
     .             -
ch   New Finite term coming from ep^2 term of the Born
     .   ggQQv_2/2d0*(
     .   +ca+ca)-
ch
     .   M0sqeps/2d0*(
     .   +cf+cf
     .  -2d0*B1g)
     .  -M0T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T13eps/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T23eps/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0sqeps/2d0*(xn+xn)*dlog(scale**2/q2)

      H1ggdelta=FinVirtgg(sn,tn)


      H1ggdelta=gsq**2/4d0*H1ggdelta
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*PoleEp1


      H1ggdelta=(H1ggdelta-Facdif)*avegg
csk      H1ggdelta=(0.d0-saveFacdif)*avegg

      return
      end

      double precision function T34H1gg(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
c      include 'breit.f'
      include 'rescoeff.f'
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision ggQQv_0,ggQQv_1,ggQQv_2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision M0T34eps,M0T13eps,M0T23eps,M0sqeps
      double precision sn,tn,FinVirtgg,FinvirtT34gg
      double precision M0T34T34eps,M0T34T13eps,M0T34T23eps
      double precision saveFacDif


      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps), _2 for O(eps^2)
      ggQQv_0=gsq**2/16d0/pi**2*M0sq

      ggQQv_1=-16d0/xn*V*(t1**2+t1*t2+t2**2)
     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
     . *gsq**2/4d0

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0


      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)*gsq**2
      M0T34eps=Tggeps(3,4)*gsq**2
      M0T34T34eps=Tgg4eps(3,4,3,4)*gsq**2
      M0T34T13eps=Tgg4eps(3,4,1,3)*gsq**2
      M0T34T23eps=Tgg4eps(3,4,2,3)*gsq**2

      
      PoleEp2=-ggQQv_0/2d0*Tgg(3,4)*(Tgg(1,1)+Tgg(2,2))
      PoleEp1=
     .   -ggQQv_0/2d0*(
     .   +Tgg4(3,4,3,3)+Tgg4(3,4,4,4)
     .  -2d0*B1g*Tgg(3,4)
     .  -(Tgg(1,1)+Tgg(2,2))*Tgg(3,4)*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(3,4,3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tgg4(3,4,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tgg4(3,4,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tgg4(3,4,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tgg4(3,4,2,4))      
ch    this is the additional 1/ep pole
     .  -ggQQv_1/2d0*(xn+xn)*Tgg(3,4)

ch    Begin the IR finite terms

      Facdif=
     .       -ggQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tgg4(3,4,3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tgg4(3,4,3,3)+Tgg4(3,4,4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +ggQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tgg4(3,4,i1,j1)
      enddo
      enddo

ch   End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif

ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-ggQQv_0/2d0*(
     .  (-2d0*B1g*Tgg(3,4)
     .  +Tgg4(3,4,3,3)+Tgg4(3,4,4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(3,4,3,4))
     .  *dlog(scale**2/q2)
     .  +(Tgg(1,1)+Tgg(2,2))*Tgg(3,4)
     .  *0.5d0*(dlog(scale**2/q2))**2)

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -ggQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tgg4(3,4,i1,j1)
      enddo
      enddo


ch    pi**2/12 term

      Facdif=Facdif+ggQQv_0/2d0
     .        *pisq/12d0*Tgg(3,4)*(Tgg(1,1)+Tgg(2,2))


      Facdif=Facdif
ch   New Finite term comiing from ep^2 term of the Born
     .   -ggQQv_2/2d0*Tgg(3,4)*(
     .   +ca+ca)
ch
     .   -M0T34eps/2d0*(
     .   +cf+cf
     .   -2d0*B1g)
     .  -M0T34T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T34T13eps/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T34T23eps/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0T34eps/2d0*(xn+xn)*dlog(scale**2/q2)

      T34H1gg=FinVirtT34gg(sn,tn)


      T34H1gg=gsq**2/4d0*T34H1gg
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*PoleEp1

      T34H1gg=(T34H1gg-Facdif)*avegg
csk      T34H1gg=(0.d0-saveFacdif)*avegg

      return
      end

      double precision function T13H1gg(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
c      include 'breit.f'
      include 'rescoeff.f'
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision ggQQv_0,ggQQv_1,ggQQv_2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision M0T34eps,M0T13eps,M0T23eps,M0sqeps
      double precision sn,tn,FinVirtgg,FinvirtT13gg
      double precision M0T13T34eps,M0T13T13eps,M0T13T23eps
      double precision saveFacDif


      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)

      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps), _2 for O(eps^2)
      ggQQv_0=gsq**2/16d0/pi**2*M0sq

      ggQQv_1=-16d0/xn*V*(t1**2+t1*t2+t2**2)
     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
     . *gsq**2/4d0

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0


      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)*gsq**2
      M0T13eps=Tggeps(1,3)*gsq**2
      M0T13T34eps=Tgg4eps(1,3,3,4)*gsq**2
      M0T13T13eps=Tgg4eps(1,3,1,3)*gsq**2
      M0T13T23eps=Tgg4eps(1,3,2,3)*gsq**2
      
      PoleEp2=-ggQQv_0/2d0*Tgg(1,3)*(Tgg(1,1)+Tgg(2,2))
      PoleEp1=
     .   -ggQQv_0/2d0*(
     .   +Tgg4(1,3,3,3)+Tgg4(1,3,4,4)
     .  -2d0*B1g*Tgg(1,3)
     .  -(Tgg(1,1)+Tgg(2,2))*Tgg(1,3)*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(1,3,3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tgg4(1,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tgg4(1,3,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tgg4(1,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tgg4(1,3,2,4)) 
ch    this is the additional 1/ep pole
     .  -ggQQv_1/2d0*(xn+xn)*Tgg(1,3)

ch    Begin the IR finite terms

      Facdif=
     .       -ggQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tgg4(1,3,3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tgg4(1,3,3,3)+Tgg4(1,3,4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +ggQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tgg4(1,3,i1,j1)
      enddo
      enddo

ch   End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif

ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-ggQQv_0/2d0*(
     .  (-2d0*B1g*Tgg(1,3)
     .  +Tgg4(1,3,3,3)+Tgg4(1,3,4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(1,3,3,4))
     .  *dlog(scale**2/q2)
     .  +(Tgg(1,1)+Tgg(2,2))*Tgg(1,3)
     .  *0.5d0*(dlog(scale**2/q2))**2)

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -ggQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tgg4(1,3,i1,j1)
      enddo
      enddo


ch    pi**2/12 term

      Facdif=Facdif+ggQQv_0/2d0
     .        *pisq/12d0*Tgg(1,3)*(Tgg(1,1)+Tgg(2,2))


      Facdif=Facdif
ch   New Finite term comiing from ep^2 term of the Born
     .   -ggQQv_2/2d0*Tgg(1,3)*(
     .   +ca+ca)
ch
     .   -M0T13eps/2d0*(
     .   +cf+cf
     .   -2d0*B1g)
     .  -M0T13T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T13T13eps/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T13T23eps/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0T13eps/2d0*(xn+xn)*dlog(scale**2/q2)

      T13H1gg=FinVirtT13gg(sn,tn)


      T13H1gg=gsq**2/4d0*T13H1gg
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*PoleEp1

      T13H1gg=(T13H1gg-Facdif)*avegg
csk      T13H1gg=(0.d0-saveFacdif)*avegg

      return
      end

      double precision function T14H1gg(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
c      include 'breit.f'
      include 'rescoeff.f'
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro
      double precision ggQQv_0,ggQQv_1,ggQQv_2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision M0T34eps,M0T14eps,M0T23eps,M0sqeps
      double precision sn,tn,FinVirtgg,FinvirtT14gg
      double precision M0T14T34eps,M0T14T13eps,M0T14T23eps
      double precision saveFacDif


      sn=2d0*dot(p,1,2)
      tn=2d0*dot(p,1,3)+mt**2
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)

      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))
      t2=1d0-t1

c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps), _2 for O(eps^2)
      ggQQv_0=gsq**2/16d0/pi**2*M0sq
      
      ggQQv_1=-16d0/xn*V*(t1**2+t1*t2+t2**2)
     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
     . *gsq**2/4d0
      
      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0

      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)*gsq**2
      M0T14eps=Tggeps(1,4)*gsq**2
      M0T14T34eps=Tgg4eps(1,4,3,4)*gsq**2
      M0T14T13eps=Tgg4eps(1,4,1,3)*gsq**2
      M0T14T23eps=Tgg4eps(1,4,2,3)*gsq**2
      
      PoleEp2=-ggQQv_0/2d0*Tgg(1,4)*(Tgg(1,1)+Tgg(2,2))
      PoleEp1=
     .   -ggQQv_0/2d0*(
     .   +Tgg4(1,4,3,3)+Tgg4(1,4,4,4)
     .  -2d0*B1g*Tgg(1,4)
     .  -(Tgg(1,1)+Tgg(2,2))*Tgg(1,4)*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(1,4,3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tgg4(1,4,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tgg4(1,4,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tgg4(1,4,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tgg4(1,4,2,4)) 
ch    this is the additional 1/ep pole
     .  -ggQQv_1/2d0*(xn+xn)*Tgg(1,4)

ch    Begin the IR finite terms

      Facdif=
     .       -ggQQv_0/2d0*(
     .      1d0/v34*L34(p)*Tgg4(1,4,3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tgg4(1,4,3,3)+Tgg4(1,4,4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +ggQQv_0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tgg4(1,4,i1,j1)
      enddo
      enddo

ch   End the IR finite terms
csk   to switch off Ft1
      Facdif = 0.d0
csk      saveFacDif = Facdif


ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-ggQQv_0/2d0*(
     .  (-2d0*B1g*Tgg(1,4)
     .  +Tgg4(1,4,3,3)+Tgg4(1,4,4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg4(1,4,3,4))
     .  *dlog(scale**2/q2)
     .  +(Tgg(1,1)+Tgg(2,2))*Tgg(1,4)
     .  *0.5d0*(dlog(scale**2/q2))**2)

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -ggQQv_0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tgg4(1,4,i1,j1)
      enddo
      enddo


ch    pi**2/12 term

      Facdif=Facdif+ggQQv_0/2d0
     .        *pisq/12d0*Tgg(1,4)*(Tgg(1,1)+Tgg(2,2))


      Facdif=Facdif
     .             -
ch   New Finite term comiing from ep^2 term of the Born
     .   ggQQv_2/2d0*Tgg(1,4)*(
     .   +ca+ca)
ch
     .   -M0T14eps/2d0*(
     .   +cf+cf
     .   -2d0*B1g)
     .  -M0T14T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T14T13eps/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T14T23eps/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0T14eps/2d0*(xn+xn)*dlog(scale**2/q2)

      T14H1gg=FinVirtT14gg(sn,tn)


      T14H1gg=gsq**2/4d0*T14H1gg
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*PoleEp1
      
      T14H1gg=(T14H1gg-Facdif)*avegg
csk      T14H1gg=(0.d0-saveFacdif)*avegg

      return
      end


ch contributions from the azimuthal average
cch  G*G interference
      double precision function msqGGav(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'rescoeff.f'
            integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,u1,factor

      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
      v34=dsqrt(1d0-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)

      msqGGav=(-4d0*gsq**2*(-1d0 + xn**2)*
     .    ((-1d0 + xn**2)*s12**2 + 2d0*xn**2*s12*u1 + 2d0*xn**2*u1**2)*
     .    (2d0*mt**4*s12**2 + 
     .     2d0*mt**2*s12*(2d0*pt**2*s12 + u1*(s12 + u1)) + 
     .    (pt**2*s12 + u1*(s12 + u1))*(3d0*pt**2*s12 + u1*(s12 + u1))))/
     .    (xn*s12**2*u1**2*(s12 + u1)**2)

      msqGGav=msqGGav*avegg

      return
      end

ch  G*D interference
      double precision function msqDGav(p)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'rescoeff.f'
            integer j,k,naem,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,u1,I1,I2,I3,I3c,I31,I32,a
      double precision bbb
      s12=2d0*dot(p,1,2)
      s34=2d0*dot(p,3,4)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)
      a=pt/mt
      if(1d0-4d0*mtrans**2/q2.le.0d0)then
        I3c=0d0
      else
        I3c=dsqrt(1d0-4d0*mtrans**2/q2)
      endif
      

      I1=(1d0 - 1d0/a**2*dlog(1d0 + a**2))/4d0

      I2=(-1d0 + (a**2 + 1d0)/a**2*dlog(1d0 + a**2))/4d0

      I31=dlog((1d0 + v34)/(1d0 - v34))/v34/s34

      I32=((1d0 - dsqrt(1d0 - v34**2))/2d0/v34*
     .      dlog((1d0 + v34)/(1d0 - v34))-
     .     dlog(mtrans**2/mt**2)-
     .     I3c*dlog((1d0 + I3c)/(1d0 - I3c)))/pt**2

      I3=s34*(I31 - I32)/4d0

      msqDGav=(-8d0*gsq**2*mt**2*(-1d0 + xn**2)*pt**2*
     .    ((I3 + 2d0*I2*xn**2*(-2d0 + xn**2) - 
     .           2d0*I1*(-1d0 + xn**2)**2)*s12**2 -
     .      2d0*xn**2*(-2d0*I1 + I3 + 2d0*(I1 - I2)*xn**2)*s12*u1 - 
     .      2d0*xn**2*(-2d0*I1 + I3 + 2d0*(I1 - I2)*xn**2)*u1**2))/
     .  (xn**2*u1**2*(s12 + u1)**2)

      msqDGav=-msqDGav*avegg

      return
      end
      



