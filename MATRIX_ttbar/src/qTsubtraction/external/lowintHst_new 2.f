C     Version that allows to separate the channels

      double precision function lowintHst(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'phasemin.f'

C
      include 'qcdcouple.f'
      include 'rescoeff.f'


c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,flgq
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx10(-nf:nf),fx20(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart
      double precision wgt,msqc(-nf:nf,-nf:nf)
      double precision xx(2),flux,BrnRat
      logical bin,first,includedipole
CC
      logical cuts
      double precision x1p,x2p,fx1p(-nf:nf),fx2p(-nf:nf)
      double precision asopi,z1,z2,alfa,beta,cut,diff
      double precision tdelta,tH1st,tH1stF,xx10,xx20,tH2st
      double precision tgaga,tcga,tgamma2
      double precision diff10,diff20,diffc10,diffc20,diffg10,diffg20
      double precision diff1f,diff2f,diffg1f,diffg2f,diffc1f,diffc2f
      double precision Pggreg,D0int,D1int,Cgq,Pgq,LF,LR
      double precision dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision Pqqint,Cqq,Cqg,Pqq,Pqg
      double precision C2qqreg,C2qqp,C2qqb,C2qg

C
      double precision beta1,H2qqdelta,H2qqD0
      common/Hstcoeff/beta1,H2qqdelta,H2qqD0

      integer order,a,b
      common/nnlo/order
CC
      integer jets,ndec,nproc
      common/parts_int/jets
      common/nproc/nproc
CC
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
     

      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      lowintHst=0d0

C     The number of jets is zero for this peice

      jets=0      

C

      W=sqrts**2



      npart=2
      call gen2(r,p,pswt,*999)




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      nvec=npart+2
      call dotem(nvec,p,s)

      call masscuts(s,*999)
      
                                                

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out

      if(cuts(p,0) .eqv. .true.) then
        goto 999
      endif
      
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts


c--- Calculate the required matrix element      


      if(nproc.eq.3) then
       call qqb_z(p,msqc)
      else
       call qqb_w(p,msqc)
      endif

            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0


      asopi=ason2pi*2


C     Compute Q2

      q2=2*dot(p,3,4)


      LF=dlog(q2/facscale**2)
      LR=dlog(q2/scale**2)


C Scaled momentum fractions

      cut=1d-8
   
C ndim here is 6 as for H->2gamma


      beta=cut+(1-cut)*r(ndim-1)
      alfa=cut+(1-cut)*r(ndim)


      xx10=xx(1)
      xx20=xx(2)

      z1=xx10**beta
      z2=xx20**alfa



c--- calculate PDF's  


      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)

       if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
       endif

       if(ggonly) then
        do j=1,nf
        fx10(j)=0d0
        fx10(-j)=0d0
        fx20(j)=0d0
        fx20(-j)=0d0
        fx1p(j)=0d0
        fx1p(-j)=0d0
        fx2p(j)=0d0
        fx2p(-j)=0d0   
        enddo
       endif

        flgq=1
        if(gqonly)flgq=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
      tH1st=0d0
      tH1stF=0d0
      tH2st=0d0

      do j=-nf,nf
      do k=-nf,nf

      if(msqc(j,k).eq.0d0) goto 75


C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq

      if(order.eq.0) goto 75

C     Start H1st: to be used later

C     H1st delta term

      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq

C     H1st: non delta terms, first leg


      tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg


      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg      


      diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq

c     gammaqq and gammaqg: second leg   


      diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq

      if(order.eq.1) goto 75

CC    End of H1st

CC    Start H2 contribution

CC    H2st gg contribution

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flgq

CC    H2st qqbar contribution from C1*C1 (without delta term)

C     regular*regular

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flgq

C     regular-delta

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq       

      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq       


CC    H2st qg contribution from C1*C1

C     regular*regular

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


C     regular-delta

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)       

      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)    

CC    H2st qqbar channel: D0(z), first leg

      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*H2qqD0/(1-z1)

      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    H2st, qqbar channel: D0(z), second leg
      
      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*H2qqD0/(1-z2)

      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, first leg

      tH2st=tH2st+fx1p(j)*C2qqreg(z1)
     &                   *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, second leg

      tH2st=tH2st+fx2p(k)*C2qqreg(z2)
     &                   *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

CC    C2qg, first leg

      tH2st=tH2st+fx1p(0)*C2qg(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)

CC    C2qg, second leg

      tH2st=tH2st+fx2p(0)*C2qg(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)

CC    Cqqbar contribution: first leg

      tH2st=tH2st+fx1p(-j)*C2qqb(z1)
     &                    *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    Cqqbar contribution: second leg
  
      tH2st=tH2st+fx2p(-k)*C2qqb(z2)
     &                    *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq  

      do a=1,nf

CC    Cqqp contribution: first leg

      if(a.ne.abs(j)) then      
       tH2st=tH2st+(fx1p(a)+fx1p(-a))*
     &       C2qqp(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      endif

CC    Cqqp contribution: second leg

      if(a.ne.abs(k)) then      
       tH2st=tH2st+(fx2p(a)+fx2p(-a))*
     &       C2qqp(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif

      enddo


 75   continue

      enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 72   xmsq=tdelta


      if(order.eq.1)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF)
      elseif(order.eq.2)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF)
     &          +asopi**2*(tdelta*H2qqdelta+tH2st)
      endif     


      lowintHst=flux*pswt*xmsq/BrnRat


      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

    

      val=lowintHst*wgt



      if (bin) then
        val=val/dfloat(itmx)
CC      call plotter(pjet,val,0)
        call plotter(p,val,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqreg(z)
      implicit none
      real *8 Pi,Z2,Z3,myli2,myli3,z,CA,CF
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      z2=Pi**2/6

      Z3=1.20205690316d0

      CF=4d0/3
      CA=3d0
      nf=5

      C2qqreg=   
     & (CF*(-344+24*Pi**2+974*z-1600*CA*z+2052*CF*z+148*nf*z-60*Pi**2*z+
     & 108*CF*Pi**2*z-1188*z**2+1584*CA*z**2-4104*CF*z**2-72*nf*z**2+
     & 72*Pi**2*z**2-216*CF*Pi**2*z**2+830*z**3+16*CA*z**3+2052*CF*z**3-
     & 76*nf*z**3-60*Pi**2*z**3+108*CF*Pi**2*z**3-
     & 272*z**4+24*Pi**2*z**4+
     & 324*CA*z*z2-1728*CF*z*z2 - 648*CA*z**2*z2 + 3456*CF*z**2*z2+
     & 324*CA*z**3*z2-1728*CF*z**3*z2 + 1188*CA*z*Z3 + 864*CF*z*Z3-
     & 324*CA*z**3*Z3-864*CF*z**3*Z3 - 108*CA*z**2*dlog(1-z) +
     & 108*CF*z**2*dlog(1-z)+108*CA*z**3*dlog(1-z)-
     & 108*CF*z**3*dlog(1-z)-
     & 216*CA*z*z2*dlog(1-z) + 216*CF*z*z2*dlog(1-z) -
     & 216*CA*z**3*z2*dlog(1-z)+216*CF*z**3*z2*dlog(1-z)-252*z*dlog(z)+
     & 348*CA*z*dlog(z)-540*CF*z*dlog(z)-
     & 60*nf*z*dlog(z)+612*z**2*dlog(z)-
     & 432*CA*z**2*dlog(z)+1404*CF*z**2*dlog(z)-744*z**3*dlog(z)+
     & 996*CA*z**3*dlog(z)-1728*CF*z**3*dlog(z)-60*nf*z**3*dlog(z)+
     & 384*z**4*dlog(z)-144*dlog(1-z)*dlog(z)+360*z*dlog(1-z)*dlog(z)-
     & 216*CA*z*dlog(1-z)*dlog(z) + 648*CF*z*dlog(1-z)*dlog(z) -
     & 432*z**2*dlog(1-z)*dlog(z) + 432*CA*z**2*dlog(1-z)*dlog(z) -
     & 1296*CF*z**2*dlog(1-z)*dlog(z) + 360*z**3*dlog(1-z)*dlog(z) -
     & 216*CA*z**3*dlog(1-z)*dlog(z) + 648*CF*z**3*dlog(1-z)*dlog(z) -
     & 144*z**4*dlog(1-z)*dlog(z)+216*CA*z*dlog(1-z)**2*dlog(z) -
     & 324*CF*z*dlog(1-z)**2*dlog(z)+216*CA*z**3*dlog(1-z)**2*dlog(z) -
     & 324*CF*z**3*dlog(1-z)**2*dlog(z)+27*z*dlog(z)**2+
     & 99*CA*z*dlog(z)**2-
     & 162*CF*z*dlog(z)**2-18*nf*z*dlog(z)**2+108*CA*z**2*dlog(z)**2 -
     & 108*CF*z**2*dlog(z)**2+45*z**3*dlog(z)**2-9*CA*z**3*dlog(z)**2+
     & 108*CF*z**3*dlog(z)**2-18*nf*z**3*dlog(z)**2-72*z**4*dlog(z)**2-
     & 108*CF*z*dlog(1-z)*dlog(z)**2-108*CF*z**3*dlog(1-z)*dlog(z)**2-
     & 18*z*dlog(z)**3 + 18*CA*z*dlog(z)**3-
     & 18*CF*z*dlog(z)**3+18*z**3*dlog(z)**3 +
     & 18*CA*z**3*dlog(z)**3+18*CF*z**3*dlog(z)**3 -
     & 72*((-1 + z)**2*(2 + (-1+3*CA-6*CF)*z + 2*z**2) -
     & 3*(CA-CF)*z*(1 + z**2)*dlog(1-z)-3*(CA-3*CF)*z*(1+z**2)*dlog(z))*
     & myli2(z) + 216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z) +
     & 216*CA*z**3*myli3(1-z)-216*CF*z**3*myli3(1-z) -
     & 432*CA*z*myli3(z) + 1080*CF*z*myli3(z) -
     & 432*CA*z**3*myli3(z)+1080*CF*z**3*myli3(z)-1944*CF*z*Z3-
     & 216*CF*z**3*Z3))/(864*(-1 + z)*z)




        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqp(z)
      implicit none
      real *8 Pi,myli2,z,CF
      integer nf

      external myli2

      Pi=3.14159265358979d0


      CF=4d0/3
     

      C2qqp=(CF*(2*(-1+z)*(-172+143*z-136*z**2+6*Pi**2*(2-z+2*z**2))-
     &   12*(z*(-21 + 30*z - 32*z**2)+
     &  6*(-2+3*z-3*z**2+2*z**3)*dlog(1-z))*
     &  dlog(z)-9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1 + z)*dlog(z)**3-
     &  72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z)))/(864d0*z)

 
      return
      end   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqb(z)
      implicit none
      real *8 Pi,Z3,myli2,myli3,z,CF,CA,C2qqp
    

      external myli2,myli3,C2qqp

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      CF=4d0/3   
      CA=3d0


       C2qqb=C2qqp(z)+
     &  (CF*(-CA+2*CF)*(45-3*Pi**2-2*Pi**2*z-45*z**2+Pi**2*z**2+9*
     &  dlog(z)+42*z*dlog(z)+33*z**2*dlog(z)+12*dlog(1-z)*dlog(z)-
     &  12*z**2*dlog(1-z)*dlog(z)-dlog(z)**3-z**2*
     &  dlog(z)**3+2*Pi**2*dlog(1+z) +
     &  2*Pi**2*z**2*dlog(1+z)-12*dlog(z)*
     &  dlog(1+z)-24*z*dlog(z)*dlog(1+z)-
     &  12*z**2*dlog(z)*dlog(1+z)+6*dlog(z)**2*dlog(1+z)+
     &  6*z**2*dlog(z)**2*dlog(1+z)-4*dlog(1+z)**3-4*z**2*dlog(1+z)**3-
     &  12*((1+z)**2+(1+z**2)*dlog(z))*myli2(-z)-
     &  12*(-1+z**2+dlog(z)+z**2*dlog(z))*myli2(z)+36*myli3(-z)+
     &  36*z**2*myli3(-z)+24*myli3(z)+24*z**2*myli3(z)+
     &  24*myli3(1d0/(1+z))+24*z**2*myli3(1d0/(1+z))-
     &  18*Z3-18*z**2*Z3))/(48*(1+z))


      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function C2qg(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3
      external myli2,myli3


      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0




        C2qg= (-35*CA)/48. - (13*CF)/32. +
     -  (43*CA)/(108.*z) + (43*CA*z)/48. +
     -   (43*CF*z)/32. + (CA*Pi**2*z)/24. + (CF*Pi**2*z)/12. -
     -   (149*CA*z**2)/216. - (5*CF*z**2)/4. - (CF*Pi**2*z**2)/12. -
     -   (3*CA*Z3)/8. + CF*Z3 - 2*CF*z*Z3 - (3*CA*z**2*Z3)/4. +
     -   2*CF*z**2*Z3 - (CA*Pi**2*Log(1-z))/24.
     -   + (3*CA*z*Log(1-z))/16. -
     -   (3*CF*z*Log(1-z))/16. - (CA*Pi**2*z*Log(1-z))/12. -
     -   (CA*z**2*Log(1-z))/4. + (CF*z**2*Log(1-z))/4. -
     -   (CA*Pi**2*z**2*Log(1-z))/12. - (CA*z*Log(1-z)**2)/8. +
     -   (CF*z*Log(1-z)**2)/8. + (CA*z**2*Log(1-z)**2)/8. -
     -   (CF*z**2*Log(1-z)**2)/8. - (CA*Log(1-z)**3)/48. -
     -   (CF*Log(1-z)**3)/48. - (CA*z*Log(1-z)**3)/8. +
     -   (CF*z*Log(1-z)**3)/24. - (CA*z**2*Log(1-z)**3)/24. -
     -   (CF*z**2*Log(1-z)**3)/24. +(7*CA*Log(z))/24.+(CF*Log(z))/4.+
     -  (CF*Pi**2*Log(z))/48.-(5*CA*z*Log(z))/12.+(15*CF*z*Log(z))/32.-
     -   (CA*Pi**2*z*Log(z))/12. - (CF*Pi**2*z*Log(z))/24. +
     -   (17*CA*z**2*Log(z))/18. - (CF*z**2*Log(z))/4. +
     -   (CF*Pi**2*z**2*Log(z))/24. - (CF*z*Log(1-z)*Log(z))/4. +
     -   (CF*z**2*Log(1-z)*Log(z))/4. + (CF*Log(1-z)**2*Log(z))/16.-
     -   (CF*z*Log(1-z)**2*Log(z))/8.
     -   + (CF*z**2*Log(1-z)**2*Log(z))/8. -
     -   (CA*Log(z)**2)/32. + (CF*Log(z)**2)/64.+(CA*z*Log(z)**2)/8.+
     -   (3*CF*z*Log(z)**2)/16. - (11*CA*z**2*Log(z)**2)/24. -
     -   (CF*z**2*Log(z)**2)/8. - (CF*Log(1-z)*Log(z)**2)/16. +
     -   (CA*z*Log(1-z)*Log(z)**2)/2.+(CF*z*Log(1-z)*Log(z)**2)/8.-
     -   (CF*z**2*Log(1-z)*Log(z)**2)/8. + (CA*Log(z)**3)/48. -
     -   (CF*Log(z)**3)/96. + (CA*z*Log(z)**3)/24.
     -   + (CF*z*Log(z)**3)/48.-
     -   (CF*z**2*Log(z)**3)/24. + (CA*Pi**2*Log(1 + z))/16. +
     -   (CA*Pi**2*z*Log(1 + z))/8. + (CA*Pi**2*z**2*Log(1 + z))/8. +
     -   (CA*Log(1-z)**2*Log(1 + z))/8. +
     -   (CA*z*Log(1-z)**2*Log(1 + z))/4. +
     -   (CA*z**2*Log(1-z)**2*Log(1 + z))/4. +
     -   (CA*z*Log(z)*Log(1 + z))/4. + (CA*z**2*Log(z)*Log(1 +z))/4. +
     -   (CA*Log(z)**2*Log(1 + z))/16. + (CA*z*Log(z)**2*Log(1+z))/8.+
     -   (CA*z**2*Log(z)**2*Log(1 + z))/8. -
     -   (CA*Log(1-z)*Log(1 + z)**2)/8. -
     -   (CA*z*Log(1-z)*Log(1 + z)**2)/4. -
     -   (CA*z**2*Log(1-z)*Log(1 + z)**2)/4. + (CA*myli2(1-z))/4d0-
     -   (CA*myli2(1-z))/(6.*z) - CA*z*myli2(1-z) +
     -   (11*CA*z**2*myli2(1-z))/12.
     -    - (CA*Log(1-z)*myli2(1-z))/8. +
     -   (CF*Log(1-z)*myli2(1-z))/8. +
     -   (CA*z*Log(1-z)*myli2(1-z))/4. -
     -   (CF*z*Log(1-z)*myli2(1-z))/4. -
     -   (CA*z**2*Log(1-z)*myli2(1-z))/4. +
     -   (CF*z**2*Log(1-z)*myli2(1-z))/4.
     -    - (CF*Log(z)*myli2(1-z))/8. +
     -   (CA*z*Log(z)*myli2(1-z))/2. + (CF*z*Log(z)*myli2(1-z))/4d0-
     -   (CF*z**2*Log(z)*myli2(1-z))/4. + (CA*z*myli2(-z))/4. +
     -   (CA*z**2*myli2(-z))/4. - (CA*Log(z)*myli2(-z))/8. -
     -   (CA*z*Log(z)*myli2(-z))/4. - (CA*z**2*Log(z)*myli2(-z))/4. +
     -   (CA*myli3(1-z))/8. - (CF*myli3(1-z))/8. -
     -   (CA*z*myli3(1-z))/4. + (CF*z*myli3(1-z))/4. +
     -   (CA*z**2*myli3(1-z))/4. - (CF*z**2*myli3(1-z))/4d0 +
     -   (3*CA*myli3(-z))/8. + (3*CA*z*myli3(-z))/4. +
     -   (3*CA*z**2*myli3(-z))/4. - (CF*myli3(z))/8. + CA*z*myli3(z)+
     -   (CF*z*myli3(z))/4. - (CF*z**2*myli3(z))/4. +
     -   (CA*myli3(1/(1 + z)))/4. + (CA*z*myli3(1/(1 + z)))/2d0 +
     -   (CA*z**2*myli3(1/(1 + z)))/2.
     -    - (CA*myli3((-1 + z)/(1 + z)))/4. -
     -   (CA*z*myli3((-1 + z)/(1 + z)))/2. -
     -   (CA*z**2*myli3((-1 + z)/(1 + z)))/2. +
     -   (CA*myli3((1 + z)/(-1 + z)))/4. +
     -   (CA*z*myli3((1 + z)/(-1 + z)))/2. +
     -   (CA*z**2*myli3((1 + z)/(-1 + z)))/2.




      return
      end


    
