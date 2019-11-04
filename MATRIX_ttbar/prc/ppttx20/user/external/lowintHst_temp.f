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
      double precision H2streg
      double precision dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision Pqqint,Cqq,Cqg,Pqq,Pqg

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

      cut=1d-7
   
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
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)

C     Start H1st: to be used later

C     H1st delta term

      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)

C     H1st: non delta terms, first leg


      tH1st=tH1st+(fx1p(j)*Cqq(z1)+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg


      tH1st=tH1st+(fx2p(k)*Cqq(z2)+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg      


      diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)

c     gammaqq and gammaqg: second leg   


      diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)

CC    End of H1st

CC    Start H2 contribution

CC    H2st gg contribution

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)

c     H2st qqbar channel: D0(z), first leg

      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*H2qqD0/(1-z1)

      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)*fx10(j)*fx20(k)*msqc(j,k)

c     H2st, qqbar channel: D0(z), second leg
      
      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*H2qqD0/(1-z2)

      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)*fx10(j)*fx20(k)*msqc(j,k)



CC    Regular contribution: first leg

      do a=-nf,nf
      do b=-nf,nf

      tH2st=tH2st+0.5d0*fx1p(a)*
     &   (-dlog(xx10))*fx20(b)*H2streg(j,k,a,b,z1)*msqc(j,k)

CC    Regular contribution: second leg

      tH2st=tH2st+0.5d0*fx2p(b)*
     &   (-dlog(xx20))*fx10(a)*H2streg(j,k,a,b,z2)*msqc(j,k)


      enddo
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



      double precision function H2streg(i,j,a,b,z)
      implicit none
      real *8 z,H2qqbqqbreg,H2qqbqq,H2qqbqqp,H2qqbqg
      integer i,j,a,b
      external H2qqbqqbreg,H2qqbqq,H2qqbqqp,H2qqbqg

      H2streg=0d0

      if(a.eq.i.and.b.eq.j) then
       H2streg=H2qqbqqbreg(z)
      elseif(a.eq.i.and.b.eq.-j) then
       H2streg=H2qqbqq(z)
      elseif(a.eq.i.and.b.eq.0) then
       H2streg=H2qqbqg(z)
      else
       H2streg=H2qqbqqp(z)
      endif

      if(b.eq.j.and.a.eq.i) then
       H2streg=H2qqbqqbreg(z)
      elseif(b.eq.j.and.a.eq.-i) then
       H2streg=H2qqbqq(z)
      elseif(b.eq.j.and.a.eq.0) then
       H2streg=H2qqbqg(z)
      else
       H2streg=H2qqbqqp(z)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function H2qqbqqbreg(z)
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

         H2qqbqqbreg=
     -       (CF*((-2*nf*(-74 + 36*z + 38*z**2 + 30*(1 + z**2)*dlog(z)+
     -            9*(1 + z**2)*dlog(z)**2))/(-1 + z) +
     -       (2*(-1 + z)*(-172 + 143*z - 136*z**2 +
     -             6*Pi**2*(2 - z + 2*z**2)) -
     -          12*(z*(-21 + 30*z - 32*z**2) +
     -             6*(-2 + 3*z - 3*z**2 + 2*z**3)*dlog(1-z))*dlog(z)-
     -          9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1+z)*dlog(z)**3-
     -          72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z))/z +
     -       (CA*(-1600 + 1584*z + 16*z**2 + 324*z2 - 648*z*z2 +
     -           324*z**2*z2+1188*z3-324*z**2*z3-108*z*dlog(1-z)+
     -           108*z**2*dlog(1-z) - 216*z2*dlog(1-z) -
     -           216*z**2*z2*dlog(1-z) + 348*dlog(z)-432*z*dlog(z)+
     -           996*z**2*dlog(z) - 216*dlog(1-z)*dlog(z) +
     -           432*z*dlog(1-z)*dlog(z) - 216*z**2*dlog(1-z)*dlog(z)+
     -          216*dlog(1-z)**2*dlog(z)+216*z**2*dlog(1-z)**2*dlog(z)+
     -           99*dlog(z)**2 + 108*z*dlog(z)**2-9*z**2*dlog(z)**2 +
     -           18*dlog(z)**3 + 18*z**2*dlog(z)**3 +
     -           216*(-(-1 + z)**2 + (1 + z**2)*dlog(1-z) +
     -           (1 + z**2)*dlog(z))*myli2(z) +
     -           216*(1 + z**2)*myli3(1-z) - 432*myli3(z) -
     -           432*z**2*myli3(z)))/(-1 + z) -
     -       (18*CF*(-174 + 348*z - 174*z**2 + 96*z2 - 192*z*z2 +
     -            96*z**2*z2 - 48*z3 + 48*z**2*z3 - 6*z*dlog(1-z) +
     -            6*z**2*dlog(1-z) - 12*z2*dlog(1-z) -
     -            12*z**2*z2*dlog(1-z) + 24*dlog(z) - 78*z*dlog(z) +
     -            102*z**2*dlog(z) - 36*dlog(1-z)*dlog(z) +
     -            72*z*dlog(1-z)*dlog(z) - 36*z**2*dlog(1-z)*dlog(z) +
     -            18*dlog(1-z)**2*dlog(z)+18*z**2*dlog(1-z)**2*dlog(z)+
     -            9*dlog(z)**2 + 6*z*dlog(z)**2 - 6*z**2*dlog(z)**2+
     -            6*dlog(1-z)*dlog(z)**2 + 6*z**2*dlog(1-z)*dlog(z)**2+
     -            dlog(z)**3 - z**2*dlog(z)**3 +
     -            12*(-2*(-1 + z)**2 + (1 + z**2)*dlog(1-z)+
     -               3*(1 + z**2)*dlog(z))*myli2(z) +
     -            12*(1 + z**2)*myli3(1-z) - 60*myli3(z) -
     -            60*z**2*myli3(z) + 108*Z3 + 12*z**2*Z3))/
     -        (-1 + z)))/432d0


        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function H2qqbqqp(z)
      implicit none
      real *8 Pi,myli2,z,CF
      integer nf

      external myli2

      Pi=3.14159265358979d0


      CF=4d0/3
     

      H2qqbqqp=(CF*(2*(-1+z)*(-172+143*z-136*z**2+6*Pi**2*(2-z+2*z**2))-
     &   12*(z*(-21 + 30*z - 32*z**2)+
     &  6*(-2+3*z-3*z**2+2*z**3)*dlog(1-z))*
     &  dlog(z)-9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1 + z)*dlog(z)**3-
     &  72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z)))/(864d0*z)

 
      return
      end   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function H2qqbqq(z)
      implicit none
      real *8 Pi,Z3,myli2,myli3,z,CF,CA,H2qqbqqp
    

      external myli2,myli3,H2qqbqqp

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      CF=4d0/3   
      CA=3d0


       H2qqbqq=H2qqbqqp(z)+
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

      double precision function H2qqbqg(z)
      implicit none
      real *8 z

      H2qqbqg=0d0

      return
      end
