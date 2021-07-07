C Originally lowint_incldip.f
C Version modified for W production to include factor (1+(as/Pi)Hst)
C Valid only for qq->W(Z) production


      double precision function lowintHst(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
C
      include 'qcdcouple.f'
C
c --- DSW. To store flavour information :
      include 'nflav.f'

c --- DSW.
      integer pflav,pbarflav
c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx10(-nf:nf),fx20(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart
      double precision wgt,msqc(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      double precision xx(2),flux,vol,vol_mass,vol3_mass,BrnRat
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      logical bin,first,includedipole
      logical creatent,dswhisto
CC
      double precision x1p,x2p,fx1p(-nf:nf),fx2p(-nf:nf)
      double precision C1qqdelta,asopi,z1,z2,alfa,beta,cut,diff
      double precision tdelta,tH1st,tH1stF,xx10,xx20,beta0
      double precision Pqqint,Cqq,Cqg,Pqq,Pqg,LF,LR,q2,dot
      external Pqqint,Cqq,Cqg,Pqq,Pqg

      integer order
      common/nnlo/order
CC
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/bypart/lord_bypart
      common/outputflags/creatent,dswhisto
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

      W=sqrts**2

c--- processes that use "gen2"     
      if((case .eq. 'W_1jet').or.(case .eq. 'Z_1jet')) then   
       npart=2
       call gen2(r,p,pswt,*999)
      else
       write(*,*)'Wrong process !'
       stop
      endif

CC   Compute Q2
 
      q2=2*dot(p,3,4)    
        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      nvec=npart+2
      call dotem(nvec,p,s)

      call masscuts(s,*999)
      
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                                 

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out

      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

c--- Calculate the required matrix elements    
  
      if     (case .eq. 'W_1jet') then
        call qqb_w(p,msqc)
      elseif (case .eq. 'Z_1jet') then
        call qqb_z(p,msqc)
      else
        write(6,*) 'Unimplemented process: case=',case
        stop 
      endif

      
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo 


      currentPDF=0
            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif

C     Resummation coefficients

      C1qqdelta=(pi**2-8)/3d0

      beta0=(33-2*nf)/12d0

      asopi=ason2pi*2
      
CC   LR,LF
    
      LR=dlog(q2/scale**2)
      LF=dlog(q2/facscale**2) 


C Scaled momentum fractions

      cut=1d-7
   
      beta=cut+(1-cut)*r(5)
      alfa=cut+(1-cut)*r(6)


      xx10=xx(1)
      xx20=xx(2)

      z1=xx10**beta
      z2=xx20**alfa


c--- calculate PDF's  

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
      tH1st=0d0
      tH1stF=0d0


      do j=-nflav,nflav
      do k=-nflav,nflav

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

      enddo
      enddo


      xmsq=tdelta+asopi*(tH1st+LF*tH1stF)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      lowintHst=flux*pswt*xmsq/BrnRat


      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

    

      val=lowintHst*wgt
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
        wtmax=dabs(val)
      endif

      if (bin) then
        val=val/dfloat(itmx)
c ---   DSW. If the user has not selected to generate
c ---   events, still call nplotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
          call nplotter(pjet,val,0)
        endif
      endif

c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = dabs(val)
        if (ran2() .lt. (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs.lt.wtmax) then
            newwt = 1d0
          else
            newwt = wtabs/wtmax
          endif
          if (newwt .gt. 1.0d0) then
            write(6,*) 'WARNING : lowintHst : event with |weight| > 1.',
     +            ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*dsign(1d0,val)
          call nplotter(pjet,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
          call decide_flavour(pflav,pbarflav)
          call storeevent(pjet,newwt,pflav,pbarflav)
        endif
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


