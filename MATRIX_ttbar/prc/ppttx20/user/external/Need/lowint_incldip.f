CC    Used to compute Higgs or W(Z) cross section at NLO only

      double precision function lowint(r,wgt)
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
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'dynamicscale.f'

CC
CC    Variables to be passed to the counterterm
CC 
      double precision qt2,qq2,shat,dot
      common/count/qt2,qq2,shat

c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
c      double precision msqa(-nf:nf,-nf:nf),n(4)
      double precision xx(2),flux,vol,vol_mass,vol3_mass,BrnRat
      logical bin,first,includedipole
      logical creatent,dswhisto
      integer nproc
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/nproc/nproc
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
      lowint=0d0

      W=sqrts**2
  

    
      npart=3     
c      call gen3(r,p,pswt,*999)
c      qq2=2*dot(p,3,4)
      call gen3m(r,p,mt,mt,0d0,pswt,*999)
      qq2=2*mt**2+2*dot(p,3,4)
      qt2=p(5,1)**2+p(5,2)**2   


      shat=2*dot(p,1,2)

      nvec=npart+2


C     Dynamic scale

      if(dynamicscale) call scaleset(qq2)


      call dotem(nvec,p,s)
      
      call masscuts(s,*999)
      
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                                 

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts


c--- Calculate the required matrix elements      

c       if(nproc.eq.3) then
c        call qqb_z_g(p,msq)
c        else
c        call qqb_w_g(p,msq)
c       endif     
       call qqb_QQb_g(p,msq) 


            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      

  777 continue    
      xmsq=0d0

      
c--- calculate PDF's  

      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)



CC     TIENI SOLO uubar
c        do j=-nf,1
c        fx1(j)=0d0
c        enddo
c        do j=3,nf
c        fx1(j)=0d0
c        enddo
c        do j=-nf,-3
c        fx2(j)=0d0
c        enddo
c        do j=-1,nf
c        fx2(j)=0d0
c        enddo
CC


      do j=-nf,nf
      do k=-nf,nf    
      
      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      xmsq=xmsq+xmsqjk
      
      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif


      
 20   continue
      enddo
      enddo


      lowint=flux*pswt*xmsq/BrnRat

      

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)


      val=lowint*wgt
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
        wtmax=dabs(val)
      endif

      if (bin) then
        val=val/dfloat(itmx)
c ---   DSW. If the user has not selected to generate
c ---   events, still call plotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
          call plotter(pjet,val,0)
        endif
      endif


      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


