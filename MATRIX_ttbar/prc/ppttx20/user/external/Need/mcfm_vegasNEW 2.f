CC    Modified to do full NNLO calculation
C
C     lowintHst is lowest order* Hst
C     
C     realvirt is real (V + 2 partons) + virtual (V + 1 parton)

C     Setting order=1 the calculation is done at NLO

      subroutine mcfm_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
************************************************************************
*                                                                      *
*  This routine should perform the sweeps of vegasnr                     *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                       *
*       myitmx  :  the number of vegasnr sweeps                          *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error
*                                                                      *
************************************************************************
      implicit none
      include 'gridinfo.f'
      include 'realwt.f'
      include 'scale.f'
      include 'facscale.f'
      include 'vegas_common.f'
      include 'PDFerrors.f'
      integer myitmx,myncall,myinit,i,j,k,nproc
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     . xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale
      character*4 part,mypart
      common/nproc/nproc
      common/part/part
      common/mypart/mypart
      common/bin/bin
      common/xreal/xreal,xreal2
      common/reset/reset,scalereset
      double precision lowint,virtint,realint
      double precision region(2*mxdim),lord_bypart(-1:1,-1:1)
C
      integer order
      common/nnlo/order
      double precision realvirt,lowintHst
C
      logical first,reset,scalereset,myreadin
      common/bypart/lord_bypart
      external lowint,virtint,realint
C      
      external realvirt,lowintHst
C
      data first/.true./
      save first

CC    HERE CHOOSE ORDER OF CALCULATION

      order=2

CC

           
c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for 'tota'
      sig=0d0
      sigr=0d0
      sigdk=0d0
      sd=0d0
      sdr=0d0
      sddk=0d0
      xreal=0d0
      xreal2=0d0
      
      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=0d0
      enddo
      enddo
      if (PDFerrors) then
        do i=0,maxPDFsets
          PDFxsec(i)=0d0
        enddo
      endif

c--- Controls behaviour of gen_njets: need to reset phase-space
c--- boundaries when going from virt to real (using tota)
c--- need to reset scale also, for special scalestart values
      reset=.false.
      scalereset=.false.

c--- Put the vegasnr parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Basic lowest-order integration
      if (part .eq. 'lord') then
       call boundregion(ndim,region)
       call vegasnr(region,ndim,lowint,myinit,myncall,myitmx,
     .               0,sig,sd,chi)
      endif

c--- Store value of part in mypart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mypart=part
      myscale=scale
      myfacscale=facscale
      
c--- If we're doing the tota integration, then set up the grid info
      if (mypart .eq. 'tota') then        
        if (first .and. (myinit .eq. 1)) then
c-- special input name for virtual grid
            ingridfile='dvegas_virt_'//ingridfile
            myreadin=readin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_virt.grid'          
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_virt.grid'
          endif
        endif
      endif        
      
CC   Virtint must contain lowest order * Hst factor

c--- Virtual integration should have one extra dimension
c--- (added and then taken away)
      if ((mypart .eq. 'virt') .or. (mypart .eq. 'tota'))  then
        part='virt'
        reset=.true.
        scalereset=.true.
        ndim=ndim-1
C        ndim=ndim+1
        call boundregion(ndim,region)
        call vegasnr(region,ndim,lowintHst,myinit,myncall,myitmx,
     .              0,sig,sd,chi)
C        ndim=ndim-1
        ndim=ndim+1
      endif

CC   Prepare the grid also for real only            

c--- If we're doing the tota integration, then set up the grid info
        if ((mypart .eq. 'tota') .or. (mypart .eq. 'real')) then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for real grid
          ingridfile(8:11)='real'
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_real.grid'          
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_real.grid'
          endif
        endif        
      endif 
 
CC Here real and virtual together

     
c--- Real integration should have three extra dimensions
c--- 'realwt' is a special option that in general should be false
c--- ('realwt' true samples the integral according to the
c---   unsubtracted real emission weight)
      if ((mypart .eq. 'real').or.(mypart.eq.'tota')) then
        part='real'
        scalereset=.true.
        ndim=ndim+3
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realvirt,myinit,myncall,myitmx,
     .              0,sigr,sdr,chi)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        if (realwt) then
          sigr=xreal
          sdr=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif
      

c--- calculate integration variables to be returned
      xinteg=sig+sigr+sigdk
      xerr=dsqrt(sd**2+sdr**2+sddk**2)      
      
c--- return part and scale to their real values
      part=mypart
      scale=myscale
      first=.false.
      
      return
      end
      
      
      subroutine boundregion(idim,region)
c--- Initializes integration region [0,1] for each variable
c--- in the idim-dimensional integration range
      implicit none
      include 'mxdim.f'
      integer i,idim
      double precision region(2*mxdim)
      
      do i=1,idim
      region(i)=0d0
      region(i+idim)=1d0
      enddo
      
      return
      end
      
      
