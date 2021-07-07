      subroutine dyinit
C
C     Main initialization routine
C
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      include 'efficiency.f'
      include 'limits.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'phasemin.f'

C
      include 'masses.f'
      include 'zerowidth.f'
   
      character*30 runstring
      common/runstring/runstring
    
C
      logical creatent,dswhisto
      double precision rtsmin,sqrts,p1ext(4),p2ext(4),
     . p(mxpart,4),val
      integer j,k,nproc
      common/rtsmin/rtsmin
      common/energy/sqrts
      common/pext/p1ext,p2ext
      common/nproc/nproc
      data p/mxpart*3d0,mxpart*4d0,mxpart*0d0,mxpart*5d0/


      write(*,*) 'CCCCCCCCCCCCCC   DYNNLO, Version 1.4   CCCCCCCCCCCCC'
      write(*,*) 'C                                                  C'
      write(*,*) 'C  Written by M.Grazzini          August 2013      C'
      write(*,*) 'C                                                  C'
      write(*,*) 'C  Please refer to:                                C'
      write(*,*) 'C  S.Catani et al., PRL 103 (2009) 082001          C'
      write(*,*) 'C  S.Catani, M.Grazzini, PRL 98 (2007) 222002      C'      
      write(*,*) 'C                                                  C'

      call setup

C Initialize efficiency variables      
      njetzero=0
      ncutzero=0
      ntotzero=0
      ntotshot=0
      

C Set-up incoming beams and PS integration cut-offs
      rtsmin=min(rtsmin,dsqrt(wsqmin+cutoff))
   

c      if(zerowidth)then
c       rtsmin=wmass
c       if(nproc.eq.3)rtsmin=zmass
c      endif

      taumin=(rtsmin/sqrts)**2
      xmin=taumin

      p1ext(4)=-half*sqrts
      p1ext(1)=0d0
      p1ext(2)=0d0
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0d0
      p2ext(2)=0d0
      p2ext(3)=+half*sqrts

* Initialize all histograms
* npart=6 is a dummy value, to ensure that all histograms are included
      npart=6
      val=1d-15   
      call plotter(p,val,1)
       
      do j=1,mxpart
      do k=1,4
      p(j,k)=0d0
      enddo
      enddo 

           

      return
      end
            
