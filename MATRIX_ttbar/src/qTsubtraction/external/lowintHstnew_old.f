C     Version that allows to separate the channels
C     Scale dependence included

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
      include 'dynamicscale.f'

c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,l,nvec,flgq
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx10(-nf:nf),fx20(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart,pb(mxpart,4)
      double precision wgt,msqc(-nf:nf,-nf:nf)
      double precision xx(2),flux,BrnRat
      logical bin,first,includedipole
CC
      logical cuts
      double precision x1p,x2p,fx1p(-nf:nf),fx2p(-nf:nf)
      double precision asopi,z1,z2,alfa,beta,cut,diff
      double precision tdelta,tH1st,tH1stF,xx10,xx20,tH2st
      double precision tH1stFg, tH1stFq
ch      double precision tgaga,tcga,tgamma2,tdeltagg
      double precision tgagag,tcgag,tgamma2g,tdeltagg
      double precision tgagaq,tcgaq,tgamma2q,tdeltaqq
ch      double precision diff10,diff20,diffc10,diffc20,diffg10,diffg20
ch      double precision diff1f,diff2f,diffg1f,diffg2f,diffc1f,diffc2f
      double precision diff10g,diff20g,diffc10g,diffc20g,
     .                 diffg10g,diffg20g
      double precision diff1fg,diff2fg,diffg1fg,diffg2fg,
     .                 diffc1fg,diffc2fg
      double precision diff10q,diff20q,diffc10q,diffc20q,
     .                 diffg10q,diffg20q
      double precision diff1fq,diff2fq,diffg1fq,diffg2fq,
     .                 diffc1fq,diffc2fq
      double precision Pggreg,D0int,D1int,Cgq,Pgq,LF,LR,H2st,H2stgq
      double precision H2ggREG,dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision Pggggreg,Pgggq,Pgqqg,Pgqqq
      double precision CgqPqq,CgqPqg,P2gg,P2gq
      double precision Pqqint,Cqq,Cqg,Pqq,Pqg
      double precision C2qqreg,C2qqp,C2qqb,C2qg

      double precision H1st

      double precision H1qdelta,H1qbdelta,H1gdelta,aaa

      common/aa/aaa

      external Pggreg,D0int,Cgq,Pgq,H2stgq,H2ggREG,P2gg,P2gq,
     &         Pggggreg,Pgqqg,Pgqqq,Pgggq,CgqPqq,CgqPqg


ch variables and functions needed for the final state radiation

      double precision s34,beta34,L34,pt,bj,tH1stfs,mtrans
      common/colour/Tqq,Tgg,Tsinglet,Toktsym,Toktasym
      double precision Tqq(1:4,1:4),myli2,H1qqdelta
      double precision H1ggdelta
      double precision Tsinglet
      double precision Toktsym
      double precision Toktasym
      double precision Tgg(1:4,1:4)

ch    kinematical factor needed for the matrix element in the gg channel

      double precision kinasym

      double precision msqcgg

ch
      integer i1,j1,i1b
      double precision qqQQv_1
      double precision costheta,scpr13,p3m,p1m,betat,xs
      double precision H0S1,IReal,loop1

ch

ch
      double precision msqGGav,msqDGav

ch
C
      double precision diff1,diff2     
      double precision Pqqqq,Pqqqg,Pqggq,Pqggg
      double precision CqqPqq,CqqPqg,CqgPgq,CqgPgg
      double precision P2qg,P2qqV,P2qqbV,P2qqS
C
      double precision beta1,H2qqdelta,H2qqD0,
     &                 H2ggdelta,H2ggd0
      double precision beta04f
      common/Hstcoeff/beta1,H2qqdelta,H2qqD0,
     &                 H2ggdelta,H2ggd0

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

      beta04f=(33-2*4)/12d0

      npart=2

      call gen2m(r,p,pswt,*999)

      do i1=1,mxpart
      do j1=1,4
      if(i1.gt.2)then
      pb(i1,j1)=p(i1,j1)
      elseif(i1.eq.1)then
      pb(1,j1)=p(2,j1)
      elseif(i1.eq.2)then
      pb(2,j1)=p(1,j1)
      endif
      enddo
      enddo




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
ch      write(*,*)sqrts,p(1,4)




c--- Calculate the required matrix element      


c      if(nproc.eq.3) then
c       call qqb_z(p,msqc)
c      else
c       call qqb_w(p,msqc)
c      endif
c      p(1,1)=0d0
c      p(2,1)=0d0
c      p(3,1)=-2564.70376259d0
c      p(4,1)=2564.70376259d0
c      p(1,2)=0d0
c      p(2,2)=0d0
c      p(3,2)=2438.473059d0
c      p(4,2)=-2438.473059d0
c      p(1,3)=-3594.7309512d0
c      p(2,3)=3502.2177108d0
c      p(3,3)=-138.1232293d0
c      p(4,3)=230.6364697d0
c      p(1,4)=-3594.7309512d0
c      p(2,4)=-3502.2177108d0
c      p(3,4)=3546.07082245d0
c      p(4,4)=3550.877839707d0

c      p(1,1)=0d0
c      p(2,1)=0d0
c      p(3,1)=0d0
c      p(4,1)=0d0
c      p(1,2)=0d0
c      p(2,2)=0d0
c      p(3,2)=0d0
c      p(4,2)=0d0
c      p(1,3)=0d0
c      p(2,3)=0d0
c      p(3,3)=0d0
c      p(4,3)=0d0
c      p(1,4)=0d0
c      p(2,4)=0d0
c      p(3,4)=0d0
c      p(4,4)=0d0
      s34=2*dot(p,3,4)
      q2=s34+dot(p,3,3)+dot(p,4,4)
ch      write(*,*)dsqrt(q2)
c      q2=130565.511818644d0
      if(dynamicscale) call scaleset(q2)
      call qqb_QQb(p,msqc,0)
ch      write(*,*)msqc(1,-1),msqc(0,0),'l'
c      write(*,*)msqc(1,-1)
c      do j=-nf,nf
c      do k=-nf,nf
c      if(j.eq.-k)then
c      msqc(j,k)=2*dot(p,3,4)/sqrts**2
c      endif
c      enddo
c      enddo

ch    Colour operators for gg channel

      kinasym=(2d0*dot(p,1,3)-2d0*dot(p,2,3))
     .        /(2d0*dot(p,1,2))

      msqcgg=msqc(0,0)

      msqc(0,0)=msqc(0,0)*(Tsinglet+Toktsym+Toktasym*kinasym**2)

      Tgg(3,3)=cf*msqc(0,0)
      Tgg(4,4)=cf*msqc(0,0)
      Tgg(1,1)=ca*msqc(0,0)
      Tgg(2,2)=ca*msqc(0,0)
      Tgg(3,4)=((ca/2d0-cf)*(Toktsym+Toktasym*kinasym**2)
     .         -cf*Tsinglet)*msqcgg
      Tgg(1,3)=-ca/4d0*(Toktsym+Toktasym*kinasym**2
     .                 +Toktasym*kinasym+2d0*Tsinglet*kinasym
     .                 +Toktsym*kinasym)*msqcgg
      Tgg(1,4)=-ca/4d0*(Toktsym+Toktasym*kinasym**2
     .                 -Toktasym*kinasym-2d0*Tsinglet*kinasym
     .                 -Toktsym*kinasym)*msqcgg
      Tgg(2,3)=Tgg(1,4)
      Tgg(2,4)=Tgg(1,3)

ch
            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
c      flux=1d0/(2d0*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0


    


C     Compute Q2
c      s34=2*dot(p,3,4)
c      q2=s34+dot(p,3,3)+dot(p,4,4)
      if(q2.le.wsqmin.or.q2.gt.wsqmax) goto 999
      
ch
ch   variables needed for the final state radiation 
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
c      write(*,*)pt,dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)
      beta34=dsqrt(1d0-4d0*mt**4/s34**2)
      betat=dsqrt(1d0-4d0*mt**2/q2)
      xs=(1d0-betat)/(1d0+betat)
ch


C     Dynamic scale

c      if(dynamicscale) call scaleset(q2)
c       write(*,*)ason2pi*8*pi**2,'c',dsqrt(q2)

      asopi=ason2pi*2

      LF=dlog(q2/facscale**2)
      LR=dlog(q2/scale**2)
ch      write(*,*)dsqrt(q2),scale,facscale


C Scaled momentum fractions

      cut=1d-8
   
C ndim here is 6 as for H->2gamma


      beta=cut+(1-cut)*r(ndim-1)
      alfa=cut+(1-cut)*r(ndim)


      xx10=xx(1)
      xx20=xx(2)
ch      xx10=0.004614313709741852d0
ch      xx20=0.702034039721575d0
ch      alfa=0.565332618728161d0
ch      beta=0.931441737138846d0
ch      msqc(0,0)=0.828074800688665d0
ch      msqc(0,0)=0d0
ch     facscale=346.6d0
ch      scale=346.6d0
ch      asopi=0.03107362017927538d0
ch      LR=0.545685369560359d0
ch      LF=0.545685369560359d0


      z1=xx10**beta
      z2=xx20**alfa

      qqQQv_1=-64d0*gsq**2/4d0



c--- calculate PDF's  
ch      write(*,*)xx10,xx20
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

CC     TIENI SOLO uubar
c        do j=-nf,1
c        fx10(j)=0d0
c        fx1p(j)=0d0
c        enddo
c        do j=3,nf
c        fx10(j)=0d0
c        fx1p(j)=0d0
c        enddo
c        do j=-nf,-3
c        fx20(j)=0d0
c        fx2p(j)=0d0
c        enddo
c        do j=-1,nf
c        fx20(j)=0d0
c        fx2p(j)=0d0
c        enddo
CC


        flgq=1
        if(gqonly)flgq=0
ch        fx10(0)=0d0
ch        fx20(0)=0d0
ch        fx1p(0)=0d0
ch        fx2p(0)=0d0
ch        do j=1,nf
ch        fx10(j)=1d0
ch        fx10(-j)=0d0
ch        fx20(j)=1d0
ch        fx20(-j)=0d0
ch        fx1p(j)=1d0
ch        fx1p(-j)=0d0
ch        fx2p(j)=1d0
ch        fx2p(-j)=0d0
ch        enddo


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
      tdeltagg=0d0
      tdeltaqq=0d0
      tH1st=0d0
      H1st=0d0
      tH1stF=0d0
      tH1stFg=0d0
      tH1stFq=0d0
      tH1stfs=0d0
      tH2st=0d0
      diff10q=0d0
      diff20q=0d0
      diff1fq=0d0
      diff2fq=0d0

      diffc10q=0d0
      diffc1fq=0d0
      diffc20q=0d0
      diffc2fq=0d0

      diffg10q=0d0
      diffg1fq=0d0
      diffg20q=0d0
      diffg2fq=0d0

      diff10q=0d0
      diff20q=0d0
      diff1fq=0d0
      diff2fq=0d0

      diffc10g=0d0
      diffc1fg=0d0
      diffc20g=0d0
      diffc2fg=0d0

      diffg10g=0d0
      diffg1fg=0d0
      diffg20g=0d0
      diffg2fg=0d0
ch      tcga=0d0
      tcgag=0d0
      tcgaq=0d0
ch      tgamma2=0d0
      tgamma2g=0d0
      tgamma2q=0d0
ch      tgaga=0d0
      tgagag=0d0
      tgagaq=0d0
      spgq1=0d0
      spgq2=0d0

      tH2sp=0d0
      H0S1=0d0
      IReal=0d0
ch      write(*,*)msqc(1,-1),'rrrr'
      H1qdelta=H1qqdelta(p)/msqc(1,-1)
ch     .         +2d0*beta0*LR 
      H1qbdelta=H1qqdelta(pb)/msqc(1,-1)
ch     .          +2d0*beta0*LR
      H1gdelta=H1ggdelta(p)/msqc(0,0)
ch     .         +2d0*beta0*LR
ch      H1qdelta=1d0
ch      H1qbdelta=1d0
ch      H1gdelta=1d0
ch      H1gdelta=0d0
ch      tdelta=msqc(0,0)

c      tdelta=msqc(1,-1)*flgq*fx10(2)*fx20(-2)+
c     .       msqc(1,-1)*flgq*fx10(-2)*fx20(2)
ch      write(*,*)msqc(2,-2),'uuuuu'

ch    check against hnnlo and dynnlo
ch      do j=1,nf
ch       msqc(j,-j)=0.668647327076049d0
ch       msqc(-j,j)=0.668647327076049d0
ch      msqc(j,-j)=0d0
ch      msqc(-j,j)=0d0
ch      enddo
ch      C1ggdelta=0d0
ch      C1qqdelta=0d0
ch      xx10=0.005876431369874011d0
ch      xx20=0.608085455323661d0
ch      alfa=0.399067733190328d0
ch      beta=0.01925242711509566d0
ch      msqc(0,0)=1.10937459535996d0
ch      facscale=173.3d0
ch      scale=173.3d0
ch      asopi=0.03400021652486929
ch      LR=1.74434433962324d0
ch      LF=1.74434433962324d0
ch      C1ggdelta=0d0
ch3.400021652486929E-002 asopi  5.876431369874011E-003 xx10
ch  0.608085455323661      xx20  0.399067733190328      alfa
ch  1.925242711509566E-002 beta   173.300000000000        173.300000000000      
ch muf,mur   1.10937459535996       0.787492532154837      msqc
ch   2.03009978130627      LF   2.03009978130627      LR
ch-5.481736956984965.876431369874011      H1gdelta=2d0*C1ggdelta+2d0*beta0*LR
ch      H1qdelta=2d0*C1qqdelta
ch      H1qbdelta=2d0*C1qqdelta
ch      H1gdelta=2d0*C1ggdelta-2d0*beta0*LR

      do j=-nf,nf
      do k=-nf,nf
ch      write(*,*)j,k,'www'
ch      write(*,*)msqc(j,k),'ooooo'
ch      if(msqc(j,k).eq.0d0) goto 75

c      if (j.eq.0) goto 75

ch       write(*,*)'aaaaaaaaaaa'
C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      write(*,*)tdelta,'a',msqc(j,k)

      if(j.eq.-k.and.j.ne.0)then
        tdeltaqq=tdeltaqq+fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.eq.0.and.k.eq.0)then
        tdeltagg=tdeltagg+fx10(j)*fx20(k)*msqc(j,k)*flgq
      endif


ch      if(j.eq.0.and.k.eq.0)then
ch       tdeltagg=fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      endif

      if(order.eq.0) goto 75

C     Start H1st: to be used later

C     H1st delta term
c       if(j.eq.-k.and.j.ne.0)then
        if(j.gt.0.and.k.lt.0)then
         tH1st=tH1st+H1qdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch         tH1st=tH1st+H1qdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        elseif(j.lt.0.and.k.gt.0)then
         tH1st=tH1st+H1qbdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch         tH1st=tH1st+H1qbdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        elseif(j.eq.0.and.k.eq.0)then
ch     Remove the lR dependence to leave the hnnlo implementation on the
ch     scale-dependence unaltered
         tH1st=tH1st+H1gdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch         tH1st=tH1st+H1gdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        endif
chhhh    Check the 1-loop amplitude used in the counterterm

c       tH1st=tH1st+2d0*H1qqdelta(p)*flgq
c       endif

c       write(*,*)2*H1qqdelta(p)/gsq**2*4d0*3d0,'ddddd'

c       write(*,*)dsqrt(-dot(p,1,3)),'a'


c      write(*,*)H1qqdelta(p)
c      write(*,*)tdelta,'a'

ch    H1st: qq channel

C     H1st: non delta terms, first leg

      if(j.eq.-k.and.j.ne.0)then
       tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg

      
      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      endif
ch      write(*,*)tH1st,'tH1st0',j,'j',k,'k',msqc(j,k),'msqc'

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg

ch    gammagg: first leg is included      


c      diff=-dlog(xx10)
c     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
ch
ch      if(j.ne.0.and.k.ne.0)then
ch       diff=-dlog(xx10)
ch     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
ch      elseif(j.eq.0.and.k.eq.0)then
ch       diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*3/(1-z1)
ch     &    +fx1p(0)*Pggreg(z1))*flgq
ch      endif

ch
ch      if(j.ne.0.and.k.ne.0)then
ch       tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
ch       tH1stFq=tH1stFq+diff*fx20(k)
ch       tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      elseif(j.lt.0.and.k.gt.0)then
ch       tH1stFqb=tH1stFqb+diff*fx20(k)*msqc(j,k)
ch       tH1stFq=tH1stFq+diff*fx20(k)
ch       tH1stFqb=tH1stFqb-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      elseif(j.eq.0.and.k.eq.0)then
ch       tH1stF=tH1stF+diff*fx20(0)*msqc(j,k)
ch       tH1stF=tH1stF-3*D0int(xx10)*fx10(0)*fx20(0)*msqc(j,k)*flgq
ch      endif

ch

c     gammaqq and gammaqg: second leg   

ch    gammagg: second leg is included


c      diff=-dlog(xx20)
c     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
ch
ch      if(j.ne.0.and.k.ne.0)then
ch       diff=-dlog(xx20)
ch     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
ch      elseif(j.eq.0.and.k.eq.0)then
ch       diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*3/(1-z2)
ch     &    +fx2p(0)*Pggreg(z2))*flgq
ch      endif

ch
ch      if(j.ne.0.and.k.ne.0)then
ch       tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
ch       tH1stFq=tH1stFq+diff*fx10(j)
ch       tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      elseif(j.lt.0.and.k.gt.0)then
ch       tH1stFqb=tH1stFqb+diff*fx10(j)*msqc(j,k)
ch       tH1stFq=tH1stFq+diff*fx10(j)
ch       tH1stFqb=tH1stFqb-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
ch      elseif(j.eq.0.and.k.eq.0)then
ch       tH1stF=tH1stF+diff*fx10(0)*msqc(j,k)
ch       tH1stF=tH1stF-3*D0int(xx20)*fx10(0)*fx20(0)*msqc(j,k)*flgq
ch      endif

ch

c     gammagg: delta term, both legs

ch      if(j.eq.0.and.k.eq.0)then
ch       tH1stF=tH1stF+2*beta0*fx10(0)*fx20(0)*msqc(j,k)*flgq
ch      endif

CC    End of H1st
       

C     H1st: muf dependence (LF factor to be added at the end)

ch    qq channel

c     gammaqq and gammaqg: first leg      

      if(j.eq.-k.and.j.ne.0)then
      diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stFq=tH1stFq+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      tH1stFq=tH1stFq-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq



c     gammaqq and gammaqg: second leg   


      diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stFq=tH1stFq+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      tH1stFq=tH1stFq-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      
      elseif(j.eq.0.and.k.eq.0)then

ch    gg channel
c     gammagg: non delta terms, first leg    


      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*3/(1-z1)
     &    +fx1p(0)*Pggreg(z1))

      tH1stF=tH1stF+diff*fx20(0)*msqc(0,0)*flgq
      tH1stF=tH1stF-3*D0int(xx10)*fx10(0)*fx20(0)*msqc(0,0)*flgq
      tH1stFg=tH1stFg+diff*fx20(0)*msqc(0,0)*flgq
      tH1stFg=tH1stFg-3*D0int(xx10)*fx10(0)*fx20(0)*msqc(0,0)*flgq

c     gammagg: non delta terms, second leg    

      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*3/(1-z2)
     &    +fx2p(0)*Pggreg(z2))

      tH1stF=tH1stF+diff*fx10(0)*msqc(0,0)*flgq
      tH1stF=tH1stF-3*D0int(xx20)*fx10(0)*fx20(0)*msqc(0,0)*flgq
      tH1stFg=tH1stFg+diff*fx10(0)*msqc(0,0)*flgq
      tH1stFg=tH1stFg-3*D0int(xx20)*fx10(0)*fx20(0)*msqc(0,0)*flgq

c     gammagg: delta term, both legs


      tH1stF=tH1stF+2*beta0*fx10(0)*fx20(0)*msqc(0,0)*flgq
      tH1stFg=tH1stFg+2*beta0*fx10(0)*fx20(0)*msqc(0,0)*flgq
      endif


ch    Final state pieces
c      tH1stfs=tH1stfs+(2d0*L34(q2,beta34,mtrans)*Tqq(3,4)
c     .       +2*dlog(scale**2/mt**2)+4d0*dlog(1+pt**2/mt**2)*Tqq(3,3))
c     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq

c      do i1=1,2
c      do j1=3,4
c      tH1stfs=tH1stfs+
c     .       (dlog(scale**2/(-2*dot(p,i1,j1)))*
c     .       dlog(mt**2/(-2*dot(p,i1,j1)))
c     .       +1/2*dlog(scale**2/(-2*dot(p,i1,j1)))**2
c     .       -1/2*LR**2
c     .       -myli2(-1*pt**2/mt**2))*Tqq(i1,j1) 
c     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq 
c      enddo
c      enddo

c       write(*,*)1d0-4*dot(p,1,3)*dot(p,2,3)/q2/mt**2,'ch'
c       write(*,*)-pt**2/mt**2,'a'

c       write(*,*)4*dot(p,1,3)*dot(p,2,3)/q2/mt**2,'ch'
c       write(*,*)1d0+pt**2/mt**2,'a'
c       write(*,*)p(1,1),p(1,2),p(1,3),p(1,4),'t'
c       write(*,*)p(2,1),p(2,2),p(2,3),p(2,4),'tbar'
       
       p3m=dsqrt(p(3,3)**2+p(3,2)**2+p(3,1)**2)
       p1m=dsqrt(p(1,3)**2+p(1,2)**2+p(1,1)**2)
       scpr13=p(1,1)*p(3,1)+p(1,2)*p(3,2)+p(1,3)*p(3,3)
c       costheta=-scpr13/p3m/p1m
       costheta=-1d0/betat*(1d0+4*dot(p,1,3)/q2)
c       write(*,*)costheta
c       write(*,*)costheta
c       write(*,*)-p3m**2+p(3,4)**2-mt**2
c       write(*,*)-(1d0+betat**2)/betat*(
c     .           -myli2(-1d0*xs*(1d0-costheta)/(1d0+costheta))
c     .           +myli2(-1d0*(1d0-costheta)/(1d0+costheta)/xs)
c     .           +4d0*dlog(xs)*dlog(dsqrt(0.5d0*(1d0+costheta)))),'ch'
c       write(*,*)1d0/beta34*L34(p),'a'
c       write(*,*)L34(p),'a'

c      write(*,*)-myli2(-1d0*pt**2/mt**2),'a'
c      write(*,*)-myli2(1d0-4d0*dot(p,1,3)*dot(p,2,3)/mt**2/q2),'ch'

c      write(*,*)dlog(1d0+pt**2/mt**2),'a'
c      write(*,*)dlog(4d0*dot(p,1,3)*dot(p,2,3)/mt**2/q2),'ch'




c       write(*,*)(1d0+betat**2)/betat/2d0,'ch'
c       write(*,*)1d0/beta34,'a'


ch    Checks with the H0*S1 type delta terms of arXiv:1307.2464v2
      H0S1=H0S1+(-2d0*(1d0+betat**2)/betat*Tqq(3,4)*(
     .           -myli2(-1d0*xs*(1d0-costheta)/(1d0+costheta))
     .           +myli2(-1d0*(1d0-costheta)/(1d0+costheta)/xs)
     .           +4d0*dlog(xs)*dlog(dsqrt(0.5d0*(1d0+costheta))))
     .-4d0*myli2(1d0-4d0*dot(p,1,3)*dot(p,2,3)/mt**2/q2)
     . *(Tqq(1,3)+Tqq(2,3))
     .+4d0*dlog(4d0*dot(p,1,3)*dot(p,2,3)/mt**2/q2)*Tqq(3,3))*
     .  ((4d0*dot(p,1,3)**2+4d0*dot(p,2,3)**2)/q2**2
     .  +2d0*mt**2/q2)*cf

c      H0S1=H0S1*
c     .  ((4d0*dot(p,1,3)**2+4d0*dot(p,2,3)**2)/q2**2
c     .  +2d0*mt**2/q2)*cf
c      H0S1=H0S1*betat/4d0/W/q2*gsq**3/(4d0*pi)**3/3d0*fbGeV2


      IReal=IReal+(1d0/beta34*L34(p)*Tqq(3,4)
     .      +2d0*dlog(1d0+pt**2/mt**2)*Tqq(3,3))
     .      *msqc(j,k)/2d0

      do i1=1,2
      do j1=3,4
       IReal=IReal-(myli2(-1d0*pt**2/mt**2)*Tqq(i1,j1))*
     .msqc(j,k)/2d0
      enddo
      enddo

c        IReal=IReal*msqc(j,k)/2d0

c      write(*,*)H0S1,'ch',j,k
c      write(*,*)tH1stfs,'a'

ch These are the relevant normalization factors,
ch In this type our Ireal agrees with H0S1

ch Delta terms from the real and virtual subtraction operator

ch    qqbar channel

ch      if(j.gt.0.and.k.lt.0)then
ch      tH1stfs=tH1stfs+(
ch     .       +1d0/beta34*L34(p)*Tqq(3,4)
ch     .       -dlog(scale**2/(q2-2d0*mt**2))*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       +dlog(scale**2/q2)*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       +(dlog(scale**2/mt**2)
ch     .       +2d0*dlog(1d0+pt**2/mt**2)
ch     .         -2d0*LR)*Tqq(3,3))
c     .        *msqc(j,k)
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq
c      write(*,*)L34(q2,beta34,mtrans)

ch      do i1=1,2
ch      do j1=3,4
ch      tH1stfs=tH1stfs+(
ch     .       (+dlog(scale**2/(-2d0*dot(p,i1,j1)))*
ch     .       dlog(mt**2/(-2d0*dot(p,i1,j1)))
ch     .       +0.5d0*dlog(scale**2/(-2*dot(p,i1,j1)))**2
ch     .       -0.5d0*LR**2
ch     .       -myli2(-1d0*pt**2/mt**2)
ch     .       -dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))*LR
ch     .        )*Tqq(i1,j1) 
c     .*msqc(j,k))
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq) 
ch      enddo
ch      enddo  

ch    qbarq channel

ch      elseif(j.lt.0.and.k.gt.0)then 
ch      tH1stfs=tH1stfs+(
ch     .       +1d0/beta34*L34(pb)*Tqq(3,4)
ch     .       -dlog(scale**2/(q2-2d0*mt**2))*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       +dlog(scale**2/q2)*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       +(dlog(scale**2/mt**2)
ch     .       +2d0*dlog(1d0+pt**2/mt**2)
ch     .         -2d0*LR)*Tqq(3,3))
c     .        *msqc(j,k)
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq
c      write(*,*)L34(q2,beta34,mtrans)

ch      do i1=1,2
ch      do j1=3,4
ch      tH1stfs=tH1stfs+(
ch     .       (+dlog(scale**2/(-2d0*dot(pb,i1,j1)))*
ch     .       dlog(mt**2/(-2d0*dot(pb,i1,j1)))
ch     .       +0.5d0*dlog(scale**2/(-2*dot(pb,i1,j1)))**2
ch     .       -0.5d0*LR**2
ch     .       -myli2(-1d0*pt**2/mt**2)
ch     .       -dlog(4d0*dot(pb,i1,j1)**2/(q2*mt**2))*LR
ch     .        )*Tqq(i1,j1) 
c     .*msqc(j,k))
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq) 
ch      enddo
ch      enddo

ch    gg channel

ch      if(j.eq.0.and.k.eq.0)then

ch      tH1stfs=tH1stfs+(
ch     .       +1d0/beta34*L34(p)*Tgg(3,4)
ch     .       -dlog(scale**2/(q2-2d0*mt**2))*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tgg(3,4)
ch     .       +dlog(scale**2/q2)*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tgg(3,4)
ch     .       +(dlog(scale**2/mt**2)
ch     .       +2d0*dlog(1d0+pt**2/mt**2)
ch     .         -2d0*LR)*Tgg(3,3))
c     .        *msqc(j,k)
ch     .       *fx10(j)*fx20(k)*flgq
c      write(*,*)L34(q2,beta34,mtrans)

ch      do i1=1,2
ch      do j1=3,4
ch      tH1stfs=tH1stfs+(
ch     .       (+dlog(scale**2/(-2d0*dot(p,i1,j1)))*
ch     .       dlog(mt**2/(-2d0*dot(p,i1,j1)))
ch     .       +0.5d0*dlog(scale**2/(-2*dot(p,i1,j1)))**2
ch     .       -0.5d0*LR**2
ch     .       -myli2(-1d0*pt**2/mt**2)
ch     .       -dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))*LR
ch     .        )*Tgg(i1,j1) 
c     .*msqc(j,k))
ch     .       *fx10(j)*fx20(k)*flgq) 
ch      enddo
ch      enddo
ch      endif
       

ch   scheme dependent delta terms
c      tH1stfs=tH1stfs+(
c     . Tqq(3,3)+Tqq(4,4)
c     . +1d0/beta34*dlog((1d0+beta34)/(1d0-beta34))
c     . *Tqq(3,4)
c     . +dlog(4d0*dot(p,1,3)**2/(q2*mt**2))*Tqq(1,3)
c     . +dlog(4d0*dot(p,1,4)**2/(q2*mt**2))*Tqq(1,4)
c     . +dlog(4d0*dot(p,2,3)**2/(q2*mt**2))*Tqq(2,3)
c     . +dlog(4d0*dot(p,2,4)**2/(q2*mt**2))*Tqq(2,4)
c     . +0.5d0*A1q)
c     . *fx10(j)*fx20(k)*qqQQv_1*aveqq*flgq



c      write(*,*)dsqrt(q2),'a'

c      tH1stfs=tH1stfs+fx10(j)*fx20(k)*msqc(j,k)*flgq
c      if(j.eq.-k.and.j.ne.0)then
c      write(*,*)msqc(j,k),'a'
c      endif


c     This term is not in the Alessandro's notes. The point is that
c     in principle this term cancels agains the one coming from Ising(FS)
c     But since we included this term in the virtual, we must have it also here
c      tH1stfs=tH1stfs-1/beta34*dlog((1+beta34)/(1-beta34))*LR*Tqq(3,4)
c     .         *fx10(j)*fx20(k)*msqc(j,k)*flgq
c      tH1stfs=tH1stfs+(LR**2*Tqq(1,1)+LR*2*B1q)
c     .         *fx10(j)*fx20(k)*msqc(j,k)*flgq

c      tH1stfs=tH1stfs+(Tqq(3,3)+Tqq(4,4)
c     .  +1d0/beta34*dlog((1d0+beta34)/(1d0-beta34)))*
c     .  Tqq(3,4)*
c     .  fx10(j)*fx20(k)*msqc(j,k)*flgq

c      do i1=1,2
c      do j1=3,4
c      tH1stfs=tH1stfs+
c     .  +dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))
c     .        *Tqq(i1,j1) 
c     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq 
c      enddo
c      enddo




c      tH1stfs=0d0


c      tH1stfs=tH1stfs+
c     .        (dlog(scale**2/mt**2)+2*dlog(1+pt**2/mt**2))*Tqq(3,3)**2
c     .        *fx10(j)*fx20(k)*msqc(j,k)*flgq

c      write(*,*)(tH1st+tH1stfs/2d0)/gsq**2*4d0*3d0,'tttt'

      if(order.eq.1) goto 75

      if(j.eq.-k.and.j.ne.0)then

ch    Contribution from the qqb->QQb hard process      

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
ch     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq 
     &            fx20(k)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)*flgq       

      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq 
     &            fx10(j)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)*flgq        


CC    H2st qg contribution from C1*C1

C     regular*regular

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


C     regular-delta

ch     C1qqdelta=0 in the hard scheme!
ch        if(j.gt.0.and.k.lt.0)then
      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k) 
ch     &            fx20(k)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)
     &            fx20(k)*H1qdelta*msqc(j,k)       

      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)  
ch     &            fx10(j)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)
     &            fx10(j)*H1qdelta*msqc(j,k)

ch        else
ch      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k) 
ch     &            fx20(k)*H1qbdelta*msqc(j,k)

ch      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)  
ch     &            fx10(j)*H1qbdelta*msqc(j,k)
ch      endif

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

CCCC Terms needed for NNLO scale dependence  CCCCCC


CC    (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


      diffg1fq=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)
     &  - Pqqint(xx10)*fx10(j)


      diffg10q=-dlog(xx10)*fx1p(0)*Pqg(z1)

      diffg2fq=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)
     &  - Pqqint(xx20)*fx20(k)


      diffg20q=-dlog(xx20)*fx2p(0)*Pqg(z2)


      tgagaq=tgagaq+2*
     #   (flgq*diffg10q*diffg20q+flgq*diffg1fq*diffg2fq
     #   +diffg10q*diffg2fq+diffg1fq*diffg20q)*msqc(j,k)


CC     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

C     First leg

      
      diff1=-dlog(xx10)*(flgq*(fx1p(j)-fx10(j)*xx10**beta)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq+fx1p(0)*(Pqqqg(z1)+Pqggg(z1)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq


C    Second leg

      
      diff2=-dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq


C     Include Pqggq

      do l=1,nf
      diff1=diff1-dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
      diff2=diff2-dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
      enddo

      tgagaq=tgagaq+diff1*fx20(k)*msqc(j,k)
      tgagaq=tgagaq+diff2*fx10(j)*msqc(j,k)



C    End of (gamma+gamma)*(gamma+gamma) term

C    Start  (C+C)*(gamma+gamma) term

c    gamma first leg, C second leg


ch      diffc2f=-dlog(xx20)*fx2p(k)*Cqq(z2)+C1qqdelta*fx20(k)
ch      diffc2fq=-dlog(xx20)*fx2p(k)*Cqq(z2)+0.5d0*H1qdelta*fx20(k)
      diffc2fq=-dlog(xx20)*fx2p(k)*Cqq(z2)

      diffc20q=-dlog(xx20)*fx2p(0)*Cqg(z2)


      tcgaq=tcgaq+msqc(j,k)*
     # (flgq*diffg10q*diffc20q+flgq*diffg1fq*diffc2fq
     #          +diffg10q*diffc2fq+diffg1fq*diffc20q)


c    C first leg, gamma second leg

ch      diffc1f=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)
ch      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)+0.5d0*H1qdelta*fx10(j)
      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)

      diffc10q=-dlog(xx10)*fx1p(0)*Cqg(z1)

      tcgaq=tcgaq+msqc(j,k)*
     # (flgq*diffc10q*diffg20q+flgq*diffc1fq*diffg2fq
     #          +diffc10q*diffg2fq+diffc1fq*diffg20q)
    

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
     &     *(-dlog(xx10))*fx20(k)*msqc(j,k) 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
     &     *(-dlog(xx20))*fx10(j)*msqc(j,k) 

c    Add Cqg*Pgq contribution

      do l=1,nf
      tcgaq=tcgaq+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcgaq=tcgaq+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq 
      enddo

CC  Start 2-loop AP

C   Gluon + pure singlet


      do l=-nf,nf
      if(l.eq.0) then
      tgamma2q=tgamma2q+fx1p(0)*P2qg(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
      tgamma2q=tgamma2q+fx2p(0)*P2qg(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      else
      tgamma2q=tgamma2q+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part

      tgamma2q=tgamma2q+fx1p(j)*P2qqV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+fx2p(k)*P2qqV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq


C   P2qq non-singlet: 1/(1-z)_+


      diff=-dlog(xx10)
     &  *(fx1p(j)-fx10(j)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq


      diff=-dlog(xx20)
     &  *(fx2p(k)-fx20(k)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq

      

C   P2qqb non singlet

      tgamma2q=tgamma2q+fx1p(-j)*P2qqbV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

      tgamma2q=tgamma2q+fx2p(-k)*P2qqbV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq



CCCCCCCCCCCC   End of NNLO scale dependence CCCCCCCCCCCCCCCCC


ch  Contribution from the gg->QQb hard process

      elseif(j.eq.0.and.k.eq.0)then

c     H2st, gg channel: D0(z), first leg

      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*H2ggD0/(1-z1)
     &    +fx1p(0)*H2ggreg(z1))

      tH2st=tH2st+0.5d0*diff*fx20(0)*msqc(0,0)*flgq
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx10)*
     .             fx10(0)*fx20(0)*msqc(0,0)*flgq

c     H2st, gg channel: D0(z), second leg
      
      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*H2ggD0/(1-z2)
     &    +fx2p(0)*H2ggreg(z2))

      tH2st=tH2st+0.5d0*diff*fx10(0)*msqc(0,0)*flgq
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx20)*
     .               fx10(0)*fx20(0)*msqc(0,0)*flgq


CCCCCC  Terms required for NNLO Scale dependence



CC    (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg: gluon channel

      diff10g=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)

      diff10g=3*diff10g+beta0*fx10(0)
     &     -dlog(xx10)*fx1p(0)*Pggreg(z1)

      diff20g=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)

      diff20g=3*diff20g+beta0*fx20(0)
     &     -dlog(xx20)*fx2p(0)*Pggreg(z2)


C    Second: gamma*gamma: gluon channel

c    First leg

      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)
     &    *(D0gggg/(1-z1)+D1gggg*dlog(1-z1)/(1-z1))
     &      +fx1p(0)*(Pggggreg(z1)+Pgqqg(z1)))
     &    +(Deltagggg-D0gggg*D0int(xx10)-D1gggg*D1int(xx10))*fx10(0)


      tgagag=tgagag+diff*flgq*fx20(0)*msqc(0,0)


c    Second leg

      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)
     &    *(D0gggg/(1-z2)+D1gggg*dlog(1-z2)/(1-z2))
     &      +fx2p(0)*(Pggggreg(z2)+Pgqqg(z2)))
     &    +(Deltagggg-D0gggg*D0int(xx20)-D1gggg*D1int(xx20))*fx20(0)

      tgagag=tgagag+diff*flgq*fx10(0)*msqc(0,0)



CC    Start  (C+C)*(gamma+gamma) term: diagonal part


c    gamma first leg, C second

      diffg10g=-dlog(xx10)
     &  *((fx1p(0)-fx10(0)*xx10**beta)*3d0/(1-z1)+Pggreg(z1)*fx1p(0))
     &    +fx10(0)*(beta0-3*D0int(xx10))

ch      diffc20=C1ggdelta*fx20(0)
ch    remove the alphas running, included at the end !
ch      diffc20g=0.5d0*H1gdelta*fx20(0)
      diffc20g=0d0

c    gamma second leg, C first

      diffg20g=-dlog(xx20)
     &  *((fx2p(0)-fx20(0)*xx20**alfa)*3d0/(1-z2)+Pggreg(z2)*fx2p(0))
     &    +fx20(0)*(beta0-3*D0int(xx20))


ch      diffc10=C1ggdelta*fx10(0)

ch      diffc10g=0.5d0*H1gdelta*fx10(0)
      diffc10g=0d0

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)

      tcgag=tcgag+CgqPqg(z1)*(-dlog(xx10))*flgq*
     &            fx1p(0)*fx20(0)*msqc(0,0)

c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)

      tcgag=tcgag+CgqPqg(z2)*(-dlog(xx20))*flgq*
     &            fx2p(0)*fx10(0)*msqc(0,0)

c    End of (C+C)*(gamma+gamma)


CC    gamma2: diagonal part

c     First leg

      diff=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diff-dlog(xx10)*P2gg(z1)*fx1p(0))
     &                *flgq*fx20(0)*msqc(0,0)


c     Second leg

      diff=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diff-dlog(xx20)*P2gg(z2)*fx2p(0))
     &                *flgq*fx10(0)*msqc(0,0)

      endif


 75   continue

      enddo
      enddo
ch      write(*,*)tH1st,'tH1st1'
ch      tH1st=0d0
ch      tH1stF=tH1stFg
      do j=1,nf

ch    gg channel

ch    H1st: gg channel
     
C     H1st: Cgq, first leg

      tH1st=tH1st+(fx1p(j)+fx1p(-j))*Cgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqc(0,0)


C     H1st: Cgq, second leg

c      elseif(j.ne.0.and.k.eq.0)then

      
      tH1st=tH1st+(fx2p(j)+fx2p(-j))*Cgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqc(0,0)

C     H1st: muf dependence: Pgq, first leg


      tH1stF=tH1stF+(-dlog(xx10))
     & *(fx1p(j)+fx1p(-j))*Pgq(z1)*fx20(0)*msqc(0,0)
      tH1stFg=tH1stFg+(-dlog(xx10))
     & *(fx1p(j)+fx1p(-j))*Pgq(z1)*fx20(0)*msqc(0,0)


C     H1st: muf dependence: Pgq, second leg


      tH1stF=tH1stF+(-dlog(xx20))
     & *(fx2p(j)+fx2p(-j))*Pgq(z2)*fx10(0)*msqc(0,0)
      tH1stFg=tH1stFg+(-dlog(xx20))
     & *(fx2p(j)+fx2p(-j))*Pgq(z2)*fx10(0)*msqc(0,0)

C NEW: Add H2 qg contribution

c First leg


      tH2st=tH2st+(fx1p(j)+fx1p(-j))*H2stgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqc(0,0)

c Second leg


      tH2st=tH2st+(fx2p(j)+fx2p(-j))*H2stgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqc(0,0)


ch    In the following we subtract the C1*C1 regular-delta terms and add
ch    H1*C1 terms
ch     regular-delta

ch      tH2st=tH2st+fx1p(0)*Cgg(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq 
ch     &            fx20(0)*(H1gdelta-C1ggdelta)*msqc(j,k)*flgq

ch      tH2st=tH2st+fx2p(0)*Cgg(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq 
ch     &            fx10(0)*(H1gdelta-C1ggdelta)*msqc(j,k)*flgq


CC    H2st qg contribution from C1*C1

C     regular*regular

ch      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*


ch      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
ch     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


C     regular-delta

ch     C1qqdelta=0 in the hard scheme!
ch      C1ggdelta=0d0
      tH2st=tH2st+(fx1p(j)+fx1p(-j))*Cgq(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k) 
     &            fx20(0)*(H1gdelta
ch     &                      -C1ggdelta
ch     &                     -beta0*LR
     &                       )*msqc(0,0)

      tH2st=tH2st+(fx2p(j)+fx2p(-j))*Cgq(z2)*(-dlog(xx20))*
ch     &            fx10(0)*0.5d0*H1gdelta*msqc(0,0)  
     &            fx10(0)*(H1gdelta
ch     &                     -C1ggdelta
ch     &                    -beta0*LR
     &                      )*msqc(0,0)
ch      write(*,*)C1ggdelta,j,'a'


CCCCC Terms needed for NNLO scale dependence


CC    Now (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


      diff1fg=diff1fg-dlog(xx10)*Pgq(z1)*(fx1p(j)+fx1p(-j))

      diff2fg=diff2fg-dlog(xx20)*Pgq(z2)*(fx2p(j)+fx2p(-j))


C     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

      tgagag=tgagag-dlog(xx10)*(Pgqqq(z1)+Pgggq(z1))*(fx1p(j)+fx1p(-j))
     &            *fx20(0)*msqc(0,0)

      tgagag=tgagag-dlog(xx20)*(Pgqqq(z2)+Pgggq(z2))*(fx2p(j)+fx2p(-j))
     &            *fx10(0)*msqc(0,0)
c      write(*,*)Pgqqq(z1),'Pgqqqz1',Pgggq(z1),'Pgggq(z1)',
ch     .         tgagag,'tgagag'


C    End of (gamma+gamma)*(gamma+gamma) term


C    Start (C+C)*(gamma+gamma)


c    Gamma first leg

      diffg1fg=diffg1fg-dlog(xx10)*(fx1p(j)+fx1p(-j))*Pgq(z1)

c    C second leg

      diffc2fg=diffc2fg-dlog(xx20)*(fx2p(j)+fx2p(-j))*Cgq(z2)

c    Gamma second leg

      diffg2fg=diffg2fg-dlog(xx20)*(fx2p(j)+fx2p(-j))*Pgq(z2)

c    C first leg

      diffc1fg=diffc1fg-dlog(xx10)*(fx1p(j)+fx1p(-j))*Cgq(z1)

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)


      tcgag=tcgag+CgqPqq(z1)*
     .   (-dlog(xx10))*(fx1p(j)+fx1p(-j))*fx20(0)*msqc(0,0) 


c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)


      tcgag=tcgag+CgqPqq(z2)*
     .     (-dlog(xx20))*(fx2p(j)+fx2p(-j))*fx10(0)*msqc(0,0) 


CC    gamma2: qg channel


c    First leg

      tgamma2g=tgamma2g
     &       -dlog(xx10)*P2gq(z1)*(fx1p(j)+fx1p(-j))*fx20(0)*msqc(0,0) 

c    Second leg

      tgamma2g=tgamma2g
     &       -dlog(xx20)*P2gq(z2)*(fx2p(j)+fx2p(-j))*fx10(0)*msqc(0,0) 



CC   Effect of spin correlations

c     First leg

      spgq1=spgq1-dlog(xx10)*(fx1p(j)+fx1p(-j))*Ggq(z1)

c    Second leg

      spgq2=spgq2-dlog(xx20)*(fx2p(j)+fx2p(-j))*Ggq(z2)


      enddo
ch      write(*,*)spgq1,spgq2,'a'


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    Check it !      

      tgagag=tgagag+2*msqc(0,0)
     # *(flgq*diff10g*diff20g+flgq*diff1fg*diff2fg
     #   +diff10g*diff2fg+diff1fg*diff20g) 
ch      write(*,*)tgagag,'tgagag2',
ch     .      diff1fg,'diff1fg',diff20g,'diff20g'
c

c    gamma first leg, C second leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg10g*diffc20g+flgq*diffg1fg*diffc2fg
     #          +diffg10g*diffc2fg+diffg1fg*diffc20g)

c    gamma second leg, C first leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg20g*diffc10g+flgq*diffg2fg*diffc1fg
     #          +diffg20g*diffc1fg+diffg2fg*diffc10g)



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    qq contribution to H2 (C1C1)


      tH2st=tH2st+msqc(0,0)*diffc1fg*diffc2fg*flgq


chc    spin correlations

ch   spin + anglular correlations

c    gg channel

      
      tH2sp=msqc(0,0)*dlog(xx10)*dlog(xx20)*
     . fx1p(0)*fx2p(0)*Ggg(z1)*Ggg(z2)*flgq

c    gq+qg channel

ch      tH2sp=tH2sp-msqc(0,0)*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
ch      tH2sp=tH2sp-msqc(0,0)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1
      tH2sp=tH2sp-msqGGav(p)*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
      tH2sp=tH2sp-msqGGav(p)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1

ch    add also the angular correlations due to D*G interference

      tH2sp=tH2sp+msqDGav(p)*fx10(0)*spgq2
      tH2sp=tH2sp+msqDGav(p)*fx20(0)*spgq1
ch      tH2sp=tH2sp+msqDGav(p)*fx10(0)*spgq2
ch      tH2sp=tH2sp+msqDGav(p)*fx20(0)*spgq1
ch      write(*,*)msqDGav(p),msqGGav(p),msqc(0,0),'a'
ch      write(*,*)msqGGav(p)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1,'b'
ch      write(*,*)msqDGav(p)*fx2p(0)*spgq1,'c'
ch      write(*,*)msqDGav(p)*fx20(0)*spgq1,'d'
ch      write(*,*)msqc(0,0)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1,'d'

c    qq channel


ch      tH2sp=tH2sp+msqc(0,0)*spgq1*spgq2*flgq
      tH2sp=tH2sp+msqGGav(p)*spgq1*spgq2*flgq


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 72   xmsq=tdelta
ch      write(*,*)tH1st,'tH1st2',tdelta,'tdelta'
ch      write(*,*)tdelta,'b'
ch      tH1stfs=tH1stfs/2d0
c      tH1stF=0d0
c      tH1stfs=0d0
ch      xmsq=0d0
c      tH1st=0d0
c      write(*,*)asopi*tH1st,xmsq,asopi*LF*tH1stF
      if(order.eq.1)then
       xmsq=xmsq
     .          +asopi*
     .          (
     .           tH1st+
     .           LF*tH1stF
ch     .         -2*beta0*LR*tdelta
     .                )
ch     .         +tH1stfs)
c       xmsq=asopi*IReal/2d0
c       xmsq=asopi*(tH1st+LF*tH1stF+tH1stfs)
      elseif(order.eq.2)then
ch       xmsq=xmsq
ch     &          +asopi*(tH1st+LF*tH1stF)
ch     &          +asopi**2*(tdeltaqq*H2qqdelta
ch     &                     +tdeltagg*H2ggdelta
ch     &                     +tH2st+tH2sp
ch     &          )

CC     add scale dependence at NNLO

ch    Checked!
ch       xmsq=xmsq+asopi**2*(
ch     &            0.5d0*beta0*LF**2*tH1stF
ch     &            +tgamma2q*LF+tgamma2g*LF
ch     &            -2*beta0*LR*(tH1st+LF*tH1stF-2*beta0*LR*tdelta)
ch     &            -beta0*LR*LF*tH1stF
ch     &            +LF*tcgaq+0.5d0*LF**2*tgagaq
ch     &            +LF*tcgag
ch     &            +0.5d0*LF**2*tgagag
ch     &            -beta0*LR*tH1st
ch     &            -2*(0.5d0*(beta0*LR)**2+beta1*LR)*tdelta
ch     &             )



c     Include missing delta term from C*gamma (no factor 2 here !)

ch      xmsq=xmsq+asopi**2*(LF*C1ggdelta*tH1stF)
       xmsq=xmsq
ch     &          +asopi*(tH1st+LF*tH1stF)
     &         +asopi**2*(tdeltaqq*H2qqdelta
     &                     +tdeltagg*H2ggdelta
     &                     +tH2st+tH2sp)

CC     add scale dependence at NNLO

       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF
     &                   +tgamma2q*LF+tgamma2g*LF
     &                   -beta0*LR*(tH1st+LF*tH1stF)
     &                   +LF*tcgaq+0.5d0*LF**2*tgagaq
     &                   +LF*tcgag+0.5d0*LF**2*tgagag)


      xmsq=xmsq+asopi**2*(LF*H1gdelta*tH1stFg)
      xmsq=xmsq+asopi**2*(LF*H1qdelta*tH1stFq)


ch      xmsq=xmsq+asopi**2*(LF*0.5d0*H1qdelta*tH1stFq)

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2gg*tdeltagg)*LF



      endif    
ch      write(*,*)asopi,'asopi',xx10,'xx10',
ch     .          xx20,'xx20',alfa,'alfa',beta,'beta',
ch     .          facscale,scale,'muf,mur',
ch     ,          msqc(0,0),msqc(1,-1),'msqc',
ch     .          LF,'LF',LR,'LR',xmsq,'xmsq',
ch     .          z1,z2,fx10(0),fx20(0),fx10(0),fx2p(0),
ch     .           H1qdelta,H1gdelta,'H1'
ch      write(*,*)xmsq,'a'
c      write(*,*)xmsq/asopi/gsq**2*4d0*3d0/10d0,'a'
ch      write(*,*)H0S1/10d0,'ch'
c      write(*,*)tH1st/asopi/gsq**2*4d0*3d0/10d0,'aa'
c      write(*,*)IReal/gsq**2*4d0*3d0/10d0,'aaa'
c      write(*,*)
c     .         PoleEp1/gsq**2*4d0
c     .         +dlog(scale**2/mt**2)*PoleEp2/gsq**2*4d0,
c     .         PoleEp1/gsq**2*4d0,
c     .            2d0*dot(p,1,2)/1000000d0
c     .           ,2d0*dot(p,1,3)/1000000d0,scale


ch      lowintHst=flux*pswt*xmsq/BrnRat


ch      Checked against SCET result!

ch      write(*,*)asopi*(H1qdelta*msqc(1,-1)
ch     .             +2*beta0*LR*msqc(1,-1)),'1loopqq'
ch      write(*,*)asopi*(H1gdelta*msqc(0,0)
ch     .              +2*beta0*LR*msqc(0,0)),'1loopgg'
ch      write(*,*)dsqrt(2*dot(p,1,2)),'2'
ch      write(*,*)2*dot(p,1,3),'3'
ch      write(*,*)2*dot(p,2,3),'4'
ch      write(*,*)as,'6'
ch      write(*,*)scale,'7' 
ch      write(*,*)aaa,'8'
      lowintHst=flux*xmsq*pswt/BrnRat
ch      write(*,*)xmsq,lowintHst,'c'

c      write(*,*)lowintHst

c      write(*,*)pswt,'ps',betat,'betat',q2


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
C
C     NEW version August 2013
C

      double precision function C2qg(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3
      external myli2,myli3


      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0


      C2qg=(688*CA-1260*CA*z-702*CF*z+1548*CA*z**2+2322*CF*z**2+ 
     -    72*CA*Pi**2*z**2+144*CF*Pi**2*z**2-1192*CA*z**3-2160*CF*z**3- 
     -    144*CF*Pi**2*z**3+324*CA*z**2*Log(1-z)-324*CF*z**2*Log(1-z)- 
     -    432*CA*z**3*Log(1-z)+432*CF*z**3*Log(1-z)- 
     -    216*CA*z**2*Log(1-z)**2+216*CF*z**2*Log(1-z)**2+ 
     -    216*CA*z**3*Log(1-z)**2-216*CF*z**3*Log(1-z)**2 + 
     -    36*CA*z*Log(1 - z)**3 - 36*CF*z*Log(1 - z)**3 - 
     -    72*CA*z**2*Log(1-z)**3 + 72*CF*z**2*Log(1-z)**3 + 
     -    72*CA*z**3*Log(1-z)**3-72*CF*z**3*Log(1-z)**3 + 
     -    504*CA*z*Log(z)+432*CF*z*Log(z)-720*CA*z**2*Log(z) + 
     -    810*CF*z**2*Log(z)+144*CA*Pi**2*z**2*Log(z)
     -    +1632*CA*z**3*Log(z)- 
     -    432*CF*z**3*Log(z) - 432*CF*z**2*Log(1-z)*Log(z) + 
     -    432*CF*z**3*Log(1-z)*Log(z)+108*CF*z*Log(1-z)**2*Log(z)- 
     -    216*CF*z**2*Log(1-z)**2*Log(z)+
     -    216*CF*z**3*Log(1-z)**2*Log(z)- 
     -    54*CA*z*Log(z)**2+27*CF*z*Log(z)**2+216*CA*z**2*Log(z)**2+ 
     -    324*CF*z**2*Log(z)**2 - 792*CA*z**3*Log(z)**2 - 
     -    216*CF*z**3*Log(z)**2 + 108*CF*z*Log(1 - z)*Log(z)**2- 
     -    864*CA*z**2*Log(1-z)*Log(z)**2-
     -    216*CF*z**2*Log(1-z)*Log(z)**2+ 
     -    216*CF*z**3*Log(1-z)*Log(z)**2+36*CA*z*Log(z)**3 - 
     -    18*CF*z*Log(z)**3+72*CA*z**2*Log(z)**3+36*CF*z**2*Log(z)**3- 
     -    72*CF*z**3*Log(z)**3 + 36*CA*Pi**2*z*Log(1 + z) + 
     -    72*CA*Pi**2*z**2*Log(1 + z) + 72*CA*Pi**2*z**3*Log(1+z)+ 
     -    432*CA*z**2*Log(z)*Log(1 + z) + 432*CA*z**3*Log(z)*Log(1+z)+ 
     -    108*CA*z*Log(z)**2*Log(1+z)+216*CA*z**2*Log(z)**2*Log(1+z)+ 
     -    216*CA*z**3*Log(z)**2*Log(1+z)-72*CA*z*Log(1+z)**3 - 
     -    144*CA*z**2*Log(1 + z)**3 - 144*CA*z**3*Log(1 + z)**3 - 
     -    72*(3*(CA - CF)*z*(1 - 2*z + 2*z**2)*Log(1 - z) + 
     -       2*CA*(2 - 3*z + 12*z**2 - 11*z**3 + 6*z**2*Log(z)))*
     -     myLi2(1-z) - 216*CA*z*
     -     (-2*z*(1 + z) + (1 + 2*z + 2*z**2)*Log(z))*myli2(-z)+ 
     -    216*CF*z*Log(z)*myli2(z)-
     -    1728*CA*z**2*Log(z)*myli2(z)- 
     -    432*CF*z**2*Log(z)*myli2(z)+432*CF*z**3*Log(z)*myli2(z)+ 
     -    216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z)- 
     -    432*CA*z**2*myli3(1-z)+432*CF*z**2*myli3(1-z) + 
     -    432*CA*z**3*myli3(1-z) - 432*CF*z**3*myli3(1-z) + 
     -    648*CA*z*myli3(-z) + 1296*CA*z**2*myli3(-z)+ 
     -    1296*CA*z**3*myli3(-z) - 216*CF*z*myli3(z) + 
     -    1728*CA*z**2*myli3(z) + 432*CF*z**2*myli3(z)- 
     -    432*CF*z**3*myli3(z) + 432*CA*z*myli3(1/(1+z))+ 
     -    864*CA*z**2*myli3(1/(1+z)) +
     -    864*CA*z**3*myli3(1/(1+z)) - 
     -    648*CA*z*Z3 + 1728*CF*z*Z3-3456*CF*z**2*Z3 - 
     -    1296*CA*z**3*Z3 + 3456*CF*z**3*Z3)/(1728d0*z)
ch     -    -CF/4d0*(z*Log(z)+0.5d0*(1-z**2)+(Pi**2-8)*z*(1-z))
     -    -1d0/4d0*z*(1-z)*CF*(PI**2/2d0-4d0)


      return
      end

      function H2stgq(z)
      implicit none
      real *8 H2stgq,H2stgq1,H2stgq2,H2stgq3,H2stgqfull
      real *8 Pi,Z3,myli2,myli3,z,CF,CA
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3d0
      CA=3d0

      nf=5


c    Full result: Check it !

            H2stgqfull=
     #      -(12072-224*nf-396*Pi**2-12444*z+224*nf*z+432*Pi**2*z+
     #      2064*z**2-52*nf*z**2-432*Pi**2*z**2-1824*z**3+
     #      72*Pi**2*z**3+1584*dlog(1-z)-120*nf*dlog(1-z)-
     #      162*Pi**2*dlog(1-z)-1584*z*dlog(1-z)+120*nf*z*dlog(1-z)-
     #      324*Pi**2*z*dlog(1-z)+414*z**2*dlog(1-z) -
     #      24*nf*z**2*dlog(1-z)-162*Pi**2*z**2*dlog(1-z) +
     #      378*dlog(1-z)**2-36*nf*dlog(1-z)**2-378*z*dlog(1-z)**2 +
     #      36*nf*z*dlog(1-z)**2+99*z**2*dlog(1-z)**2 -
     #      18*nf*z**2*dlog(1-z)**2+60*dlog(1-z)**3 -
     #      60*z*dlog(1-z)**3+30*z**2*dlog(1-z)**3+1296*dlog(z) +
     #      6318*z*dlog(z)-288*z**2*dlog(z)+1584*z**3*dlog(z) +
     #      2376*dlog(1-z)*dlog(z)-2592*z*dlog(1-z)*dlog(z) +
     #      972*z**2*dlog(1-z)*dlog(z)-432*z**3*dlog(1-z)*dlog(z) +
     #      324*dlog(1-z)**2*dlog(z)+1620*z*dlog(1-z)**2*dlog(z) +
     #      486*z**2*dlog(1-z)**2*dlog(z)-900*z*dlog(z)**2 -
     #      189*z**2*dlog(z)**2-216*z**3*dlog(z)**2 -
     #      324*dlog(1-z)*dlog(z)**2+324*z*dlog(1-z)*dlog(z)**2 -
     #      162*z**2*dlog(1-z)*dlog(z)**2+84*z*dlog(z)**3 +
     #      66*z**2*dlog(z)**3+270*Pi**2*dlog(1+z)+
     #      432*Pi**2*z*dlog(1+z)+216*Pi**2*z**2*dlog(1+z)-
     #      324*z**2*dlog(z)*dlog(1+z)-1296*dlog(1-z)*dlog(z)*dlog(1+z)-
     #      2592*z*dlog(1-z)*dlog(z)*dlog(1+z) -
     #      1296*z**2*dlog(1-z)*dlog(z)*dlog(1+z) +
     #      324*dlog(z)**2*dlog(1+z)+324*z*dlog(z)**2*dlog(1+z) +
     #      162*z**2*dlog(z)**2*dlog(1+z)+648*dlog(z)*dlog(1+z)**2+
     #      1296*z*dlog(z)*dlog(1+z)**2+648*z**2*dlog(z)*dlog(1+z)**2 -
     #      216*dlog(1+z)**3-216*z*dlog(1+z)**3 -
     #      108*z**2*dlog(1+z)**3 -
     #      324*(z**2+2*(1+z)**2*dlog(1-z)+(2+2*z+z**2)*dlog(z) -
     #      2*dlog(1+z)-4*z*dlog(1+z)-2*z**2*dlog(1+z))*myli2(-z)-108*
     #       (-22+24*z-6*z**2+4*z**3-6*(1+z)**2*dlog(1-z) +
     #         3*(6-2*z+3*z**2)*dlog(z)+6*dlog(1+z)+12*z*dlog(1+z)+
     #         6*z**2*dlog(1+z))*myli2(z) +
     #      648*dlog(1-z)*myli2((1-z)/(1+z)) +
     #      1296*z*dlog(1-z)*myli2((1-z)/(1+z)) +
     #      648*z**2*dlog(1-z)*myli2((1-z)/(1+z)) -
     #      648*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      1296*z*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      648*z**2*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      648*dlog(1-z)*myli2((-1+z)/(1+z)) -
     #      1296*z*dlog(1-z)*myli2((-1+z)/(1+z)) -
     #      648*z**2*dlog(1-z)*myli2((-1+z)/(1+z)) +
     #      648*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      1296*z*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      648*z**2*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      1944*myli3(-z)+1944*z*myli3(-z) +
     #      972*z**2*myli3(-z)+3240*myli3(z)-648*z*myli3(z) +
     #      1620*z**2*myli3(z)+1296*myli3(1/(1+z)) +
     #      1296*z*myli3(1/(1+z))+648*z**2*myli3(1/(1+z)) -
     #      5184*Z3+3240*z*Z3-2592*z**2*Z3)/(324d0*z)
ch     #      +CF**2*3d0*z/4d0+CF*CA/z*((1+z)*Log(z)+2*(1-z)-
ch     #      (5+Pi**2)/4d0*z**2)
     #      -1d0/2d0*z*CF*(11d0+3d0*Pi**2)/2d0


C     Subtract Ggg*Ggq term
    
      H2stgq=H2stgqfull-(-4)*(2*(1-z)+(1+z)*dlog(z))/4d0

      return
      end


C     Regular part of H2st in gg channel

      function H2ggREG(z)
      implicit none
      real *8 H2ggREG
      real *8 Pi,Z3,myli2,myli3,z
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      nf=5

c      H2ggREG=(16164-226*nf+396*Pi**2-17790*z-52*nf*z+1512*Pi**2*z- 
c     -  5064*z**2 + 390*nf*z**2 + 684*Pi**2*z**2 + 5784*z**3 + 
c     -  314*nf*z**3 - 972*Pi**2*z**3 - 11100*z**4 - 164*nf*z**4 - 
c     -  1080*Pi**2*z**4 + 12006*z**5 - 262*nf*z**5 - 540*Pi**2*z**5- 
c     -  9648*dlog(1-z) + 864*Pi**2*dlog(1-z) + 7992*z*dlog(1-z) - 
c     -  1728*Pi**2*z*dlog(1-z) + 7434*z**2*dlog(1-z) + 
c     -  18*nf*z**2*dlog(1-z) - 4176*z**3*dlog(1-z) + 
c     -  864*Pi**2*z**3*dlog(1-z) + 2214*z**4*dlog(1-z) - 
c     -  18*nf*z**4*dlog(1-z) - 864*Pi**2*z**4*dlog(1-z) - 
c     -  3816*z**5*dlog(1-z) + 864*Pi**2*z**5*dlog(1-z) - 
c     -  432*dlog(1-z)**3-2592*z*dlog(1-z)**3-864*z**2*dlog(1-z)**3+ 
c     -  2160*z**3*dlog(1-z)**3 + 1296*z**4*dlog(1-z)**3 + 
c     -  432*z**5*dlog(1-z)**3 + 3060*dlog(z) - 216*Pi**2*dlog(z)+ 
c     -  693*z*dlog(z) - 222*nf*z*dlog(z) + 432*Pi**2*z*dlog(z) - 
c     -  4203*z**2*dlog(z)-204*nf*z**2*dlog(z)+2727*z**3*dlog(z) + 
c     -  222*nf*z**3*dlog(z)-216*Pi**2*z**3*dlog(z)+1143*z**4*dlog(z)+ 
c     -  204*nf*z**4*dlog(z)+216*Pi**2*z**4*dlog(z)-3420*z**5*dlog(z)- 
c     -  216*Pi**2*z**5*dlog(z) + 2376*dlog(1-z)*dlog(z)- 
c     -  2592*z*dlog(1-z)*dlog(z) + 216*z**2*dlog(1-z)*dlog(z)+ 
c     -  216*z**3*dlog(1-z)*dlog(z) - 2592*z**4*dlog(1-z)*dlog(z)+ 
c     -  2376*z**5*dlog(1-z)*dlog(z) - 1620*dlog(1-z)**2*dlog(z) - 
c     -  7452*z*dlog(1-z)**2*dlog(z) - 2916*z**2*dlog(1-z)**2*dlog(z)+ 
c     -  6156*z**3*dlog(1-z)**2*dlog(z)+4212*z**4*dlog(1-z)**2*dlog(z)+ 
c     -  972*z**5*dlog(1-z)**2*dlog(z)
c     -  -1188*dlog(z)**2-2295*z*dlog(z)**2- 
c     -  54*nf*z*dlog(z)**2 - 783*z**2*dlog(z)**2-30*nf*z**2*dlog(z)**2+ 
c     -  891*z**3*dlog(z)**2+54*nf*z**3*dlog(z)**2+1971*z**4*dlog(z)**2+ 
c     -  30*nf*z**4*dlog(z)**2+1404*z**5*dlog(z)**2 + 
c     -  1620*dlog(1-z)*dlog(z)**2 + 6804*z*dlog(1-z)*dlog(z)**2 + 
c     -  2268*z**2*dlog(1-z)*dlog(z)**2-6156*z**3*dlog(1-z)*dlog(z)**2- 
c     -  4212*z**4*dlog(1-z)*dlog(z)**2 - 972*z**5*dlog(1-z)*dlog(z)**2- 
c     -  432*dlog(z)**3 - 756*z*dlog(z)**3 - 8*nf*z*dlog(z)**3 - 
c     -  216*z**2*dlog(z)**3 - 8*nf*z**2*dlog(z)**3+756*z**3*dlog(z)**3+ 
c     -  8*nf*z**3*dlog(z)**3+648*z**4*dlog(z)**3+8*nf*z**4*dlog(z)**3+ 
c     -  108*z**5*dlog(z)**3 - 108*Pi**2*dlog(1+z) - 
c     -  108*Pi**2*z*dlog(1+z) - 108*Pi**2*z**2*dlog(1+z) + 
c     -  108*Pi**2*z**3*dlog(1+z) + 108*Pi**2*z**4*dlog(1+z) + 
c     -  108*Pi**2*z**5*dlog(1+z) + 324*dlog(z)**2*dlog(1+z) + 
c     -  324*z*dlog(z)**2*dlog(1+z) + 324*z**2*dlog(z)**2*dlog(1+z) - 
c     -  324*z**3*dlog(z)**2*dlog(1+z) - 324*z**4*dlog(z)**2*dlog(1+z) - 
c     -  324*z**5*dlog(z)**2*dlog(1+z) - 648*dlog(z)*dlog(1+z)**2 - 
c     -  648*z*dlog(z)*dlog(1+z)**2 - 648*z**2*dlog(z)*dlog(1+z)**2 + 
c     -  648*z**3*dlog(z)*dlog(1+z)**2 + 648*z**4*dlog(z)*dlog(1+z)**2 + 
c     -  648*z**5*dlog(z)*dlog(1+z)**2 + 216*dlog(1+z)**3 + 
c     -  216*z*dlog(1+z)**3 + 216*z**2*dlog(1+z)**3 - 
c     -  216*z**3*dlog(1+z)**3 - 216*z**4*dlog(1+z)**3 - 
c     -  216*z**5*dlog(1+z)**3 + 
c     -  648*(-1 + z)*(1 + z + z**2)**2*dlog(z)*myli2(-z) + 
c     -  216*(-11 - 42*z - 19*z**2 + 27*z**3 + 30*z**4 + 15*z**5 + 
c     -  12*(-1 - 6*z - 2*z**2 + 5*z**3 + 3*z**4 + z**5)*dlog(1-z) - 
c     -  3*(-1 - 9*z - z**2 + 9*z**3 + 5*z**4 + z**5)*dlog(z))*myli2(z)
c     -  + 648*myli3(-z) + 648*z*myli3(-z) + 
c     -  648*z**2*myli3(-z) - 648*z**3*myli3(-z) - 
c     -  648*z**4*myli3(-z) - 648*z**5*myli3(-z) + 
c     -  1944*myli3(z) + 12312*z*myli3(z) + 5832*z**2*myli3(z) - 
c     -  8424*z**3*myli3(z) - 4536*z**4*myli3(z) - 
c     -  3240*z**5*myli3(z) + 2592*myli3(z/(-1 + z)) + 
c     -  15552*z*myli3(z/(-1 + z)) + 5184*z**2*myli3(z/(-1 + z)) - 
c     -  12960*z**3*myli3(z/(-1 + z)) - 
c     -  7776*z**4*myli3(z/(-1 + z)) - 2592*z**5*myli3(z/(-1 + z))- 
c     -  1296*myli3(z/(1 + z)) - 1296*z*myli3(z/(1 + z)) - 
c     -  1296*z**2*myli3(z/(1 + z)) + 1296*z**3*myli3(z/(1 + z)) + 
c     -  1296*z**4*myli3(z/(1 + z)) + 1296*z**5*myli3(z/(1 + z)) - 
c     -  7776*Z3 + 4212*z*Z3 - 4212*z**2*Z3 - 
c     -  648*z**3*Z3 + 10368*z**4*Z3 - 5832*z**5*Z3)/
c     -  (72d0*z*(-1 + z**2))



        H2ggREG=-(-10776+226*nf+396*Pi**2+12408*z+52*nf*z-432*Pi**2*z+ 
     -   1656*z**2 - 390*nf*z**2 + 36*Pi**2*z**2 - 2388*z**3 - 
     -   314*nf*z**3 + 36*Pi**2*z**3 + 9120*z**4 + 164*nf*z**4 - 
     -   432*Pi**2*z**4 - 10020*z**5 + 262*nf*z**5 + 396*Pi**2*z**5 + 
     -   54*z**2*dlog(1-z) - 18*nf*z**2*dlog(1-z) - 54*z**4*dlog(1-z)+ 
     -   18*nf*z**4*dlog(1-z)-648*dlog(z)-6957*z*dlog(z)+ 
     -   222*nf*z*dlog(z)-693*z**2*dlog(z)+204*nf*z**2*dlog(z)+ 
     -   2133*z**3*dlog(z)-222*nf*z**3*dlog(z)+1341*z**4*dlog(z) - 
     -   204*nf*z**4*dlog(z)+4824*z**5*dlog(z)-2376*dlog(1-z)*dlog(z)+ 
     -   2592*z*dlog(1-z)*dlog(z)-216*z**2*dlog(1-z)*dlog(z)- 
     -   216*z**3*dlog(1-z)*dlog(z)+2592*z**4*dlog(1-z)*dlog(z)- 
     -   2376*z**5*dlog(1-z)*dlog(z) + 324*dlog(1-z)**2*dlog(z)- 
     -   324*z*dlog(1-z)**2*dlog(z)+324*z**2*dlog(1-z)**2*dlog(z) + 
     -   324*z**3*dlog(1-z)**2*dlog(z)-324*z**4*dlog(1-z)**2*dlog(z)+ 
     -   324*z**5*dlog(1-z)**2*dlog(z) + 675*z*dlog(z)**2+ 
     -   54*nf*z*dlog(z)**2-297*z**2*dlog(z)**2+30*nf*z**2*dlog(z)**2+ 
     -   513*z**3*dlog(z)**2-54*nf*z**3*dlog(z)**2+297*z**4*dlog(z)**2- 
     -   30*nf*z**4*dlog(z)**2 - 1188*z**5*dlog(z)**2 + 
     -   324*dlog(1-z)*dlog(z)**2 - 324*z*dlog(1-z)*dlog(z)**2 + 
     -   324*z**2*dlog(1-z)*dlog(z)**2 + 324*z**3*dlog(1-z)*dlog(z)**2- 
     -   324*z**4*dlog(1-z)*dlog(z)**2 + 324*z**5*dlog(1-z)*dlog(z)**2- 
     -   108*z*dlog(z)**3 + 8*nf*z*dlog(z)**3-216*z**2*dlog(z)**3 + 
     -   8*nf*z**2*dlog(z)**3+108*z**3*dlog(z)**3-8*nf*z**3*dlog(z)**3+ 
     -   216*z**4*dlog(z)**3-8*nf*z**4*dlog(z)**3-108*z**5*dlog(z)**3+ 
     -   108*Pi**2*dlog(1+z) + 108*Pi**2*z*dlog(1+z) + 
     -   108*Pi**2*z**2*dlog(1+z) - 108*Pi**2*z**3*dlog(1+z) - 
     -   108*Pi**2*z**4*dlog(1+z) - 108*Pi**2*z**5*dlog(1+z) - 
     -   324*dlog(z)**2*dlog(1+z) - 324*z*dlog(z)**2*dlog(1+z) - 
     -   324*z**2*dlog(z)**2*dlog(1+z) + 324*z**3*dlog(z)**2*dlog(1+z)+ 
     -   324*z**4*dlog(z)**2*dlog(1+z) + 324*z**5*dlog(z)**2*dlog(1+z)+ 
     -   648*dlog(z)*dlog(1+z)**2 + 648*z*dlog(z)*dlog(1+z)**2 + 
     -   648*z**2*dlog(z)*dlog(1+z)**2 - 648*z**3*dlog(z)*dlog(1+z)**2- 
     -   648*z**4*dlog(z)*dlog(1+z)**2 - 648*z**5*dlog(z)*dlog(1+z)**2- 
     -   216*dlog(1+z)**3 - 216*z*dlog(1+z)**3 - 216*z**2*dlog(1+z)**3+ 
     -   216*z**3*dlog(1+z)**3 + 216*z**4*dlog(1+z)**3 + 
     -   216*z**5*dlog(1+z)**3 - 
     -   648*(-1 + z)*(1 + z + z**2)**2*dlog(z)*myli2(-z) + 
     -   216*(-((-1 + z)**2*(11 + 10*z + 10*z**2 + 11*z**3)) + 
     -   3*(3 - z + 3*z**2 + z**3 - 3*z**4 + z**5)*dlog(z))*myli2(z)
     -    - 648*myli3(-z) - 648*z*myli3(-z) - 
     -   648*z**2*myli3(-z) + 648*z**3*myli3(-z) + 
     -   648*z**4*myli3(-z) + 648*z**5*myli3(-z) - 
     -   3240*myli3(z) + 648*z*myli3(z) - 3240*z**2*myli3(z) - 
     -   648*z**3*myli3(z) + 3240*z**4*myli3(z) - 
     -   648*z**5*myli3(z) + 1296*myli3(z/(1 + z)) + 
     -   1296*z*myli3(z/(1 + z)) + 1296*z**2*myli3(z/(1 + z)) - 
     -   1296*z**3*myli3(z/(1 + z)) - 1296*z**4*myli3(z/(1 + z)) - 
     -   1296*z**5*myli3(z/(1 + z)) + 3888*Z3 - 6804*z*Z3 + 
     -   1620*z**2*Z3 + 4536*z**3*Z3 - 3888*z**4*Z3 + 
     -   4536*z**5*Z3)/(72d0*z*(-1 + z**2))

C     Subtract Ggg*Ggg

      H2ggREG=H2ggREG-(-9)*(2*(1-z)+(1+z)*dlog(z))/4d0

      return
      end

CC Spin correlations

      function Ggq(z)
      real *8 z,Ggq,CF

      CF=4d0/3
      
      Ggq=CF*(1-z)/z
    

      return
      end

      function Ggg(z)
      real *8 z,Ggg,CA

      CA=3d0
      
      Ggg=CA*(1-z)/z
    

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
      double precision p(mxpart,4),pt
      double precision myli2,dot
      double precision r34(1:2)

      beta34t=dsqrt(1d0-(2*mtrans**2/(q2-2*mtrans**2))**2)
      c=dsqrt((1d0-beta34)/(1d0+beta34))
      ct=dsqrt((1d0-beta34t)/(1d0+beta34t))
      do i=1,2
       r34(i)=dot(p,i,3)/dot(p,i,4)
      enddo
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)
      beta34=dsqrt(1d0-mt**4/(dot(p,3,4)**2))
    
c      L34=dlog((1+beta34)/(1-beta34))*dlog(mtrans**2/mt**2)-
c     .    2*myli2(2*beta34/(1+beta34))-1d0/4d0*
c     .    dlog((1+beta34)/(1-beta34))**2+
c     .    2*(myli2(1-c*ct)+myli2(1-c/ct)+1/2*dlog(ct)**2)

      L34=0.5d0*dlog((1d0+beta34)/(1d0-beta34))*dlog(mtrans**4/mt**4)
     .    -2d0*myli2(2d0*beta34/(1d0+beta34))
     .    -0.25d0*dlog((1d0+beta34)/(1d0-beta34))**2
      do i=1,2
       L34=L34+myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*r34(i))
     .        +myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))/r34(i))
     .        +0.5d0*dlog(r34(i))**2
      enddo


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
c      include 'breit.f'
      include 'rescoeff.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
      common/colour/Tqq
      double precision Tqq(1:4,1:4)



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtqq

c      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
c      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )


*
        sn=2d0*dot(p,1,2)
        tn=2d0*dot(p,1,3)+mt**2
	un = 2d0*mt**2 - sn - tn
*
        x  = 4d0*mt**2/(Dsqrt(sn-4d0*mt**2)+Dsqrt(sn))**2 
        y  = - tn/mt**2
        z  = - un/mt**2
*
        upx = 1d0+x
	umx = 1d0-x
	upxp = 1d0+xp
	umxp = 1d0-xp
	upy = 1d0+y
	umy = 1d0-y
	upz = 1d0+z
	umz = 1d0-z
*
c        Nc = 3.d0
c	Nl = 5.d0
	Nh = 1d0

      pt=dsqrt(p(4,1)**2+p(4,2)**2)
c      write(*,*)pt,dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


c      xlf=dfloat(nflav)
c      xmu=scale
      rmuom2=dlog(scale**2/mt**2)
      naem=0
      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))


c      call dotem(4,p,s)
c      call virteval(t1,ro,qqsym,qqasy,ggsym)  
      




      xlf=dfloat(nflav)
      t2=1d0-t1
      tbar=t1-0.25d0*ro
      ubar=t2-0.25d0*ro
      b=dsqrt(1d0-ro)   
      xlp=0.5d0*(1d0+b)
      xlm=0.5d0*(1d0-b)
      vlpm=dlog(xlp/xlm)     
      vlsm=dlog(4d0/ro)      
      vltm=dlog(4d0*t1/ro)      
      vlwm=dlog(4d0*t2/ro)      
      vlbl=dlog(b/xlm)       
      vdw=myli2(1d0-ro/(4d0*t2))-0.5d0*vlwm**2      
      vdt=myli2(1d0-ro/(4d0*t1))-0.5d0*vltm**2
      vdmp=myli2(-xlm/xlp)    
      vdmb=myli2(-xlm/b)+0.5d0*vlbl**2  

c--- Q-Qbar and Qbar-Q contributions

      f1=(vlpm**2/2d0-2d0*vdmb-pisq/3d0)/b
      f2=(-b*vlsm+vlpm**2/4d0+vdmp+pisq/12d0)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75d0*vlpm**2
     . +3d0*vdmp+pisq/4d0+2d0*b**3)/b**5
c      f4 = (vlpm**2/4d0+vdmp+pisq/12d0)/b
c      f5t1 = (vltm**2+vdt+pisq/6d0)/t1**3
c      f5t2 = (vlwm**2+vdw+pisq/6d0)/t2**3

      
c---  Singular parts
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
c      write(*,*)gsq
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
c      qqQQv_0=msqc/aveqq
c      write(*,*)gsq,'b',dsqrt(q2)

c      write(*,*)qqQQv_0*aveqq,'b'
c      if (scheme .eq. 'tH-V') then
      qqQQv_1=-64d0*gsq**2/4d0
ch       qqQQv_1=0d0
c      else
c      qqQQv_1=0d0
c      endif 
      
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

ch      if(j.gt.0.and.k.lt.0)then
ch      tH1stfs=tH1stfs+(
ch     .      
ch     .       -dlog(scale**2/(q2-2d0*mt**2))*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       +dlog(scale**2/q2)*
ch     .  dlog((1d0+beta34)/(1d0-beta34))/beta34*Tqq(3,4)
ch     .       
ch     .       
ch     .         -2d0*LR)*Tqq(3,3))
c     .        *msqc(j,k)
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq
c      write(*,*)L34(q2,beta34,mtrans)

ch      do i1=1,2
ch      do j1=3,4
ch      tH1stfs=tH1stfs+(
ch     .       
ch
ch     .       -0.5d0*LR**2
ch     .       -myli2(-1d0*pt**2/mt**2)
ch     .       -dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))*LR
ch     .        )*Tqq(i1,j1) 
c     .*msqc(j,k))
ch     .       *fx10(j)*fx20(k)*msqc(j,k)*flgq) 
ch      enddo
ch      enddo  


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


      Facdif=Facdif-qqQQv_1/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
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

c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
ch      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
ch      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
ch     .           /gsq**2*4d0
ch      write(*,*)H1qqdelta/gsq**2*4d0
      H1qqdelta=(H1qqdelta-Facdif)*aveqq
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
c      include 'breit.f'
      include 'rescoeff.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1,ggQQv_1,
     . ggQQv_2,ggQQov_1,ggQQsv_1,ggQQsv_1t1,ggQQsv_1t2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
c      common/colour/Tqq
      common/colour/Tqq,Tgg,Tsinglet,Toktsym,Toktasym
      double precision Tsinglet
      double precision Toktsym
      double precision Toktasym
      double precision Tqq(1:4,1:4)
      double precision Tgg(1:4,1:4)



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtgg

c      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
c      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )



        sn=2d0*dot(p,1,2)
        tn=2d0*dot(p,1,3)+mt**2
	un = 2d0*mt**2 - sn - tn

        x  = 4d0*mt**2/(Dsqrt(sn-4d0*mt**2)+Dsqrt(sn))**2 
        y  = - tn/mt**2
        z  = - un/mt**2

        upx = 1d0+x
	umx = 1d0-x
	upxp = 1d0+xp
	umxp = 1d0-xp
	upy = 1d0+y
	umy = 1d0-y
	upz = 1d0+z
	umz = 1d0-z
*
c        Nc = 3.d0
c	Nl = 5.d0
	Nh = 1d0

      pt=dsqrt(p(4,1)**2+p(4,2)**2)
c      write(*,*)pt,dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)


c      xlf=dfloat(nflav)
c      xmu=scale
      rmuom2=dlog(scale**2/mt**2)
      naem=0
      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      ro=4d0*mt**2/s12
      t1=-s13/s12
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))


c      call dotem(4,p,s)
c      call virteval(t1,ro,qqsym,qqasy,ggsym)  
      




      xlf=dfloat(nflav)
      t2=1d0-t1
      tbar=t1-0.25d0*ro
      ubar=t2-0.25d0*ro
      b=dsqrt(1d0-ro)   
      xlp=0.5d0*(1d0+b)
      xlm=0.5d0*(1d0-b)
      vlpm=dlog(xlp/xlm)     
      vlsm=dlog(4d0/ro)      
      vltm=dlog(4d0*t1/ro)      
      vlwm=dlog(4d0*t2/ro)      
      vlbl=dlog(b/xlm)       
      vdw=myli2(1d0-ro/(4d0*t2))-0.5d0*vlwm**2      
      vdt=myli2(1d0-ro/(4d0*t1))-0.5d0*vltm**2
      vdmp=myli2(-xlm/xlp)    
      vdmb=myli2(-xlm/b)+0.5d0*vlbl**2  

c--- Q-Qbar and Qbar-Q contributions

      f1=(vlpm**2/2d0-2d0*vdmb-pisq/3d0)/b
      f2=(-b*vlsm+vlpm**2/4d0+vdmp+pisq/12d0)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75d0*vlpm**2
     . +3d0*vdmp+pisq/4d0+2d0*b**3)/b**5
c      f4 = (vlpm**2/4d0+vdmp+pisq/12d0)/b
c      f5t1 = (vltm**2+vdt+pisq/6d0)/t1**3
c      f5t2 = (vlwm**2+vdw+pisq/6d0)/t2**3

      
c---  Singular parts
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
c      write(*,*)gsq
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
c      qqQQv_0=msqc/aveqq
c      write(*,*)gsq,'b',dsqrt(q2)

c      write(*,*)qqQQv_0*aveqq,'b'
c      if (scheme .eq. 'tH-V') then
c      qqQQv_1=-64d0*gsq**2/4d0
c      ggQQv_1=8d0/xn*V*(V/t1/t2-2d0*xnsq)
c     . *(-(t1**2+t2**2)-1d0)*gsq**2/4d0*avegg
c      write(*,*)
c     . +16d0/xn*V*(t1**2+t1*t2+t2**2)
c     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
c     . *gsq**2/4d0*avegg,'a'
c      write(*,*)ggQQv_1,'r'
c      write(*,*)Tgg(1,1),'t'
      ggQQv_1=-16d0/xn*V*(t1**2+t1*t2+t2**2)
     . *(V*(t1**2+t2**2)-2d0*t1*t2)/t1/t2
     . *gsq**2/4d0*avegg

ch      write(*,*)ggQQv_1,8d0/xn*V*(V/t1/t2-2d0*xnsq)
ch     . *(-(t1**2+t2**2)-1d0)
ch     . *gsq**2/4d0*avegg

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0*avegg
      ggQQov_1=8d0/xn*V*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0*avegg
      ggQQsv_1=8d0/xn*V*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0*avegg
      ggQQsv_1t1=16d0*xnsq*V
     . *(-2d0*t2/t1+2d0*t2**2)
     . *gsq**2/4d0*avegg
      ggQQsv_1t2=16d0*xnsq*V
     . *(-2d0*t1/t2+2d0*t1**2)
     . *gsq**2/4d0*avegg

ch      write(*,*)ggQQv_1/gsq**2*4d0/avegg,
ch     .          ggQQv_2/gsq**2*4d0/avegg,
ch     .      1d0/xn*(ggQQov_1-V*ggQQsv_1)/2d0/gsq**2*4d0/avegg,
ch     .      1d0/xn*(ggQQov_1)/2d0/gsq**2*4d0/avegg,
ch     .   (ggQQsv_1*xn/2d0-ggQQsv_1t1/4d0)/gsq**2*4d0/avegg,
ch     .   (ggQQsv_1-ggQQsv_1t1/4d0)/gsq**2*4d0/avegg,
ch     .   (ggQQsv_1*xn/2d0-ggQQsv_1t2/4d0)/gsq**2*4d0/avegg,
ch     .          s12,s13+mt**2,s23+mt**2

ch      ggQQv_1=0d0

ch      ggQQv_2=1d0
ch      ggQQov_1=0d0
ch      ggQQsv_1=0d0
ch      ggQQsv_1t1=0d0
ch      ggQQsv_1t2=0d0
    
c      else
c      qqQQv_1=0d0
c      endif      
      
c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth

c--- overall factor in the virtual terms is (fourpi*mass^2)^(epsilon)      
c--- to match with our usual definition, we should thus multiply by
c--- (musq/mass^2)^(epsilon)
c--- Note that the overall factor in this definition of the virtual
c--- functions also differs from ours in that we have 1/Gamma(1-ep).
c--- This gives an extra (-pisqo6) compared to Gamma(1-ep)/Gamma(1-2*ep)
c      epin=epinv+rmuom2
c      epin2=epinv**2+epinv*rmuom2+half*rmuom2**2-pisqo6
      
      PoleEp2=-1d0/2d0*(Tgg(1,1)+Tgg(2,2))
      PoleEp1=
     .   -1d0/2d0*(
     .   +Tgg(3,3)+Tgg(4,4)
     .  -2d0*B1g*Tgg(1,1)/ca
     .  -(Tgg(1,1)+Tgg(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tgg(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tgg(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tgg(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tgg(2,4))
      
ch    this is the additional 1/ep pole
c     .  -ggQQv_1/2d0*(Tgg(1,1)+Tgg(2,2))
     .  -ggQQv_1/2d0*(xn+xn)



ch    Begin the IR finite terms

      Facdif=
     .       -1d0/2d0*(
     .      1d0/v34*L34(p)*Tgg(3,4)
     .     +dlog(1d0+pt**2/mt**2)*(Tgg(3,3)+Tgg(4,4))
     .                    )
      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       +1d0/2d0*
     .        myli2(-1d0*pt**2/mt**2)*Tgg(i1,j1)
      enddo
      enddo

ch   End the IR finite terms

ch      Facdif=
ch     .       1d0/2d0*(
ch     .  0.5d0*dlog(scale**2/mt**2)*(Tgg(3,3)+Tgg(4,4)))

ch      do i1=1,2
ch      do j1=3,4
ch      Facdif=Facdif+1d0/2d0*((
ch     .  +dlog(scale**2/(-2d0*dot(p,i1,j1)))*
ch     .   dlog(mt**2/(-2d0*dot(p,i1,j1)))
ch     .  +0.5d0*(dlog(scale**2/(-2d0*dot(p,i1,j1))))**2
ch     .  -0.5d0*(dlog(scale**2/q2))**2)*Tgg(i1,j1))
ch      enddo
ch      enddo

ch Finite terms appearing from the expansion of (scale**2/q2)^ep-
ch interfered with the poles.
      Facdif=Facdif-1d0/2d0*(
     .  (-2d0*B1g*Tgg(1,1)/ca
     .  +Tgg(3,3)+Tgg(4,4)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tgg(3,4))
     .  *dlog(scale**2/q2)
     .  +(Tgg(1,1)+Tgg(2,2))
     .  *0.5d0*(dlog(scale**2/q2))**2)

      do i1=1,2
      do j1=3,4
      Facdif=Facdif
     .       -1d0/2d0*
     .   dlog(4d0*dot(p,i1,j1)**2/(q2*mt**2))           
     .  *dlog(scale**2/q2)*Tgg(i1,j1)
      enddo
      enddo


ch    pi**2/12 term

      Facdif=Facdif+1d0/2d0
     .        *pisq/12d0*(Tgg(1,1)+Tgg(2,2))




      Facdif=Facdif-
ch   New Finite term comiing from ep^2 term of the Born
     .   ggQQv_2/2d0*(
     .   +ca+ca)-
ch
     .   ggQQv_1/2d0*(
     .   +cf+cf
     .  -2d0*B1g)
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  -1d0/xn*(ggQQov_1-V*ggQQsv_1)/4d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   (ggQQsv_1*xn/2d0-ggQQsv_1t1/4d0)/2d0*
c     .   (ggQQsv_1-ggQQsv_1t1/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   (ggQQsv_1*xn/2d0-ggQQsv_1t2/4d0)/2d0*
c     .   (ggQQsv_1-ggQQsv_1t2/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -ggQQv_1/2d0*(xn+xn)*dlog(scale**2/q2)

      H1ggdelta=FinVirtgg(sn,tn)


      H1ggdelta=gsq**2/4d0*H1ggdelta*avegg
     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .         +pisq/12d0*PoleEp2
     .         +dlog(scale**2/mt**2)*(
     .         PoleEp1
ch     .         +1d0/3d0*PoleEp1/2d0/ca
     .                            )
ch     .         -0.5d0*dlog(q2/mt**2)**2*PoleEp2
ch     .          )

ch
ch     .         +2d0*(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
ch     .         +3d0*PoleEp2

ch this is checked and agrees numerically
c      H1ggdelta=gsq**2/4d0*H1ggdelta
c     .         +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
c     .         /avegg
c     .         +pisq/12d0*PoleEp2
c     .         +dlog(scale**2/mt**2)*PoleEp1
c     .         /avegg

ch

c      write(*,*)H1ggdelta/gsq**2*4d0,'f'
c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
ch      write(*,*)PoleEp2/gsq**2*4d0/avegg,sn/1d6,tn/1d6,mt,scale
ch      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
ch     .           /gsq**2*4d0/avegg
ch      write(*,*)H1ggdelta/gsq**2*4d0/avegg
      H1ggdelta=H1ggdelta-Facdif
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


      return
      end

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
      double precision t1,u1

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

      s12=2d0*dot(p,1,2)
      s34=2d0*dot(p,3,4)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
ch      write(*,*)t1+u1+s12,'e'
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
ch      write(*,*)mt**2-u1*t1/s12,pt**2,'t'
      mtrans=dsqrt(mt**2+pt**2)
      a=pt/mt
      I3c=dsqrt(1d0-4d0*mtrans**2/q2)

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

      msqDGav=msqDGav*avegg
ch      write(*,*)avegg,gsq,xn,'l'

      return
      end
      








      
