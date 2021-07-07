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
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'

c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,l,nvec,flgq,flqq
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
      double precision diff1fgq,diff2fgq,diffg1fgq,diffg2fgq,
     .                 diffc1fgq,diffc2fgq
      double precision diff1fgqb,diff2fgqb,diffg1fgqb,diffg2fgqb,
     .                 diffc1fgqb,diffc2fgqb
      double precision diff10q,diff20q,diffc10q,diffc20q,
     .                 diffg10q,diffg20q
      double precision diff1fq,diff2fq,diffg1fq,diffg2fq,
     .                 diffc1fq,diffc2fq
      double precision Pggreg,D0int,D1int,Cgq,Pgq,LF,LR,H2st,H2stgq
      double precision H2ggREG,dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision spgq1q,spgq2q,spgq1qb,spgq2qb
      double precision Pggggreg,Pgggq,Pgqqg,Pgqqq
      double precision CgqPqq,CgqPqg,P2gg,P2gq
      double precision Pqqint,Cqq,Cqg,Pqq,Pqg
      double precision C2qqreg,C2qqp,C2qqb,C2qg,C2gq,C2qgnew

      double precision H1st

      double precision H1qdelta,H1qbdelta,H1gdelta,aaa
      double precision H1ggdelta,H1qqdelta

      common/aa/aaa

      external Pggreg,D0int,Cgq,Pgq,H2stgq,H2ggREG,P2gg,P2gq,
     &         Pggggreg,Pgqqg,Pgqqq,Pgggq,CgqPqq,CgqPqg


ch variables and functions needed for the final state radiation

      double precision s34,beta34,L34,pt,bj,tH1stfs,mtrans

ch    kinematical factor needed for the matrix element in the gg channel

      double precision kinasym

      double precision msqcgg

ch
      integer i1,j1,i1b,i
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
      s34=2*dot(p,3,4)
      q2=s34+dot(p,3,3)+dot(p,4,4)
      if(dynamicscale) call scaleset(q2)
      call qqb_QQb(p,msqc,1)
      call col_operators(p)
            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
c      flux=1d0/(2d0*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0

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
ch      xx10=0.03806738946629785d0
ch      xx20=0.536062971556476d0
ch      alfa=0.03469564397716076d0
ch      beta=0.870166713756036d0
ch      msqc(0,0)=0.589741270721332d0
ch      facscale=173.3d0
ch      scale=173.3d0
ch      asopi=0.03400021652486929d0
ch      LR=2.38614996331854d0
ch      LF=2.38614996331854d0
chxmsq=  1.185379790314491E-003 xx10,xx20=  3.806738946629785E-002
ch  0.536062971556476      alfa,beta=  3.469564397716076E-002
ch  0.870166713756036      asopi=  3.400021652486929E-002 scale,facscale=
ch   173.300000000000        173.300000000000      lr,lf=
ch2.38614996331854     
ch   2.38614996331854      msqb,msqc=  0.589741270721332     
ch  0.610803290115765     


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
ch        fx10(j)=0d0
ch        fx10(-j)=0d0
ch        fx20(j)=0d0
ch        fx20(-j)=0d0
ch        fx1p(j)=0d0
ch        fx1p(-j)=0d0
ch        fx2p(j)=0d0
ch        fx2p(-j)=0d0
ch        enddo

        flqq=1
        if(qqonly)flqq=0


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

      diff10g=0d0
      diff20g=0d0
      diff1fg=0d0
      diff1fgq=0d0
      diff1fgqb=0d0
      diff2fg=0d0
      diff2fgq=0d0
      diff2fgqb=0d0

      diffc10g=0d0
      diffc1fg=0d0
      diffc1fgq=0d0
      diffc1fgqb=0d0
      diffc20g=0d0
      diffc2fg=0d0
      diffc2fgq=0d0
      diffc2fgqb=0d0

      diffg10g=0d0
      diffg1fg=0d0
      diffg1fgq=0d0
      diffg1fgqb=0d0
      diffg20g=0d0
      diffg2fg=0d0
      diffg2fgq=0d0
      diffg2fgqb=0d0
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
ch      H1qdelta=2*C1qqdelta
ch      H1qbdelta=2*C1qqdelta
ch      H1gdelta=2*C1ggdelta-2*beta0*LR
ch      H1gdelta=0d0
ch      tdelta=msqc(0,0)

c      tdelta=msqc(1,-1)*flgq*fx10(2)*fx20(-2)+
c     .       msqc(1,-1)*flgq*fx10(-2)*fx20(2)
ch      write(*,*)msqc(2,-2),'uuuuu'

ch    check against hnnlo and dynnlo
ch      do j=1,nf
ch       msqc(j,-j)=0.610803290115765d0
ch       msqc(-j,j)=0.610803290115765d0
ch      msqc(j,-j)=0d0
ch      msqc(-j,j)=0d0
ch      enddo
ch      msqc(0,0)=0d0
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
C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq*flqq

      if(j.eq.-k.and.j.ne.0)then
        tdeltaqq=tdeltaqq+fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.eq.0.and.k.eq.0)then
        tdeltagg=tdeltagg+fx10(j)*fx20(k)*msqc(j,k)*flgq
      endif


      if(order.eq.0) goto 75

C     Start H1st: to be used later

ch     H1st delta term
c       if(j.eq.-k.and.j.ne.0)then
        if(j.gt.0.and.k.lt.0)then
         tH1st=tH1st+H1qdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        elseif(j.lt.0.and.k.gt.0)then
         tH1st=tH1st+H1qbdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        elseif(j.eq.0.and.k.eq.0)then
         tH1st=tH1st+H1gdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        endif

ch    H1st quark channel
C     H1st: non delta terms, first leg

      if(j.eq.-k.and.j.ne.0)then
       tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg

      
      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      endif

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


      if(order.eq.1) goto 75

      if(j.eq.-k.and.j.ne.0)then

ch    Contribution from the qqb->QQb hard process      

CC    Start H2 contribution

CC    H2st gg contribution

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flgq
     &                                                   *flqq

CC    H2st qqbar contribution from C1*C1 (without delta term)

C     regular*regular

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flgq
     &                                                   *flqq

C     regular-delta

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq 
     &            fx20(k)*H1qdelta*msqc(j,k)*flgq*flqq       

      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq 
     &            fx10(j)*H1qdelta*msqc(j,k)*flgq*flqq        


CC    H2st qg contribution from C1*C1

C     regular*regular

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flqq

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flqq


C     regular-delta

ch     C1qqdelta=0 in the hard scheme!
      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k) 
ch     &            fx20(k)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)
     &            fx20(k)*H1qdelta*msqc(j,k)*flqq       

      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
ch     &            fx10(j)*C1qqdelta*msqc(j,k)  
ch     &            fx10(j)*(0.5d0*H1qdelta-beta0*LR)*msqc(j,k)
     &            fx10(j)*H1qdelta*msqc(j,k)*flqq

CC    H2st qqbar channel: D0(z), first leg

      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*H2qqD0/(1-z1)

      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)*flgq*flqq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq*flqq

CC    H2st, qqbar channel: D0(z), second leg
      
      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*H2qqD0/(1-z2)

      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)*flgq*flqq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq*flqq

CC    C2qq, regular part, first leg

      tH2st=tH2st+fx1p(j)*C2qqreg(z1)
     &                   *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
     &                                                   *flqq

CC    C2qq, regular part, second leg

      tH2st=tH2st+fx2p(k)*C2qqreg(z2)
     &                   *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
     &                                                   *flqq

CC    C2qg, first leg

      tH2st=tH2st+fx1p(0)*C2qg(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)
     &                                                  *flqq

CC    C2qg, second leg

      tH2st=tH2st+fx2p(0)*C2qg(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)
     &                                                  *flqq

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
     #   +diffg10q*diffg2fq+diffg1fq*diffg20q)*msqc(j,k)*flqq


CC     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

C     First leg

      
      diff1=(-dlog(xx10)*(flgq*(fx1p(j)-fx10(j)*xx10**beta)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq+fx1p(0)*(Pqqqg(z1)+Pqggg(z1)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq)*flqq


C    Second leg

      
      diff2=(-dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq)*flqq


C     Include Pqggq

      do l=1,nf
      if(l.eq.k)then
       diff1=diff1-dlog(xx10)*(fx1p(l)+fx1p(-l)*flqq)
     &                       *Pqggq(z1)*flgq
       diff2=diff2-dlog(xx20)*(fx2p(l)*flqq+fx2p(-l))
     &                       *Pqggq(z2)*flgq
      elseif(l.eq.-k)then
       diff1=diff1-dlog(xx10)*(fx1p(l)*flqq+fx1p(-l))
     &                                *Pqggq(z1)*flgq
       diff2=diff2-dlog(xx20)*(fx2p(l)+fx2p(-l)*flqq)
     &                       *Pqggq(z2)*flgq
      else
       diff1=diff1-dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
       diff2=diff2-dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
      endif
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
     #          +diffg10q*diffc2fq+diffg1fq*diffc20q)*flqq


c    C first leg, gamma second leg

ch      diffc1f=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)
ch      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)+0.5d0*H1qdelta*fx10(j)
      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)

      diffc10q=-dlog(xx10)*fx1p(0)*Cqg(z1)

      tcgaq=tcgaq+msqc(j,k)*
     # (flgq*diffc10q*diffg20q+flgq*diffc1fq*diffg2fq
     #          +diffc10q*diffg2fq+diffc1fq*diffg20q)*flqq
    

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
     &     *(-dlog(xx10))*fx20(k)*msqc(j,k)*flqq 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
     &     *(-dlog(xx20))*fx10(j)*msqc(j,k)*flqq 

c    Add Cqg*Pgq contribution

      do l=1,nf
      if(l.eq.k)then
      tcgaq=tcgaq+(fx1p(l)+fx1p(-l)*flqq)*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcgaq=tcgaq+(fx2p(l)*flqq+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      elseif(l.eq.-k)then
      tcgaq=tcgaq+(fx1p(l)*flqq+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcgaq=tcgaq+(fx2p(l)+fx2p(-l)*flqq)*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      else
      tcgaq=tcgaq+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcgaq=tcgaq+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif 
      enddo

CC  Start 2-loop AP

C   Gluon + pure singlet


      do l=-nf,nf
      if(l.eq.0) then
      tgamma2q=tgamma2q+fx1p(0)*P2qg(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flqq
      tgamma2q=tgamma2q+fx2p(0)*P2qg(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flqq
      elseif(l.eq.-k) then
      tgamma2q=tgamma2q+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq*flqq
      tgamma2q=tgamma2q+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      elseif(l.eq.k) then
      tgamma2q=tgamma2q+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq*flqq
      else
      tgamma2q=tgamma2q+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part

      tgamma2q=tgamma2q+fx1p(j)*P2qqV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq*flqq
      tgamma2q=tgamma2q+fx2p(k)*P2qqV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq*flqq


C   P2qq non-singlet: 1/(1-z)_+


      diff=-dlog(xx10)
     &  *(fx1p(j)-fx10(j)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq
     &                                                    *flqq


      diff=-dlog(xx20)
     &  *(fx2p(k)-fx20(k)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq
     &                                                    *flqq

      

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

      tH2st=tH2st+0.5d0*diff*fx20(0)*msqc(0,0)*flgq*flqq
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx10)*
     .             fx10(0)*fx20(0)*msqc(0,0)*flgq*flqq

c     H2st, gg channel: D0(z), second leg
      
      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*H2ggD0/(1-z2)
     &    +fx2p(0)*H2ggreg(z2))*flqq

      tH2st=tH2st+0.5d0*diff*fx10(0)*msqc(0,0)*flgq*flqq
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx20)*
     .               fx10(0)*fx20(0)*msqc(0,0)*flgq*flqq


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


      tgagag=tgagag+diff*flgq*fx20(0)*msqc(0,0)*flqq


c    Second leg

      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)
     &    *(D0gggg/(1-z2)+D1gggg*dlog(1-z2)/(1-z2))
     &      +fx2p(0)*(Pggggreg(z2)+Pgqqg(z2)))
     &    +(Deltagggg-D0gggg*D0int(xx20)-D1gggg*D1int(xx20))*fx20(0)

      tgagag=tgagag+diff*flgq*fx10(0)*msqc(0,0)*flqq



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

      tcgag=tcgag+CgqPqg(z1)*(-dlog(xx10))*flgq*flqq*
     &            fx1p(0)*fx20(0)*msqc(0,0)

c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)

      tcgag=tcgag+CgqPqg(z2)*(-dlog(xx20))*flgq*flqq*
     &            fx2p(0)*fx10(0)*msqc(0,0)

c    End of (C+C)*(gamma+gamma)


CC    gamma2: diagonal part

c     First leg

      diff=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diff-dlog(xx10)*P2gg(z1)*fx1p(0))
     &                *flgq*flqq*fx20(0)*msqc(0,0)


c     Second leg

      diff=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diff-dlog(xx20)*P2gg(z2)*fx2p(0))
     &                *flgq*flqq*fx10(0)*msqc(0,0)

      endif


 75   continue

      enddo
      enddo

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


ch      tH2st=tH2st+(fx1p(j)+fx1p(-j))*H2stgq(z1)
ch     & *(-dlog(xx10))*fx20(0)*msqc(0,0)

c Second leg


ch      tH2st=tH2st+(fx2p(j)+fx2p(-j))*H2stgq(z2)
ch     & *(-dlog(xx20))*fx10(0)*msqc(0,0)


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
ch    There is no regular-regular contribution in the gluon channel,
ch    since Cgg(z)~delta(1-z)
ch      tH2st=tH2st+(fx1p(j)+fx1p(-j))*Cgq(z1)*(-dlog(xx10))*
ch     &            fx2p(0)*Cgg(z2)*(-dlog(xx20))*msqc(j,k)

ch      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
ch     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)

CC    C2gq, first leg

cch      tH2st=tH2st+(fx1p(j)+fx1p(-j))*C2gq(z1)*(-dlog(xx10))
cch     ,                 *fx20(0)*msqc(0,0)

      tH2st=tH2st+(fx1p(j)+fx1p(-j))*H2stgq(z1)*(-dlog(xx10))
     ,                 *fx20(0)*msqc(0,0)*flqq

ch      write(*,*)H2stgq(z1)-C2gq(z1),'a'
CC    C2gq, second leg

cch      tH2st=tH2st+(fx2p(j)+fx2p(-j))*C2gq(z2)*(-dlog(xx20))
cch     ,                 *fx10(0)*msqc(0,0)
      tH2st=tH2st+(fx2p(j)+fx2p(-j))*H2stgq(z2)*(-dlog(xx20))
     ,                 *fx10(0)*msqc(0,0)*flqq



C     regular-delta

ch     C1qqdelta=0 in the hard scheme!
ch      C1ggdelta=0d0
      tH2st=tH2st+(fx1p(j)+fx1p(-j))*Cgq(z1)*(-dlog(xx10))*
ch     &            fx20(k)*C1qqdelta*msqc(j,k) 
     &            fx20(0)*(H1gdelta
ch     &                      -C1ggdelta
ch     &                     -beta0*LR
     &                       )*msqc(0,0)*flqq

      tH2st=tH2st+(fx2p(j)+fx2p(-j))*Cgq(z2)*(-dlog(xx20))*
ch     &            fx10(0)*0.5d0*H1gdelta*msqc(0,0)  
     &            fx10(0)*(H1gdelta
ch     &                     -C1ggdelta
ch     &                    -beta0*LR
     &                      )*msqc(0,0)*flqq



CCCCC Terms needed for NNLO scale dependence


CC    Now (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


      diff1fg=diff1fg-dlog(xx10)*Pgq(z1)*(fx1p(j)+fx1p(-j))
ch      diff1fgq=diff1fgq-dlog(xx10)*Pgq(z1)*fx1p(j)
ch      diff1fgqb=diff1fgqb-dlog(xx10)*Pgq(z1)*fx1p(-j)

      diff2fg=diff2fg-dlog(xx20)*Pgq(z2)*(fx2p(j)+fx2p(-j))
ch      diff2fgq=diff2fgq-dlog(xx20)*Pgq(z2)*fx2p(j)
ch      diff2fgqb=diff2fgqb-dlog(xx20)*Pgq(z2)*fx2p(-j)


C     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

      tgagag=tgagag-dlog(xx10)*(Pgqqq(z1)+Pgggq(z1))*(fx1p(j)+fx1p(-j))
     &            *fx20(0)*msqc(0,0)*flqq

      tgagag=tgagag-dlog(xx20)*(Pgqqq(z2)+Pgggq(z2))*(fx2p(j)+fx2p(-j))
     &            *fx10(0)*msqc(0,0)*flqq
c      write(*,*)Pgqqq(z1),'Pgqqqz1',Pgggq(z1),'Pgggq(z1)',
ch     .         tgagag,'tgagag'


C    End of (gamma+gamma)*(gamma+gamma) term


C    Start (C+C)*(gamma+gamma)


c    Gamma first leg

      diffg1fg=diffg1fg-dlog(xx10)*(fx1p(j)+fx1p(-j))*Pgq(z1)
ch      diffg1fgq=diffg1fgq-dlog(xx10)*fx1p(j)*Pgq(z1)
ch      diffg1fgqb=diffg1fgqb-dlog(xx10)*fx1p(-j)*Pgq(z1)

c    C second leg

      diffc2fg=diffc2fg-dlog(xx20)*(fx2p(j)+fx2p(-j))*Cgq(z2)
ch      diffc2fgq=diffc2fgq-dlog(xx20)*fx2p(j)*Cgq(z2)
ch      diffc2fgqb=diffc2fgqb-dlog(xx20)*fx2p(-j)*Cgq(z2)

c    Gamma second leg

      diffg2fg=diffg2fg-dlog(xx20)*(fx2p(j)+fx2p(-j))*Pgq(z2)
ch      diffg2fgq=diffg2fgq-dlog(xx20)*fx2p(j)*Pgq(z2)
ch      diffg2fgqb=diffg2fgqb-dlog(xx20)*fx2p(-j)*Pgq(z2)

c    C first leg

      diffc1fg=diffc1fg-dlog(xx10)*(fx1p(j)+fx1p(-j))*Cgq(z1)
ch      diffc1fgq=diffc1fgq-dlog(xx10)*fx1p(j)*Cgq(z1)
ch      diffc1fgqb=diffc1fgqb-dlog(xx10)*fx1p(-j)*Cgq(z1)

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)


      tcgag=tcgag+CgqPqq(z1)*
     .   (-dlog(xx10))*(fx1p(j)+fx1p(-j))*fx20(0)*msqc(0,0)*flqq 


c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)


      tcgag=tcgag+CgqPqq(z2)*
     .     (-dlog(xx20))*(fx2p(j)+fx2p(-j))*fx10(0)*msqc(0,0)*flqq 


CC    gamma2: qg channel


c    First leg

      tgamma2g=tgamma2g
     &       -dlog(xx10)*P2gq(z1)*(fx1p(j)+fx1p(-j))*fx20(0)*msqc(0,0)
     &                                                      *flqq 

c    Second leg

      tgamma2g=tgamma2g
     &       -dlog(xx20)*P2gq(z2)*(fx2p(j)+fx2p(-j))*fx10(0)*msqc(0,0)
     &                                                      *flqq 



CC   Effect of spin correlations

c     First leg

      spgq1=spgq1-dlog(xx10)*(fx1p(j)+fx1p(-j))*Ggq(z1)

c    Second leg

      spgq2=spgq2-dlog(xx20)*(fx2p(j)+fx2p(-j))*Ggq(z2)


      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    Check it !      

      tgagag=tgagag+2*msqc(0,0)
     # *(flgq*flqq*diff10g*diff20g
     #   +diff10g*diff2fg*flqq+diff1fg*diff20g*flqq) 

c

c    gamma first leg, C second leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*flqq*diffg10g*diffc20g
     #  +diffg10g*diffc2fg*flqq
     #  +diffg1fg*diffc20g*flqq)

c    gamma second leg, C first leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*flqq*diffg20g*diffc10g
     #   +diffg20g*diffc1fg*flqq
     #   +diffg2fg*diffc10g*flqq)

      do i=1,nf
      do j=1,nf
      diffg1fgq=-dlog(xx10)*fx1p(i)*Pgq(z1)
      diffg1fgqb=-dlog(xx10)*fx1p(-i)*Pgq(z1)
      diffg2fgq=-dlog(xx20)*fx2p(j)*Pgq(z2)
      diffg2fgqb=-dlog(xx20)*fx2p(-j)*Pgq(z2)
      diffc1fgq=-dlog(xx10)*fx1p(i)*Cgq(z1)
      diffc1fgqb=-dlog(xx10)*fx1p(-i)*Cgq(z1)
      diffc2fgq=-dlog(xx20)*fx2p(j)*Cgq(z2)
      diffc2fgqb=-dlog(xx20)*fx2p(-j)*Cgq(z2)
      diff1fgq=-dlog(xx10)*Pgq(z1)*fx1p(i)
      diff1fgqb=-dlog(xx10)*Pgq(z1)*fx1p(-i)
      diff2fgq=-dlog(xx20)*Pgq(z2)*fx2p(j)
      diff2fgqb=-dlog(xx20)*Pgq(z2)*fx2p(-j)
      if(i.eq.j)then
       tgagag=tgagag+2*msqc(0,0)
     # *(flgq*diff1fgq*diff2fgq
     #   +flgq*diff1fgqb*diff2fgqb
     #   +flgq*flqq*diff1fgq*diff2fgqb
     #   +flgq*flqq*diff1fgqb*diff2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg1fgq*diffc2fgq
     #  +flgq*diffg1fgqb*diffc2fgqb
     #  +flgq*flqq*diffg1fgq*diffc2fgqb
     #  +flgq*flqq*diffg1fgqb*diffc2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg2fgq*diffc1fgq
     #  +flgq*diffg2fgqb*diffc1fgqb
     #  +flgq*flqq*diffg2fgq*diffc1fgqb
     #  +flgq*flqq*diffg2fgqb*diffc1fgq)
      else
       tgagag=tgagag+2*msqc(0,0)
     # *(flgq*diff1fgq*diff2fgq
     #   +flgq*diff1fgqb*diff2fgqb
     #   +flgq*diff1fgq*diff2fgqb
     #   +flgq*diff1fgqb*diff2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg1fgq*diffc2fgq
     #  +flgq*diffg1fgqb*diffc2fgqb
     #  +flgq*diffg1fgq*diffc2fgqb
     #  +flgq*diffg1fgqb*diffc2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg2fgq*diffc1fgq
     #  +flgq*diffg2fgqb*diffc1fgqb
     #  +flgq*diffg2fgq*diffc1fgqb
     #  +flgq*diffg2fgqb*diffc1fgq)
      endif
      enddo
      enddo



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    qq contribution to H2 (C1C1)

      do i=1,nf
      do j=1,nf
      diffc1fgq=-dlog(xx10)*fx1p(i)*Cgq(z1)
      diffc1fgqb=-dlog(xx10)*fx1p(-i)*Cgq(z1)
      diffc2fgq=-dlog(xx20)*fx2p(j)*Cgq(z2)
      diffc2fgqb=-dlog(xx20)*fx2p(-j)*Cgq(z2)
      if(i.eq.j)then
       tH2st=tH2st+msqc(0,0)*flgq*(
     .        diffc1fgq*diffc2fgq+
     .        diffc1fgqb*diffc2fgqb+
     .        diffc1fgq*diffc2fgqb*flqq+
     .        diffc1fgqb*diffc2fgq*flqq)
      else
       tH2st=tH2st+msqc(0,0)*flgq*(
     .        diffc1fgq*diffc2fgq+
     .        diffc1fgqb*diffc2fgqb+
     .        diffc1fgq*diffc2fgqb+
     .        diffc1fgqb*diffc2fgq)
      endif
      enddo
      enddo


chc    spin correlations

ch   spin + anglular correlations

c    gg channel

      
      tH2sp=msqc(0,0)*dlog(xx10)*dlog(xx20)*
     . fx1p(0)*fx2p(0)*Ggg(z1)*Ggg(z2)*flgq*flqq

c    gq+qg channel

ch      tH2sp=tH2sp-msqc(0,0)*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
ch      tH2sp=tH2sp-msqc(0,0)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1
      tH2sp=tH2sp-msqGGav(p)*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
     .                                         *flqq
      tH2sp=tH2sp-msqGGav(p)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1
     .                                         *flqq
ch      write(*,*)msqGGav(p)*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2,'a'
ch      write(*,*)msqGGav(p)*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1,'b'

ch    add also the angular correlations due to D*G interference

      tH2sp=tH2sp+msqDGav(p)*fx10(0)*spgq2*flqq
      tH2sp=tH2sp+msqDGav(p)*fx20(0)*spgq1*flqq
ch      write(*,*)'old',H2stgq(z1),'new',C2gq(z1)
ch      write(*,*)'oldqg',C2qg(z1),'newqg',C2qgnew(z1)
ch      write(*,*)'diff',H2stgq(z1)-C2gq(z1)

c    qq channel

ch      tH2sp=tH2sp+msqc(0,0)*spgq1*spgq2*flgq
      do i=1,nf
      do j=1,nf
      spgq1q=-dlog(xx10)*fx1p(i)*Ggq(z1)
      spgq1qb=-dlog(xx10)*fx1p(-i)*Ggq(z1)
      spgq2q=-dlog(xx20)*fx2p(j)*Ggq(z2)
      spgq2qb=-dlog(xx20)*fx2p(-j)*Ggq(z2)
      if(i.eq.j)then
       tH2sp=tH2sp+msqGGav(p)*flgq*(
     .        spgq1q*spgq2q+
     .        spgq1qb*spgq2qb+
     .        spgq1q*spgq2qb*flqq+
     .        spgq1qb*spgq2q*flqq)
      else
       tH2sp=tH2sp+msqGGav(p)*flgq*(
     .        spgq1q*spgq2q+
     .        spgq1qb*spgq2qb+
     .        spgq1q*spgq2qb+
     .        spgq1qb*spgq2q)
      endif
      enddo
      enddo
ch      tH2sp=tH2sp+msqGGav(p)*spgq1*spgq2*flgq


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 72   xmsq=tdelta
      if(order.eq.1)then
       xmsq=xmsq
     .          +asopi*
     .          (
     .           tH1st*flqq+
     .           LF*tH1stF*flqq
ch     .         -2*beta0*LR*tdelta
     .                )
      elseif(order.eq.2)then
c     Include missing delta term from C*gamma (no factor 2 here !)

ch      xmsq=xmsq+asopi**2*(LF*C1ggdelta*tH1stF)
       xmsq=xmsq
ch     &          +asopi*(tH1st+LF*tH1stF)
     &         +asopi**2*(tdeltaqq*H2qqdelta*flqq
     &                     +tdeltagg*H2ggdelta*flqq
     &                     +tH2st+tH2sp)

CC     add scale dependence at NNLO

       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF*flqq
     &                   +tgamma2q*LF+tgamma2g*LF
     &                   -beta0*LR*(tH1st+LF*tH1stF)*flqq
     &                   +LF*tcgaq+0.5d0*LF**2*tgagaq
     &                   +LF*tcgag+0.5d0*LF**2*tgagag)


      xmsq=xmsq+asopi**2*(LF*H1gdelta*tH1stFg)*flqq
      xmsq=xmsq+asopi**2*(LF*H1qdelta*tH1stFq)*flqq


ch      xmsq=xmsq+asopi**2*(LF*0.5d0*H1qdelta*tH1stFq)

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2gg*tdeltagg)*LF*flqq



      endif

ch      write(*,*)'xmsq=',xmsq,'xx10,xx20=',xx10,xx20,
ch     .          'alfa,beta=',alfa,beta,
ch     .          'asopi=',asopi,'scale,facscale=',
ch     .          scale,facscale,'lr,lf=',LR,LF,
ch     .          'msqb,msqc=',msqc(0,0),msqc(1,-1)   

ch     Checked against SCET result!

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
      if(lowintHst.ne.lowintHst)then
        stop
      endif

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
ch   hard scheme!
     & -1/4d0*(1-z)*CF*(Pi**2/2d0-4)


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


       C2qqb=
     &        C2qqp(z)+
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
ch    hard scheme
     -    -1d0/4d0*z*(1-z)*CF*(PI**2/2d0-4d0)


      return
      end

      double precision function C2qgnew(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3
      external myli2,myli3


      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0

ch  Cgecked against C2qg !
      C2qgnew=CA*(-1/12d0/z*(1-z)*(11*z**2-z+2)*myli2(1-z)
     -            +(2*z**2-2*z+1)*(
     -               myli3(1-z)/8d0-myli2(1-z)*Log(1-z)/8d0
     -               +1/48d0*Log(1-z)**3)
     -            +(2*z**2+2*z+1)*(
     -               3*myli3(-z)/8d0+myli3(1/(1+z))/4d0-
     -               myli2(-z)*Log(z)/8d0-1/24d0*Log(1+z)**3+
     -               1/16d0*Log(z)**2*Log(1+z)+1/48d0*Pi**2*Log(1+z))
     -            +z/4d0*(1+z)*myli2(-z)+z*myli3(z)
     -            -1/2d0*z*myli2(1-z)*Log(z)-z*myli2(z)*Log(z)
     -            -3/8d0*(2*z**2+1)*Z3-149/216d0*z**2
     -            -1/96d0*(44*z**2-12*z+3)*Log(z)**2
     -            +1/72d0*(68*z**2+6*Pi**2*z-30*z+21)*Log(z)
     -            +Pi**2*z/24d0+43*z/48d0+43/108d0/z
     -            +1/48d0*(2*z+1)*Log(z)**3-1/2d0*z*Log(1-z)*Log(z)**2
     -            -1/8d0*(1-z)*z*Log(1-z)**2
     -            +1/4d0*z*(1+z)*Log(z)*Log(1+z)
     -            +1/16d0*(3-4*z)*z*Log(1-z)-35/48d0)
     -       +CF*((2*z**2-2*z+1)*(
     -               Z3-myli3(1-z)/8d0-myli3(z)/8d0
     -               +myli2(1-z)*Log(1-z)/8d0+myli2(z)*Log(z)/8d0
     -               -1/48d0*Log(1-z)**3+1/16d0*Log(z)*Log(1-z)**2
     -               +1/16d0*Log(z)**2*Log(1-z))
     -            -3*z**2/8d0-1/96d0*(4*z**2-2*z+1)*Log(z)**3
     -            +1/64d0*(-8*z**2+12*z+1)*Log(z)**2
     -            +1/32d0*(-8*z**2+23*z+8)*Log(z)+5/24d0*Pi**2*z*(1-z)
     -            +11*z/32d0+1/8d0*z*(1-z)*Log(1-z)**2
     -            -1/4d0*z*(1-z)*Log(1-z)*Log(z)
     -            -1/16d0*z*(3-4*z)*Log(1-z)
     -            -9/32d0)
     -       -CF/4d0*(z*Log(z)+1/2d0*(1-z**2)
     -                +(Pi**2/2d0-4)*z*(1-z))
ch     hard scheme
     -    -1d0/4d0*z*(1-z)*CF*(PI**2/2d0-4d0)


      return
      end

      double precision function C2gq(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3,H2gq
      external myli2,myli3
      integer nf

      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0
      nf=5

ch    H2gq Does not agree with H2stgqfull
ch    Eq. (23) of arxiv:1106.4652v2
      H2gq=CF**2*(1/48d0*(2-z)*Log(z)**3-1/32d0*(3*z+4)*Log(z)**2+
     -            5/16d0*(z-3)*Log(z)+
     -            1/12d0*(1/z+z/2d0-1)*Log(1-z)**3+
     -            1/16d0*(z+6/z-6)*Log(1-z)**2+
     -            (5*z/8d0+2/z-2)*Log(1-z)+5/8d0-13/16d0*z)+
     -     CF*nf*(1/24d0/z*(1+(1-z)**2)*Log(1-z)**2+
     -            1/18d0*(z+5/z-5)*Log(1-z)-14/27d0+14/27d0/z+
     -            13/108d0*z)+
     -     CF*CA*(-(1+(1+z)**2)/2d0/z*myli3(1/(1+z))+
     -            (1/2d0-5/2d0/z-5/4d0*z)*myli3(z)-
     -            3/4d0/z*(1+(1+z)**2)*myli3(-z)+
     -            (2-11/6d0/z-z/2d0+z**2/3d0+
     -              (-1/2d0+3/2d0/z+3*z/4d0)*Log(z))*myli2(z)+
     -            (z/4d0+(1+(1+z)**2)*Log(z)/4d0/z)*myli2(-z)+
     -            (1+(1+z)**2)*Log(1+z)**3/12d0/z-
     -            1/24d0/z*((1+(1+z)**2)*(3*Log(z)**2+Pi**2)
     -              -6*z**2*Log(z))*Log(1+z)-
     -            (1+(1-z)**2)*Log(1-z)**3/24d0/z+
     -            1/48d0/z*(6*(1+(1-z)**2)*Log(z)-5*z**2
     -              -22*(1-z))*Log(1-z)**2+
     -            1/72d0/z*(-152+152*z-43*z**2+
     -              6*(-22+24*z-9*z**2+4*z**3)*Log(z)+
     -              9*(1+(1-z)**2)*Log(z)**2)*Log(1-z)-
     -            1/12d0*(1+z/2d0)*Log(z)**3+
     -            1/48d0*(36+9*z+8*z**2)*Log(z)**2+
     -            (-107/24d0-1/z+z/12d0-11/9d0*z**2)*Log(z)+
     -            1/z*(4*Z3-503/54d0+11/36d0*Pi**2)+
     -            1007/108d0-Pi**2/3d0-5/2d0*Z3+
     -            z*(Pi**2/3d0+2*Z3-133/108d0)+
     -            z**2*(38/27d0-Pi**2/18d0))

ch Eq. (30) of arxiv:1106.4652v2
      C2gq=H2gq+CF**2*3/8d0*z
     .         +CF*CA/z*((1+z)*Log(z)+2*(1-z)-(5+Pi**2)/8d0*z**2)
ch      hard scheme
     .         -1d0/4d0*z*CF*((5+Pi**2)*CA-3d0*CF)/2d0
ch       C2gq=H2gq

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
ch     #      -1d0/2d0*z*CF*(11d0+3d0*Pi**2)/2d0


C     Subtract Ggg*Ggq term
    
ch      H2stgq=H2stgqfull-(-4)*(2*(1-z)+(1+z)*dlog(z))/z
      H2stgq=H2stgqfull

      H2stgq=H2stgqfull+CF**2*3/8d0*z
     .         +CF*CA/z*((1+z)*Log(z)+2*(1-z)-(5+Pi**2)/8d0*z**2)
ch      hard scheme
     .         -1d0/4d0*z*CF*((5+Pi**2)*CA-3d0*CF)/2d0


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

      L34new=0.5d0*dlog((1d0+beta34)/(1d0-beta34))*dlog(mtrans**4/mt**4)
     .    -2d0*myli2(2d0*beta34/(1d0+beta34))
     .    -0.25d0*dlog((1d0+beta34)/(1d0-beta34))**2
       L34new=L34new+
     .         2d0*myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*dexp(y34))
     .       +2d0*myli2(1d0-dsqrt((1d0-beta34)/(1d0+beta34))*dexp(-y34))
     .        +y34**2
ch      write(*,*)L34new/L34,'aaa' L34=L34new

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
c      include 'breit.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin



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


      Facdif=Facdif
     .     -qqQQv_1/2d0*(
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

c      write(*,*)H1qqdelta/gsq**2*4d0,'fqq'

c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
c      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
c      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
c     .           /gsq**2*4d0,'dqq'
ch      write(*,*)H1qqdelta/gsq**2*4d0
      H1qqdelta=(H1qqdelta-Facdif)*aveqq

ccch
ch      H1qqdelta=qqQQv_0/2d0*aveqq*(
ch     .   +Tqq(3,3)+Tqq(4,4)
ch     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
ch     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
ch     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
ch     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
ch     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))






ccch
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
c      include 'breit.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtT34qq

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

cch      write(*,*)qqQQv_0/2d0*aveqq*(
cch     .   +Tqq(3,3)+Tqq(4,4)
cch     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
cch     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
cch     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
cch     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
cch     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4)),'low'

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


      Facdif=Facdif
     .     -qqQQv_1/2d0*(
     .   +Tqq(3,3)+Tqq(4,4)
     .  -2d0*B1q
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))
     .  -qqQQv_1/2d0*(Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T34H1qq=FinVirtT34qq(sn,tn)
ch      write(*,*)Tqq(3,4)*PoleEp2*4d0/gsq**2,
ch     .   Tqq(3,4)*PoleEp1*4d0/gsq**2,'bbb'


      T34H1qq=gsq**2/4d0*T34H1qq
     .        +Tqq(3,4)*(
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1)

c      write(*,*)T34H1qq/gsq**2*4d0,'fqq34'

c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
c      write(*,*)Tqq(3,4)*PoleEp2/gsq**2*4d0,
c     .           sn/1d6,tn/1d6,mt,scale
c      write(*,*)Tqq(3,4)*((PoleEp1+dlog(scale**2/mt**2)
c     .          *PoleEp2)
c     .           /gsq**2*4d0),'dqq34'
ch      write(*,*)H1qqdelta/gsq**2*4d0
      T34H1qq=(T34H1qq-Tqq(3,4)*Facdif)*aveqq

ccch

ch      T34H1qq=qqQQv_0/2d0*aveqq*Tqq(3,4)*(
ch     .   +Tqq(3,3)+Tqq(4,4)
ch     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4)
ch     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq(1,3)
ch     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq(1,4)
ch     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq(2,3)
ch     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq(2,4))





ccch
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
c      include 'breit.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtT13qq,FinVirtqq

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
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(1,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(1,3,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(1,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(1,3,2,4))
     .  -qqQQv_1/2d0*Tqq(1,3)*
     .     (Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T13H1qq=FinVirtT13qq(sn,tn)

ch      write(*,*)PoleEp2*4d0/gsq**2,
ch     .          PoleEp1*4d0/gsq**2,'ddd'

ch      write(*,*)scale,'aa',dsqrt(q2),'aaa'
ch      scale=mt
ch      write(*,*)gsq**2/4d0*(T13H1qq-FinVirtT13qq(sn,tn))
ch     .             *aveqq,'a'
ch      scale=dsqrt(q2)


      T13H1qq=gsq**2/4d0*T13H1qq
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1

c      write(*,*)T13H1qq/gsq**2*4d0,'fqq13'

ch      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
ch     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
c      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
c      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
c     .           /gsq**2*4d0,'dqq13'
ch      write(*,*)H1qqdelta/gsq**2*4d0
      T13H1qq=(T13H1qq-Facdif)*aveqq

ccch
ch      T13H1qq=qqQQv_0/2d0*aveqq*(
ch     .   +Tqq(1,3)*(Tqq(3,3)+Tqq(4,4)
ch     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
ch     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(1,3,1,3)
ch     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(1,4,1,3)
ch     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(2,3,1,3)
ch     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(2,4,1,3))






ccch
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
c      include 'breit.f'
      include 'rescoeff.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . qqQQv_0,qqQQv_1
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtT23qq,FinVirtqq,
     &                 FinVirtT34qq,FinvirtT13qq

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
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  +dlog((1d0+v34)/(1d0-v34))/v34*Tqq(3,4))
     .  +dlog((2d0*dot(p,1,3))**2/(q2*mt**2))*Tqq4(2,3,1,3)
     .  +dlog((2d0*dot(p,1,4))**2/(q2*mt**2))*Tqq4(2,3,1,4)
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2))*Tqq4(2,3,2,3)
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2))*Tqq4(2,3,2,4))
     .  -qqQQv_1/2d0*Tqq(2,3)*
     .     (Tqq(1,1)+Tqq(2,2))*dlog(scale**2/q2)

      T23H1qq=FinVirtT23qq(sn,tn)
ch      T23H1qq=-(FinVirtT13qq(sn,tn)+
ch     .          FinVirtT34qq(sn,tn)+
ch     .           cf*Finvirtqq(sn,tn))
ch      write(*,*)PoleEp2*4d0/gsq**2,
ch     .          PoleEp1*4d0/gsq**2,'ddd'
      write(*,*)FinVirtT23qq(sn,tn)+
     . FinVirtT13qq(sn,tn)+
     . FinVirtT34qq(sn,tn)+
     . cf*FinVirtqq(sn,tn),'aaaa'
      write(*,*)FinVirtT23qq(sn,tn),'bb'


      T23H1qq=gsq**2/4d0*T23H1qq
     .        +0.5d0*PoleEp2*dlog(scale**2/mt**2)**2
     .        +pisq/12d0*PoleEp2
     .        +dlog(scale**2/mt**2)*PoleEp1

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
      T23H1qq=(T23H1qq-Facdif)*aveqq
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
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
      integer j,k,naem,i1,j1
      double precision q2,v34,mtrans,beta34t,c,ct,dot
      double precision myli2,p(mxpart,4),msqc
      double precision s34,L34,pt,s12,s13,s23
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . ggQQv_0,ggQQv_1,
     . ggQQv_2,ggQQov_1,ggQQsv_1,ggQQsv_1t1,ggQQsv_1t2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
      double precision M0T34eps,M0T13eps,M0T23eps,
     .                 M0sqeps



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
      ggQQv_0=gsq**2/16d0/pi**2*M0sq
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
     . *gsq**2/4d0

ch      write(*,*)ggQQv_1,8d0/xn*V*(V/t1/t2-2d0*xnsq)
ch     . *(-(t1**2+t2**2)-1d0)
ch     . *gsq**2/4d0*avegg

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0
      ggQQov_1=8d0/xn*V*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1=8d0/xn*V*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1t1=16d0*xnsq*V
     . *(-2d0*t2/t1+2d0*t2**2)
     . *gsq**2/4d0
      ggQQsv_1t2=16d0*xnsq*V
     . *(-2d0*t1/t2+2d0*t1**2)
     . *gsq**2/4d0


      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)
     . *gsq**2
      M0T34eps=Tggeps(3,4)*gsq**2
      M0T13eps=Tggeps(1,3)*gsq**2
      M0T23eps=Tggeps(2,3)*gsq**2

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
c     .  -ggQQv_1/2d0*(Tgg(1,1)+Tgg(2,2))
     .  -ggQQv_1/2d0*(xn+xn)

ch
ch     Checked !
cch      write(*,*)ggQQv_1,'b',
cch     .                M0sqeps,'a',
cch     .          Mgs0eps(1,1)*c1gs+
cch     .          Mgs0eps(2,2)*c2gs+
cch     .          Mgs0eps(3,3)*c3gs,'aa',
cch     .  1d0/xn*(ggQQov_1-V*ggQQsv_1)
cch     .    /2d0,'b34',
cch     .  M0T34eps,'a34',
cch     .  (ggQQsv_1*xn/2d0-ggQQsv_1t1/4d0)
cch     .   ,'b13',
cch     .   M0T13eps,'a13',
cch     .   (ggQQsv_1*xn/2d0-ggQQsv_1t2/4d0)
cch     .   ,'b23',
cch     .   M0T23eps,'a23'

ch



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
ch   New Finite term comiing from ep^2 term of the Born
     .   ggQQv_2/2d0*(
     .   +ca+ca)-
ch
     .   M0sqeps/2d0*(
     .   +cf+cf
     .  -2d0*B1g)
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  -M0T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T13eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t1/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T23eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t2/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0sqeps/2d0*(xn+xn)*dlog(scale**2/q2)

      H1ggdelta=FinVirtgg(sn,tn)


      H1ggdelta=gsq**2/4d0*H1ggdelta
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
c      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
c      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
c     .           /gsq**2*4d0,'d'
ch      write(*,*)H1ggdelta/gsq**2*4d0/avegg
      H1ggdelta=(H1ggdelta-Facdif)*avegg
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . ggQQv_0,ggQQv_1,
     . ggQQv_2,ggQQov_1,ggQQsv_1,ggQQsv_1t1,ggQQsv_1t2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
      double precision M0T34eps,M0T13eps,M0T23eps,
     .                 M0sqeps,M0T34T34eps,M0T34T13eps,M0T34T23eps



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtgg,FinvirtT34gg

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
      ggQQv_0=gsq**2/16d0/pi**2*M0sq
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
     . *gsq**2/4d0

ch      write(*,*)ggQQv_1,8d0/xn*V*(V/t1/t2-2d0*xnsq)
ch     . *(-(t1**2+t2**2)-1d0)
ch     . *gsq**2/4d0*avegg

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0
      ggQQov_1=8d0/xn*V*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1=8d0/xn*V*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1t1=16d0*xnsq*V
     . *(-2d0*t2/t1+2d0*t2**2)
     . *gsq**2/4d0
      ggQQsv_1t2=16d0*xnsq*V
     . *(-2d0*t1/t2+2d0*t1**2)
     . *gsq**2/4d0


      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)
     . *gsq**2
      M0T34eps=Tggeps(3,4)*gsq**2
      M0T34T34eps=Tgg4eps(3,4,3,4)*gsq**2
      M0T34T13eps=Tgg4eps(3,4,1,3)*gsq**2
      M0T34T23eps=Tgg4eps(3,4,2,3)*gsq**2

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
c     .  -ggQQv_1/2d0*(Tgg(1,1)+Tgg(2,2))
     .  -ggQQv_1/2d0*(xn+xn)*Tgg(3,4)

ch
ch     Checked !
cch      write(*,*)ggQQv_1,'b',
cch     .                M0sqeps,'a',
cch     .          Mgs0eps(1,1)*c1gs+
cch     .          Mgs0eps(2,2)*c2gs+
cch     .          Mgs0eps(3,3)*c3gs,'aa',
cch     .  1d0/xn*(ggQQov_1-V*ggQQsv_1)
cch     .    /2d0,'b34',
cch     .  M0T34eps,'a34',
cch     .  (ggQQsv_1*xn/2d0-ggQQsv_1t1/4d0)
cch     .   ,'b13',
cch     .   M0T13eps,'a13',
cch     .   (ggQQsv_1*xn/2d0-ggQQsv_1t2/4d0)
cch     .   ,'b23',
cch     .   M0T23eps,'a23'

ch



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
     .             -
ch   New Finite term comiing from ep^2 term of the Born
     .   ggQQv_2/2d0*Tgg(3,4)*(
     .   +ca+ca)-
ch
     .   M0T34eps/2d0*(
     .   +cf+cf
     .  -2d0*B1g)
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  -M0T34T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T34T13eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t1/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T34T23eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t2/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0T34eps/2d0*(xn+xn)*dlog(scale**2/q2)

      T34H1gg=FinVirtT34gg(sn,tn)


      T34H1gg=gsq**2/4d0*T34H1gg
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

c      write(*,*)T34H1gg/gsq**2*4d0,'f34'
c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
c      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
c      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
c     .           /gsq**2*4d0,'d34'
ch      write(*,*)H1ggdelta/gsq**2*4d0/avegg
      T34H1gg=(T34H1gg-Facdif)*avegg
c      write(*,*)aveqq

c      H1qqdelta=(-Facdif)*aveqq/2d0
c      write(*,*)x,y,z
     


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
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,
     . ggQQv_0,ggQQv_1,
     . ggQQv_2,ggQQov_1,ggQQsv_1,ggQQsv_1t1,ggQQsv_1t2
      double precision FacDif,Add,PoleEp1,PoleEp2
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
      double precision M0T34eps,M0T13eps,M0T23eps,
     .                 M0sqeps,M0T13T34eps,M0T13T13eps,M0T13T23eps



      double precision en,the,sn,tn,un,x,xp,y,z,z2,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,FinVirtgg,FinvirtT13gg

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
      ggQQv_0=gsq**2/16d0/pi**2*M0sq
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
     . *gsq**2/4d0

ch      write(*,*)ggQQv_1,8d0/xn*V*(V/t1/t2-2d0*xnsq)
ch     . *(-(t1**2+t2**2)-1d0)
ch     . *gsq**2/4d0*avegg

      ggQQv_2=8d0/xn*V*(V/t1/t2-2d0*xnsq)
     . *gsq**2/4d0
      ggQQov_1=8d0/xn*V*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1=8d0/xn*V*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
     . *gsq**2/4d0
      ggQQsv_1t1=16d0*xnsq*V
     . *(-2d0*t2/t1+2d0*t2**2)
     . *gsq**2/4d0
      ggQQsv_1t2=16d0*xnsq*V
     . *(-2d0*t1/t2+2d0*t1**2)
     . *gsq**2/4d0


      M0sqeps=(Mgs0eps(1,1)*c1gs+
     .          Mgs0eps(2,2)*c2gs+
     .          Mgs0eps(3,3)*c3gs)
     . *gsq**2
      M0T13eps=Tggeps(1,3)*gsq**2
      M0T13T34eps=Tgg4eps(1,3,3,4)*gsq**2
      M0T13T13eps=Tgg4eps(1,3,1,3)*gsq**2
      M0T13T23eps=Tgg4eps(1,3,2,3)*gsq**2

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
c     .  -ggQQv_1/2d0*(Tgg(1,1)+Tgg(2,2))
     .  -ggQQv_1/2d0*(xn+xn)*Tgg(1,3)

ch
ch     Checked !
cch      write(*,*)ggQQv_1,'b',
cch     .                M0sqeps,'a',
cch     .          Mgs0eps(1,1)*c1gs+
cch     .          Mgs0eps(2,2)*c2gs+
cch     .          Mgs0eps(3,3)*c3gs,'aa',
cch     .  1d0/xn*(ggQQov_1-V*ggQQsv_1)
cch     .    /2d0,'b34',
cch     .  M0T34eps,'a34',
cch     .  (ggQQsv_1*xn/2d0-ggQQsv_1t1/4d0)
cch     .   ,'b13',
cch     .   M0T13eps,'a13',
cch     .   (ggQQsv_1*xn/2d0-ggQQsv_1t2/4d0)
cch     .   ,'b23',
cch     .   M0T23eps,'a23'

ch



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
     .             -
ch   New Finite term comiing from ep^2 term of the Born
     .   ggQQv_2/2d0*Tgg(1,3)*(
     .   +ca+ca)-
ch
     .   M0T13eps/2d0*(
     .   +cf+cf
     .  -2d0*B1g)
c     .  -(Tqq(1,1)+Tqq(2,2))*dlog(q2/mt**2)
     .  -M0T13T34eps/2d0*
     .   dlog((1d0+v34)/(1d0-v34))/v34-
     .   M0T13T13eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t1/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,3))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,4))**2/(q2*mt**2)))-
     .   M0T13T23eps/2d0*
c     .   (ggQQsv_1-ggQQsv_1t2/4d0)/2d0*
     .  (dlog((2d0*dot(p,1,4))**2/(q2*mt**2))
     .  +dlog((2d0*dot(p,2,3))**2/(q2*mt**2)))
     .  -M0T13eps/2d0*(xn+xn)*dlog(scale**2/q2)

      T13H1gg=FinVirtT13gg(sn,tn)


      T13H1gg=gsq**2/4d0*T13H1gg
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

c      write(*,*)T13H1gg/gsq**2*4d0,'f13'
c      write(*,*)(gsq**2/4d0*FinVirt(sn,tn))
c     .         /gsq**2*4d0*3d0*aveqq,'qqdelta'
c      write(*,*)(H1qqdelta-Facdif)/gsq**2*4d0
c      H1qqdelta=0d0
c      Add=0d0
c      Facdif=0d0
c      write(*,*)PoleEp2/gsq**2*4d0,sn/1d6,tn/1d6,mt,scale
c      write(*,*)(PoleEp1+dlog(scale**2/mt**2)*PoleEp2)
c     .           /gsq**2*4d0,'d13'
ch      write(*,*)H1ggdelta/gsq**2*4d0/avegg
      T13H1gg=(T13H1gg-Facdif)*avegg
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
      double precision t1,u1,factor

      s12=2d0*dot(p,1,2)
      s13=2d0*dot(p,1,3)
      s23=2d0*dot(p,2,3)
      q2=2d0*dot(p,3,4)+2d0*mt**2
      t1=s13
      u1=s23
ch      write(*,*)t1,u1,'eeee'
      v34=dsqrt(1d0-4d0*mt**4/((2d0*dot(p,3,4))**2))      
      pt=dsqrt(p(4,1)**2+p(4,2)**2)
      mtrans=dsqrt(mt**2+pt**2)

      msqGGav=(-4d0*gsq**2*(-1d0 + xn**2)*
     .    ((-1d0 + xn**2)*s12**2 + 2d0*xn**2*s12*u1 + 2d0*xn**2*u1**2)*
     .    (2d0*mt**4*s12**2 + 
     .     2d0*mt**2*s12*(2d0*pt**2*s12 + u1*(s12 + u1)) + 
     .    (pt**2*s12 + u1*(s12 + u1))*(3d0*pt**2*s12 + u1*(s12 + u1))))/
     .    (xn*s12**2*u1**2*(s12 + u1)**2)
ch      write(*,*)pt**2-u1*t1/s12+mt**2,'eeee'
ch      factor=(-mt**4*s12**2+8d0*mt**2*s12*t1*u1+
ch     .        2d0*t1*u1*(t1**2-4d0*t1*u1+u1**2))/mt**4/s12**2

      msqGGav=msqGGav*avegg
ch      write(*,*)factor*avegg-msqGGav,'eeee'
ch      msqGGav=factor*avegg

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
      double precision bbb
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

      msqDGav=-msqDGav*avegg
ch      bbb=(2*Cf*I1-2*Ca*I2+(Ca-Cf)*I3)*mt**2*pt**2/xn/t1**2/u1**2*
ch     .    (-16*gsq**2*(xn**2-1d0)*((xn**2-1d0)*s12**2-2*xn**2*t1*u1
ch     .    ))+
ch     .    (I3-2*I2)*mt**2*pt**2/t1**2/u1**2*
ch     .    (8*gsq**2*xn**2*(xn**2-1)*(u1**2+t1**2))

ch      bbb=bbb*avegg

ch      write(*,*)msqDGav+bbb,'bbb'

      return
      end
      








      
