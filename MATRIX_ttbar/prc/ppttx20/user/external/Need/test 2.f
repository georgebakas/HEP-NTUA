CC    Counterterm to be subtracted from real+virt to get a finite
CC    cross section at qt->0

C     Version that allows to separate also qg channel

C     Scale dependence included up to NNLO

      double precision function countint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'efficiency.f'
      include 'masses.f'
      include 'limits.f'
C
      include 'jetlabel.f'
      include 'qcdcouple.f'
      include 'phasemin.f'
      include 'rescoeff.f'
      include 'dynamicscale.f'
      include 'gammacusp.f'
      include 'projected_amplitudes.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'
C
      integer ih1,ih2,j,k,l,nd,nmax,nmin,nvec,order
      integer i1,j1,i,i1p,j1p
      integer nproc
      common/nproc/nproc
      double precision vector(mxdim),W,val,xint
      double precision sqrts,qtmax
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4)
      double precision m3,m4,m5,qtcut,xqtcut
CC
      logical cuts
      double precision x1,x2,dot,ptrans(mxpart,4)

ch
      double precision ptransb(mxpart,4)

ch
      double precision q2,qt2,shat,Itilde
      double precision fx10(-nf:nf),fx20(-nf:nf)
      double precision fx1p(-nf:nf),fx2p(-nf:nf)
      double precision alfa,beta,diff,Pqq,Pqg,Pqqint,Cqq,Cqg
ch
      double precision Pggreg,Pgq
ch
      double precision xjacq2,xjacqt2,xth,x3,almin,almax
      double precision xmio,fluxborn,pswt0
      double precision shad,yq,zmax,tauh,Vol,y3
      double precision xx0(2),xx10,xx20
      double precision sig1,sig2,LR,LF
      double precision sig11,sig12,sig11q,sig12q,sig11g,sig12g
      double precision sig21,sig22,sig23,sig24
ch
      double precision sig21q,sig22q,sig23q,sig24q
      double precision sig21g,sig22g,sig23g,sig24g
ch
      double precision tdelta,tH1st,tH1stF
ch     .        ,tgaga,tcga,tgamma2
ch
      double precision tdeltaq,tdeltaqb,tH1stq,tH1stFq,tdeltag,
     .                 tH1stFg,tH1stg,tH1stFqb,tH1stqb
      double precision tgagaq,tcgaq,tgamma2q
      double precision tgagag,tcgag,tgamma2g
ch
      double precision LL1,LL2,LL3,LL4
      double precision z1,z2,diff1,diff2,cut
ch
      double precision diff1q,diff2q,diffq
      double precision diff1g,diff2g,diffg
      double precision diff1fg,diff2fg
      double precision diff1fgq,diff2fgq
      double precision diff1fgqb,diff2fgqb
ch
      double precision D0int,D1int
      double precision Pqqqq,Pqqqg,Pqggq,Pqggg
      double precision Pgggq,Pgqqq,CgqPqq,P2gq
      double precision CqqPqq,CqqPqg,CqgPgq,CqgPgg
      double precision P2qg,P2qqV,P2qqbV,P2qqS
      double precision Pggggreg,Pgqqg,
     .                 CgqPqg,P2gg,Cgq
ch      double precision diffg10,diffg20,diffc10,diffc20
ch      double precision diffg1f,diffg2f,diffc1f,diffc2f
ch
      double precision diffg10q,diffg20q,diffc10q,diffc20q
      double precision diffg1fq,diffg2fq,diffc1fq,diffc2fq
      double precision diffg10g,diffg20g,diffc10g,diffc20g
      double precision diffg1fg,diffg2fg,diffc1fg,diffc2fg
      double precision diffg1fgq,diffg2fgq,diffc1fgq,diffc2fgq
      double precision diffg1fgqb,diffg2fgqb,diffc1fgqb,diffc2fgqb
      double precision diff10g,diff20g
ch
      external Itilde,Pqq,Pqg,Cqq,Cqg,Pqqint,D0int,D1int
      external Pqqqq,Pqqqg,Pqggq,Pqggg,CqqPqq,CqqPqg,CqgPgq,CqgPgg
      external P2qqV,P2qqbV,P2qg,P2qqS,Pggreg,Pgq
      external Pggggreg,Pgqqg
      external CgqPqg,P2gg,Cgq
      external L34,myLi2,myli3
      double precision L34,pt,myLi2,myli3
      double precision beta34,gammacuspprime2q,F1q,F1qb,F1g

ch
      double precision H1qqdelta,H1ggdelta,
     .                 H1qdelta,H1qbdelta,H1gdelta,
     .                 T34H1qq,T13H1qq,T23H1qq,
     .                 Gamma1H1q,Gamma1H1qb,
     .                 T34H1q,T13H1q,T33H1q,T44H1q,
     .                 T14H1q,T23H1q,T24H1q,
     .                 T34H1qb,T13H1qb,T33H1qb,T44H1qb,
     .                 T14H1qb,T23H1qb,T24H1qb,
     .                 Gamma1H1g,
     .                 T34H1g,T13H1g,T33H1g,T44H1g,
     .                 T14H1g,T23H1g,T24H1g,
     .                 T34H1gg,T13H1gg,T14H1gg

ch
      common/xmio/xmio
      common/xx0/xx0
      common/qtcut/xqtcut
      common/nnlo/order 

CC
CC    Variables passed from virtint or lowint
CC
      common/count/qt2,q2,shat


ch
      double precision y34,v34,logfin,betat,vs,vsqmin,vsqmax,
     .               Gamma1qM0M0,Gamma1gM0M0,Gamma1q,Gamma1qb,
     .               Gamma1g,Gamma2q,Gamma2qb,Gamma2g,
     .               Gamma1sq_q,Gamma1sq_qb,Gamma1sq_g,
     .               Gamma1F1g


CC
      integer n2,n3,sgnj,sgnk,flgq,flqq,flqqb
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      logical creatent,dswhisto
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/incldip/incldip
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
      pswt=0d0
      countint=0d0 

      do nd=0,1
      xmsq(nd)=0d0
      enddo     

c      if(zerowidth) then

CC Check if q2 is the proper interval
      
ch      if(q2.lt.wsqmin.or.q2.gt.wsqmax) goto 999
c      xjacq2=pi*mass3*width3
c      else

CC   Generate q2 again, up to wsqmax

      x3=vector(6)
c      vsqmax=1d0/wsqmin
c      vsqmin=1d0/wsqmax
c      vs=(vsqmax-vsqmin)*x3+vsqmin
c      q2=1/vs
c      xjacq2=(vsqmax-vsqmin)*q2**2
      q2=wsqmin+x3*(wsqmax-wsqmin)
      xjacq2=wsqmax-wsqmin
c      xjacq2=shat-4*mt**2

c      endif

CC   Generate qt2 up to qtmax

      xth=vector(3)

CC  Now compute qtcut from xqtcut

       qtcut=xqtcut*dsqrt(q2)

       if(xth.lt.0.02d0) goto 999
       qt2=qtcut**2*dexp(1d0/xth-1)


CC    Jacobian for qt2

      xjacqt2=1d0/xth**2*qt2



      shad=sqrts**2

      xmio=dsqrt(qt2/q2)

      
      npart=3      
      nvec=npart+2

      Vol=1d0
   

CC   Dynamic scale

      if(dynamicscale) call scaleset(q2)
ch      call scaleset(q2)
ch      write(*,*)ason2pi*2d0*pi,'asq2'
ch      call scaleset(mt**2)
CC   LR,LF

ch      xmio=dsqrt(qt2/scale**2)

ch      write(*,*)ason2pi*2d0*pi,'as'

    
      LR=dlog(q2/scale**2)
      LF=dlog(q2/facscale**2) 

   
   

CC   LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
CC   In this way normalization is fixed to dsigma/dqt2


      LL1=Itilde(1)/q2**2
ch      LL1=-1d0/qt2
      LL2=Itilde(2)/q2**2
ch      LL2=4d0/qt2*dlog(xmio)
      LL3=Itilde(3)/q2**2
ch      LL3=-12d0/qt2*dlog(xmio)**2
      LL4=Itilde(4)/q2**2
ch      LL4=32d0/qt2*dlog(xmio)**3

  

CC Generate BORN momenta for counterterm
      
c      call genBORN2(q2,shat,vector,ptrans,pswt0,*999)
ch
ch      write(*,*)shat,'aaa'
      call genBORN2m(q2,shat,vector,ptrans,pswt0,*999)
      v34=dsqrt(1-4d0*mt**4/((2d0*dot(ptrans,3,4))**2))
ch      beta34=-0.5d0*dlog((1d0+v34)/(1d0-v34))
ch      v34=0.66988751038089656d0
      beta34=0.5d0*dlog((1d0+v34)/(1d0-v34))
      pt=dsqrt(ptrans(4,1)**2+ptrans(4,2)**2)

ch
      call storeptilde(1,ptrans)

CC Here we have to check if the counterevent passes the cuts

       jets=0
       incldip(1)=cuts(ptrans,0)
c       if (incldip(1)) goto 999

CC Compute Born matrix element


c      if(nproc.eq.3)then
c      call qqb_z(ptrans,msqc)
c      else
c      call qqb_w(ptrans,msqc)
c      endif
ch

      do i1=1,mxpart
      do j1=1,4
      if(i1.gt.2)then
      ptransb(i1,j1)=ptrans(i1,j1)
      elseif(i1.eq.1)then
      ptransb(1,j1)=ptrans(2,j1)
      elseif(i1.eq.2)then
      ptransb(2,j1)=ptrans(1,j1)
      endif
      enddo
      enddo
      call qqb_QQb(ptrans,msqc,1)
ch      write(*,*)msqc(1,-1),ason2pi*2*pi,'msqc'
ch      write(*,*)msqc(0,0)/avegg/gsq**2,'a'

!     Checked! The new implementation of color correlations are in agreement with the previous one

      call col_operators(ptrans)

      betat=dsqrt(1d0-4*mt**2/q2)
ch  


C Scaled momentum fractions

      cut=1d-7
   

      beta=cut+(1-cut)*vector(8)
      alfa=cut+(1-cut)*vector(9)

      xx10=xx0(1)
      xx20=xx0(2)

ch      xx10=0.426883676710647d0
ch      xx20=0.832506199168213d0
ch      alfa=0.491091579871728d0
ch      beta=0.196308205115838d0
ch      msqc(0,0)=5.24695775775360d0
ch      msqc(0,0)=0d0
ch      LL1=-0.000000004267013419514816d0
ch      LL2=-0.0000001009997713546848d0
ch      LL3=-0.000001792969733154202d0
ch      LL4=-0.00002821048441927710d0
ch      facscale=173.3d0
ch      scale=173.3d0
ch      ason2pi=0.01700411926779227d0
ch      LR=5.24695775775360d0
ch      LF=5.24695775775360d0
ch      do j=1,nf
ch       msqc(j,-j)=6.24349255330604d0
ch       msqc(-j,j)=6.24349255330604d0
ch      msqc(j,-j)=0d0
ch      msqc(-j,j)=0d0
ch      enddo
      

      z1=xx10**beta
      z2=xx20**alfa
ch      write(*,*)z1*z2,xx10*xx20,'a'

      y34=0.5d0*dlog(xx10/xx20)

c      call genCT2(q2,qt2,y34,vector,p,pswt,*999)
ch
       call genCT2m(q2,qt2,y34,vector,p,pswt,*999)
ch

ch      ptrans(1,1)=0.0000000000000000d0
ch      ptrans(1,2)=0.0000000000000000d0
ch      ptrans(1,3)=-11.895367708481288d0
ch      ptrans(1,4)=-11.895367708481288d0
ch      ptrans(2,1)=0.0000000000000000d0
ch      ptrans(2,2)=0.0000000000000000d0
ch      ptrans(2,3)=2962.6349316411265d0
ch      ptrans(2,4)=-2962.6349316411265d0
ch      ptrans(3,1)=-38.027766888667841d0
ch      ptrans(3,2)=-1.5260942462916707d0
ch      ptrans(3,3)=-1961.1857127875460d0
ch      ptrans(3,4)=1969.1954524780372d0
ch      ptrans(4,1)=38.027766888667841d0
ch      ptrans(4,2)=1.5260942462916707d0
ch      ptrans(4,3)=-989.55385114509909d0
ch      ptrans(4,4)=1005.3348468715708d0
ch      q2=140966.52759145221d0
ch      pt=38.058376452332105d0

ch    -2*Re(Gamma1) with as/pi normalization

      Gamma1q=0.5d0*(
     .       Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tqq(2,4))

ch      write(*,*)Gamma1q,'g1q'

      Gamma1qb=0.5d0*(
     .       Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptransb,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptransb,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptransb,2,4)**2/(q2*mt**2))*Tqq(2,4))


ch      write(*,*)Gamma1qb*msqc(1,-1),'g1qb'

      Gamma1g=0.5d0*(
     .       Tgg(3,3)+Tgg(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tgg(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tgg(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tgg(2,4))

      if(order.eq.2)then
      Gamma1sq_q=
     .                0.5d0*(Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4))*Gamma1q
     .               +0.5d0*(Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4))*0.5d0*
     .      (dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tqq(2,4))

      do i1=1,2
      do j1=1,2
       Gamma1sq_q=Gamma1sq_q+1d0/4d0*(
     . dlog(4d0*dot(ptrans,i1,3)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptrans,j1,3)**2/(q2*mt**2))*Tqq4(i1,3,j1,3)+
     . dlog(4d0*dot(ptrans,i1,4)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptrans,j1,3)**2/(q2*mt**2))*Tqq4(i1,4,j1,3)+
     . dlog(4d0*dot(ptrans,i1,3)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptrans,j1,4)**2/(q2*mt**2))*Tqq4(i1,3,j1,4)+
     . dlog(4d0*dot(ptrans,i1,4)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptrans,j1,4)**2/(q2*mt**2))*Tqq4(i1,4,j1,4))
      enddo
      enddo

ch      write(*,*)Gamma1sq_q,'a'

ch      write(*,*)Gamma1sq_q*msqc(1,-1),'count'

      Gamma1sq_qb=
     .                0.5d0*(Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4))*Gamma1qb
     .               +0.5d0*(Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4))*0.5d0*
     .      (dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptransb,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptransb,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptransb,2,4)**2/(q2*mt**2))*Tqq(2,4))

      do i1=1,2
      do j1=1,2
       Gamma1sq_qb=Gamma1sq_qb+1d0/4d0*(
     . dlog(4d0*dot(ptransb,i1,3)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptransb,j1,3)**2/(q2*mt**2))*Tqq4(i1,3,j1,3)+
     . dlog(4d0*dot(ptransb,i1,4)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptransb,j1,3)**2/(q2*mt**2))*Tqq4(i1,4,j1,3)+
     . dlog(4d0*dot(ptransb,i1,3)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptransb,j1,4)**2/(q2*mt**2))*Tqq4(i1,3,j1,4)+
     . dlog(4d0*dot(ptransb,i1,4)**2/(q2*mt**2))*
     . dlog(4d0*dot(ptransb,j1,4)**2/(q2*mt**2))*Tqq4(i1,4,j1,4))
      enddo
      enddo

ch      write(*,*)Gamma1sq_qb*msqc(1,-1),'countqb'

      Gamma1sq_g=0.25d0*(
     .       Tgg(3,3)**2+Tgg(4,4)**2
     .      +2d0*Tgg(3,3)*Tgg(4,4)
     .      +1d0/v34**2*dlog((1d0+v34)/(1d0-v34))**2*Tgg4(3,4,3,4)
     .      +2d0/v34*dlog((1d0+v34)/(1d0-v34))
     .              *Tgg(3,4)*(Tgg(3,3)+Tgg(4,4))
     .      +4d0*dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))**2
     .          *Tgg4(1,3,1,3)
     .      +4d0*dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))**2
     .          *Tgg4(2,3,2,3)
     .      +4d0*dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))
     .         *dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))
     .         *(Tgg4(1,3,2,3)+Tgg4(2,3,1,3))
     .      +4d0*(Tgg(3,3)+Tgg(4,4))*
     .        (dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg(1,3)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg(2,3))
     .      +2d0/v34*dlog((1d0+v34)/(1d0-v34))
     .              *dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))
     .              *Tgg4(3,4,1,3)
     .      +2d0/v34*dlog((1d0+v34)/(1d0-v34))
     .              *dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))
     .              *Tgg4(3,4,2,3)
     .      +2d0/v34*dlog((1d0+v34)/(1d0-v34))
     .              *dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))
     .              *Tgg4(1,3,3,4)
     .      +2d0/v34*dlog((1d0+v34)/(1d0-v34))
     .              *dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))
     .              *Tgg4(2,3,3,4))




ch    Eq.(31) of 1408:4564v1
      F1q=(Tqq(3,3)+Tqq(4,4))*dlog((pt**2+mt**2)/mt**2)
     .    +(Tqq(3,3)+Tqq(4,4)+2d0*Tqq(3,4))*myli2(-pt**2/mt**2)
     .    +1d0/v34*L34(ptrans)*Tqq(3,4)

ch      write(*,*)F1q,'F1q'

      F1qb=(Tqq(3,3)+Tqq(4,4))*dlog((pt**2+mt**2)/mt**2)
     .    +(Tqq(3,3)+Tqq(4,4)+2d0*Tqq(3,4))*myli2(-pt**2/mt**2)
     .    +1d0/v34*L34(ptransb)*Tqq(3,4)

      F1g=(Tgg(3,3)+Tgg(4,4))*dlog((pt**2+mt**2)/mt**2)
     .    +(Tgg(3,3)+Tgg(4,4)+2d0*Tgg(3,4))*myli2(-pt**2/mt**2)
     .    +1d0/v34*L34(ptrans)*Tgg(3,4)

      Gamma1F1g=-0.5d0*(2d0**myli2(-pt**2/mt**2)+1d0/v34*L34(ptrans))
     .          *0.5d0*(
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg4(1,3,3,4)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tgg4(1,4,3,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg4(2,3,3,4)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tgg4(2,4,3,4)
     .      -dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg4(3,4,1,3)
     .      -dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tgg4(3,4,1,4)
     .      -dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg4(3,4,2,3)
     .      -dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tgg4(3,4,2,4))


ch    The part of Eq.(10) of 0908:3676v2 that does not depend on gammacusp
ch    beta34 in that Eq. is equal to beta34 used in the code -i*Pi
ch    So there are some real numbers implicitly written in the Eq., which are
ch    recovered here:  Needs to be checked !

cch   corresponds to beta34=-1/2*log[(1+v)/(1-v)]+I*Pi
ch      v34=0.2d0/msqc(0,0)
ch      beta34=-0.5d0*dlog((1+v34)/(1-v34))

ch      gammacuspprime2q=ca/2d0*
ch     .    (-5d0*Pi**2/6+Z3+beta34**2
ch     .     +1d0/v34**2*(myli3((1d0-v34)/(1d0+v34))
ch     .                  +beta34*myli2((1d0-v34)/(1d0+v34))
ch     .                  +beta34**3/3d0+19d0/6*Pi**2*beta34
ch     .                  -Z3)
ch     .     +1d0/v34*(myli2((1d0-v34)/(1d0+v34))
ch     .                -2d0*beta34*dlog(2d0*v34/(1d0-v34))
ch     .                +5/6d0*pi**2*beta34+3d0*beta34**2
ch     .                -19/6d0*pi**2
ch     .                -beta34**3/3d0))

ch      write(*,*)gammacuspprime2q*2d0/ca,'bb'

cch   corresponds to beta34=1/2*log[(1+v)/(1-v)]-I*Pi
      gammacuspprime2q=ca/2d0*
     .    (-5d0*Pi**2/6+Z3
     .    +beta34**2
     .     +1d0/v34**2*(myli3((1d0-v34)/(1d0+v34))
     .                  +beta34*myli2((1d0-v34)/(1d0+v34))
     .                  +beta34**3/3d0-5d0/6*Pi**2*beta34
     .                  -Z3)
     .     +1d0/v34*(myli2((1d0-v34)/(1d0+v34))
     .                -2d0*beta34*dlog(2d0*v34/(1d0+v34))
     .                +5/6d0*pi**2*beta34-beta34**2
     .                +5/6d0*pi**2
     .                -beta34**3/3d0))
ch      write(*,*)gammacuspprime2q,'aa'
ch     .   )

ch      write(*,*)gammacuspprime2q*2d0/ca,'cc'

ch      v34=0.5d0
ch      beta34=-0.5d0*dlog((1d0+v34)/(1d0-v34))

ch      gammacuspprime2q=ca/2d0*
ch     .    (Pi**2/6+Z3+beta34**2
ch     .     +1d0/v34**2*(myli3((1d0-v34)/(1d0+v34))
ch     .                  -beta34*myli2((1d0-v34)/(1d0+v34))
ch     .                  -beta34**3/3d0-1d0/6*Pi**2*beta34
ch     .                  -Z3)
ch     .     +1d0/v34*(myli2((1d0-v34)/(1d0+v34))
ch     .                +2d0*beta34*dlog(2d0*v34/(1d0-v34))
ch     .                +1d0/6d0*pi**2*beta34+3d0*beta34**2
ch     .                -1/6d0*pi**2
ch     .                +beta34**3/3d0))
ch      write(*,*)gammacuspprime2q,'a'
ch      write(*,*)3d0*beta34**2,
ch     .   2d0*beta34*dlog(2d0*v34/(1d0-v34)),'aa'

ch      write(*,*)dlog((1+betat)/(1-betat))/beta34,'a'
ch      write(*,*)(dexp(beta34)+dexp(-beta34))/
ch     .    (dexp(beta34)-dexp(-beta34))*v34,'b'
ch      write(*,*)dexp(-2*beta34)/((1d0-v34)/(1d0+v34)),'c'


ch    Corresponds to beta43=1/2*Log[(1+v)/(1-v)]
ch      gammacuspprime2q=ca/2d0*
ch     .    (Pi**2/6+Z3+beta34**2
ch     .     +1d0/v34**2*(myli3((1d0-v34)/(1d0+v34))
ch     .                  +beta34*myli2((1d0-v34)/(1d0+v34))
ch     .                  +beta34**3/3d0+1d0/6*Pi**2*beta34
ch     .                  -Z3)
ch     .     +1d0/v34*(myli2((1d0-v34)/(1d0+v34))
ch     .                -2d0*beta34*dlog(2d0*v34/(1d0+v34))
ch     .                -1d0/6d0*pi**2*beta34-beta34**2
ch     .                -1/6d0*pi**2
ch     .                -beta34**3/3d0))


ch    -2*Re(Gamma2) with as/pi normalization
      Gamma2q=0.5d0*gammacusp2q*(
     .      1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tqq(2,4))
     .         +gammacuspprime2q*Tqq(3,4)
     .         -2d0*gammaQ1
     .         +0.5d0*beta0*F1q
ch      write(*,*)Gamma2q,'Gamma2q'
ch      write(*,*)gammacusp2q,'eeeeeee'

ch      write(*,*)1d0/2*Gamma1q*F1q,Gamma1q,'aaa'


      Gamma2qb=0.5d0*gammacusp2q*(
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptransb,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptransb,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptransb,2,4)**2/(q2*mt**2))*Tqq(2,4))
     .         +gammacuspprime2q*Tqq(3,4)
     .         -2d0*gammaQ1
     .         +0.5d0*beta0*F1qb

      Gamma2g=0.5d0*gammacusp2q*(
     .      1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tgg(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tgg(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tgg(2,4))
     .         +gammacuspprime2q*Tgg(3,4)
     .         -2d0*gammaQ1
     .         +0.5d0*beta0*F1g
     .         +0.5d0*Gamma1F1g
      write(*,*)CA**2*msqc(0,0)*4d0,'eps^-4'
      write(*,*)msqc(0,0),'a'

ch      H1qdelta=2*H1qqdelta(ptrans)/msqc(1,-1)
ch      H1qbdelta=2*H1qqdelta(ptransb)/msqc(1,-1)
ch      H1gdelta=2*H1ggdelta(ptrans)/msqc(0,0)
ch       write(*,*)gammaQ1,'ssss'

       H1qdelta=H1qqdelta(ptrans)/msqc(1,-1)
ch     .          +2*beta0*LR
cch      H1qdelta=2*C1qqdelta
cch      H1qbdelta=2*C1qqdelta
ch       write(*,*)H1qdelta*msqc(1,-1),'h1q'


       H1qbdelta=H1qqdelta(ptransb)/msqc(1,-1)
ch       write(*,*)H1qbdelta*msqc(1,-1),'h1qb'
ch      write(*,*)H1qdelta,H1qbdelta,'a'

       H1gdelta=H1ggdelta(ptrans)/msqc(0,0)
ch      write(*,*)ason2pi*2d0*(H1qdelta*msqc(1,-1)
ch     .              +2*beta0*LR*msqc(1,-1)
ch     .    ),'1loopqq'
ch      write(*,*)ason2pi*2d0*(H1gdelta*msqc(0,0)
ch     .              +2*beta0*LR*msqc(0,0)),'1loopgg'
ch      write(*,*)dsqrt(2*dot(ptrans,1,2)),'2'
ch      write(*,*)2*dot(ptrans,1,3),'3'
ch      write(*,*)2*dot(ptrans,2,3),'4'
ch      write(*,*)ason2pi*2d0*pi,'6'
ch      write(*,*)scale,'7' 
ch      write(*,*)dsqrt(qt2),'qt'


ch     The colour-correlated 1-loop amplitudes squared

       T34H1q=T34H1qq(ptrans)/msqc(1,-1)

       T34H1qb=T34H1qq(ptransb)/msqc(1,-1)

       T34H1g=T34H1gg(ptrans)/msqc(0,0)

ch      write(*,*)scale,'bb',dsqrt(q2),'bbb'

       T13H1q=T13H1qq(ptrans)/msqc(1,-1)

ch      write(*,*)scale,'cc',dsqrt(q2),'ccc'

       T13H1qb=T13H1qq(ptransb)/msqc(1,-1)

       T13H1g=T13H1gg(ptrans)/msqc(0,0)
ch      write(*,*)ason2pi*2d0*(T34H1q*msqc(1,-1)
ch     .             +2*beta0*Tqq(3,4)*LR*msqc(1,-1)),'T341loopqq'
ch      write(*,*)ason2pi*2d0*(T13H1q*msqc(1,-1)
ch     .             +2*beta0*LR*Tqq(1,3)*msqc(1,-1)),'T131loopqq'

ch       write(*,*)T23H1qq(ptrans)/msqc(1,-1),'a'

ch     Use colour conservation to compute the rest

       T24H1q=T13H1q
       T24H1qb=T13H1qb
       T24H1g=T13H1g
       T33H1q=Cf*H1qdelta
       T33H1qb=Cf*H1qbdelta
       T33H1g=Cf*H1gdelta
       T44H1q=Cf*H1qdelta
       T44H1qb=Cf*H1qbdelta
       T44H1g=Cf*H1gdelta
       T23H1q=-(T34H1q+T13H1q+T33H1q)
ch       write(*,*)T23H1q,'aa'
       T23H1g=-(T34H1g+T13H1g+T33H1g)
ch       T14H1g=T14H1gg(ptrans)/msqc(0,0)
ch       write(*,*)T23H1qq(ptrans)/msqc(1,-1)/T23H1q,'b'
ch       scale=mt
ch       T23H1q=T23H1qq(ptrans)/msqc(1,-1)
ch       T23H1qb=T23H1qq(ptransb)/msqc(1,-1)
ch       write(*,*)(T13H1q-T23H1q)/2d0/beta0/dlog(q2/mt**2)
ch       write(*,*)T14H1q+T23H1q+
ch     .           T34H1q+T33H1q,'cccc'
ch       write(*,*)T13H1qb+T23H1qb+
ch     .           T34H1qb+T33H1qb,'dddd'
ch       write(*,*)T23H1q,'a'

       
       T23H1qb=-(T34H1qb+T13H1qb+T33H1qb)
       T14H1q=T23H1q
       T14H1qb=T23H1qb
       T14H1g=T23H1g
ch      write(*,*)dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2)
ch     .   )*T13H1qb,'a'
ch      write(*,*)dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2)
ch     .   )*T23H1q,'b'


ch

ch       T14H1g=T14H1gg(ptrans)/msqc(0,0)

ch
ch       write(*,*)T23H1q,'b'
ch       write(*,*)T34H1q/H1qdelta*2*xn,'a'
ch      write(*,*)scale,'dd',dsqrt(q2),'ddd'
ch       write(*,*)-2*beta0*msqc(1,-1)*dlog(q2/mt**2)
ch     .            *(xn**2-2d0)/2d0/xn,'b'

ch Using the colour-corrleated amplitudes compute the action of 
ch  <M0|(-2*Re(Gamma1)|M1>+c.c
      Gamma1H1q=0.5d0*(
     .       T33H1q+T44H1q
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *T34H1q
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*T13H1q
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*T14H1q
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*T23H1q
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*T24H1q)
ch      write(*,*)Gamma1H1q*msqc(1,-1),'low'

      Gamma1H1qb=0.5d0*(
     .       T33H1qb+T44H1qb
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *T34H1qb
     .      +dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2))*T13H1qb
     .      +dlog(4d0*dot(ptransb,1,4)**2/(q2*mt**2))*T14H1qb
     .      +dlog(4d0*dot(ptransb,2,3)**2/(q2*mt**2))*T23H1qb
     .      +dlog(4d0*dot(ptransb,2,4)**2/(q2*mt**2))*T24H1qb)

      Gamma1H1g=0.5d0*(
     .       T33H1g+T44H1g
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *T34H1g
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*T13H1g
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*T14H1g
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*T23H1g
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*T24H1g)
ch      write(*,*)Gamma1H1qb*msqc(1,-1),'lowqb'


ch       write(*,*)T14H1g+T13H1g+
ch     .           T34H1g+cf*H1gdelta,'ccccg'


ch      H1gdelta=2*C1ggdelta-2*beta0*LR
ch      Gamma1q=0d0
ch      Gamma1qb=0d0
ch      Gamma1g=0d0
      else
       H1qdelta=0d0
       H1qbdelta=0d0
       H1gdelta=0d0
       Gamma1H1q=0d0
       Gamma1H1qb=0d0
       Gamma1sq_q=0d0
       Gamma1sq_qb=0d0
       Gamma2q=0d0
       Gamma2qb=0d0
       Gamma2g=0d0
       Gamma1sq_g=0d0
      endif
ch       Gamma2g=0d0
ch       Gamma1sq_g=0d0
ch       Gamma1g=0d0

ch      do j=1,nf
ch       msqc(j,-j)=5.24349255330604d0
ch       msqc(-j,j)=5.24349255330604d0
ch      msqc(j,-j)=0d0
ch      msqc(-j,j)=0d0
ch      enddo
ch      H1qdelta=2d0*C1qqdelta
ch      H1qbdelta=2d0*C1qqdelta
cch      H1qdelta=0d0
cch      H1qbdelta=0d0
ch      H1gdelta=2d0*C1ggdelta-2d0*beta0*LR
ch      Gamma1q=0d0
ch      Gamma1qb=0d0
ch      msqc(0,0)=0d0
ch      Gamma1H1q=0d0
ch      Gamma1H1qb=0d0
ch      Gamma1sq_qb=0d0
ch      Gamma1sq_q=0d0
ch      Gamma2q=0d0
ch      Gamma2qb=0d0
ch      Gamma1g=0d0
ch      write(*,*)ason2pi*2d0*(H1qdelta*msqc(1,-1)
ch     .             +2*beta0*LR*msqc(1,-1)),'1loopqq'
ch      write(*,*)ason2pi*2d0*(H1gdelta*msqc(0,0)
ch     .              +2*beta0*LR*msqc(0,0)),'1loopgg'
ch      write(*,*)dsqrt(2*dot(ptrans,1,2)),'2'
ch      write(*,*)2*dot(ptrans,1,3),'3'
ch      write(*,*)2*dot(ptrans,2,3),'4'
ch      write(*,*)as,'6'
ch      write(*,*)scale,'7'

ch      write(*,*)msqc(1,-1),msqc(2,-2),msqc(0,0)
c      write(*,*)H1qdelta,H1gdelta,'a'

ch      H1qdelta=H1qdelta+cf*pi**2/6d0
ch      H1qbdelta=H1qbdelta+cf*pi**2/6d0
ch      H1gdelta=H1gdelta+ca*pi**2/6d0
ch      H1qdelta=0d0
ch      H1qbdelta=0d0
ch      H1gdelta=0d0
ch        write(*,*)H1qdelta,H1gdelta,'a'
ch        write(*,*)2*H1qqdelta(ptrans)/msqc(1,-1),
ch     .            2*H1qqdelta(ptransb)/msqc(1,-1),
ch     .            2*H1ggdelta(ptrans)/msqc(0,0),'b'



ch
       incldip(1)=cuts(ptrans,0)
       if (incldip(1)) goto 999
ch       if (dsqrt(qt2).lt. 0.0001d0) goto 999
ch       if(xmio.gt.0.0025d0.or.xmio.lt.0.002d0)goto 999

           
c--- calculate PDF's  

c      if(xx10.lt.1d-5)write(*,*)q2,xx10
c      if(xx20.lt.1d-5)write(*,*)q2,xx20

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)


CC Switch off gluon !!

      if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
      endif

CC Gluon only !

      if(ggonly) then
       do j=1,5
       fx10(j)=0d0
       fx10(-j)=0d0
       fx1p(j)=0d0
       fx1p(-j)=0d0
       fx20(j)=0d0
       fx20(-j)=0d0
       fx2p(j)=0d0
       fx2p(-j)=0d0
       enddo
      endif

       flgq=1
       if(gqonly)flgq=0
ch        do j=1,5
ch       if(j.ne.2.and.j.ne.1)then
ch       fx10(j)=0d0
ch       fx10(-j)=0d0
ch       fx1p(j)=0d0
ch       fx1p(-j)=0d0
ch       fx20(j)=0d0
ch       fx20(-j)=0d0
ch       fx2p(j)=0d0
ch       fx2p(-j)=0d0
ch       else
ch       fx10(j)=0d0
ch       fx10(-j)=0d0
ch       fx1p(j)=0d0
ch       fx1p(-j)=0d0
ch       fx20(j)=0d0
ch       fx20(-j)=0d0
ch       fx2p(j)=0d0
ch       fx2p(-j)=0d0
ch       endif
ch       enddo

        
ch        fx10(0)=0d0
ch        fx20(0)=0d0
ch        fx1p(0)=0d0
ch        fx2p(0)=0d0

ch       do j=1,5
ch       fx10(-j)=0d0
ch       fx1p(-j)=0d0
ch       fx20(-j)=0d0
ch       fx2p(-j)=0d0
ch       fx10(j)=0d0
ch       fx1p(j)=0d0
ch       fx20(j)=0d0
ch       fx2p(j)=0d0
ch       enddo
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

ch      do j=1,5
ch       fx10(-j)=0d0
ch       fx1p(-j)=0d0
ch       fx20(j)=0d0
ch       fx2p(j)=0d0
ch      enddo

cch Select only the b quarks
ch      do j=1,5
ch      if(j.ne.5)then
ch       fx10(j)=0d0
ch       fx10(-j)=0d0
ch       fx1p(j)=0d0
ch       fx1p(-j)=0d0
ch       fx20(j)=0d0
ch       fx20(-j)=0d0
ch       fx2p(j)=0d0
ch       fx2p(-j)=0d0
ch      else
ch        fx10(-5)=0d0
ch        fx20(5)=0d0
ch        fx1p(-5)=0d0
ch        fx2p(5)=0d0
ch      endif
ch      enddo

cch

        flqq=1
        if(qqonly)flqq=0


!       flqqb flag works correctly, CHecked !
        flqqb=1
        if(qqbonly)flqqb=0



C Flux for Born cross section


       fluxborn=fbGeV2/(2*q2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start construction of the counterterm

        tdelta=0d0
        tdeltaq=0d0
        tdeltaqb=0d0
        tdeltag=0d0
        tH1st=0d0
        tH1stq=0d0
        tH1stqb=0d0
        tH1stg=0d0
        tH1stF=0d0
        tH1stFq=0d0
        tH1stFqb=0d0
        tH1stFg=0d0
ch        tgaga=0d0
ch        tcga=0d0
ch        tgamma2=0d0
        tgagaq=0d0
        tcgaq=0d0
        tgamma2q=0d0
        tgagag=0d0
        tcgag=0d0
        tgamma2g=0d0

ch        diffc10=0d0
ch        diffc1f=0d0
ch        diffc20=0d0
ch        diffc2f=0d0

ch        diffg10=0d0
ch        diffg1f=0d0
ch        diffg20=0d0
ch        diffg2f=0d0

        diffc10q=0d0
        diffc1fq=0d0
        diffc20q=0d0
        diffc2fq=0d0

        diffg10q=0d0
        diffg1fq=0d0
        diffg20q=0d0
        diffg2fq=0d0

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

        diff10g=0d0
        diff20g=0d0

        diff1fg=0d0
        diff1fgq=0d0
        diff1fgqb=0d0
        diff2fg=0d0
        diff2fgq=0d0
        diff2fgqb=0d0

        sig1=0d0
        sig2=0d0

        sig11q=0d0
        sig12q=0d0
        sig11g=0d0
        sig12g=0d0
        sig11=0d0
        sig12=0d0
        sig21=0d0      
        sig22=0d0
        sig23=0d0
        sig24=0d0

        sig21q=0d0      
        sig22q=0d0
        sig23q=0d0
        sig24q=0d0

        sig21g=0d0      
        sig22g=0d0
        sig23g=0d0
        sig24g=0d0

      
ch      write(*,*)z1,z2,xx10,xx20,facscale,beta,alfa,'1'
ch      write(*,*)msqc(1,-1),msqc(0,0),'2'
      do j=-nf,nf
      do k=-nf,nf


c      if(k.ne.0.or.j.ne.0) goto 75
ch
      if(msqc(j,k).eq.0d0) goto 75

C     Simplest term without convolutions
      
      if(j.gt.0.and.k.lt.0)then
       tdeltaq=tdeltaq+fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.lt.0.and.k.gt.0)then
       tdeltaqb=tdeltaqb+fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.eq.0.and.k.eq.0)then
       tdeltag=tdeltag+fx10(j)*fx20(k)*msqc(j,k)*flgq*flqqb
      endif


C     Start H1st: to be used later

C     H1st delta term

ch      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq

        if(j.gt.0.and.k.lt.0)then
         tH1stq=tH1stq+H1qdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
        elseif(j.lt.0.and.k.gt.0)then
         tH1stqb=tH1stqb+H1qbdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
         elseif(j.eq.0.and.k.eq.0)then
          tH1stg=tH1stg+H1gdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
     .                                                   *flqqb
        endif

ch        if(j.ne.0.and.k.ne.0)then
ch         tH1stq=tH1stq+H1qdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch         elseif(j.eq.0.and.k.eq.0)then
ch          tH1stg=tH1stg+H1gdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch     .                                                   *flqqb
ch        endif


ch        if(j.gt.0.and.k.lt.0)then
ch         tH1stq=tH1stq+2*C1qqdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch        elseif(j.lt.0.and.k.gt.0)then
ch         tH1stq=tH1stq+2*C1qqdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch         elseif(j.eq.0.and.k.eq.0)then
ch          tH1stg=tH1stg+2*C1ggdelta*msqc(j,k)*fx10(j)*fx20(k)*flgq
ch        endif
c       tH1st


      if(j.ne.0.and.k.ne.0)then

C     H1st: non delta terms, first leg
        if(j.gt.0.and.k.lt.0)then
       tH1stq=tH1stq+(fx1p(j)*Cqq(z1)*flgq
     &               +fx1p(0)*Cqg(z1)*flqqb)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
        elseif(j.lt.0.and.k.gt.0)then
       tH1stqb=tH1stqb+(fx1p(j)*Cqq(z1)*flgq
     &               +fx1p(0)*Cqg(z1)*flqqb)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
        endif


C     H1st: non delta terms, second leg

        if(j.gt.0.and.k.lt.0)then
       tH1stq=tH1stq+(fx2p(k)*Cqq(z2)*flgq
     &                +fx2p(0)*Cqg(z2)*flqqb)         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
        elseif(j.lt.0.and.k.gt.0)then
       tH1stqb=tH1stqb+(fx2p(k)*Cqq(z2)*flgq
     &                +fx2p(0)*Cqg(z2)*flqqb)         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
        endif

      endif
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg

ch    gammagg: first leg is included      


c      diff=-dlog(xx10)
c     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
ch
      if(j.ne.0.and.k.ne.0)then
       diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq
     &     +fx1p(0)*Pqg(z1)*flqqb)
      elseif(j.eq.0.and.k.eq.0)then
       diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*3/(1-z1)
     &    +fx1p(0)*Pggreg(z1))*flgq*flqqb
      endif

ch
ch      if(j.ne.0.and.k.ne.0)then
       if(j.gt.0.and.k.lt.0)then
       tH1stFq=tH1stFq+diff*fx20(k)*msqc(j,k)
       tH1stFq=tH1stFq-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.lt.0.and.k.gt.0)then
       tH1stFqb=tH1stFqb+diff*fx20(k)*msqc(j,k)
       tH1stFqb=tH1stFqb-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.eq.0.and.k.eq.0)then
       tH1stFg=tH1stFg+diff*fx20(0)*msqc(j,k)
       tH1stFg=tH1stFg-3*D0int(xx10)*fx10(0)*fx20(0)*msqc(j,k)*flgq
     &                                                       *flqqb
      endif

ch

c     gammaqq and gammaqg: second leg   

ch    gammagg: second leg is included


c      diff=-dlog(xx20)
c     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
ch
      if(j.ne.0.and.k.ne.0)then
       diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq
     &     +fx2p(0)*Pqg(z2)*flqqb)
      elseif(j.eq.0.and.k.eq.0)then
       diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*3/(1-z2)
     &    +fx2p(0)*Pggreg(z2))*flgq*flqqb
      endif

ch
ch      if(j.ne.0.and.k.ne.0)then
       if(j.gt.0.and.k.lt.0)then
       tH1stFq=tH1stFq+diff*fx10(j)*msqc(j,k)
       tH1stFq=tH1stFq-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.lt.0.and.k.gt.0)then
       tH1stFqb=tH1stFqb+diff*fx10(j)*msqc(j,k)
       tH1stFqb=tH1stFqb-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq
      elseif(j.eq.0.and.k.eq.0)then
       tH1stFg=tH1stFg+diff*fx10(0)*msqc(j,k)
       tH1stFg=tH1stFg-3*D0int(xx20)*fx10(0)*fx20(0)*msqc(j,k)*flgq
     &                                                       *flqqb
      endif

ch

c     gammagg: delta term, both legs

      if(j.eq.0.and.k.eq.0)then
       tH1stFg=tH1stFg+2*beta0*fx10(0)*fx20(0)*msqc(j,k)*flgq
     &                                                 *flqqb
      endif


CC    End of H1st
ch

      if(order.eq.1) goto 75

ch    qqb channel:

      if(j.ne.0.and.k.ne.0)then

CC    Now (gamma+gamma)*(gamma+gamma) term: to be used later

C     First part: one gamma for each leg: FLGQ here is non trivial ! DONE




      diffg1fq=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)
     &  - Pqqint(xx10)*fx10(j)


      diffg10q=-dlog(xx10)*fx1p(0)*Pqg(z1)*flqqb

      diffg2fq=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)
     &  - Pqqint(xx20)*fx20(k)


      diffg20q=-dlog(xx20)*fx2p(0)*Pqg(z2)*flqqb


      tgagaq=tgagaq+2*
     #   (flgq*diffg10q*diffg20q+flgq*diffg1fq*diffg2fq
ch     #   (flgq*diffg10q*diffg20q+diffg1fq*diffg2fq
     #   +diffg10q*diffg2fq+diffg1fq*diffg20q)*msqc(j,k)*flqq



CC     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

C     First leg

      
      diff1q=(-dlog(xx10)*(flgq*(fx1p(j)-fx10(j)*xx10**beta)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq
     &    +fx1p(0)*(Pqqqg(z1)+Pqggg(z1))*flqqb)
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq)*flqq


C    Second leg

      
      diff2q=(-dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq
     &    +fx2p(0)*(Pqqqg(z2)+Pqggg(z2))*flqqb)
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq)*flqq


C     Include Pqggq

      do l=1,nf
      if(l.eq.j)then
       diff1q=diff1q-dlog(xx10)*(fx1p(l)*flqq
     &                          +fx1p(-l)*flqqb)*Pqggq(z1)*flgq
       diff2q=diff2q-dlog(xx20)*(fx2p(l)*flqqb
     &                           +flqq*fx2p(-l))*Pqggq(z2)*flgq
      elseif(l.eq.-j)then
       diff1q=diff1q-dlog(xx10)*(fx1p(l)*flqqb
     &                           +flqq*fx1p(-l))*Pqggq(z1)*flgq
       diff2q=diff2q-dlog(xx20)*(fx2p(l)*flqq
     &                          +fx2p(-l)*flqqb)*Pqggq(z2)*flgq
      else
       diff1q=diff1q-dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
     &                                            *flqqb
       diff2q=diff2q-dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
     &                                            *flqqb
      endif
      enddo

      tgagaq=tgagaq+diff1q*fx20(k)*msqc(j,k)
      tgagaq=tgagaq+diff2q*fx10(j)*msqc(j,k)




C    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE

C    Start  (C+C)*(gamma+gamma) term

c    gamma first leg, C second leg


ch      diffc2fq=-dlog(xx20)*fx2p(k)*Cqq(z2)+C1qqdelta*fx20(k)
      diffc2fq=-dlog(xx20)*fx2p(k)*Cqq(z2)
ch     #          0.5d0*(H1qdelta+H1qbdelta)/2d0*fx20(k)
ch     #          +0.5d0*H1qdelta*fx20(k)

      diffc20q=-dlog(xx20)*fx2p(0)*Cqg(z2)*flqqb


      tcgaq=tcgaq+flqq*msqc(j,k)*
     # (flgq*diffg10q*diffc20q+flgq*diffg1fq*diffc2fq
     #          +diffg10q*diffc2fq+diffg1fq*diffc20q)



c    C first leg, gamma second leg

ch      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)
      diffc1fq=-dlog(xx10)*fx1p(j)*Cqq(z1)
ch     #          0.5d0*(H1qdelta+H1qbdelta)/2d0*fx10(j)
ch     #          +0.5d0*H1qdelta*fx10(j)

      diffc10q=-dlog(xx10)*fx1p(0)*Cqg(z1)*flqqb

      tcgaq=tcgaq+flqq*msqc(j,k)*
     # (flgq*diffc10q*diffg20q+flgq*diffc1fq*diffg2fq
     #          +diffc10q*diffg2fq+diffc1fq*diffg20q)

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx1p(j)*CqqPqq(z1)*flgq
     &       +fx1p(0)*(CqqPqg(z1)+CqgPgg(z1))*flqqb)
     &     *(-dlog(xx10))*fx20(k)*msqc(j,k)*flqq 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcgaq=tcgaq
     &     +(fx2p(k)*CqqPqq(z2)*flgq
     &       +fx2p(0)*(CqqPqg(z2)+CqgPgg(z2))*flqqb)
     &     *(-dlog(xx20))*fx10(j)*msqc(j,k)*flqq 

c    Add Cqg*Pgq contribution

      do l=1,nf
      if(l.eq.j)then
       tcgaq=tcgaq+(fx1p(l)*flqq
     &              +fx1p(-l)*flqqb)*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
       tcgaq=tcgaq+(fx2p(l)*flqqb
     &              +flqq*fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq 
      elseif(l.eq.-j)then
       tcgaq=tcgaq+(fx1p(l)*flqqb
     &              +flqq*fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
       tcgaq=tcgaq+(fx2p(l)*flqq
     &              +fx2p(-l)*flqqb)*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      else
       tcgaq=tcgaq+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
     &                                          *flqqb 
       tcgaq=tcgaq+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
     &                                          *flqqb  
      endif
      enddo


CC  Start 2-loop AP

C   Gluon + pure singlet


      do l=-nf,nf
      if(l.eq.0) then
      tgamma2q=tgamma2q+fx1p(0)*P2qg(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flqq
     &                                *flqqb
      tgamma2q=tgamma2q+fx2p(0)*P2qg(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flqq
     &                                *flqqb
      elseif(l.eq.j)then
      tgamma2q=tgamma2q+flqq*fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+flqqb*fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      elseif(l.eq.-j)then
      tgamma2q=tgamma2q+flqqb*fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+flqq*fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      else
      tgamma2q=tgamma2q+flqqb*fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2q=tgamma2q+flqqb*fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part

      tgamma2q=tgamma2q+fx1p(j)*P2qqV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq*flqq
      tgamma2q=tgamma2q+fx2p(k)*P2qqV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq*flqq


C   P2qq non-singlet: 1/(1-z)_+


      diffq=-dlog(xx10)
     &  *(fx1p(j)-fx10(j)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diffq*fx20(k)*msqc(j,k)*flgq
     &                                                     *flqq


      diffq=-dlog(xx20)
     &  *(fx2p(k)-fx20(k)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2q=tgamma2q+2d0/3*Kappa*diffq*fx10(j)*msqc(j,k)*flgq
     &                                                     *flqq

      

C   P2qqb non singlet

      tgamma2q=tgamma2q+fx1p(-j)*P2qqbV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq*flqqb

      tgamma2q=tgamma2q+fx2p(-k)*P2qqbV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq*flqqb

      elseif(j.eq.0.and.k.eq.0)then

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

      diffg=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)
     &    *(D0gggg/(1-z1)+D1gggg*dlog(1-z1)/(1-z1))
     &      +fx1p(0)*(Pggggreg(z1)+Pgqqg(z1)))
     &    +(Deltagggg-D0gggg*D0int(xx10)-D1gggg*D1int(xx10))*fx10(0)


      tgagag=tgagag+diffg*flqq*flgq*flqqb*fx20(0)*msqc(0,0)


c    Second leg

      diffg=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)
     &    *(D0gggg/(1-z2)+D1gggg*dlog(1-z2)/(1-z2))
     &      +fx2p(0)*(Pggggreg(z2)+Pgqqg(z2)))
     &    +(Deltagggg-D0gggg*D0int(xx20)-D1gggg*D1int(xx20))*fx20(0)

      tgagag=tgagag+diffg*flqq*flgq*flqqb*fx10(0)*msqc(0,0)



CC    Start  (C+C)*(gamma+gamma) term: diagonal part


c    gamma first leg, C second

      diffg10g=-dlog(xx10)
     &  *((fx1p(0)-fx10(0)*xx10**beta)*3d0/(1-z1)+Pggreg(z1)*fx1p(0))
     &    +fx10(0)*(beta0-3*D0int(xx10))

ch      diffc20g=0.5d0*H1gdelta*fx20(0)
      diffc20g=0d0
ch      diffc20g=(ca*pi**2/12d0)*fx20(0)
ch      diffc20g=0d0

c    gamma second leg, C first

      diffg20g=-dlog(xx20)
     &  *((fx2p(0)-fx20(0)*xx20**alfa)*3d0/(1-z2)+Pggreg(z2)*fx2p(0))
     &    +fx20(0)*(beta0-3*D0int(xx20))


ch      diffc10g=0.5d0*H1gdelta*fx10(0)
      diffc10g=0d0
ch      diffc10g=(ca*pi**2/12d0)*fx10(0)
ch       diffc10g=0d0

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)

      tcgag=tcgag+CgqPqg(z1)*(-dlog(xx10))*flqqb*
     &                 flqq*flgq*fx1p(0)*fx20(0)*msqc(0,0) 

c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)

      tcgag=tcgag+CgqPqg(z2)*(-dlog(xx20))*flqqb*
     &                 flqq*flgq*fx2p(0)*fx10(0)*msqc(0,0) 

c    End of (C+C)*(gamma+gamma)


CC    gamma2: diagonal part

c     First leg

      diffg=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diffg-dlog(xx10)*P2gg(z1)*fx1p(0))
     &                *flqq*flgq*flqqb*fx20(0)*msqc(0,0)


c     Second leg

      diffg=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)  

      tgamma2g=tgamma2g+(1.5d0*Kappa*diffg-dlog(xx20)*P2gg(z2)*fx2p(0))
     &                *flqq*flgq*flqqb*fx10(0)*msqc(0,0)

      endif

 75   continue

      enddo
      enddo

ch    do here gq 

      do j=1,nf

C     H1st: Cgq, first leg


      tH1stg=tH1stg+(fx1p(j)+fx1p(-j))*Cgq(z1)
     & *(-dlog(xx10))*fx20(0)*flqqb*msqc(0,0)


C     H1st: Cgq, second leg

      
      tH1stg=tH1stg+(fx2p(j)+fx2p(-j))*Cgq(z2)
     & *(-dlog(xx20))*fx10(0)*flqqb*msqc(0,0)

C     gg channel H1st: muf dependence: Pgq, first leg

      tH1stFg=tH1stFg+(-dlog(xx10))
     & *(fx1p(j)+fx1p(-j))*Pgq(z1)*fx20(0)*flqqb*msqc(0,0)


C     gg channel H1st: muf dependence: Pgq, second leg


      tH1stFg=tH1stFg+(-dlog(xx20))
     & *(fx2p(j)+fx2p(-j))*Pgq(z2)*fx10(0)*flqqb*msqc(0,0)

CC    End of H1st

      if(order.eq.1) goto 78

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
     &            *fx20(0)*msqc(0,0)*flqq*flqqb

      tgagag=tgagag-dlog(xx20)*(Pgqqq(z2)+Pgggq(z2))*(fx2p(j)+fx2p(-j))
     &            *fx10(0)*msqc(0,0)*flqq*flqqb


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


      tcgag=tcgag+CgqPqq(z1)*(-dlog(xx10))*(fx1p(j)+fx1p(-j))
     .                *fx20(0)*msqc(0,0)*flqq*flqqb 


c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)


      tcgag=tcgag+CgqPqq(z2)*(-dlog(xx20))*(fx2p(j)+fx2p(-j))
     .                   *fx10(0)*msqc(0,0)*flqq*flqqb 



CC    gamma2: qg channel


c    First leg

      tgamma2g=tgamma2g
     &   -dlog(xx10)*P2gq(z1)*(fx1p(j)+fx1p(-j))*fx20(0)*msqc(0,0)
     &                       *flqq*flqqb 

c    Second leg

      tgamma2g=tgamma2g
     &   -dlog(xx20)*P2gq(z2)*(fx2p(j)+fx2p(-j))*fx10(0)*msqc(0,0)
     &                       *flqq*flqqb 


 78   continue


      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    Check it !      

      tgagag=tgagag+2*msqc(0,0)
     # *(flgq*diff10g*diff20g
     #   +diff10g*diff2fg+diff1fg*diff20g)*flqq*flqqb

c

c    gamma first leg, C second leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg10g*diffc20g
     #  +diffg10g*diffc2fg
     #  +diffg1fg*diffc20g)*flqq*flqqb

c    gamma second leg, C first leg

      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg20g*diffc10g
     #  +diffg20g*diffc1fg
     #  +diffg2fg*diffc10g)*flqq*flqqb

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
     # *(flgq*diff1fgq*diff2fgq*flqqb
     #   +flgq*diff1fgqb*diff2fgqb*flqqb
     #   +flgq*flqq*diff1fgq*diff2fgqb
     #   +flgq*flqq*diff1fgqb*diff2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg1fgq*diffc2fgq*flqqb
     #  +flgq*diffg1fgqb*diffc2fgqb*flqqb
     #  +flgq*flqq*diffg1fgq*diffc2fgqb
     #  +flgq*flqq*diffg1fgqb*diffc2fgq)
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg2fgq*diffc1fgq*flqqb
     #  +flgq*diffg2fgqb*diffc1fgqb*flqqb
     #  +flgq*flqq*diffg2fgq*diffc1fgqb
     #  +flgq*flqq*diffg2fgqb*diffc1fgq)
      else
       tgagag=tgagag+2*msqc(0,0)
     # *(flgq*diff1fgq*diff2fgq
     #   +flgq*diff1fgqb*diff2fgqb
     #   +flgq*diff1fgq*diff2fgqb
     #   +flgq*diff1fgqb*diff2fgq)*flqqb
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg1fgq*diffc2fgq
     #  +flgq*diffg1fgqb*diffc2fgqb
     #  +flgq*diffg1fgq*diffc2fgqb
     #  +flgq*diffg1fgqb*diffc2fgq)*flqqb
      tcgag=tcgag+msqc(0,0)*
     # (flgq*diffg2fgq*diffc1fgq
     #  +flgq*diffg2fgqb*diffc1fgqb
     #  +flgq*diffg2fgq*diffc1fgqb
     #  +flgq*diffg2fgqb*diffc1fgq)*flqqb
      endif
      enddo
      enddo
ch      tgagag=0d0
ch      tcgag=0d0


CC   First order

ch Include also the gg channel

c      sig12=-0.5d0*A1q*tdelta
c      sig11=-B1q*tdelta-tH1stF
      sig12q=-0.5d0*A1q*(tdeltaq+tdeltaqb)*flqq
      sig11q=(-B1q*(tdeltaq+tdeltaqb)
     .       -tH1stFq-tH1stFqb
ch     .   )*flqq
     .       +Gamma1q*tdeltaq+Gamma1qb*tdeltaqb)*flqq
ch     .       +Gamma1q*(tdeltaq+tdeltaqb))*flqq
ch      write(*,*)ason2pi*2d0*((-B1q+Gamma1q)*msqc(1,-1)),'sig11q'
ch      write(*,*)ason2pi*2d0*(Gamma1q*msqc(1,-1)),'sig11q',
ch     .    LL1*ason2pi*2d0*((Gamma1q)*msqc(1,-1))+
ch     .    LL2*ason2pi*2d0*(-0.5d0*A1q*msqc(1,-1))
ch      write(*,*)ason2pi*2d0*(-0.5d0*A1q*msqc(1,-1)),'sig12q',
ch     .    ason2pi*2d0*sig12q,z2


      sig12g=-0.5d0*A1g*tdeltag*flqq
      sig11g=(-B1g*tdeltag-tH1stFg
     .       +Gamma1g*tdeltag)*flqq

      sig12=sig12q+sig12g
      sig11=sig11q+sig11g
ch      sig11=sig11g
ch      write(*,*)sig12q,sig11q,'3'
ch      write(*,*)sig12g,sig11g,'4'

ch

CC   Second order

      sig24q=(A1q)**2/8*(tdeltaq+tdeltaqb)*flqq

      sig24g=(A1g)**2/8*tdeltag*flqq

      sig24=sig24q+sig24g
c      sig24=sig24q
       
      sig23q=-beta0*A1q/3*(tdeltaq+tdeltaqb)*flqq
     .           -0.5d0*A1q*sig11q
      
      sig23g=-beta0*A1g/3*tdeltag*flqq
     .           -0.5d0*A1g*sig11g

c      sig23=sig23q
      sig23=sig23q+sig23g

      sig22q=0.5d0*(beta0*A1q*LR-A2q)*(tdeltaq+tdeltaqb)
     &                               *flqq
     &     -0.5d0*A1q*(tH1stq+tH1stqb
     &                 +LF*tH1stFq+LF*tH1stFqb)*flqq
     &     -0.5d0*(B1q-beta0)*sig11q
     &     +0.5d0*B1q*(tH1stFq+tH1stFqb)*flqq
     &     -(Gamma1q*tH1stFq+Gamma1qb*tH1stFqb)*flqq
     &     -0.5d0*B1q*(Gamma1q*tdeltaq+Gamma1qb*tdeltaqb)*flqq
ch     &     -0.5d0*B1q*Gamma1q*(tdeltaq+tdeltaqb)*flqq
     &     +0.5d0*tgagaq
     &     +0.5d0*Gamma1sq_q*tdeltaq*flqq
     &     +0.5d0*Gamma1sq_qb*tdeltaqb*flqq

      sig22g=0.5d0*(beta0*A1g*LR-A2g)*tdeltag*flqq
     &     -0.5d0*A1g*(tH1stg+LF*tH1stFg)*flqq
     &     -0.5d0*(B1g-beta0)*sig11g
     &     +0.5d0*B1g*tH1stFg*flqq
     &     -Gamma1g*tH1stFg*flqq
     &     -0.5d0*B1g*Gamma1g*tdeltag*flqq
     &     +0.5d0*tgagag
     &     +0.5d0*Gamma1sq_g*tdeltag*flqq
C    Add mur dependence from H1st
ch     &     +beta0*A1g*LR*tdeltag

      sig22=sig22q+sig22g
c      sig22=sig22q



      sig21q=-beta0*LR*sig11q
     &       -B1q*(tH1stq+tH1stqb
     &             +LF*tH1stFq+LF*tH1stFqb)
     &                           *flqq
     &     +Gamma1q*(tH1stq
     &               -H1qdelta*tdeltaq
ch     &               -H1qdelta*(tdeltaq
ch     &               +tdeltaqb)
     &                +LF*tH1stFq)*flqq
     &     +Gamma1qb*(tH1stqb
     &               -H1qbdelta*tdeltaqb
ch     &               -H1qdelta*(tdeltaq
ch     &               +tdeltaqb)
     &                +LF*tH1stFqb)*flqq
     &     -LF*tgagaq
     &     -B2q*(tdeltaq+tdeltaqb)*flqq
     &     +beta0*(tH1stq+tH1stqb
     &               -H1qdelta*tdeltaq
     &               -H1qbdelta*tdeltaqb
ch     &               -H1qdelta*(tdeltaq
ch     &               +tdeltaqb)
     &       )*flqq
     &     -tcgaq
     &     -tgamma2q
     &     +(Gamma2q*tdeltaq+Gamma2qb*tdeltaqb)*flqq
     &     +(Gamma1H1q*tdeltaq+
ch     &       )
     &       Gamma1H1qb*
     &       tdeltaqb)*flqq
c     Include missing delta term from C*gamma (no factor 2 here !)
ch      sig21q=sig21q-2*C1qqdelta*tH1stFq
      sig21q=sig21q-H1qdelta*tH1stFq*flqq
     &             -H1qbdelta*tH1stFqb*flqq

C     Include missing term from contact term in 2 loop AP

      sig21q=sig21q-2*Delta2qq*(tdeltaq+tdeltaqb)*flqq
ch
ch      sig21q=sig21q+2*beta0*LR*(B1q*tdeltaq+tH1stFq)
ch
      sig21g=-beta0*LR*sig11g-B1g*(tH1stg+LF*tH1stFg)
     &                           *flqq
     &     +Gamma1g*(tH1stg-H1gdelta*tdeltag+LF*tH1stFg)*flqq
     &     -LF*tgagag
     &     -B2g*tdeltag*flqq
     &     +beta0*(tH1stg-H1gdelta*tdeltag)*flqq
     &     -tcgag
     &     -tgamma2g
     &     +Gamma2g*tdeltag*flqq
     &     +Gamma1H1g*tdeltag*flqq
c     Include missing delta term from C*gamma (no factor 2 here !)
ch      sig21g=sig21g-2*C1ggdelta*tH1stFg
      sig21g=sig21g-H1gdelta*tH1stFg*flqq
ch      sig21g=sig21g-H1gdelta
C     Include missing term from contact term in 2 loop AP
      sig21g=sig21g-2*Delta2gg*tdeltag*flqq
C     NEW: Include mur dependence from H1st
ch      sig21g=sig21g+2*beta0*LR*(B1g*tdeltag+tH1stFg)

      sig21=sig21q+sig21g
ch      sig21=sig21q
ch      write(*,*)LR,LF,'logs'
ch      write(*,*)sig24q,sig23q,sig22q,sig21q,'5'
ch      write(*,*)sig24g,sig23g,sig22g,sig21g,'6'


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ch      write(*,*)sig24,sig23,sig22,sig21,'2'
CC Include as/pi factors and sum O(as) and O(as^2) contributions

      sig1=sig12*LL2+sig11*LL1
ch      sig1=sig11*LL1
      sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1
ch      sig2=sig24g*LL4+sig23g*LL3+sig22g*LL2+sig21g*LL1
ch      sig2=sig24*LL4+sig23*LL3+sig22*LL2




        sig1=sig1*ason2pi*2
        sig2=sig2*(ason2pi*2)**2

        if(order.eq.1)then
c          xmsq(1)=-sig1
          xmsq(1)=-sig1
c          xmsq(1)=-tdelta
        else
ch          xmsq(1)=-(sig1+sig2)
          xmsq(1)=-sig2
        endif

ch        write(*,*)'sig2=',sig2
ch        write(*,*)'sig1=',sig1,sig11,sig12,tH1stFg
ch        write(*,*)'sig24q=',(ason2pi*2)**2*A1q**2/8d0*msqc(1,-1)
ch        write(*,*)'sig23q=',(ason2pi*2)**2*msqc(1,-1)*
ch     .    (-beta0*A1q/3-0.5d0*A1q*Gamma1q)
ch        write(*,*)'sig22q=',(ason2pi*2)**2*msqc(1,-1)*(
ch     &    0.5d0*(
ch     &           beta0*A1q*LR
ch     &           -A2q)
ch     &     -0.5d0*A1q*H1qdelta
ch     &     -0.5d0*(-beta0)*Gamma1q
ch     &     +0.5d0*Gamma1sq_q)
ch      write(*,*)'sig21q=',(ason2pi*2)**2*msqc(1,-1)*(
ch     &          -beta0*LR*sig11q
ch     &     +Gamma2q
ch     &     +Gamma1H1q
ch     &     -B2q)
ch      write(*,*)
ch     &  -(ason2pi*2)**2*msqc(1,-1)*beta0*LR*sig11q,'sig11LR',
ch     &  (ason2pi*2)**2*msqc(1,-1)*Gamma2q,'Gamma2',
ch     &  (ason2pi*2)**2*msqc(1,-1)*Gamma1H1q,'Gamma1H1',
ch     &  -(ason2pi*2)**2*msqc(1,-1)*B2q,'B2'

ch      write(*,*)(ason2pi*2)**2*A1q**2/8d0*msqc(1,-1)*LL4+
ch     &  LL3*(ason2pi*2)**2*msqc(1,-1)*
ch     .    (-beta0*A1q/3-0.5d0*A1q*Gamma1q)+
ch     &     LL2*(ason2pi*2)**2*msqc(1,-1)*(
ch     &    0.5d0*(
ch     &           beta0*A1q*LR
ch     &           -A2q)
ch     &     -0.5d0*A1q*H1qdelta
ch     &     -0.5d0*(-beta0)*Gamma1q
ch     &     +0.5d0*Gamma1sq_q)+
ch     &  LL1*(ason2pi*2)**2*msqc(1,-1)*(
ch     &          -beta0*LR*Gamma1q
ch     &     +Gamma2q
ch     &     +Gamma1H1q
ch     &     +2*beta0*LR*Gamma1q
ch     &     -B2q)
ch      write(*,*)(ason2pi*2)**2*A1q**2/8d0*msqc(1,-1)*LL4,'ll4',
ch     &  LL3*(ason2pi*2)**2*msqc(1,-1)*
ch    .    (-beta0*A1q/3-0.5d0*A1q*Gamma1q),'ll3',
ch     &     LL2*(ason2pi*2)**2*msqc(1,-1)*(
ch     &    0.5d0*(
ch     &           beta0*A1q*LR
ch     &           -A2q
ch     &        )
ch     &     -0.5d0*A1q*H1qdelta
ch     &     -0.5d0*(-beta0)*Gamma1q
ch     &     +0.5d0*Gamma1sq_q
ch     &      ),'ll2',
ch     &  LL1*(ason2pi*2)**2*msqc(1,-1)*(
ch     &          -beta0*LR*Gamma1q
ch     &     +Gamma2q
ch     &     +Gamma1H1q
ch     &     +2*beta0*LR*Gamma1q
ch     &     -B2q
ch     &    ),'ll1'


ch     .    ,'ason2pi=',ason2pi,
ch     .  'xx10,xx20=',xx10,xx20,
ch     .  'alfa,beta=',alfa,beta,
ch     .  'LL1,LL2,LL3,LL4=',LL1,LL2,LL3,LL4,scale,
ch     .  facscale,
ch     .  'Lr,LF=',LR,LF,
ch     .  'msqb,msqc=',msqc(0,0),msqc(1,-1),
ch     .  'q2,xth=',q2,xth

c           xmsq(1)=sig1+sig2

CC Include iacobians

c      xmsq(1)=xmsq(1)*xjacqt2*xjacq2*q2/shad/Vol

      xmsq(1)=xmsq(1)*xjacq2*xjacqt2*q2*pswt0/shad/Vol
ch     .      *dsqrt(qt2)/dsqrt(q2)


      countint=0d0
      xint=0d0


C Multiply by BORN phase space weight

c        xmsq(1)=xmsq(1)*fluxborn*pswt0/BrnRat

        xmsq(1)=fluxBorn*xmsq(1)/BrnRat


 77    continue



c---Add to total

        xint=xmsq(1)
        val=xmsq(1)*wgt
        

c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(1,pjet)
          call dotem(nvec,pjet,s)
          val=val/dfloat(itmx)       
ch          call plotter(ptrans,val,1)       
          call plotter(p,val,0)
        endif

     
      countint=xint
     
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  countint=0d0
      ntotzero=ntotzero+1
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC qq splitting function (with asopi normalization)

      function Pqq(z)
      implicit none
      real *8 Pqq,z
      Pqq=2d0/3*(1+z**2)/(1-z)
      return
      end

CC qg splitting function (with asopi normalization)

      function Pqg(z)
      implicit none
      real *8 Pqg,z
      Pqg=0.25d0*(1-2*z*(1-z))
      return
      end

CC gq splitting function (with asopi normalization)

      function Pgq(z)
      implicit none
      real *8 Pgq,z
      Pgq=2d0/3*(1+(1-z)**2)/z
      return
      end

CC Non delta term in Cqq coefficient (with asopi normalization)

      function Cqq(z)
      implicit none
      real *8 Cqq,z
      Cqq=2d0/3*(1-z)
      return
      end

CC gg splitting function (not normalized !)



CC Cqg coefficient (with asopi normalization)

      function Cqg(z)
      implicit none
      real *8 Cqg,z
      Cqg=0.5d0*z*(1-z)
      return
      end

CC Cgq coefficient (with asopi normalization)

      function Cgq(z)
      implicit none
      real *8 Cgq,z
      Cgq=2d0/3*z
      return
      end


CC Integral of Pqq=1/2 CF (1+x^2)/(1-x) from 0 to z

      function Pqqint(z)
      implicit none
      real *8 Pqqint,z
      Pqqint=-2d0/3*(z+z**2/2+2*dlog(1-z))
      return
      end

CC Integral of 1/(1-x) from 0 to z

      function D0int(z)
      implicit none
      real *8 D0int,z
      D0int=-dlog(1-z)
      return
      end

CC Integral of log(1-x)/(1-x) from 0 to z

      function D1int(z)
      implicit none
      real *8 D1int,z
      D1int=-0.5d0*dlog(1-z)**2
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                P*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Regular part of Pqq*Pqq (checked !)

      function Pqqqq(z)
      implicit none
      real *8 Pqqqq,z
      Pqqqq=4d0/9*(-4*dlog(z)/(1-z)-2*(1-z)
     &  +(1+z)*(3*dlog(z)-4*dlog(1-z)-3))
      return
      end


CC Pqq*Pqg (checked !)

      function Pqqqg(z) 
      implicit none
      real *8 Pqqqg,z
      Pqqqg=1d0/3*((z**2+(1-z)**2)*dlog((1-z)/z)
     &  -(z-0.5d0)*dlog(z)+z-0.25d0)
      return
      end

CC Pqg*Pgq (checked !)

      function Pqggq(z)
      implicit none
      real *8 Pqggq,z
      Pqggq=1d0/3*(2d0/3/z+(1+z)*dlog(z)-2d0/3*z**2-0.5d0*(z-1))
      return
      end


CC Full Pqg*Pgg (checked !)

      function Pqggg(z)
      implicit none
      real *8 Pqggg,z,beta0,Pqg
      integer nf
      external Pqg
      nf=5
      beta0=(33-2*nf)/12d0
      Pqggg=1.5d0*(1/3d0/z+(z**2-z+0.5d0)*dlog(1-z)
     &     +(2*z+0.5d0)*dlog(z)+0.25d0+2*z-31d0/12*z**2)

      Pqggg=Pqggg+beta0*Pqg(z)
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                C*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Cqq*Pqq (without delta term in Cqq) (checked !)

      function CqqPqq(z)
      implicit none
      real *8 CqqPqq,z
      CqqPqq=2d0/9*(1-z)*(4*dlog(1-z)-2*dlog(z)-1)
      return
      end

CC Cqq*Pqg (without delta term in Cqq) (checked !)

      function CqqPqg(z)
      implicit none
      real *8 CqqPqg,z
      CqqPqg=(-2+z+z**2-(1+2*z)*dlog(z))/6d0
      return
      end

CC Cqg*Pgq (checked !)

      function CqgPgq(z) 
      implicit none
      real *8 CqgPgq,z
      CqgPgq=(1d0/3/z-1+2*z**2/3-z*dlog(z))/3d0
      return
      end

CC Cqg*Pgg (checked !)

      function CqgPgg(z)
      implicit none
      real *8 CqgPgg,z,beta0
      integer nf
      nf=5
      beta0=(33-2*nf)/12d0
      CqgPgg=3d0/4*(2*z*(1-z)*dlog(1-z)-4*z*dlog(z)
     &      +1d0/3/z-1-5*z+17d0*z**2/3)+beta0/2*z*(1-z)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           Two loop AP:  pqq of ESW is my 3/2 Pqq
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Pqq NS: Eq. (4.107) ESW (no 1/(1-x)_+ and delta term)

      function P2qqV(x)
      implicit none
      real *8 x,P2qqV,Pqq,pi
      integer nf
      external Pqq

      pi=3.14159265358979d0
      nf=5

      P2qqV=16d0/9*(-(2*dlog(x)*dlog(1-x)+1.5d0*dlog(x))*3d0/2*Pqq(x)
     &     -(1.5d0+3.5d0*x)*dlog(x)-0.5d0*(1+x)*dlog(x)**2-5*(1-x))
     &     +4*((0.5d0*dlog(x)**2+11d0/6*dlog(x))*3d0/2*Pqq(x)
     &     -(67d0/18-pi**2/6)*(1+x)
     &     +(1+x)*dlog(x)+20d0/3*(1-x))
     &     +2d0/3d0*nf*(-dlog(x)*Pqq(x)+10d0/9*(1+x)-4d0/3*(1-x))

c     Change to as/pi normalization

      P2qqV=P2qqV/4  

      return
      end


C    Pqqb NS: Eq. (4.108) ESW

      function P2qqbV(x)
      implicit none
      real *8 x,P2qqbV,Pqq,S2
      external Pqq,S2

      P2qqbV=-2d0/9*(3d0*Pqq(-x)*S2(x)+2*(1+x)*dlog(x)+4*(1-x))
      
c     Change to as/pi normalization

      P2qqbV=P2qqbV/4 

      return
      end



C    Pqg Singlet: Eq. (4.110) ESW (ESW Pqg is 4 times my Pqg)

      function P2qg(x)
      implicit none
      real *8 x,P2qg,Pqg,pi,S2,logx,logomxsx
      external Pqg,S2

      pi=3.14159265358979d0
      logx=dlog(x)
      logomxsx=dlog((1-x)/x)

      P2qg=2d0/3*(4-9*x-(1-4*x)*logx-(1-2*x)*logx**2+4*dlog(1-x)
     &    +(2*logomxsx**2-4*logomxsx-2d0/3*pi**2+10d0)*4*Pqg(x))
     &    +1.5d0*(182d0/9+14d0/9*x+40d0/9/x+(136d0/3*x-38d0/3)*logx
     &    -4*dlog(1-x)-(2+8*x)*logx**2+8*Pqg(-x)*S2(x)
     &    +(-logx**2+44d0/3*logx-2*dlog(1-x)**2+4*dlog(1-x)+pi**2/3
     &    -218d0/9)*4*Pqg(x))

c     Change to as/pi normalization

      P2qg=P2qg/4d0
  
c     Divide by 2 to eliminate 2nf factor

      P2qg=P2qg/2d0

      return
      end

C     Pqq Pure Singlet appearing in ESW Eq. (4.95)
C     PSqq=PSqqb
C     Obtained through Eq.(4.101)
C     PSqq=1/2/nf (P2qq-P2qqbV-P2qqV) (contains only CF TR=2/3)
CC gg splitting function (not normalized !)

      function Pgg(x)
      implicit none
      real *8 Pgg,x

      Pgg=1d0/(1-x)+1d0/x-2+x*(1-x)

      return
      end


CC gg splitting function: regular part (with asopi normalization)


      function Pggreg(z)
      implicit none
      real *8 Pggreg,z
      Pggreg=3*((1-2*z)/z+z*(1-z))
      return
      end
      function P2qqS(x)
      implicit none
      real *8 P2qqS,x

      P2qqS=2d0/3*(20 - 18*x + 54*x**2 - 56*x**3
     &    +3*x*(3 + 15*x + 8*x**2)*dlog(x) 
     &    - 9*x*(1 + x)*dlog(x)**2)/(9*x)
      
      P2qqS=P2qqS/4

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                P*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function Pggggreg(x) ! checked !
      implicit none
      real *8 Pggggreg,x
      real *8 Pggreg,beta0
      integer nf
      external Pggreg

      nf=5
      beta0=(33-2*nf)/12d0

      Pggggreg=-9*dlog(x)/(1-x)+2*beta0*Pggreg(x)
     &   +9*(3*(1-x)+11d0/3/x*(x**3-1d0)+2d0/3*dlog(1-x)*Pggreg(x)
     &   +dlog(x)*(x**2-3*x-1d0/x))

      return
      end


      function Pgggq(x) ! checked !
      implicit none
      real *8 Pgggq,x
      real *8 Pgq,beta0
      integer nf
      external Pgq

      nf=5
      beta0=(33-2*nf)/12d0

      Pgggq=2*((1+(1-x)**2)/x*dlog(1-x)-2*(1+x+1d0/x)*dlog(x)
     &    +4-31d0/6/x+x/2+2d0/3*x**2)+beta0*Pgq(x)

      return
      end

      function Pgqqq(x) ! checked !
      implicit none
      real *8 Pgqqq,x
      
      Pgqqq=4d0/9*((2-x)*dlog(x)+dlog(1-x)*(2*x+4/x-4)+2-x/2)
      
      return
      end

      function Pgqqg(x) ! checked !
      implicit none
      real *8 Pgqqg,x
      integer nf

      nf=5 
      Pgqqg=1d0/6*(1+4d0/3/x-x-4*x**2/3+2*(1+x)*dlog(x))

      Pgqqg=2*nf*Pgqqg

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                C*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      function CgqPqq(x) ! checked !
      implicit none
      real *8 CgqPqq,x
      
      CgqPqq=2d0/9*(2+x+4*x*dlog(1-x)-2*x*dlog(x))

      return
      end



      function CgqPqg(x) ! checked !
      implicit none
      real *8 CgqPqg,x
      integer nf

      nf=5       
      CgqPqg=1d0/6*(1+x-2*x**2+2*x*dlog(x))

      CgqPqg=2*nf*CgqPqg

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           Two loop AP:  pqq of ESW is my 3/2 Pqq
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      function P2gg(x)
      implicit none
      real *8 P2gg,S2,pgg,x,pi
      integer nf
      external pgg,S2

      nf=5

      pi=3.14159265358979d0

      P2gg=2d0/3*nf*(-16+8*x+20d0/3*x**2+4d0/3/x-(6+10*x)*dlog(x)
     &    -(2+2*x)*dlog(x)**2)+1.5d0*nf*(2-2*x+26d0/9*(x**2-1d0/x)
     &    -4d0/3*(1+x)*dlog(x)-20d0/9*(1/x-2+x*(1-x)))
     &    +9*(27d0/2*(1-x)+67d0/9*(x**2-1d0/x)
     &        -(25d0/3-11d0/3*x+44d0/3*x**2)*dlog(x)
     &        +4*(1+x)*dlog(x)**2+2*pgg(-x)*S2(x)
     &    +(67d0/9-pi**2/3)*(1/x-2+x*(1-x))
     &     +(-4*dlog(x)*dlog(1-x)+dlog(x)**2)*pgg(x))
     

      P2gg=P2gg/4d0

      return
      end


      function P2gq(x)
      implicit none
      real *8 P2gq,Pgq,S2,x,logx,logomx,pi
      external Pgq,S2
      integer nf
     
      pi=3.14159265358979d0

      nf=5

      logx=dlog(x)
      logomx=dlog(1-x)

      P2gq=16d0/9*(-2.5d0-3.5d0*x+(2+3.5d0*x)*logx
     &     -(1-0.5d0*x)*logx**2-2*x*logomx
     &     -(3*logomx+logomx**2)*1.5d0*Pgq(x))
     &     +4d0*(28d0/9+65d0/18*x+44d0/9*x**2-(12+5*x+8d0/3*x**2)*logx
     &     +(4+x)*logx**2+2*x*logomx+S2(x)*1.5d0*Pgq(-x)
     &     +(0.5d0-2*logx*logomx+0.5d0*logx**2+11d0/3*logomx+logomx**2
     &     -pi**2/6)*1.5d0*Pgq(x))
     &     +2d0/3*nf*(-4d0*x/3-(20d0/9+4d0/3*logomx)*1.5d0*Pgq(x))


      P2gq=P2gq/4

      return
      end


C    S2: Eq. (4.114) ESW

      function S2(x)
      implicit none
      real *8 x,pi,S2,myli2
      external myli2      
      pi=3.14159265358979d0

      S2=-2*myli2(-x)+0.5d0*dlog(x)**2-2*dlog(x)*dlog(1+x)-pi**2/6
      return
      end
