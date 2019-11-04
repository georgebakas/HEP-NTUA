      subroutine gamma1squared(pin, mu_Q, m_t, value, channel)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
c      double precision pin(mxpart*4)
      double precision pin(20)
c      double precision pout(5,4)
      double precision pout(mxpart,4)
c      double precision H1qqdelta, T23H1qq, T34H1qq, T13H1qq
c      double precision H1ggdelta, T23H1gg, T34H1gg, T13H1gg
c      double precision H1, H1_T34, H1_T13, H1_T23
      double precision mu_Q
c mt could be set only once !!!
      double precision m_t
      double precision value
      double precision countint
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
c      call col_operators(pout)
c      if (channel.eq.1) then
c         H1 = H1ggdelta(pout)
c         H1_T34 = T34H1gg(pout)
c         H1_T13 = T13H1gg(pout)
cc         H1_T23 = T23H1gg(pout)
c      elseif (channel.eq.2) then
c         H1 = H1qqdelta(pout)
c         H1_T34 = T34H1qq(pout)
c         H1_T13 = T13H1qq(pout)
c         H1_T23 = T23H1qq(pout)
c      endif

      value = countint(pout)

      return
      end



CC    Counterterm to be subtracted from real+virt to get a finite
CC    cross section at qt->0

C     Version that allows to separate also qg channel

C     Scale dependence included up to NNLO

      double precision function countint(ptrans)
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
      double precision r_cut
      common/xmio/xmio
      common/xx0/xx0
      common/qtcut/xqtcut
      common/nnlo/order 
      common/r_cut/r_cut

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

      q2=(ptrans(3,4)+ptrans(4,4))**2
     &     -(ptrans(3,1)+ptrans(4,1))**2
     &     -(ptrans(3,2)+ptrans(4,2))**2
     &     -(ptrans(3,3)+ptrans(4,3))**2

      qt2=(ptrans(3,1)+ptrans(4,1))**2
     &     +(ptrans(3,2)+ptrans(4,2))**2




      v34=dsqrt(1-4d0*mt**4/((2d0*dot(ptrans,3,4))**2))
      beta34=0.5d0*dlog((1d0+v34)/(1d0-v34))
      pt=dsqrt(ptrans(4,1)**2+ptrans(4,2)**2)


      call qqb_QQb(ptrans,msqc,1)

      call col_operators(ptrans)

      betat=dsqrt(1d0-4*mt**2/q2)



      Gamma1q=0.5d0*(
     .       Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tqq(2,4))


      Gamma1qb=0.5d0*(
     .       Tqq(3,3)+Tqq(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tqq(3,4)
     .      +dlog(4d0*dot(ptransb,1,3)**2/(q2*mt**2))*Tqq(1,3)
     .      +dlog(4d0*dot(ptransb,1,4)**2/(q2*mt**2))*Tqq(1,4)
     .      +dlog(4d0*dot(ptransb,2,3)**2/(q2*mt**2))*Tqq(2,3)
     .      +dlog(4d0*dot(ptransb,2,4)**2/(q2*mt**2))*Tqq(2,4))



      Gamma1g=0.5d0*(
     .       Tgg(3,3)+Tgg(4,4)
     .      +1d0/v34*dlog((1d0+v34)/(1d0-v34))
     .      *Tgg(3,4)
     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tgg(1,3)
     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tgg(1,4)
     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tgg(2,3)
     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tgg(2,4))

c      if(order.eq.2)then
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



       H1qdelta=H1qqdelta(ptrans)/msqc(1,-1)
       H1qbdelta=H1qqdelta(ptransb)/msqc(1,-1)
       H1gdelta=H1ggdelta(ptrans)/msqc(0,0)


ch     The colour-correlated 1-loop amplitudes squared

       T34H1q=T34H1qq(ptrans)
c/msqc(1,-1)

       T34H1qb=T34H1qq(ptransb)
c/msqc(1,-1)

       T34H1g=T34H1gg(ptrans)
c/msqc(0,0)

       T13H1q=T13H1qq(ptrans)
c/msqc(1,-1)

       T13H1qb=T13H1qq(ptransb)
c/msqc(1,-1)

       T13H1g=T13H1gg(ptrans)
c/msqc(0,0)

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
       T23H1g=-(T34H1g+T13H1g+T33H1g)

       
       T23H1qb=-(T34H1qb+T13H1qb+T33H1qb)
       T14H1q=T23H1q
       T14H1qb=T23H1qb
       T14H1g=T23H1g


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


c      else
c       H1qdelta=0d0
c       H1qbdelta=0d0
c       H1gdelta=0d0
c 
c      Gamma1H1q=0d0
c       Gamma1H1qb=0d0
c       Gamma1sq_q=0d0
c       Gamma1sq_qb=0d0
c       Gamma2q=0d0
c       Gamma2qb=0d0
c       Gamma2g=0d0
c       Gamma1sq_g=0d0
c      endif

      countint = Gamma1sq_q
c      countint = Gamma1H1q
c      countint = Gamma2q
c      countint = Gamma2g
c      countint = Gamma1sq_g

      return
      end

