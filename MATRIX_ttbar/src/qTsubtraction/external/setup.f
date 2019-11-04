      subroutine setup
      implicit none
      include 'constants.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'alfacut.f'
      include 'pdfiset.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'nlooprun.f'
      include 'rescoeff.f'
      include 'jetcuts.f'
      include 'flags.f'
      include 'dipolescale.f'
      include 'gammacusp.f'
      include 'born_col_correl.f'
CC
      include 'vegas_common.f'
      include 'nwz.f'
      include 'ewcharge.f'
      include 'dynamicscale.f'
      include 'lhapdf.f'
      include 'inner_prod.f'
      character *2 plabel(mxpart)
      common/plabel/plabel
      integer order,notag,nqcdjets,nqcdstart,isub,nproc,ndec,nd
      common/nnlo/order
      common/notag/notag
      common/nqcdjets/nqcdjets,nqcdstart
      common/isub/isub
      double precision BrnRat,gamgambr,wwbr,zzbr,br0
      common/BrnRat/BrnRat
      common/nproc/nproc

      double precision xqtcut
      common/qtcut/xqtcut      

      double precision beta1,H2qqdelta,H2qqD0,
     &                 H2ggdelta,H2ggD0
      common/Hstcoeff/beta1,H2qqdelta,H2qqD0,
     &                 H2ggdelta,H2ggD0

      logical isol
      common/isol/isol

C     Labels that identify the charged leptons

      integer i1,i2
      common/isolabel/i1,i2

      logical int
      common/int/int

      logical dorebin
      common/dorebin/dorebin

      character *50 prefix
      character *36 pdfstring
      integer nset
      common/prefix/nset,prefix
      common/pdfstring/pdfstring

      character*4 part
      character*30 runstring
      integer j,i,ii,jj
      logical makecuts
      integer nmin,nmax,n2,n3,n30
      integer ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,rseed
      double precision rtsmin,sroot,LT
      double precision Mwmin,Mwmax
      double precision Rcut
      double precision ran2,randummy
      double precision cmass,bmass
      double precision mass2,width2,mass3,width3,vmass
      double precision amz,alphas
      double precision brwen,brzee,brtau,brtop
      double precision gammaq2test      
      character *3 str1
      character *10 str2
      character *38 str3
      character *50 string


      common/couple/amz
      
      common/breit/n2,n3,mass2,width2,mass3,width3
      
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin
      common/mwminmax/Mwmin,Mwmax
 

      common/part/part
      common/runstring/runstring
      common/energy/sroot
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/ranno/idum

      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/qmass/cmass,bmass

      common/rseed/rseed
      save /ranno/

      logical lhapdfs
      common/lhapdfs/lhapdfs

c      mt = 173.3d0
c      scale = mt
c      facscale = mt

c        as=alphas(scale,amz,nlooprun)
      as=1
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CH    Colour matrices

CH    qq channel

      do j=1,4
       Tqq(j,j)=cf
      enddo
      Tqq(1,2)=ca/2d0-cf
      Tqq(3,4)=Tqq(1,2)
      Tqq(1,3)=ca/2d0-2d0*cf
      Tqq(1,4)=-ca+2d0*cf
      Tqq(2,3)=Tqq(1,4)
      Tqq(2,4)=Tqq(1,3)
      do j=1,4
      do i=j+1,4
       Tqq(i,j)=Tqq(j,i)
      enddo
      enddo
ch    non-trivial 4-colour matrices
      
      do i=1,4
        do j=i,4
           do ii=1,4
              do jj=1,4
                 Tqq4(i,j,ii,jj)=Tqq(i,j)*Tqq(ii,jj)
              enddo
           enddo
        enddo
      enddo
      
c      do i=1,4
c         Tqq4(i,i,i,i)=Tqq(i,i)*Tqq(i,i)
c         do ii=1,4
c            Tqq4(i,i,ii,ii)=Tqq(i,i)*Tqq(ii,ii)
c         enddo
c         do j=i+1,4
c            do ii=1,4
c               Tqq4(i,j,ii,ii)=Tqq(i,j)*Tqq(ii,ii)
c               do jj=ii+1,4
c                  Tqq4(i,j,ii,jj)=Tqq(i,j)*Tqq(ii,jj)
c                  Tqq4(i,i,ii,jj)=Tqq(i,i)*Tqq(ii,jj)
c               enddo
c            enddo
c         enddo
c      enddo

      Tqq4(1,3,1,3)=cf/2d0/xn+Tqq(1,3)**2
      Tqq4(2,3,2,3)=cf/2d0/xn+Tqq(2,3)**2
      Tqq4(1,3,2,3)=-cf/2d0/xn+Tqq(1,3)*Tqq(2,3)
      Tqq4(2,3,1,3)=Tqq4(1,3,2,3)
      Tqq4(1,4,1,4)=Tqq4(2,3,2,3)
      Tqq4(2,4,2,4)=Tqq4(1,3,1,3)
      Tqq4(1,4,2,4)=Tqq4(2,3,1,3)
      Tqq4(2,4,1,4)=Tqq4(1,3,2,3)
      Tqq4(1,3,2,4)=Tqq4(1,3,1,3)
      Tqq4(1,3,1,4)=Tqq4(1,3,2,3)
      Tqq4(1,4,1,3)=Tqq4(2,3,1,3)
      Tqq4(1,4,2,3)=Tqq4(2,3,2,3)
      Tqq4(2,3,1,4)=Tqq4(2,3,2,3)
      Tqq4(2,3,2,4)=Tqq4(2,3,1,3)
      Tqq4(2,4,1,3)=Tqq4(1,3,1,3)
      Tqq4(2,4,2,3)=Tqq4(1,3,2,3)

ch    The rest of 4-colour operators are just a product of corresponing 2-colour operators

ch    Set the inner product of orthogonal basis

      c1qs=xn**2
      c2qs=(xn**2-1d0)/4d0
      c1gs=xn*(xn**2-1d0)
      c2gs=xn*(xn**2-1d0)/2d0
      c3gs=(xn**2-4d0)*(xn**2-1d0)/2d0/xn

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC    Resummation coefficients

      beta0=(33-2*nf)/12d0


      beta1=(153d0-19*nf)/24d0

      Kappa=67/6d0-(pi**2)/2d0-5d0/9d0*nf

      A1q=4d0/3
      A2q=0.5d0*A1q*Kappa
      B1q=-2d0
      A1g=3d0
      A2g=0.5d0*A1g*Kappa
      B1g=-2*beta0


ch      B2q=4d0/9*(pi**2-3d0/4-12*Z3)+(11d0/9*pi**2-193d0/12+6*Z3)
ch     & +nf/6d0*(17d0/3-4d0/9*pi**2)
ch      write(*,*)B2q,'b'

ch   B2q in the hard scheme !

      B2q=1d0/16*
     &    (cf**2*(-3d0+4d0*pi**2-48d0*Z3)
     &     +cf*ca*(-17d0/3-44d0/9*pi**2+24d0*Z3)
     &     +cf*nf*(2d0/3+8d0/9*pi**2))
     &     +beta0*cf*pi**2/6d0

      gammaq2test=1d0/16*
     &    (cf**2*(-3d0/2+2*pi**2-24*Z3)
     &     +cf*ca*(-961d0/54-11*Pi**2/6d0+26d0*Z3)
     &     +cf*nf*(130d0/27+2*pi**2/3))
C     Delta term in c1qq coefficient
ch      B2q=B2q-2*gammaq2test

      C1qqdelta=(pi**2-8)/3d0

C     Delta term in P2qq splitting function (as/pi normalization)

      Delta2qq=16d0/9*(3d0/8-pi**2/2+6*Z3)
     &   +4*(17d0/24+11d0*pi**2/18-3*Z3)-2d0/3*nf*(1d0/6+2*pi**2/9d0)

      Delta2qq=Delta2qq/4d0

CC    Coefficients of D0 and D1 in P*P (as/pi normalization)

      D0qqqq=8d0/3
      D1qqqq=32d0/9


CC    Coefficients of delta(1-z) in P*P

      Deltaqqqq=4d0/9*(9d0/4-2*pi**2/3d0)

C     H2qq contribution: coefficient of delta(1-z)

ch    H2qqdelta=76.82+-0.25
      H2qqdelta=76.82d0
ch      H2qqdelta=0d0

C     H2qq contribution: coefficient of D0(z)

      H2qqD0=-404d0/27+(56d0*nf)/81+14*Z3

      gammacusp2q=(67/36d0-Pi**2/12)*ca-5/18d0*nf
      gammaQ1=1d0/16*(Cf*ca*(2*Pi**2/3-98d0/9-4*Z3)
     &                +20d0/9*Cf*nf)

ch    gg channel

C     B2g from qt code (checked !)

ch      B2g=9*(23/24d0+(11*pi**2)/18d0-3*Z3/2d0)+
ch     /    2*nf/3d0-3*nf*(1/12d0+pi**2/9d0)-11/2d0

ch   B2g in the hard scheme !

      B2g=1d0/16*
     &    (ca**2*(-64d0/3d0-24d0*Z3)
     &     +16d0/3*ca*nf
     &     +4d0*cf*nf)
     &     +beta0*ca*pi**2/6d0

      

C     Delta term in c1gg coefficient
CH
      C1ggdelta=(11+3*pi**2)/4d0

CH
C     Delta term in P2gg splitting function (as/pi normalization)

      Delta2gg=9*(8d0/3+3*Z3)-2d0/3*nf-2*nf

      Delta2gg=Delta2gg/4d0

CC    Coefficients of D0 and D1 in P*P (as/pi normalization)

      D0gggg=6*beta0
      D1gggg=18d0

CC    Coefficients of delta(1-z) in P*P

      Deltagggg=beta0**2-3d0/2*pi**2

C     H2gg contribution: coefficient of delta(1-z)

ch      LT=2*dlog(hmass/mtop)

ch      H2ggdelta=11399d0/144+19d0/8*LT+133*Pi**2/8+13d0/16*Pi**4
ch     #          -55d0/2*Z3+nf*(-1189d0/144+2d0*LT/3-5*Pi**2/12)

ch      LT=2d0*LT

C     H2gg contribution: coefficient of D0(z)

!     Use the fitted value from 8TeV run
ch    H2ggdelta=103.8+-0.9
      H2ggdelta=103.8d0
ch      H2ggdelta=0d0

      H2ggD0=-101d0/3+14d0*nf/9d0+63d0/2*Z3


      return

      end
      
