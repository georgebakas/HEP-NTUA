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
      integer j,i
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

      lhapdfs=.false.
      isol=.false.

      virtonly=.false.
      realonly=.false.

      noglue=.false.
      ggonly=.false.
      gqonly=.false. 
      qqonly=.false.
      qqbonly=.false.

      nmin=1
      nmax=2

      clustering=.true.
      colourchoice=0
c      rtsmin=40d0

ch
      rtsmin=2*mt
ch
      cutoff=0.0005d0

CC    

      Qflag=.true.
      Gflag=.true.

CC
      aii=1d0
      aif=1d0
      afi=1d0
      aff=1d0

      inclusive=.true.
 
CC   Parameters used to define jets 
CC   Logical variable 'algorithm' can be taken to be 'ktal', 'ankt' or 'cone'

      algorithm='ankt'      

      ptjetmin=0d0
      etajetmin=0d0
      etajetmax=20d0

      Rcut=0.4d0

CC    Dynamic scale (if true muf=mur=q)

      dynamicscale=.false.

      removebr=.false.
      makecuts=.true.     


CC    Adjust the grid at each iteration

      dorebin=.true.


CC    Read a previously saved grid

      readin=.false.


      writeout=.false.
      ingridfile=''
      outgridfile=''

CC    Read inputfile


      read(*,*) sroot
      read(*,*) ih1,ih2
      read(*,*) nproc !decay mode
      read(*,*) scale,facscale  ! mur,muf
      read(*,*) order        
      read(*,*) part
      read(*,*) zerowidth
      read(*,*) Mwmin,Mwmax
      read(*,*) itmx1,ncall1
      read(*,*) itmx2,ncall2
      read(*,*) rseed          ! seed
      read(*,*) iset,nset     
      read(*,*) PDFname,PDFmember
      read(*,*) runstring

      if(Mwmin.lt.2d0*mt)then
       Mwmin=2d0*mt
      endif


CC    Set all factorization scales to facscale
CC    to avoid problems when dynamicscale=.false.

      do nd=0,40
       dipscale(nd)=facscale
      enddo

CC    Cut on qt/Q

      xqtcut=0.0005d0
ch      xqtcut=0.0005d0
ch      xqtcut=1d0/8000d0


C     Limits on invariant mass of vector boson decay products
C     (irrelevant if zerowidth=.true.)


      wsqmin=Mwmin**2
      wsqmax=Mwmax**2
      
       

C     Check if the limits are compatible with sroot

      if(wsqmax.gt.(sroot**2)) wsqmax=sroot**2


      do j=1,mxpart
      plabel(j)=''
      enddo


      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio

      BrnRat=1d0

      call coupling

      call cstring(pdfstring)

      if(lhapdfs.eqv.(.false.)) then
       write(6,*)'CCCCCCCCCC      Parton Distributions    CCCCCCCCCCCC'
       write(6,*)'C                                                  C'
       write(6,*)'C       ', pdfstring,'       C'
       write(6,*)'C                                                  C'
      endif
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

      notag=0
      nqcdjets=0
      isub=0


      
      mb=0

    
    
     
      if((ih1.eq.1).and.(ih2.eq.-1)) then
       str1='ppb'
      elseif((ih1.eq.1).and.(ih2.eq.1)) then
       str1='pp'
      else
       write(*,*)'Initial state not allowed'
      endif

ch
        nqcdjets=0
        ndim=4
        mass2=mt
        ndec=2
        n2=0
        n3=0
        i1=3
        i2=4
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        str3='-> Q-(p3)+Qb-(p4)'


ch


CC    Here ndim is the number of dimensions for Born

c      if(nproc.eq.1) then
C
C     W+->l+nubar
C
c       str2=' -> W+ -> ' 
c       ndec=2
c       ndim=4    
c       plabel(3)='nl'
c       plabel(4)='ea'
c       plabel(5)='pp'
c       plabel(6)='pp'
c       n2=0
c       i1=4
c       i2=4
c       nwz=1
c       mass3=wmass
c       width3=wwidth

c       if (removebr) then
c        call branch(brwen,brzee,brtau,brtop)
c        BrnRat=brwen
c       endif

c      str3=' nu(p3)+e+(p4)'

c      elseif(nproc.eq.2) then
C
C      W-=>l-nu
C
c       str2=' -> W- -> ' 
c       ndec=2
c       ndim=4
c       plabel(3)='el'
c       plabel(4)='na'
c       plabel(5)='pp'
c       plabel(6)='pp'
c       n2=0
c       i1=3
c       i2=3
c       nwz=-1
c       mass3=wmass
c       width3=wwidth
     
       
     

c       if (removebr) then
c        call branch(brwen,brzee,brtau,brtop)
c        BrnRat=brwen
c       endif

c      str3=' e^-(p3)+nu~(p4)'

c      elseif(nproc.eq.3) then
C
c      Z->e+e-
C
c       str2=' -> Z -> ' 
c       ndec=2
c       ndim=4
c       plabel(3)='el'
c       plabel(4)='ea'
c       plabel(5)='pp'
c       plabel(6)='pp'
c       n2=0
c       i1=3
c       i2=4
c       nwz=0
c       mass3=zmass
c       width3=zwidth

c       l1=le     
c       r1=re

c      q1=0 switch off the photon

c       q1=-1

c       int=.false.
c       str3=' e-(p3)+e+(p4)'
       
     

c       if (removebr) then
c        call branch(brwen,brzee,brtau,brtop)
c        BrnRat=brzee
c       endif




c      else
c       write(*,*)'Wrong decay channel'
c       stop
c      endif


C     New: decide if using Breit-Wigner or not

c      n3=0
c      n3=1
c      vmass=zmass
c      if(nproc.eq.3) vmass=zmass

c      if(mwmin.lt.vmass.and.mwmax.gt.vmass) n3=1



CCCCCCCCCCCCC

      call strcat(str1,str2,string)
      call strcat(string,str3,string)

     
      call cstring(string)      
     

      write(6,*)'C                                                  C'

      if(order.eq.0) then
      write(6,*)'C         Computing LO cross section for           C'
      elseif(order.eq.1) then
      write(6,*)'C        Computing NLO cross section for           C'
      elseif(order.eq.2) then
      write(6,*)'C        Computing NNLO cross section for          C'
      else
      write(*,*)'Order can be 0,1 or 2 !'
       stop
      endif
  
      write(6,*)'C                                                  C' 
      write(*,*)'C',string,'C'
      write(6,*)'C                                                  C' 
      write(6,96)sroot
      write(6,*)'C                                                  C'
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 


 96   format(' C           at  Sqrt(s)=',f8.2,' GeV               C')

      nqcdjets=0
    

      call ckmfill(nwz)

    

CCCCCCCCCCCCCCCCCC


c--- set-up the random number generator with a negative seed
      idum=-abs(rseed)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)


c--- check that we have a valid value of 'part'
      if ( (part .ne. 'lord') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
          write(6,*) 'part=',part,' is not a valid option'
          stop     
      endif      



        as=alphas(scale,amz,nlooprun)
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
      
