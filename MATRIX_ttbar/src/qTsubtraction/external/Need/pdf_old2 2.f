      SUBROUTINE fdist(ih,x,q,fx)
      implicit none

      include 'pdfiset.f'
      include 'constants.f'
      real * 8 FX(-nf:nf)
      REAL * 8 X,Q,Q2,DUV,DDV,DDEL,DUDB,DSB,DGL,CHR,BOT,UB,DB,DQ1
      REAL * 8 BBAR,CBAR,DSBAR,PHOT
      REAL * 8 ctq4fn,ctq5pdf,ctq6pdf

C     For Alekhin pdfs

      real*8 pdfs09(0:8),dpdfs09(0:8,25)
      real*8 pdfs06(0:9),dpdfs06(0:9,23)
C
C     For Reya pdfs
       double precision JR09VFNNLOxuv,JR09VFNNLOxdv,JR09VFNNLOxgl,
     &                  JR09VFNNLOxub,JR09VFNNLOxdb,JR09VFNNLOxsb,
     &                  JR09VFNNLOxcb,JR09VFNNLOxbb,
     &                  JR09VFNNLOalphas,
     &                  xuv(-13:13),xdv(-13:13),xgl(-13:13),xub(-13:13),
     &                  xdb(-13:13),xsb(-13:13),xcb(-13:13),xbb(-13:13),
     &                  alphas(-13:13),
     &                  exuv,exdv,exgl,exub,exdb,exsb,excb,exbb,ealphas
       double precision GJR08VFNSxuv,GJR08VFNSxdv,GJR08VFNSxgl,
     &                  GJR08VFNSxub,GJR08VFNSxdb,GJR08VFNSxsb,
     &                  GJR08VFNSxcb,GJR08VFNSxbb,
     &                  GJR08VFNSalphas
C
      integer j,mode,ih
      integer NPDF,NPAR

      character *50 prefix,prefix1
      integer nset
      common/prefix/nset,prefix


      
     
      Q2=Q**2
  
   

C Fix to prevent undefined math operations for x=1.
C Assumes that all structure functions vanish for x=1.


      IF(1-X.EQ.0) THEN
         DO J=-NF,NF
            FX(J) = 0
         ENDDO
         RETURN
      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C       CTEQ4

        if (iset.eq.1) then
           mode=3
           DGL=Ctq4Fn (mode,0, x, Q)*x
           UB=Ctq4Fn (mode,-1, x, Q)*x
           DB=Ctq4Fn (mode,-2, x, Q)*x
           DSB=Ctq4Fn (mode,-3, x, Q)*x
           CHR=Ctq4Fn (mode,-4, x, Q)*x
           BOT=Ctq4Fn (mode,-5, x, Q)*x
           DUV=Ctq4Fn (mode,1, x, Q)*x - UB
           DDV=Ctq4Fn (mode,2, x, Q)*x - DB
          elseif (iset.eq.2) then
           mode=1
           DGL=Ctq4Fn (mode,0, x, Q)*x
           UB=Ctq4Fn (mode,-1, x, Q)*x
           DB=Ctq4Fn (mode,-2, x, Q)*x
           DSB=Ctq4Fn (mode,-3, x, Q)*x
           CHR=Ctq4Fn (mode,-4, x, Q)*x
           BOT=Ctq4Fn (mode,-5, x, Q)*x
           DUV=Ctq4Fn (mode,1, x, Q)*x - UB
           DDV=Ctq4Fn (mode,2, x, Q)*x - DB 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         

C MRST98
          elseif ((ISET.LT.17).AND.(ISET.GT.10)) THEN
           mode=iset-10
           call mrs98(x,q2,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C CTEQ5

          elseif ((ISET.LE.29).AND.(ISET.GT.20)) THEN          
           DGL=Ctq5Pdf (0, X, Q)*x
           UB=Ctq5Pdf (-1, X, Q)*x
           DB=Ctq5Pdf (-2, X, Q)*x
           DSB=Ctq5Pdf (-3, X, Q)*x
           CHR=Ctq5Pdf (-4, X, Q)*x
           BOT=Ctq5Pdf (-5, X, Q)*x
           DUV=Ctq5Pdf (1, X, Q)*x - UB
           DDV=Ctq5Pdf (2, X, Q)*x - DB

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MRST99
          elseif (iset.eq.30) then
           mode=1
           call mrs99(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MRST2001

      ELSEIF ((ISET.LT.45).AND.(ISET.GT.40)) THEN
      mode=iset-40
      call mrst2001(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

      ELSEIF ((ISET.LT.49).AND.(ISET.GT.44)) THEN
      mode=iset-44
      call mrstnnlo(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MRST2002 LO

      ELSEIF (ISET.eq.49) THEN
      mode=iset-48
      call mrstlo(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C CTEQ6

      elseif ((ISET.LE.54).AND.(ISET.GT.50)) THEN          
           DGL=Ctq6Pdf (0, X, Q)*x
           UB=Ctq6Pdf (-1, X, Q)*x
           DB=Ctq6Pdf (-2, X, Q)*x
           DSB=Ctq6Pdf (-3, X, Q)*x
           CHR=Ctq6Pdf (-4, X, Q)*x
           BOT=Ctq6Pdf (-5, X, Q)*x
           DUV=Ctq6Pdf (1, X, Q)*x - UB
           DDV=Ctq6Pdf (2, X, Q)*x - DB

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C CTEQ6.6

      elseif ((ISET.LE.59).AND.(ISET.GT.54)) THEN          
           DGL=Ctq6Pdf (0, X, Q)*x
           UB=Ctq6Pdf (-1, X, Q)*x
           DB=Ctq6Pdf (-2, X, Q)*x
           DSB=Ctq6Pdf (-3, X, Q)*x
           CHR=Ctq6Pdf (-4, X, Q)*x
           BOT=Ctq6Pdf (-5, X, Q)*x
           DUV=Ctq6Pdf (1, X, Q)*x - UB
           DDV=Ctq6Pdf (2, X, Q)*x - DB

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MRST2002

      elseif ((iset.eq.61).or.(iset.eq.62)) then
      mode=iset-60
      call mrst2002(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MRST2004

      elseif ((iset.eq.71).or.(iset.eq.72)) then
      mode=iset-70
      call mrst2004(x,q,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Alekhin NNLO A06
 
      elseif (iset.eq.75) then

      call a06(x,q2,pdfs06,dpdfs06,npdf,npar)
      if(nset.eq.0) then
        duv=pdfs06(1)
        ddv=pdfs06(2)
        dgl=pdfs06(3)
        ub=pdfs06(4)
        dsb=pdfs06(5)
        db=pdfs06(6)
        chr=pdfs06(7)
        bot=pdfs06(8)
       else
        duv=pdfs06(1)+dpdfs06(1,nset)
        ddv=pdfs06(2)+dpdfs06(2,nset)
        dgl=pdfs06(3)+dpdfs06(3,nset)
        ub=pdfs06(4)+dpdfs06(4,nset)
        dsb=pdfs06(5)+dpdfs06(5,nset)
        db=pdfs06(6)+dpdfs06(6,nset)
        chr=pdfs06(7)+dpdfs06(7,nset)
        bot=pdfs06(8)+dpdfs06(8,nset)
       endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Alekhin NNLO A09

      elseif (iset.eq.85) then     
       call a09(x,q2,pdfs09,dpdfs09,5,nset)
       if(nset.eq.0) then
        duv=pdfs09(1)
        ddv=pdfs09(2)
        dgl=pdfs09(3)
        ub=pdfs09(4)
        dsb=pdfs09(5)
        db=pdfs09(6)
        chr=pdfs09(7)
        bot=pdfs09(8)
       else
        duv=pdfs09(1)+dpdfs09(1,nset)
        ddv=pdfs09(2)+dpdfs09(2,nset)
        dgl=pdfs09(3)+dpdfs09(3,nset)
        ub=pdfs09(4)+dpdfs09(4,nset)
        dsb=pdfs09(5)+dpdfs09(5,nset)
        db=pdfs09(6)+dpdfs09(6,nset)
        chr=pdfs09(7)+dpdfs09(7,nset)
        bot=pdfs09(8)+dpdfs09(8,nset)
       endif



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C MSTW 2008

      elseif ((iset.gt.89).and.(iset.lt.93)) then

      prefix1 = prefix(1:len_trim(prefix))//'.68cl' ! 68% C.L. errors

      if(nset.eq.0)prefix1=prefix
    

      CALL GetAllPDFs(prefix1,nset,x,q,
     #    DUV,DDV,UB,DB,DSB,DSBAR,CHR,Cbar,BOT,bbar,DGL,phot)



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C GJR08VFLO   Dynamical Parton Distribution Functions

      elseif ((iset.gt.64).and.(iset.lt.66)) then

          DUV = GJR08VFNSxuv(x,Q2,14)
          DDV = GJR08VFNSxdv(x,Q2,14)
          DGL = GJR08VFNSxgl(x,Q2,14)
          UB = GJR08VFNSxub(x,Q2,14)
          DB = GJR08VFNSxdb(x,Q2,14)
          DSB = GJR08VFNSxsb(x,Q2,14)
          CHR = GJR08VFNSxcb(x,Q2,14)
          BOT = GJR08VFNSxbb(x,Q2,14)
!          alphas(set) = GJR08VFNSalphas(Q2,14)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C GJR08VFNLO   Dynamical Parton Distribution Functions

      elseif ((iset.gt.65).and.(iset.lt.67)) then

          DUV = GJR08VFNSxuv(x,Q2,nset)
          DDV = GJR08VFNSxdv(x,Q2,nset)
          DGL = GJR08VFNSxgl(x,Q2,nset)
          UB = GJR08VFNSxub(x,Q2,nset)
          DB = GJR08VFNSxdb(x,Q2,nset)
          DSB = GJR08VFNSxsb(x,Q2,nset)
          CHR = GJR08VFNSxcb(x,Q2,nset)
          BOT = GJR08VFNSxbb(x,Q2,nset)
!          alphas(set) = GJR08VFNSalphas(Q2,nset)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C JR09VFNNLO   Dynamical Parton Distribution Functions

      elseif ((iset.gt.66).and.(iset.lt.68)) then

          DUV = JR09VFNNLOxuv(x,Q2,nset)
          DDV = JR09VFNNLOxdv(x,Q2,nset)
          DGL = JR09VFNNLOxgl(x,Q2,nset)
          UB = JR09VFNNLOxub(x,Q2,nset)
          DB = JR09VFNNLOxdb(x,Q2,nset)
          DSB = JR09VFNNLOxsb(x,Q2,nset)
          CHR = JR09VFNNLOxcb(x,Q2,nset)
          BOT = JR09VFNNLOxbb(x,Q2,nset)
!          alphas(set) = JR09VFNNLOalphas(Q2,nset)



      ELSE
        WRITE(*,*)'NO SUCH DISTRIBUTION'
        STOP
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      if(iset.lt.90) then

      DSBAR=DSB
      CBAR=CHR
      BBAR=BOT

      endif

CC   BE CAREFUL !!

CC Different from HNNLO:
CC   u->2 d->1


c for protons
      if (ih.eq.1) then
      FX(0)=DGL/x
      FX(2)=(DUV+UB)/x
      FX(1)=(DDV+DB)/x
      FX(3)=DSB/x
      FX(4)=CHR/x
      FX(5)=BOT/x      
      FX(-2)=UB/x
      FX(-1)=DB/x
      FX(-3)=DSBAR/x
      FX(-4)=CBAR/x
      FX(-5)=BBAR/x
c for anti-protons
      elseif (ih.eq.-1) then
      FX(0)=DGL/x
      FX(-2)=(DUV+UB)/x
      FX(-1)=(DDV+DB)/x
      FX(-3)=DSB/x
      FX(-4)=CHR/x
      FX(-5)=BOT/x
      FX(2)=UB/x
      FX(1)=DB/x
      FX(3)=DSBAR/x
      FX(4)=CBAR/x
      FX(5)=BBAR/x
      endif 


       RETURN
      END



      subroutine mrs98(x,q2,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C****************************************************************C
C								 C
C     This is a package for the new MRS 1998 parton              C
C     distributions. The format is similar to the previous       C
C     (1996) MRS-R series.                                       C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV, MSbar factorization is assumed,     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C								 C
C     TEMPORARY NAMING SCHEME:                                   C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  1     FT08A  central gluon, a_s    300      0.1175   0.00561  C
C  2     FT09A  higher gluon          300      0.1175   0.00510  C
C  3     FT11A  lower gluon           300      0.1175   0.00408  C
C  4     FT24A  lower a_s             229      0.1125   0.00586  C
C  5     FT23A  higher a_s            383      0.1225   0.00410  C
C						                 C
C						                 C
C      The corresponding grid files are called ft08a.dat etc.    C
C							  	 C
C      The reference is:                                         C
C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
C      Univ. Durham preprint DTP/98/??, hep-ph/??????? (1998)    C
C                                                                C
C      Comments to : W.J.Stirling@durham.ac.uk                   C
C  
c  and for the LO sets
C     TEMPORARY NAMING SCHEME:                                   C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C								 C
C  6     LO05A  central gluon, a_s    174      0.1250   0.01518  C
C  7     LO09A  higher gluon          174      0.1250   0.01616  C
C  8     LO10A  lower gluon           174      0.1250   0.01533  C
C  9     LO01A  lower  a_s             136      0.1200   0.01652  C
C  10     LO07A  higher a_s           216      0.1300   0.01522  C
C						                 C
C						                 C
C      The corresponding grid files are called lt05a.dat etc.    C
c                                                                C
C								 C
C****************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
          if(mode.eq.1) then
        call mrs981(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrs982(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrs983(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrs984(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.5) then
        call mrs985(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
c from here LO 
      elseif(mode.eq.6) then
        call mrs986(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.7) then
        call mrs987(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.8) then
        call mrs988(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.9) then
        call mrs989(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.10) then
        call mrs9810(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ')
      return
      end

      subroutine mrs981(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/ft08a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs982(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/ft09a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs983(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/ft11a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      
      subroutine mrs984(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/ft24a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs985(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/ft23a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      

C****************************************************************C
C								 C
C     This is a package for the new MRS LO 1998 parton           C
C     distributions. The format is similar to the previous       C
C     (1996) MRS-R series.                                       C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV,                                     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C                                                                C
C     Reference Martin, Roberts, Stirling and Thorne             C
C     Durham preprint DTP/98/52 (August 1998)                    C
C								 C
C     TEMPORARY NAMING SCHEME:                                   C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  6     LO05A  central gluon, a_s    174      0.1250   0.01518  C
C  7     LO09A  higher gluon          174      0.1250   0.01616  C
C  8     LO10A  lower gluon           174      0.1250   0.01533  C
C  9     LO01A  lower a_s             136      0.1200   0.01652  C
C  10     LO07A  higher a_s           216      0.1300   0.01522  C
C						                 C
C						                 C

      subroutine mrs986(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/lo05a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs987(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/lo09a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs988(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/lo10a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      
      subroutine mrs989(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/lo01a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs9810(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/lo07a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end


      subroutine mrs99(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C****************************************************************C
C								 C
C     This is a package for the new **corrected** MRST parton    C
C     distributions. The format is similar to the previous       C
C     (1998) MRST series.                                        C
C								 C
C     NOTE: 7 new sets are added here, corresponding to shifting C
C     the small x HERA data up and down by 2.5%, and by varying  C
C     the charm and strange distributions, and by forcing a      C
C     larger d/u ratio at large x.                               C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV, MSbar factorization is assumed,     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C								 C
C     NAMING SCHEME:                                             C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
C  2     COR02  higher gluon          300      0.1175   0.00497  C
C  3     COR03  lower gluon           300      0.1175   0.00398  C
C  4     COR04  lower a_s             229      0.1125   0.00585  C
C  5     COR05  higher a_s            383      0.1225   0.00384  C
C  6     COR06  quarks up             303.3    0.1178   0.00497  C
C  7     COR07  quarks down           290.3    0.1171   0.00593  C
C  8     COR08  strange up            300      0.1175   0.00524  C
C  9     COR09  strange down          300      0.1175   0.00524  C
C  10    C0R10  charm up              300      0.1175   0.00525  C
C  11    COR11  charm down            300      0.1175   0.00524  C
C  12    COR12  larger d/u            300      0.1175   0.00515  C
C						                 C
C      The corresponding grid files are called cor01.dat etc.    C
C							  	 C
C      The reference is:                                         C
C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
C      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
C                                                                C
C      Comments to : W.J.Stirling@durham.ac.uk                   C
C                                                                C
C								 C
C****************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
c      write(6,*)q,q2
c      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99
c      if(x.lt.xmin.or.x.gt.xmax)       print 98
          if(mode.eq.1) then
        call mrs991(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrs992(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrs993(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrs994(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.5) then
        call mrs995(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.6) then
        call mrs996(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.7) then
        call mrs997(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.8) then
        call mrs998(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.9) then
        call mrs999(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.10) then
        call mrs9910(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.11) then
        call mrs9911(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.12) then
        call mrs9912(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ')
      return
      end

      subroutine mrs991(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor01.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs992(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor02.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs993(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor03.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      
      subroutine mrs994(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor04.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs995(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor05.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      

      subroutine mrs996(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor06.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs997(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor07.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      


      subroutine mrs998(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor08.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs999(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor09.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      


      subroutine mrs9910(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor10.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrs9911(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor11.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      
      subroutine mrs9912(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='Pdfdata/cor12.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      

 
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     From now 2001 MRST sets     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine mrstlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 LO parton            C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201xxx                                  C
C                                                               C
C  There is 1 pdf set corresponding to mode = 1                 C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.220 C
C  corresponding to alpha_s(M_Z) of 0.130                       C
C  This set reads a grid whose first number is 0.02868          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrstlo1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrstlo1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/lo2002.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 


 


      subroutine mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 NLO parton           C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0110215                                  C
C                                                               C
C  There are 4 pdf sets corresponding to mode = 1, 2, 3, 4      C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.323 C
C  corresponding to alpha_s(M_Z) of 0.119                       C
C  This set reads a grid whose first number is 0.00927          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.290         C
C  corresponding to alpha_s(M_Z) of 0.117                       C
C  This set reads a grid whose first number is 0.00953          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.362         C
C  corresponding to alpha_s(M_Z) of 0.121                       C
C  This set reads a grid whose first number is 0.00889          C
C                                                               C
C  Mode=4 gives the set MRST2001J which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.353(alpha_s(M_Z) = 0.121 C 
C  This set reads a grid whose first number is 0.00826          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst20011(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrst20012(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrst20013(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrst20014(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrst20011(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/alf119.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      subroutine mrst20012(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/alf117.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst20013(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/alf121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst20014(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/j121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

   


      subroutine mrstnnlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the MRST 2002 NNLO parton distributionsC
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201127                                  C
C                                                               C
C  There are 4 pdf sets corresponding to mode = 1, 2, 3, 4      C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.235 C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the `average' of the slow and fast evolutions    C 
C  This set reads a grid whose first number is 0.00725          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.235         C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the fast evolution                               C 
C  This set reads a grid whose first number is 0.00734          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.235         C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the slow evolution                               C 
C  This set reads a grid whose first number is 0.00739          C
C                                                               C
C  Mode=4 gives the set MRSTNNLOJ which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.267(alpha_s(M_Z) =0.1180 C 
C  This set reads a grid whose first number is 0.00865          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrstnnlo1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrstnnlo2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrstnnlo3(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrstnnlo4(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrstnnlo1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/vnvalf1155.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      subroutine mrstnnlo2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/vnvalf1155.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrstnnlo3(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/vnvalf1155b.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrstnnlo4(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/vnvalf1180j.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine jeppe1(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      PARAMETER(NNX=49,MMY=37)
      dimension xx(nx),yy(my),ff(nx,my),ff1(NNX,MMY),ff2(NNX,MMY),
     xff12(NNX,MMY),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe2(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx(xx,nx,x)
      m=locx(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

    


      real*8 function  polderiv(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end


      subroutine mrst2002(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2002 updated NLO and      C
C  NNLO parton distributions.                                   C 
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0211080                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with alpha_s(M_Z,NLO) = 0.1197      C  
C  This set reads a grid whose first number is 0.00949          C
C                                                               C
C  Mode=2 gives the NNLO set with alpha_s(M_Z,NNLO) = 0.1154    C
C  This set reads a grid whose first number is 0.00685          C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/mrst2002nlo.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/mrst2002nnlo.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst2004(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C                                                               C
C  This is a package for the new MRST 2004 'physical gluon' NLO C
C  and NNLO parton distributions.                               C
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0410230                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with Lambda(4) = 347 MeV            C
C  This set reads a grid called mrst2004nlo.dat                 C
C  whose first number is 0.00910                                C
C                                                               C
C  Mode=2 gives the NNLO set with Lambda(4) = 251 MeV           C
C  This set reads a grid called mrst2004nnlo.dat                C
C  whose first number is 0.00673                                C
C                                                               C
C  These fits use a new, physically motivated parametrisation   C
C  for the gluon at the starting scale, Q_0^2 = 1 GeV^2         C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst12(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrst22(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrst12(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/mrst2004nlo.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe12(nx,nq,xxl,qql,f1,cc1)
      call jeppe12(nx,nq,xxl,qql,f2,cc2)
      call jeppe12(nx,nq,xxl,qql,f3,cc3)
      call jeppe12(nx,nq,xxl,qql,f4,cc4)
      call jeppe12(nx,nq,xxl,qql,f6,cc6)
      call jeppe12(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe12(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe12(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe22(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe22(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst22(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='Pdfdata/mrst2004nnlo.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe12(nx,nq,xxl,qql,f1,cc1)
      call jeppe12(nx,nq,xxl,qql,f2,cc2)
      call jeppe12(nx,nq,xxl,qql,f3,cc3)
      call jeppe12(nx,nq,xxl,qql,f4,cc4)
      call jeppe12(nx,nq,xxl,qql,f6,cc6)
      call jeppe12(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe12(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe12(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe22(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe22(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe22(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine jeppe12(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      parameter(nnx=49,mmy=37)
      dimension xx(nx),yy(my),ff(nnx,mmy),ff1(nnx,mmy),ff2(nnx,mmy),
     xff12(nnx,mmy),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe22(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx(xx,nx,x)
      m=locx(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end


C============================================================================
C                CTEQ Parton Distribution Functions: Version 4
C                          June 21, 1996
C
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C   Ref: MSUHEP-60426, CTEQ-604, e-Print Archive: hep-ph/9606399
C
C   This package contains 9 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C   Iset   PDF      Description             Alpha_s(Mz)  Q0(GeV)  Table_File
C ---------------------------------------------------------------------------
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ for which Q starts at a lower value of 0.7 GeV.  
C   The Table_Files are assumed to be in the working directory.
C   
C   The function Ctq4Fn (Iset, Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar)
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.

C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq4Fn (Iset, Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Character Flnm(10)*11
      Common
     > / CtqPar_5_2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder
      Data (Flnm(I), I=1,10)
     > / 'cteq4m.tbl ', 'cteq4d.tbl ', 'cteq4l.tbl '
     > , 'cteq4a1.tbl', 'cteq4a2.tbl', 'cteq4m.tbl ', 'cteq4a4.tbl'
     > , 'cteq4a5.tbl', 'cteq4hj.tbl', 'cteq4lq.tbl' /
      Data Isetold, Isetmin, Isetmax / -987, 1, 10 /
      save Flnm, Isetold, Isetmin, Isetmax

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         If (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
            Print *, 'Invalid Iset number in Ctq4Fn :', Iset
            Stop
         Endif
         IU= NextUt()
         Open(IU, File='Pdfdata/'//Flnm(Iset), Status='OLD', Err=100)
         Call ReadTbl (IU)
         Close (IU)
         Isetold=Iset
      Endif

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq4Fn: ', X
        Stop
      Endif
      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq4Fn: ', Q
        Stop
      Endif
      If (Iparton .lt. -NfMx .or. Iparton .gt. NfMx) Then
        Print *, 'Iparton out of range in Ctq4Fn: ', Iparton
        Stop
      Endif

      Ctq4Fn = PartonX (Iparton, X, Q)
      if(Ctq4Fn.lt.0.D0)  Ctq4Fn = 0.D0

      Return

 100  Print *, ' Data file ', Flnm(Iset), ' cannot be opened '
     >//'in Ctq4Fn!!'
      Stop
C                             ********************
      End

      Function NextUt()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 50, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUt = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C


  

C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.

C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar_5_2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq5Pdf: ', X
        Stop
      Endif
      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq5Pdf: ', Q
        Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar_5_1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar_5_1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar_5_2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq5 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=9)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5m.tbl', 'cteq5d.tbl', 'cteq5l.tbl', 'cteq5hj.tbl'
     > , 'cteq5hq.tbl', 'cteq5f3.tbl', 'cteq5f4.tbl'
     > , 'cteq5m1.tbl', 'ctq5hq1.tbl'  /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
            Print *, 'Invalid Iset number in SetCtq5 :', Iset
            Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File='Pdfdata/'//Tablefile, Status='OLD', Err=100)
         Endif
         Call ReadTbl (IU)
         Close (IU)
         Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
C                             ********************
      End

      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar_5_1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar_5_2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

!      Function NextUn()
!C                                 Returns an unallocated FORTRAN i/o unit.
!      Logical EX
!C
!      Do 10 N = 10, 300
!         INQUIRE (UNIT=N, OPENED=EX)
!         If (.NOT. EX) then
!            NextUn = N
!            Return
!         Endif
! 10   Continue
!      Stop ' There is no available I/O unit. '
!C               *************************
!      End
!C

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


C============================================================================
C                CTEQ Parton Distribution Functions: Version 6
C                             January 24, 2002, v6.0
C                             April 10, 2002, v6.1
C
C   Ref: "New Generation of Parton Distributions with
C         Uncertainties from Global QCD Analysis"
C   By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       hep-ph/0201195
C
C   This package contains 3 standard sets of CTEQ6 PDF's and 40 up/down sets
C   with respect to CTEQ6M PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
C   4    CTEQ6L1  Leading Order           0.130     215   165    cteq6l1.tbl
C     ------------------------------
C   1xx  CTEQ6M1xx  +/- w.r.t. CTEQ6M     0.118     326   226    cteq6m1xx.tbl
C    (where xx=01--40)
C ---------------------------------------------------------------------------
C   ** ALL fits are obtained by using the same coupling strength 
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1 
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in hep-ph/0201195.
C
C   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   Pumplin@pa.msu.edu or Tung@pa.msu.edu.
C
C===========================================================================
!
!      Function Ctq6Pdf (Iparton, X, Q)
!      Implicit Double Precision (A-H,O-Z)
!      Logical Warn
!      Common
!     > / CtqPar_6_2 / Nx, Nt, NfMx
!     > / QCDtable /  Alambda, Nfl, Iorder
!
!      Data Warn /.true./
!      save Warn
!
!      If (X .lt. 0D0 .or. X .gt. 1D0) Then
!        Print *, 'X out of range in Ctq6Pdf: ', X
!        Stop
!      Endif
!      If (Q .lt. Alambda) Then
!        Print *, 'Q out of range in Ctq6Pdf: ', Q
!        Stop
!      Endif
!      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
!         If (Warn) Then
!C        put a warning for calling extra flavor.
!             Warn = .false.
!             Print *, 'Warning: Iparton out of range in Ctq6Pdf: '
!     >              , Iparton
!         Endif
!         Ctq6Pdf = 0D0
!         Return
!      Endif
!
!      Ctq6Pdf = PartonX6 (Iparton, X, Q)
!      if(Ctq6Pdf.lt.0D0)  Ctq6Pdf = 0D0
!
!      Return
!
!C                             ********************
!      End
!
!      Subroutine SetCtq6 (Iset)
!      Implicit Double Precision (A-H,O-Z)
!      Parameter (Isetmax0=4)
!      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
!      Data (Flnm(I), I=1,Isetmax0)
!     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l'/
!      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,101,140/
!      save
!
!C             If data file not initialized, do so.
!      If(Iset.ne.Isetold) then
!         IU= NextUn6()
!         If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
!            Tablefile=Flnm(Iset)//'.tbl'
!         Elseif (Iset.eq.Isetmax0) Then
!            Tablefile=Flnm(Iset)//'1.tbl'
!         Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
!            write(nn,'(I3)') Iset
!            Tablefile=Flnm(1)//nn//'.tbl'
!         Else
!            Print *, 'Invalid Iset number in SetCtq6 :', Iset
!            Stop
!         Endif
!         Open(IU, File='Pdfdata/'//Tablefile, Status='OLD', Err=100)
! 21      Call ReadTbl6 (IU)
!         Close (IU)
!         Isetold=Iset
!      Endif
!      Return
!
! 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
!     >//'in SetCtq6!!'
!      Stop
!C                             ********************
!      End
!
!      Subroutine ReadTbl6 (Nu)
!      Implicit Double Precision (A-H,O-Z)
!      Character Line*80
!      PARAMETER (MXX = 96, MXQ = 20, MXF = 5)
!      PARAMETER (MXPQX = (MXF + 3) * MXQ * MXX)
!      Common
!     > / CtqPar_6_1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
!     > / CtqPar_6_2 / Nx, Nt, NfMx
!     > / XQrange / Qini, Qmax, Xmin
!     > / QCDtable /  Alambda, Nfl, Iorder
!     > / Masstbl / Amass(6)
!
!      Read  (Nu, '(A)') Line
!      Read  (Nu, '(A)') Line
!      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
!      Iorder = Nint(Dr)
!      Nfl = Nint(Fl)
!      Alambda = Al
!
!      Read  (Nu, '(A)') Line
!      Read  (Nu, *) NX,  NT, NfMx
!
!      Read  (Nu, '(A)') Line
!      Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)
!
!      Read  (Nu, '(A)') Line
!      Read  (Nu, *) XMIN, (XV(I), I =0, NX)
!
!      Do 11 Iq = 0, NT
!         TV(Iq) = Log(Log (TV(Iq) /Al))
!   11 Continue
!C
!C                  Since quark = anti-quark for nfl>2 at this stage,
!C                  we Read  out only the non-redundent data points
!C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)
!
!      Nblk = (NX+1) * (NT+1)
!      Npts =  Nblk  * (NfMx+3)
!      Read  (Nu, '(A)') Line
!      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)
!
!      Return
!C                        ****************************
!      End
!
!      Function NextUn6()
!C                                 Returns an unallocated FORTRAN i/o unit.
!      Logical EX
!C
!      Do 10 N = 10, 300
!         INQUIRE (UNIT=N, OPENED=EX)
!         If (.NOT. EX) then
!            NextUn6 = N
!            Return
!         Endif
! 10   Continue
!      Stop ' There is no available I/O unit. '
!C               *************************
!      End
!C
!
!      SUBROUTINE POLINT6 (XA,YA,N,X,Y,DY)
!
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!C                                        Adapted from "Numerical Recipes"
!      PARAMETER (NMAX=10)
!      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
!      NS=1
!      DIF=ABS(X-XA(1))
!      DO 11 I=1,N
!        DIFT=ABS(X-XA(I))
!        IF (DIFT.LT.DIF) THEN
!          NS=I
!          DIF=DIFT
!        ENDIF
!        C(I)=YA(I)
!        D(I)=YA(I)
!11    CONTINUE
!      Y=YA(NS)
!      NS=NS-1
!      DO 13 M=1,N-1
!        DO 12 I=1,N-M
!          HO=XA(I)-X
!          HP=XA(I+M)-X
!          W=C(I+1)-D(I)
!          DEN=HO-HP
!          IF(DEN.EQ.0.)PAUSE
!          DEN=W/DEN
!          D(I)=HP*DEN
!          C(I)=HO*DEN
!12      CONTINUE
!        IF (2*NS.LT.N-M)THEN
!          DY=C(NS+1)
!        ELSE
!          DY=D(NS)
!          NS=NS-1
!        ENDIF
!        Y=Y+DY
!13    CONTINUE
!      RETURN
!      END
!
!      Function PartonX6 (IPRTN, XX, QQ)
!
!c  Given the parton distribution function in the array U in
!c  COMMON / PEVLDT / , this routine interpolates to find
!c  the parton distribution at an arbitray point in x and q.
!c
!      Implicit Double Precision (A-H,O-Z)
!
!      Parameter (MXX = 96, MXQ = 20, MXF = 5)
!      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))
!
!      Common
!     > / CtqPar_6_1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
!     > / CtqPar_6_2 / Nx, Nt, NfMx
!     > / XQrange / Qini, Qmax, Xmin
!
!      Dimension fvec(4), fij(4)
!      Dimension xvpow(0:mxx)
!      Data OneP / 1.00001d0 /
!      Data xpow / 0.3d0 /       !**** choice of interpolation variable
!      Data nqvec / 4 /
!      Data ientry / 0 /
!      Save ientry,xvpow
!
!c store the powers used for interpolation on first call...
!      if(ientry .eq. 0) then
!         ientry = 1
!
!         xvpow(0) = 0D0
!         do i = 1, nx
!            xvpow(i) = xv(i)**xpow
!         enddo
!      endif
!
!      X = XX
!      Q = QQ
!      tt = dlog(dlog(Q/Al))
!
!c      -------------    find lower end of interval containing x, i.e.,
!c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
!      JLx = -1
!      JU = Nx+1
! 11   If (JU-JLx .GT. 1) Then
!         JM = (JU+JLx) / 2
!         If (X .Ge. XV(JM)) Then
!            JLx = JM
!         Else
!            JU = JM
!         Endif
!         Goto 11
!      Endif
!C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
!C                           |---|---|---|...|---|-x-|---|...|---|---|
!C                     x     0  Xmin               x                 1
!C
!      If     (JLx .LE. -1) Then
!        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
!        Stop
!      ElseIf (JLx .Eq. 0) Then
!         Jx = 0
!      Elseif (JLx .LE. Nx-2) Then
!
!C                For interrior points, keep x in the middle, as shown above
!         Jx = JLx - 1
!      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then
!
!C                  We tolerate a slight over-shoot of one (OneP=1.00001),
!C              perhaps due to roundoff or whatever, but not more than that.
!C                                      Keep at least 4 points >= Jx
!         Jx = JLx - 2
!      Else
!        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
!        Stop
!      Endif
!C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.
!
!C                       This is the variable to be interpolated in
!      ss = x**xpow
!
!      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then
!
!c     initiation work for "interior bins": store the lattice points in s...
!      svec1 = xvpow(jx)
!      svec2 = xvpow(jx+1)
!      svec3 = xvpow(jx+2)
!      svec4 = xvpow(jx+3)
!
!      s12 = svec1 - svec2
!      s13 = svec1 - svec3
!      s23 = svec2 - svec3
!      s24 = svec2 - svec4
!      s34 = svec3 - svec4
!
!      sy2 = ss - svec2
!      sy3 = ss - svec3
!
!c constants needed for interpolating in s at fixed t lattice points...
!      const1 = s13/s23
!      const2 = s12/s23
!      const3 = s34/s23
!      const4 = s24/s23
!      s1213 = s12 + s13
!      s2434 = s24 + s34
!      sdet = s12*s34 - s1213*s2434
!      tmp = sy2*sy3/sdet
!      const5 = (s34*sy2-s2434*sy3)*tmp/s12
!      const6 = (s1213*sy2-s12*sy3)*tmp/s34
!
!      EndIf
!
!c         --------------Now find lower end of interval containing Q, i.e.,
!c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
!      JLq = -1
!      JU = NT+1
! 12   If (JU-JLq .GT. 1) Then
!         JM = (JU+JLq) / 2
!         If (tt .GE. TV(JM)) Then
!            JLq = JM
!         Else
!            JU = JM
!         Endif
!         Goto 12
!       Endif
!
!      If     (JLq .LE. 0) Then
!         Jq = 0
!      Elseif (JLq .LE. Nt-2) Then
!C                                  keep q in the middle, as shown above
!         Jq = JLq - 1
!      Else
!C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
!        Jq = Nt - 3
!
!      Endif
!C                                   This is the interpolation variable in Q
!
!      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
!c                                        store the lattice points in t...
!      tvec1 = Tv(jq)
!      tvec2 = Tv(jq+1)
!      tvec3 = Tv(jq+2)
!      tvec4 = Tv(jq+3)
!
!      t12 = tvec1 - tvec2
!      t13 = tvec1 - tvec3
!      t23 = tvec2 - tvec3
!      t24 = tvec2 - tvec4
!      t34 = tvec3 - tvec4
!
!      ty2 = tt - tvec2
!      ty3 = tt - tvec3
!
!      tmp1 = t12 + t13
!      tmp2 = t24 + t34
!
!      tdet = t12*t34 - tmp1*tmp2
!
!      EndIf
!
!
!c get the pdf function values at the lattice points...
!
!      If (Iprtn .GE. 3) Then
!         Ip = - Iprtn
!      Else
!         Ip = Iprtn
!      EndIf
!      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1
!
!      Do it = 1, nqvec
!
!         J1  = jtmp + it*(NX+1)
!
!       If (Jx .Eq. 0) Then
!C                          For the first 4 x points, interpolate x^2*f(x,Q)
!C                           This applies to the two lowest bins JLx = 0, 1
!C            We can not put the JLx.eq.1 bin into the "interrior" section
!C                           (as we do for q), since Upd(J1) is undefined.
!         fij(1) = 0
!         fij(2) = Upd(J1+1) * XV(1)**2
!         fij(3) = Upd(J1+2) * XV(2)**2
!         fij(4) = Upd(J1+3) * XV(3)**2
!C
!C                 Use Polint6 which allows x to be anywhere w.r.t. the grid
!
!         Call Polint6 (XVpow(0), Fij(1), 4, ss, Fx, Dfx)
!
!         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
!C                                              Pdf is undefined for x.eq.0
!       ElseIf  (JLx .Eq. Nx-1) Then
!C                                                This is the highest x bin:
!
!        Call Polint6 (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)
!
!        Fvec(it) = Fx
!
!       Else
!C                       for all interior points, use Jon's in-line function
!C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
!         sf2 = Upd(J1+1)
!         sf3 = Upd(J1+2)
!
!         g1 =  sf2*const1 - sf3*const2
!         g4 = -sf2*const3 + sf3*const4
!
!         Fvec(it) = (const5*(Upd(J1)-g1)
!     &               + const6*(Upd(J1+3)-g4)
!     &               + sf2*sy3 - sf3*sy2) / s23
!
!       Endif
!
!      enddo
!C                                   We now have the four values Fvec(1:4)
!c     interpolate in t...
!
!      If (JLq .LE. 0) Then
!C                         1st Q-bin, as well as extrapolation to lower Q
!        Call Polint6 (TV(0), Fvec(1), 4, tt, ff, Dfq)
!
!      ElseIf (JLq .GE. Nt-1) Then
!C                         Last Q-bin, as well as extrapolation to higher Q
!        Call Polint6 (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
!      Else
!C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
!C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
!C                         the full range QV(0:Nt)  (in contrast to XV)
!        tf2 = fvec(2)
!        tf3 = fvec(3)
!
!        g1 = ( tf2*t13 - tf3*t12) / t23
!        g4 = (-tf2*t34 + tf3*t24) / t23
!
!        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
!     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)
!
!        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
!      EndIf
!
!      PartonX6 = ff
!
!      Return
!C                                       ********************
!      End
!



C============================================================================
C                CTEQ Parton Distribution Functions: version 6.0-6.6
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C                             August 6, 2003, v6.11
C                             December 12, 2004, v6.12
C                             December 4, 2006, v6.5 (CTEQ6.5M series added)
C                             March 23, 2007, v6.51 (CTEQ6.5S/C series added)
C                             April 24, 2007, v6.52 (minor improvement)
C                             March 30, 2008, v6.6 
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
C       JHEP 0310:046(2003), hep-ph/0303013
C
C   Ref[3]: "Neutrino dimuon Production and Strangeness Asymmetry of the Nucleon"
C       By: F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky, H.L. Lai, S. Kretzer, J.F. Owens, W.K. Tung
C       Eur. Phys. J. C40:145(2005), hep-ph/0312323
C
C   Ref[4]: "CTEQ6 Parton Distributions with Heavy Quark Mass Effects"
C       By: S. Kretzer, H.L. Lai, F. Olness, W.K. Tung
C       Phys. Rev. D69:114005(2004), hep-ph/0307022
C
C   Ref[5]: "Heavy Quark Mass Effects in Deep Inelastic Scattering and Global QCD Analysis"
C       By : W.K. Tung, H.L. Lai, A. Belyaev, J. Pumplin, D. Stump, C.-P. Yuan
C       JHEP 0702:053(2007), hep-ph/0611254
C
C   Ref[6]: "The Strange Parton Distribution of Nucleon: Global Analysis and Applications"
C       By : H.L. Lai, P. Nadolsky, J. Pumplin, D. Stump, W.K. Tung, C.-P. Yuan
C       JHEP 0704:089,2007, hep-ph/0702268
C
C   Ref[7]: "The Charm Content of the Nucleon"
C       By : J. Pumplin, H.L. Lai, W.K. Tung
C       Phys.Rev.D75:054029,2007, hep-ph/0701220

C   Ref[8]: "Implications of CTEQ global analysis for collider observables"
C       By : P. M. Nadolsky, H.-L. Lai, Q.-H. Cao, J. Huston, J. Pumplin, D. R. Stump, W.-K. Tung, C.-P. Yuan
C       arXiv:0802.0007 [hep-ph], submitted to Phys. Rev. D. 
C

C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
C   (4) 5 special sets for strangeness study from Ref[3].
C   (5) 1 special set for heavy quark study from Ref[4].
C   (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref[5].
C   (7) 8 sets of PDFs resulting from the strangeness study, Ref[6].
C   (8) 7 sets of PDFs resulting from the charm study, Ref[7].
C   (9) CTEQ6.6M and its 44 up/down eigenvector sets from Ref[8].
C  (10) Fits with nonperturbative charm from the study in  Ref[8].
C  (11) Fits with alternative values of the strong coupling strength from the study in Ref[8].


C  Details about the calling convention are:
C --------------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File   Ref
C ================================================================================
C Standard, "best-fit", sets:                 
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl    [1]
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl    [1]
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl    [1]
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl   [1]
C 200    CTEQ6.1M: updated CTEQ6M (see below, under "uncertainty" section)     [2]
C 400    CTEQ6.6M; the 2008 set (see below, under "uncertainty" section)       [8]
C
C --------------------------
C  Special sets with nonperturbative charm at Q_0=1.3 GeV from Ref [8]
C --------------------------
C 450    CTEQ6.6C1   BHPS model for IC    0.118     326   226    ctq66.c1.pds
C 451    CTEQ6.6C2   BHPS model for IC    0.118     326   226    ctq66.c2.pds
C 452    CTEQ6.6C3   Sea-like model       0.118     326   226    ctq66.c3.pds
C 453    CTEQ6.6C4   Sea-like model       0.118     326   226    ctq66.c4.pds
C     Momentum Fraction carried by c+cbar=2c at Q0=1.3 GeV:
C    Iset:     451  452   453   454 
C Mom. frac:  0.01 0.035  0.01  0.035


C --------------------------
C  Special CTEQ6.6 sets with alternative values of strong coupling strength [8]
C --------------------------
C 460    CTEQ6.6A1                        0.125           328    ctq66.a1.pds
C 461    CTEQ6.6A2                        0.122           281    ctq66.a2.pds
C 462    CTEQ6.6A3                        0.114           179    ctq66.a3.pds
C 463    CTEQ6.6A4                        0.112           159    ctq66.a4.pds

C --------------------------
C Special sets for strangeness study:  Ref.[3]
C --------------------------
C  11    CTEQ6A   Class A                 0.118     326   226    cteq6sa.pds
C  12    CTEQ6B   Class B                 0.118     326   226    cteq6sb.pds
C  13    CTEQ6C   Class C                 0.118     326   226    cteq6sc.pds
C  14    CTEQ6B+  Large [S-]              0.118     326   226    cteq6sb+.pds
C  15    CTEQ6B-  Negative [S-]           0.118     326   226    cteq6sb-.pds
C --------------------------
C Special set for Heavy Quark study:   Ref.[4]
C --------------------------
C  21    CTEQ6HQ                          0.118     326   226    cteq6hq.pds
C --------------------------
C Released sets for strangeness study:  Ref.[6]
C -------------------------- s=sbr
C  30    CTEQ6.5S0   Best-fit             0.118     326   226    ctq65.s+0.pds
C  31    CTEQ6.5S1   Low s+               0.118     326   226    ctq65.s+1.pds
C  32    CTEQ6.5S2   High s+              0.118     326   226    ctq65.s+2.pds
C  33    CTEQ6.5S3   Alt Low s+           0.118     326   226    ctq65.s+3.pds
C  34    CTEQ6.5S4   Alt High s+          0.118     326   226    ctq65.s+4.pds
C -------------------------- s!=sbr
C          strangeness asymmetry <x>_s-
C  35    CTEQ6.5S-0  Best-fit    0.0014    0.118     326   226    ctq65.s-0.pds
C  36    CTEQ6.5S-1  Low        -0.0010    0.118     326   226    ctq65.s-1.pds
C  37    CTEQ6.5S-2  High        0.0050    0.118     326   226    ctq65.s-2.pds
C --------------------------
C Released sets for charm study:  Ref.[7]
C --------------------------
C  40    CTEQ6.5C0   no intrinsic charm   0.118     326   226    ctq65.c0.pds
C  41    CTEQ6.5C1   BHPS model for IC    0.118     326   226    ctq65.c1.pds
C  42    CTEQ6.5C2   BHPS model for IC    0.118     326   226    ctq65.c2.pds
C  43    CTEQ6.5C3   Meson cloud model    0.118     326   226    ctq65.c3.pds
C  44    CTEQ6.5C4   Meson cloud model    0.118     326   226    ctq65.c4.pds
C  45    CTEQ6.5C5   Sea-like model       0.118     326   226    ctq65.c5.pds
C  46    CTEQ6.5C6   Sea-like model       0.118     326   226    ctq65.c6.pds
C
C     Momentum Fraction carried by c,cbar at Q0=1.3 GeV:
C    Iset:charm  ,cbar     | Iset:charm  ,cbar     | Iset:charm  ,cbar
C    41: 0.002857,0.002857 | 43: 0.003755,0.004817 | 45: 0.005714,0.005714
C    42: 0.010000,0.010000 | 44: 0.007259,0.009312 | 46: 0.012285,0.012285
C
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects, Ref[5]:  central fit: CTEQ6.5M (=CTEQ65.00)
C                              -----------------------
C  3xx  CTEQ65.xx  +/- sets               0.118     326   226    ctq65.xx.pds
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 300      is CTEQ65.00 (=CTEQ6.5M),
C             301/302 are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects and free strangeness, Ref[8]:  
C                central fit: CTEQ6.6M (=CTEQ66.00)
C                              -----------------------
C  4xx  CTEQ66.xx  +/- sets               0.118     326   226    ctq66.xx.pds
C        where xx = 01-44: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 400      is CTEQ66.00 (=CTEQ6.6M),
C             401/402 are CTEQ66.01/02, +/- sets of 1st eigenvector, ... etc.

C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for 
C    *  10^-8 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.6 series;
C    *  10^-7 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.5S/C series;
C    *  10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV) for CTEQ6, CTEQ6.1 series;
C
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   nadolsky@pa.msu.edu, pumplin@pa.msu.edu or tung@pa.msu.edu.
C
C===========================================================================

      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0d0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Ctq6Pdf = 0D0
        Return
      Endif

      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif

      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if (Ctq6Pdf.lt.0.D0) Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=8)
!CM      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Character Flnm(Isetmax0)*6, nn*3,nset2*2, Tablefile*40
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.', 'ctq66.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
      Data Isetmin4,Isetmax4 /400,444/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetmnC08,IsetmxC08 /450,453/
      Data IsetmnAS08,IsetmxAS08 /460,463/

      Data IsetHQ /21/
      Common /Setchange/ Isetch
!CM
      character *50 prefix,prefix1
      integer nset
      common/prefix/nset,prefix
!CM
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
        fmtpds=.true.

        If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                  Iset = 1,2,3 for 6m, 6d, 6l
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'.tbl'
        Elseif (Iset.eq.4) Then
C                                                             4  (2nd LO fit)
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'1.tbl'
        Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(1)//nn//'.tbl'
        Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(5)//nn(2:3)//'.tbl'
        Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15
          If(Iset.eq.11) then
            Tablefile=Flnm(6)//'a.pds'
          Elseif(Iset.eq.12) then
            Tablefile=Flnm(6)//'b.pds'
          Elseif(Iset.eq.13) then
            Tablefile=Flnm(6)//'c.pds'
          Elseif(Iset.eq.14) then
            Tablefile=Flnm(6)//'b+.pds'
          Elseif(Iset.eq.15) then
            Tablefile=Flnm(6)//'b-.pds'
          Endif
        Elseif (Iset.eq.IsetHQ) Then
C                                                               21
          TableFile='cteq6hq.pds'
        Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'s+'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37
          Is = Iset - 5
          write(nn,'(I2)') Is
          Tablefile=Flnm(7)//'s-'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'c'//nn(2:2)//'.pds'
        Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340
          write(nn,'(I3)') Iset
          Tablefile=Flnm(7)//nn(2:3)//'.pds'
        Elseif (Iset.ge.Isetmin4 .and. Iset.le.Isetmax4) Then
C                                                    (Cteq6.6)  400 - 444   
!CM          write(nn,'(I3)') Iset
          write(nset2,'(I2.2)') nset
!CM          Tablefile=Flnm(8)//nn(2:3)//'.pds'
          Tablefile=Flnm(8)//nset2//'.pds'
        Elseif (Iset.ge.IsetmnC08 .and. Iset.le.IsetmxC08) Then
C                                                   (Cteq6.6C)  450 - 453
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.ge.IsetmnAS08 .and. Iset.le.IsetmxAS08) Then
C                                                   (Cteq6.6AS)  460 - 463
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'a'//nn(3:3)//'.pds'
        Else
          Print *, 'Invalid Iset number in SetCtq6 :', Iset
          Stop
        Endif
        IU= NextUn()
        Open(IU, File='Pdfdata/'//Tablefile, Status='OLD', Err=100)
 21     Call Readpds (IU,fmtpds)
        Close (IU)
        Isetold=Iset
        Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >  //'in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C



C----------------------------------------------------------------------
C--   Fortran interpolation code for MSTW PDFs, building on existing
C--   MRST Fortran code and Jeppe Andersen's C++ code.
C--   Three user interfaces:
C--    call GetAllPDFs(prefix,ih,x,q,upv,dnv,usea,dsea,
C--                    str,sbar,chm,cbar,bot,bbar,glu,phot)
C--    call GetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)
C--    xf = GetOnePDF(prefix,ih,x,q,f)
C--   See enclosed example.f for usage.
C--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>.
C----------------------------------------------------------------------

C----------------------------------------------------------------------

C--   Traditional MRST-like interface: return all flavours.
C--   (Note the additional "sbar", "cbar", "bbar" and "phot"
C--   compared to previous MRST releases.)
      subroutine GetAllPDFs(prefix,ih,x,q,
     &     upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot)
      implicit none
      integer ih
      double precision x,q,upv,dnv,usea,dsea,str,sbar,chm,cbar,
     &     bot,bbar,glu,phot,GetOnePDF,up,dn,sv,cv,bv
      character*(*) prefix

C--   Quarks.
      dn  = GetOnePDF(prefix,ih,x,q,1)
      up  = GetOnePDF(prefix,ih,x,q,2)
      str = GetOnePDF(prefix,ih,x,q,3)
      chm = GetOnePDF(prefix,ih,x,q,4)
      bot = GetOnePDF(prefix,ih,x,q,5)

C--   Valence quarks.
      dnv = GetOnePDF(prefix,ih,x,q,7)
      upv = GetOnePDF(prefix,ih,x,q,8)
      sv  = GetOnePDF(prefix,ih,x,q,9)
      cv  = GetOnePDF(prefix,ih,x,q,10)
      bv  = GetOnePDF(prefix,ih,x,q,11)
      
C--   Antiquarks = quarks - valence quarks.
      dsea = dn - dnv
      usea = up - upv
      sbar = str - sv
      cbar = chm - cv
      bbar = bot - bv

C--   Gluon.
      glu = GetOnePDF(prefix,ih,x,q,0)

C--   Photon (= zero unless considering QED contributions).
      phot = GetOnePDF(prefix,ih,x,q,13)

      return
      end

C----------------------------------------------------------------------

C--   Alternative LHAPDF-like interface: return PDFs in an array.
      subroutine GetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)
      implicit none
      integer ih,f
      double precision x,q,xpdf(-6:6),xphoton,xvalence,GetOnePDF
      character*(*) prefix

      do f = 1, 6
         xpdf(f) = GetOnePDF(prefix,ih,x,q,f) ! quarks
         xvalence = GetOnePDF(prefix,ih,x,q,f+6) ! valence quarks
         xpdf(-f) = xpdf(f) - xvalence ! antiquarks
      end do
      xpdf(0) = GetOnePDF(prefix,ih,x,q,0) ! gluon
      xphoton = GetOnePDF(prefix,ih,x,q,13) ! photon
      
      return
      end

C----------------------------------------------------------------------

C--   Get only one parton flavour 'f', using PDG notation:
C--    f =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
C--      = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
C--   Can also get valence quarks directly:
C--    f =  7, 8, 9,10,11,12.
C--      = dv,uv,sv,cv,bv,tv.
C--   Photon: f = 13.
      double precision function GetOnePDF(prefix,ih,x,q,f)
      implicit none
      logical warn,fatal
      parameter(warn=.false.,fatal=.true.)
C--   Set warn=.true. to turn on warnings when extrapolating.
C--   Set fatal=.false. to return zero instead of terminating when
C--    invalid input values of x and q are used.
      integer ih,f,nhess,nx,nq,np,nqc0,nqb0,nqc,nqb,n,m,ip,io,
     &     alphaSorder,nExtraFlavours
      double precision x,q,xmin,xmax,qsqmin,qsqmax,mc2,mb2,eps,
     &     dummy,qsq,xlog,qsqlog,res,res1,anom,ExtrapolatePDF,
     &     InterpolatePDF,distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      parameter(nx=64,nq=48,np=12,nqc0=4,nqb0=14,
     &     nqc=nq-nqc0,nqb=nq-nqb0)
      parameter(xmin=1d-6,xmax=1d0,qsqmin=1d0,qsqmax=1d9,eps=1d-6)
      parameter(nhess=2*20)
      character set*2,prefix*(*),filename*60,oldprefix(0:nhess)*50
      character dummyChar,dummyWord*50
      double precision ff(np,nx,nq)
      double precision qq(nq),xx(nx),cc(np,0:nhess,nx,nq,4,4)
      double precision xxl(nx),qql(nq)
C--   Store distance along each eigenvector, tolerance,
C--   heavy quark masses and alphaS parameters in COMMON block.
      common/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder
      save
      data xx/1d-6,2d-6,4d-6,6d-6,8d-6,
     &     1d-5,2d-5,4d-5,6d-5,8d-5,
     &     1d-4,2d-4,4d-4,6d-4,8d-4,
     &     1d-3,2d-3,4d-3,6d-3,8d-3,
     &     1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     &	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     &	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     &	   .5d0,.525d0,.55d0,.575d0,.6d0,.625d0,.65d0,.675d0,
     &     .7d0,.725d0,.75d0,.775d0,.8d0,.825d0,.85d0,.875d0,
     &     .9d0,.925d0,.95d0,.975d0,1d0/
      data qq/1.d0,
     &     1.25d0,1.5d0,0.d0,0.d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,
     &     1d1,1.2d1,0.d0,0.d0,2.6d1,4d1,6.4d1,1d2,
     &     1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     &     1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     &     1.8d6,3.2d6,5.6d6,1d7,1.8d7,3.2d7,5.6d7,1d8,
     &     1.8d8,3.2d8,5.6d8,1d9/

      if (f.lt.-6.or.f.gt.13) then
         print *,"Error: invalid parton flavour = ",f
         stop
      end if

      if (ih.lt.0.or.ih.gt.nhess) then
         print *,"Error: invalid eigenvector number = ",ih
         stop
      end if

C--   Check if the requested parton set is already in memory.
      if (oldprefix(ih).ne.prefix) then

C--   Start of initialisation for eigenvector set "i" ...
C--   Do this only the first time the set "i" is called,
C--   OR if the prefix has changed from the last time.

C--   Check that the character arrays "oldprefix" and "filename"
C--   are large enough.
         if (len_trim(prefix).gt.len(oldprefix(ih))) then
            print *,"Error in GetOnePDF: increase size of oldprefix"
            stop
         else if (len_trim(prefix)+7.gt.len(filename)) then
            print *,"Error in GetOnePDF: increase size of filename"
            stop
         end if

         write(set,'(I2.2)') ih  ! convert integer to string
C--   Remove trailing blanks from prefix before assigning filename.
         filename = prefix(1:len_trim(prefix))//'.'//set//'.dat'
C--   Line below can be commented out if you don't want this message.
         print *,"Reading PDF grid from ",filename(1:len_trim(filename))
         open(unit=33,file=filename,iostat=io,status='old')
         if (io.ne.0) then
            print *,"Error in GetOnePDF: can't open ",
     &           filename(1:len_trim(filename))
            stop
         end if

C--   Read header containing heavy quark masses and alphaS values.
         read(33,*) 
         read(33,*)
         read(33,*) dummyChar,dummyWord,dummyWord,dummyChar,
     &        distance,tolerance
         read(33,*) dummyChar,dummyWord,dummyChar,mCharm
         read(33,*) dummyChar,dummyWord,dummyChar,mBottom
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSQ0
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSMZ
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSorder
         read(33,*) dummyChar,dummyWord,dummyChar,nExtraFlavours
         read(33,*)
         read(33,*)
         mc2=mCharm**2
         mb2=mBottom**2
         qq(4)=mc2
         qq(5)=mc2+eps
         qq(14)=mb2
         qq(15)=mb2+eps

C--   Check that the heavy quark masses are sensible.
         if (mc2.lt.qq(3).or.mc2.gt.qq(6)) then
            print *,"Error in GetOnePDF: invalid mCharm = ",mCharm
            stop
         end if
         if (mb2.lt.qq(13).or.mb2.gt.qq(16)) then
            print *,"Error in GetOnePDF: invalid mBottom = ",mBottom
            stop
         end if
         
C--   The nExtraFlavours variable is provided to aid compatibility
C--   with future grids where, for example, a photon distribution
C--   might be provided (cf. the MRST2004QED PDFs).
         if (nExtraFlavours.lt.0.or.nExtraFlavours.gt.1) then
            print *,"Error in GetOnePDF: invalid nExtraFlavours = ",
     &           nExtraFlavours
            stop
         end if

C--   Now read in the grids from the grid file.
         do n=1,nx-1
            do m=1,nq
               if (nExtraFlavours.gt.0) then
                  if (alphaSorder.eq.2) then ! NNLO
                     read(33,'(12(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,12)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     read(33,'(10(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9),ff(12,n,m)
                  end if
               else             ! nExtraFlavours = 0
                  if (alphaSorder.eq.2) then ! NNLO
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(11(1pe12.4))',iostat=io)
     &                 (ff(ip,n,m),ip=1,11)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(9(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9)
                  end if
               end if
               if (io.ne.0) then
                  print *,"Error in GetOnePDF reading ",filename
                  stop
               end if
            enddo
         enddo

C--   Check that ALL the file contents have been read in.
         read(33,*,iostat=io) dummy
         if (io.eq.0) then
            print *,"Error in GetOnePDF: not at end of ",filename
            stop
         end if
         close(unit=33)

C--   PDFs are identically zero at x = 1.
         do m=1,nq
            do ip=1,np
               ff(ip,nx,m)=0d0
            enddo
         enddo

         do n=1,nx
            xxl(n)=log10(xx(n))
         enddo
         do m=1,nq
            qql(m)=log10(qq(m))
         enddo

C--   Initialise all parton flavours.
         do ip=1,np
            call InitialisePDF(ip,np,ih,nhess,nx,nq,nqc0,nqb0,
     &           xxl,qql,ff,cc)
         enddo

         oldprefix(ih) = prefix

C--   ... End of initialisation for eigenvector set "ih".

      end if                    ! oldprefix(ih).ne.prefix

C----------------------------------------------------------------------

      qsq=q*q
C--   If mc2 < qsq < mc2+eps, then qsq = mc2+eps.
      if (qsq.gt.qq(nqc0).and.qsq.lt.qq(nqc0+1)) qsq = qq(nqc0+1)
C--   If mb2 < qsq < mb2+eps, then qsq = mb2+eps.
      if (qsq.gt.qq(nqb0).and.qsq.lt.qq(nqb0+1)) qsq = qq(nqb0+1)
      
      xlog=log10(x)
      qsqlog=log10(qsq)

      res = 0.d0

      if (f.eq.0) then          ! gluon
         ip = 1
      else if (f.ge.1.and.f.le.5) then ! quarks
         ip = f+1
      else if (f.le.-1.and.f.ge.-5) then ! antiquarks
         ip = -f+1
      else if (f.ge.7.and.f.le.11) then ! valence quarks
         ip = f
      else if (f.eq.13) then    ! photon
         ip = 12
      else if (abs(f).ne.6.and.f.ne.12) then
         if (warn.or.fatal) print *,"Error in GetOnePDF: f = ",f
         if (fatal) stop
      end if
      
      if (x.le.0.d0.or.x.gt.xmax.or.q.le.0.d0) then

         if (warn.or.fatal) print *,"Error in GetOnePDF: x,qsq = ",
     &        x,qsq
         if (fatal) stop

      else if (abs(f).eq.6.or.f.eq.12) then ! set top quarks to zero
         
         res = 0.d0

      else if (qsq.lt.qsqmin) then ! extrapolate to low Q^2

         if (warn) then
            print *, "Warning in GetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         if (x.lt.xmin) then    ! extrapolate to low x

            res = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         else                   ! do usual interpolation
            
            res = InterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = InterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         end if

C--   Calculate the anomalous dimension, dlog(xf)/dlog(qsq),
C--   evaluated at qsqmin.  Then extrapolate the PDFs to low
C--   qsq < qsqmin by interpolating the anomalous dimenion between
C--   the value at qsqmin and a value of 1 for qsq << qsqmin.
C--   If value of PDF at qsqmin is very small, just set
C--   anomalous dimension to 1 to prevent rounding errors.
         if (abs(res).ge.1.D-5) then
            anom = (res1-res)/res/0.01D0
         else
            anom = 1.D0
         end if
         res = res*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin)

      else if (x.lt.xmin.or.qsq.gt.qsqmax) then ! extrapolate

         if (warn) then
            print *, "Warning in GetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         res = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)
         
         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if

      else                      ! do usual interpolation
         
         res = InterpolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)

         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if
            
      end if
      
      GetOnePDF = res

      return
      end

C----------------------------------------------------------------------

      subroutine InitialisePDF(ip,np,ih,nhess,nx,my,myc0,myb0,
     &     xx,yy,ff,cc)
      implicit none
      integer nhess,ih,nx,my,myc0,myb0,j,k,l,m,n,ip,np
      double precision xx(nx),yy(my),ff(np,nx,my),
     &     ff1(nx,my),ff2(nx,my),ff12(nx,my),ff21(nx,my),
     &     yy0(4),yy1(4),yy2(4),yy12(4),z(16),
     &     cl(16),cc(np,0:nhess,nx,my,4,4),iwt(16,16),
     &     polderiv1,polderiv2,polderiv3,d1,d2,d1d2,xxd

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     &     -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     &     2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     &     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     &     0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     &     0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     &     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     &     9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     &     -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     &     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     &     -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     &     4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/

      do m=1,my
         ff1(1,m)=polderiv1(xx(1),xx(2),xx(3),
     &        ff(ip,1,m),ff(ip,2,m),ff(ip,3,m))
         ff1(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff(ip,nx-2,m),ff(ip,nx-1,m),ff(ip,nx,m))
         do n=2,nx-1
            ff1(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff(ip,n-1,m),ff(ip,n,m),ff(ip,n+1,m))
         enddo
      enddo

C--   Calculate the derivatives at qsq=mc2,mc2+eps,mb2,mb2+eps
C--   in a similar way as at the endpoints qsqmin and qsqmax.
      do n=1,nx
         do m=1,my
            if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff2(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff(ip,n,m),ff(ip,n,m+1),ff(ip,n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff2(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff(ip,n,m-2),ff(ip,n,m-1),ff(ip,n,m))
            else
               ff2(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff(ip,n,m-1),ff(ip,n,m),ff(ip,n,m+1))
            end if
         end do
      end do

C--   Calculate the cross derivatives (d/dx)(d/dy).
      do m=1,my
         ff12(1,m)=polderiv1(xx(1),xx(2),xx(3),
     &        ff2(1,m),ff2(2,m),ff2(3,m))
         ff12(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff2(nx-2,m),ff2(nx-1,m),ff2(nx,m))
         do n=2,nx-1
            ff12(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff2(n-1,m),ff2(n,m),ff2(n+1,m))
         enddo
      enddo

C--   Calculate the cross derivatives (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff21(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff1(n,m),ff1(n,m+1),ff1(n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff21(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff1(n,m-2),ff1(n,m-1),ff1(n,m))
            else
               ff21(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff1(n,m-1),ff1(n,m),ff1(n,m+1))
            end if
         end do
      end do

C--   Take the average of (d/dx)(d/dy) and (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            ff12(n,m)=0.5*(ff12(n,m)+ff21(n,m))
         end do
      end do

      do n=1,nx-1
         do m=1,my-1
            d1=xx(n+1)-xx(n)
            d2=yy(m+1)-yy(m)
            d1d2=d1*d2
            
            yy0(1)=ff(ip,n,m)
            yy0(2)=ff(ip,n+1,m)
            yy0(3)=ff(ip,n+1,m+1)
            yy0(4)=ff(ip,n,m+1)
            
            yy1(1)=ff1(n,m)
            yy1(2)=ff1(n+1,m)
            yy1(3)=ff1(n+1,m+1)
            yy1(4)=ff1(n,m+1)
            
            yy2(1)=ff2(n,m)
            yy2(2)=ff2(n+1,m)
            yy2(3)=ff2(n+1,m+1)
            yy2(4)=ff2(n,m+1)
            
            yy12(1)=ff12(n,m)
            yy12(2)=ff12(n+1,m)
            yy12(3)=ff12(n+1,m+1)
            yy12(4)=ff12(n,m+1)
            
            do k=1,4
               z(k)=yy0(k)
               z(k+4)=yy1(k)*d1
               z(k+8)=yy2(k)*d2
               z(k+12)=yy12(k)*d1d2
            enddo
            
            do l=1,16
               xxd=0.d0
               do k=1,16
                  xxd=xxd+iwt(k,l)*z(k)
               enddo
               cl(l)=xxd
            enddo
            l=0
            do k=1,4
               do j=1,4
                  l=l+1
                  cc(ip,ih,n,m,k,j)=cl(l)
               enddo
            enddo
         enddo
      enddo
      return
      end

C----------------------------------------------------------------------

      double precision function InterpolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,locx,l,m,n,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,t,u

      n=locx(xx,nx,x)
      m=locx(yy,my,y)
      
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(ip,ih,n,m,l,4)*u+cc(ip,ih,n,m,l,3))*u
     .        +cc(ip,ih,n,m,l,2))*u+cc(ip,ih,n,m,l,1)
      enddo

      InterpolatePDF = z

      return
      end

C----------------------------------------------------------------------

      double precision function ExtrapolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,locx,n,m,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,f0,f1,z0,z1,InterpolatePDF
      
      n=locx(xx,nx,x)           ! 0: below xmin, nx: above xmax
      m=locx(yy,my,y)           ! 0: below qsqmin, my: above qsqmax
      
C--   If extrapolation in small x only:
      if (n.eq.0.and.m.gt.0.and.m.lt.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(1),y,nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(2),y,nx,my,xx,yy,cc)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
C--   If extrapolation into large q only:
      else if (n.gt.0.and.m.eq.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,x,yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,x,yy(my-1),nx,my,xx,yy,cc)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
C--   If extrapolation into large q AND small x:
      else if (n.eq.0.and.m.eq.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(1),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(1),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(2),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(2),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.0.d0.and.z1.gt.0.d0) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if
      else
         print *,"Error in ExtrapolatePDF"
         stop
      end if

      ExtrapolatePDF = z      

      return
      end

C----------------------------------------------------------------------

      integer function locx(xx,nx,x)
C--   returns an integer j such that x lies inbetween xx(j) and xx(j+1).
C--   nx is the length of the array with xx(nx) the highest element.
      implicit none
      integer nx,jl,ju,jm
      double precision x,xx(nx)
      if(x.eq.xx(1)) then
         locx=1
         return
      endif
      if(x.eq.xx(nx)) then
         locx=nx-1  
         return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
         jl=jm
      else
         ju=jm
      endif
      go to 1
    2 locx=jl
      return
      end

C----------------------------------------------------------------------

      double precision function polderiv1(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x1 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv1=(x3*x3*(y1-y2)+2.d0*x1*(x3*(-y1+y2)+x2*(y1-y3))
     &     +x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv2(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x2 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv2=(x3*x3*(y1-y2)-2.d0*x2*(x3*(y1-y2)+x1*(y2-y3))
     &     +x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv3(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x3 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv3=(x3*x3*(-y1+y2)+2.d0*x2*x3*(y1-y3)+x1*x1*(y2-y3)
     &     +x2*x2*(-y1+y3)+2.d0*x1*x3*(-y2+y3))/
     &     ((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

C----------------------------------------------------------------------



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLO (To be published)
!!          (See also: Phys. Rev. D79 (2009) 074023 and arXiv:0810.4274)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This package contains the JR09 VFNS NNLO(MSbar) dynamical parton
!! distributions of the nucleon and their associated exact iterative
!! solutions for alpha_s.
!!
!! The sets resulting from displacements in the parameter space with
!! respect to the central value of the NNLO(MSbar) fit along the plus
!! (minus) directions of the (rescaled) eigenvectors of the hessian
!! matrix are included. The tolerance parameter for these displacements
!! was chosen to be T = 4.535 for a total of 1568 fitted points. Since
!! alpha_s was included as a free parameter in the error estimation the
!! use of the provided alpha_s solution for each set is mandatory for
!! uncertainty studies (the difference on alpha_s for different
!! eigenvector sets can be as large as 10% at low scales).
!!
!! The grids are generated for 10^-9 <= x <= 1 and Qo^2 <= Q^2 <= 10^8
!! (GeV^2) where Qo^2 = 0.55 GeV^2 for the NNLO distributions. Outside
!! these ranges the output is either set to NaN or obtained through
!! extrapolation and should NOT be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The routines use a modification of the standard multidimensional
!! linear interpolation routine FINT (CERNLIB E104) distributed as the
!! file 'dfint.f'.
!! The file './JR09VFNNLO.grd', where ./ means the path from the working
!! directory to the file, is read.
!! For questions, comments, problems etc please contact:
!! pjimenez@het.physik.uni-dortmund.de
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLOinit:
!!  Initialization routine of the package to be called (only once)
!!  before using any of the other subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR09VFNNLOx'parton'(x,Q2,set) with 'parton' = uv,dv,gl,ub,db,sb,cb,bb:
!!      Parton distribution 'parton' (times x).
!!        x == Bjorken-x.
!!        Q2 == Q**2 (GeV**2) == Factorization scale
!!                            == Renormalization scale.
!!        set == set to be used.
!! JR09VFNNLOalphas(Q2,set):
!!      Value of alpha_s (no additional 2pi or 4pi factors).
!!        Q2 == Q**2 (GeV**2) == Renormalization scale.
!!        set == set to be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set == Index specifying the set to be used:
!!        0 NNLO(MSbar).
!!        1,2,...,13 set corresponding to a displacement +T with respect
!!                 to the set 0 in the direction of the ith eigenvector.
!!       -1,-2,...,-13 the same for displacements -T.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      block data JR09VFNNLO
       implicit none
       integer shape(2)
       double precision grid(217)
       common /JR09VFNNLOgrid/ grid,shape
       data shape /118,99/
       data grid
     &   /1d-9,1.25d-9,1.6d-9,2d-9,2.5d-9,3.16d-9,4d-9,5d-9,6.3d-9,8d-9,
     &    1d-8,1.25d-8,1.6d-8,2d-8,2.5d-8,3.16d-8,4d-8,5d-8,6.3d-8,8d-8,
     &    1d-7,1.25d-7,1.6d-7,2d-7,2.5d-7,3.16d-7,4d-7,5d-7,6.3d-7,8d-7,
     &    1d-6,1.25d-6,1.6d-6,2d-6,2.5d-6,3.16d-6,4d-6,5d-6,6.3d-6,8d-6,
     &    1d-5,1.25d-5,1.6d-5,2d-5,2.5d-5,3.16d-5,4d-5,5d-5,6.3d-5,8d-5,
     &    1d-4,1.25d-4,1.6d-4,2d-4,2.5d-4,3.16d-4,4d-4,5d-4,6.3d-4,8d-4,
     &    1d-3,1.25d-3,1.6d-3,2d-3,2.5d-3,3.16d-3,4d-3,5d-3,6.3d-3,8d-3,
     &    1d-2,1.25d-2,1.6d-2,2d-2,2.5d-2,3.16d-2,4d-2,5d-2,6.3d-2,8d-2,
     &    0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0,
     &    0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0,
     &    0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0,
     &    0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0,
     &    0.9d0,0.920d0,0.94d0,0.960d0,0.98d0,1d0,
     &    0.3d0,0.31d0,0.35d0,0.375d0,0.4d0,0.45d0,0.5d0,0.51d0,0.525d0,
     &    0.55d0,0.575d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0,
     &    1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0,
     &    1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1,
     &    1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2,
     &    1d3,1.25d3,1.6d3,2d3,2.5d3,3.16d3,4d3,5d3,6.3d3,8d3,
     &    1d4,1.25d4,1.6d4,2d4,2.5d4,3.16d4,4d4,5d4,6.3d4,8d4,
     &    1d5,1.25d5,1.6d5,2d5,2.5d5,3.16d5,4d5,5d5,6.3d5,8d5,
     &    1d6,1.25d6,1.6d6,2d6,2.5d6,3.16d6,4d6,5d6,6.3d6,8d6,
     &    1d7,1.25d7,1.6d7,2d7,2.5d7,3.16d7,4d7,5d7,6.3d7,8d7,1d8/
      end block data JR09VFNNLO

      subroutine JR09VFNNLOinit
       implicit none
       integer shape(2),i,j,k
       double precision NaN,grid(217),
     &                  alphas(99,-13:13),
     &                  xuv(118,99,-13:13),
     &                  xdv(118,99,-13:13),
     &                  xgl(118,99,-13:13),
     &                  xub(118,99,-13:13),
     &                  xdb(118,99,-13:13),
     &                  xsb(118,99,-13:13),
     &                  xcb(118,99,-13:13),
     &                  xbb(118,99,-13:13)
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOalphasc/ alphas
       common /JR09VFNNLOxuvc/ xuv
       common /JR09VFNNLOxdvc/ xdv
       common /JR09VFNNLOxglc/ xgl
       common /JR09VFNNLOxubc/ xub
       common /JR09VFNNLOxdbc/ xdb
       common /JR09VFNNLOxsbc/ xsb
       common /JR09VFNNLOxcbc/ xcb
       common /JR09VFNNLOxbbc/ xbb
       NaN=0d0
       NaN=0d0/NaN
       open(10,file='Pdfdata/JR09VFNNLO.grd')
       do 1 k=-13,13
          do 2 j=1,99
             read(10,*) alphas(j,k)
    2     continue
    1  continue
       do 10 k=-13,13
          do 11 j=1,99
             do 12 i=1,118
                read(10,*) xuv(i,j,k)
   12        continue
   11     continue
   10  continue
       do 20 k=-13,13
          do 21 j=1,99
             do 22 i=1,118
                read(10,*) xdv(i,j,k)
   22        continue
   21     continue
   20  continue
       do 30 k=-13,13
          do 31 j=1,99
             do 32 i=1,118
                read(10,*) xgl(i,j,k)
   32        continue
   31     continue
   30  continue
       do 40 k=-13,13
          do 41 j=1,99
             do 42 i=1,118
                read(10,*) xub(i,j,k)
   42        continue
   41     continue
   40  continue
       do 50 k=-13,13
          do 51 j=1,99
             do 52 i=1,118
                read(10,*) xdb(i,j,k)
   52        continue
   51     continue
   50  continue
       do 60 k=-13,13
          do 61 j=1,99
             do 62 i=1,118
                read(10,*) xsb(i,j,k)
   62        continue
   61     continue
   60  continue
       do 70 k=-13,13
          do 71 j=1,99
             do 72 i=1,118
                read(10,*) xcb(i,j,k)
   72        continue
   71     continue
   70  continue
       do 80 k=-13,13
          do 81 j=1,99
             do 82 i=1,118
                read(10,*) xbb(i,j,k)
   82        continue
   81     continue
   80  continue
       close(10)
       do 1000 k=-13,13
          do 1001 j=1,9
             do 1002 i=1,118
                xuv(i,j,k)=NaN
                xdv(i,j,k)=NaN
                xgl(i,j,k)=NaN
                xub(i,j,k)=NaN
                xdb(i,j,k)=NaN
                xsb(i,j,k)=NaN
                xcb(i,j,k)=NaN
                xbb(i,j,k)=NaN
 1002        continue
 1001     continue
 1000  continue
       return
      end subroutine JR09VFNNLOinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOalphas(Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),alphas(99,-13:13),arg,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOalphasc/ alphas
       arg = Q2
       JR09VFNNLOalphas = dfint(1,arg,shape(2),grid(119),alphas(1,nset))
       return
      end function JR09VFNNLOalphas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxuv(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xuv(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxuvc/ xuv
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxuv = dfint(2,arg,shape,grid,xuv(1,1,nset))
       return
      end function JR09VFNNLOxuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxdv(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xdv(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxdvc/ xdv
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxdv = dfint(2,arg,shape,grid,xdv(1,1,nset))
       return
      end function JR09VFNNLOxdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxgl(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xgl(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxglc/ xgl
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxgl = dfint(2,arg,shape,grid,xgl(1,1,nset))
       return
      end function JR09VFNNLOxgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxub(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xub(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxubc/ xub
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxub = dfint(2,arg,shape,grid,xub(1,1,nset))
       return
      end function JR09VFNNLOxub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxdb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xdb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxdbc/ xdb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxdb = dfint(2,arg,shape,grid,xdb(1,1,nset))
       return
      end function JR09VFNNLOxdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxsb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xsb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxsbc/ xsb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxsb = dfint(2,arg,shape,grid,xsb(1,1,nset))
       return
      end function JR09VFNNLOxsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxcb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xcb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxcbc/ xcb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxcb = dfint(2,arg,shape,grid,xcb(1,1,nset))
       return
      end function JR09VFNNLOxcb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR09VFNNLOxbb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xbb(118,99,-13:13),arg(2),x,Q2,dfint
       common /JR09VFNNLOgrid/ grid,shape
       common /JR09VFNNLOxbbc/ xbb
       arg(1) = x
       arg(2) = Q2
       JR09VFNNLOxbb = dfint(2,arg,shape,grid,xbb(1,1,nset))
       return
      end function JR09VFNNLOxbb


      block data GJR08VFNS
       implicit none
       integer shape(2)
       double precision grid(217)
       common /GJR08VFNSgrid/ grid,shape
       data shape /118,99/
       data grid
     &   /1d-9,1.25d-9,1.6d-9,2d-9,2.5d-9,3.16d-9,4d-9,5d-9,6.3d-9,8d-9,
     &    1d-8,1.25d-8,1.6d-8,2d-8,2.5d-8,3.16d-8,4d-8,5d-8,6.3d-8,8d-8,
     &    1d-7,1.25d-7,1.6d-7,2d-7,2.5d-7,3.16d-7,4d-7,5d-7,6.3d-7,8d-7,
     &    1d-6,1.25d-6,1.6d-6,2d-6,2.5d-6,3.16d-6,4d-6,5d-6,6.3d-6,8d-6,
     &    1d-5,1.25d-5,1.6d-5,2d-5,2.5d-5,3.16d-5,4d-5,5d-5,6.3d-5,8d-5,
     &    1d-4,1.25d-4,1.6d-4,2d-4,2.5d-4,3.16d-4,4d-4,5d-4,6.3d-4,8d-4,
     &    1d-3,1.25d-3,1.6d-3,2d-3,2.5d-3,3.16d-3,4d-3,5d-3,6.3d-3,8d-3,
     &    1d-2,1.25d-2,1.6d-2,2d-2,2.5d-2,3.16d-2,4d-2,5d-2,6.3d-2,8d-2,
     &    0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0,
     &    0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0,
     &    0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0,
     &    0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0,
     &    0.9d0,0.920d0,0.94d0,0.960d0,0.98d0,1d0,
     &    0.3d0,0.31d0,0.35d0,0.375d0,0.4d0,0.45d0,0.5d0,0.51d0,0.525d0,
     &    0.55d0,0.575d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0,
     &    1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0,
     &    1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1,
     &    1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2,
     &    1d3,1.25d3,1.6d3,2d3,2.5d3,3.16d3,4d3,5d3,6.3d3,8d3,
     &    1d4,1.25d4,1.6d4,2d4,2.5d4,3.16d4,4d4,5d4,6.3d4,8d4,
     &    1d5,1.25d5,1.6d5,2d5,2.5d5,3.16d5,4d5,5d5,6.3d5,8d5,
     &    1d6,1.25d6,1.6d6,2d6,2.5d6,3.16d6,4d6,5d6,6.3d6,8d6,
     &    1d7,1.25d7,1.6d7,2d7,2.5d7,3.16d7,4d7,5d7,6.3d7,8d7,1d8/
      end block data GJR08VFNS

      subroutine GJR08VFNSinit
       implicit none
       integer shape(2),i,j,k
       double precision NaN,grid(217),
     &                  alphas(99,-13:14),
     &                  xuv(118,99,-13:14),
     &                  xdv(118,99,-13:14),
     &                  xgl(118,99,-13:14),
     &                  xub(118,99,-13:14),
     &                  xdb(118,99,-13:14),
     &                  xsb(118,99,-13:14),
     &                  xcb(118,99,-13:14),
     &                  xbb(118,99,-13:14)
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSalphasc/ alphas
       common /GJR08VFNSxuvc/ xuv
       common /GJR08VFNSxdvc/ xdv
       common /GJR08VFNSxglc/ xgl
       common /GJR08VFNSxubc/ xub
       common /GJR08VFNSxdbc/ xdb
       common /GJR08VFNSxsbc/ xsb
       common /GJR08VFNSxcbc/ xcb
       common /GJR08VFNSxbbc/ xbb
       NaN=0d0
       NaN=0d0/NaN
       open(10,file='Pdfdata/GJR08VFNS.grd')
       do 1 k=-13,14
          do 2 j=1,99
             read(10,*) alphas(j,k)
    2     continue
    1  continue
       do 10 k=-13,14
          do 11 j=1,99
             do 12 i=1,118
                read(10,*) xuv(i,j,k)
   12        continue
   11     continue
   10  continue
       do 20 k=-13,14
          do 21 j=1,99
             do 22 i=1,118
                read(10,*) xdv(i,j,k)
   22        continue
   21     continue
   20  continue
       do 30 k=-13,14
          do 31 j=1,99
             do 32 i=1,118
                read(10,*) xgl(i,j,k)
   32        continue
   31     continue
   30  continue
       do 40 k=-13,14
          do 41 j=1,99
             do 42 i=1,118
                read(10,*) xub(i,j,k)
   42        continue
   41     continue
   40  continue
       do 50 k=-13,14
          do 51 j=1,99
             do 52 i=1,118
                read(10,*) xdb(i,j,k)
   52        continue
   51     continue
   50  continue
       do 60 k=-13,14
          do 61 j=1,99
             do 62 i=1,118
                read(10,*) xsb(i,j,k)
   62        continue
   61     continue
   60  continue
       do 70 k=-13,14
          do 71 j=1,99
             do 72 i=1,118
                read(10,*) xcb(i,j,k)
   72        continue
   71     continue
   70  continue
       do 80 k=-13,14
          do 81 j=1,99
             do 82 i=1,118
                read(10,*) xbb(i,j,k)
   82        continue
   81     continue
   80  continue
       close(10)
       do 1000 k=-13,13
          do 1001 j=1,6
             do 1002 i=1,118
                xuv(i,j,k)=NaN
                xdv(i,j,k)=NaN
                xgl(i,j,k)=NaN
                xub(i,j,k)=NaN
                xdb(i,j,k)=NaN
                xsb(i,j,k)=NaN
                xcb(i,j,k)=NaN
                xbb(i,j,k)=NaN
 1002        continue
 1001     continue
 1000  continue
       return
      end subroutine GJR08VFNSinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSalphas(Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),alphas(99,-13:14),arg,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSalphasc/ alphas
       arg = Q2
       GJR08VFNSalphas = dfint(1,arg,shape(2),grid(119),alphas(1,nset))
       return
      end function GJR08VFNSalphas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxuv(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xuv(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxuvc/ xuv
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxuv = dfint(2,arg,shape,grid,xuv(1,1,nset))
       return
      end function GJR08VFNSxuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxdv(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xdv(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxdvc/ xdv
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxdv = dfint(2,arg,shape,grid,xdv(1,1,nset))
       return
      end function GJR08VFNSxdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxgl(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xgl(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxglc/ xgl
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxgl = dfint(2,arg,shape,grid,xgl(1,1,nset))
       return
      end function GJR08VFNSxgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxub(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xub(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxubc/ xub
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxub = dfint(2,arg,shape,grid,xub(1,1,nset))
       return
      end function GJR08VFNSxub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxdb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xdb(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxdbc/ xdb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxdb = dfint(2,arg,shape,grid,xdb(1,1,nset))
       return
      end function GJR08VFNSxdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxsb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xsb(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxsbc/ xsb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxsb = dfint(2,arg,shape,grid,xsb(1,1,nset))
       return
      end function GJR08VFNSxsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxcb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xcb(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxcbc/ xcb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxcb = dfint(2,arg,shape,grid,xcb(1,1,nset))
       return
      end function GJR08VFNSxcb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08VFNSxbb(x,Q2,nset)
       implicit none
       integer shape(2),nset
       double precision grid(217),xbb(118,99,-13:14),arg(2),x,Q2,dfint
       common /GJR08VFNSgrid/ grid,shape
       common /GJR08VFNSxbbc/ xbb
       arg(1) = x
       arg(2) = Q2
       GJR08VFNSxbb = dfint(2,arg,shape,grid,xbb(1,1,nset))
       return
      end function GJR08VFNSxbb



!! CERNLIB E104 modified to be used with (G)JR GRIDS:
!! Name changed from fint to dfint.
!! Real variables changed to double precision.
!! External references to CERNLIB (error handling) routines removed.
          DOUBLE PRECISION FUNCTION DFINT(NARG,ARG,NENT,ENT,TABLE)
          INTEGER   NENT(9), INDEX(32)
          DOUBLE PRECISION ARG(9),   ENT(9),   TABLE(9), WEIGHT(32)
          DFINT  =  0d0
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1d0
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0d0)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DFINT  =  DFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      WRITE(*,1000) NARG
          STOP
1000      FORMAT( 7X, 24HFUNCTION DFINT... NARG =,I6,
     +              17H NOT WITHIN RANGE)
          END

