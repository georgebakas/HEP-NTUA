************************************************************************
*
* to be compiled as 
* g77 qqQQ_NLO.for -o qqQQ_NLO.exe 
* 
************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function FinVirtqq(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 CSoneepm2qq,CSoneepm1qq,CSoneep0qq,
     $                 CSoneepm2Kqq,
     $                 CSoneepm1Kqq,CSoneep0Kqq,
     $                 PoleEPm2qq,PoleEPm1qq,PoleEP0qq
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtqq  = PoleEP0qq(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w) 

ch       write(*,*)PoleEPm2qq(s,t,tm,w)
ch       write(*,*)PoleEPm1qq(s,t,tm,w) 
ch     $             - DLog(tm**2/w**2)*PoleEPm2qq(s,t,tm,w)
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end

      double precision function FinVirtgg(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,ypzp2,ypzm2,sqrypzp2,
     $                 sqrypzm2,CSoneepm2gg,CSoneepm1gg,CSoneep0gg,
     $                 CSoneepm1Kgg,CSoneep0Kgg,CSoneepm2Kgg,
     $                 PoleEPm2gg,PoleEPm1gg,PoleEP0gg
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w

      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='ggQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c       s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1gg(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtgg  = PoleEP0gg(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1gg(s,t,tm,w) 
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end

      double precision function FinVirtT13qq(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T13CSoneepm2qq,T13CSoneepm1qq,T13CSoneep0qq,
     $                 T13CSoneepm2Kqq,
     $                 T13CSoneepm1Kqq,T13CSoneep0Kqq,
     $                 T13PoleEPm2qq,T13PoleEPm1qq,T13PoleEP0qq
     $                 ,T13PoleEP0qqc
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT13qq  = T13PoleEP0qq(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w) 
ch       write(*,*)T13PoleEP0qq(s,t,tm,w)/T13PoleEP0qqc(s,t,tm,w),'a'

ch       write(*,*)T13PoleEPm2qq(s,t,tm,w),
ch     $    T13PoleEPm1qq(s,t,tm,w),'ccc'
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end

      double precision function FinVirtT23qq(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T13CSoneepm2qq,T13CSoneepm1qq,T13CSoneep0qq,
     $                 T13CSoneepm2Kqq,
     $                 T13CSoneepm1Kqq,T13CSoneep0Kqq,
     $                 T23PoleEPm2qq,T23PoleEPm1qq,T23PoleEP0qq
     $                 ,T13PoleEP0qqc
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT23qq  = T23PoleEP0qq(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w) 
ch       write(*,*)T13PoleEP0qq(s,t,tm,w)/T13PoleEP0qqc(s,t,tm,w),'a'

ch       write(*,*)T23PoleEPm2qq(s,t,tm,w),
ch     $    T23PoleEPm1qq(s,t,tm,w),'ccc'
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end


      double precision function FinVirtT34qq(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T34CSoneepm2qq,T34CSoneepm1qq,T34CSoneep0qq,
     $                 T34CSoneepm2Kqq,
     $                 T34CSoneepm1Kqq,T34CSoneep0Kqq,
     $                 T34PoleEPm2qq,T34PoleEPm1qq,T34PoleEP0qq
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT34qq  = T34PoleEP0qq(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w)
ch       write(*,*)T34PoleEPm2qq(s,t,tm,w),
ch     $    T34PoleEPm1qq(s,t,tm,w),'aaa' 
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end 

      double precision function FinVirtT34gg(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T34CSoneepm2qq,T34CSoneepm1qq,T34CSoneep0qq,
     $                 T34CSoneepm2Kqq,
     $                 T34CSoneepm1Kqq,T34CSoneep0Kqq,
     $                 T34PoleEPm2qq,T34PoleEPm1qq,PoleEP0T34
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT34gg  = PoleEP0T34(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w)
ch       write(*,*)T34PoleEPm2qq(s,t,tm,w),
ch     $    T34PoleEPm1qq(s,t,tm,w),'aaa' 
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end

      double precision function FinVirtT13gg(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T34CSoneepm2qq,T34CSoneepm1qq,T34CSoneep0qq,
     $                 T34CSoneepm2Kqq,
     $                 T34CSoneepm1Kqq,T34CSoneep0Kqq,
     $                 T34PoleEPm2qq,T34PoleEPm1qq,PoleEP0T13
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT13gg  = PoleEP0T13(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w)
ch       write(*,*)T34PoleEPm2qq(s,t,tm,w),
ch     $    T34PoleEPm1qq(s,t,tm,w),'aaa' 
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end

      double precision function FinVirtT14gg(s,t)
      implicit none
      include 'scale.f'
      include 'masses.f'
      double precision en,w,tm,pi,hM,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,Nc,Nl,Nh,
     $                 stev,ttev,tmin,tmax,rnorm,
     $                 T34CSoneepm2qq,T34CSoneepm1qq,T34CSoneep0qq,
     $                 T34CSoneepm2Kqq,
     $                 T34CSoneepm1Kqq,T34CSoneep0Kqq,
     $                 T34PoleEPm2qq,T34PoleEPm1qq,PoleEP0T14
*
      Integer I
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
*
c      write(6,'('' '')')
c      write(6,'(''s (in TeV^2)              = '',$)')
c      read(5,*) stev
c      write(6,'(''mt (top mass in GeV)      = '',$)')
c      read(5,*) tm
c      write(6,'(''mu (renorm. scale in GeV) = '',$)')
c      read(5,*) w
      tm=mt
      w=scale
*
******* Allowed range for t: *****
*
c      tmin = - stev/2.d0 * ( 1.d0 - Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
c      tmax = - stev/2.d0 * ( 1.d0 + Dsqrt(1.d0-4.d-6*tm**2/stev) ) 
c     $       + 1.d-6*tm**2
*
**********************************
*
c      write(6,'('' '')')
c      write(6,'(''t (in TeV^2) belongs to the range '')')
c      write(6,'(1pe22.5,'' < t < '',1pe22.5)')
c     $         tmax, tmin
c      write(6,'(''t (in TeV^2)              = '',$)')
c      read(5,*) ttev
*
c       if (ttev.le.tmax.OR.ttev.ge.tmin) then
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'(''           t is outside the range '')')
c       write(6,'(''###############################################'')')
c        stop
c       endif
*
c      write(6,'('' '')')
c      write(6,'(''# Choose the normalization factor:'')')
c      write(6,'(''# 0  ->  Exp[-EulerGamma*ep]'')')
c      write(6,'(''# 1  ->  Gamma[1+ep]'')')
c      read(5,*) rnorm
*
c      write(6,'('' '')')
c      write(6,'(''# The results will be also written in qqQQN.res'')')
c       write(6,'(''#'')')
c       write(6,'(''#'')')
c       open(3,FILE='qqQQ_NLO.res')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
c       write(3,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(3,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(3,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(3,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(3,'(''#'')')
c       write(3,'(''###############################################'')')
c       write(3,'(''#'')')
*
************************************************************************
* 
c        s = 1.d+6*stev
c	t = 1.d+6*ttev
*
c       CSoneepm2 = PoleEPm2(s,t,tm,w)
c       CSoneepm1 = PoleEPm1(s,t,tm,w) 
c     $             - DLog(tm**2/w**2)*PoleEPm2(s,t,tm,w)
       FinVirtT14gg  = PoleEP0T14(s,t,tm,w) 
c     $             + 1.d0/2.d0*DLog(tm**2/w**2)**2*PoleEPm2(s,t,tm,w)
c     $             - DLog(tm**2/w**2)*PoleEPm1(s,t,tm,w)
ch       write(*,*)T34PoleEPm2qq(s,t,tm,w),
ch     $    T34PoleEPm1qq(s,t,tm,w),'aaa' 
*
* Normalization factor Exp[-EulerGamma*ep]:
*
c       CSoneepm2K = CSoneepm2
c       CSoneepm1K = CSoneepm1
c       CSoneep0K  = CSoneep0 + z2/2.d0*CSoneepm2
*
*
c       if (rnorm.eq.0) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Exp[-EulerGamma*ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2K
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1K
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0K
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
c       if (rnorm.eq.1) then
c       write(6,'(''###############################################'')')
c       write(6,'(''#'')')
c       write(6,'(''# s (in TeV^2)              = '',F22.9)') stev
c       write(6,'(''# t (in TeV^2)              = '',F22.9)') ttev
c       write(6,'(''# mt (top mass in GeV)      = '',F22.9)') tm
c       write(6,'(''# mu (renorm. scale in GeV) = '',F22.9)') w
c       write(6,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(6,'(''#'')')
c       write(6,'(''###############################################'')')
c       write(6,'('' '')')
c       write(6,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(6,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(6,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(6,'('' '')')
c       write(6,'(''###############################################'')')
*
*
c       write(3,'(''# Normalization factor  Gamma[1+ep]'')') 
c       write(3,'(''#'')') 
c       write(3,'('' '')')
c       write(3,'(''  |M|^2 = 1/ep^2*('',1pe22.8,'')'')') CSoneepm2
c       write(3,'(''        + 1/ep  *('',1pe22.8,'')'')') CSoneepm1
c       write(3,'(''        +        ('',1pe22.8,'')'')') CSoneep0
c       write(3,'('' '')')
c       write(3,'(''###############################################'')')
c       endif
*
*
c15    continue
c       close(3)
      return
      end


*
*
*
************************************************************************
************************************************************************
*
*             SUBROUTINES for HPLs
*
************************************************************************
************************************************************************
**  hplog: a subroutine for the evaluation of harmonic polylogarithms
**         Version 1.1         26/10/2004
**  described in: 
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of the Harmonic 
**                            Polylogarithms up to Weight 4
**                            (hep-ph/0107173; Comp.Phys.Comm. 141 (2001) 296)
**  the harmonic polylogarithms are defined in: 
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
**  changes with respect to version 1.0 (12/07/2001):
**  improved Chebyshev expansions by factoring out leading behaviour
**  to improve the relative acurracy for very small values of the arguments
**************************************************************************
      subroutine hplog(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                      Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 4; 
** Hc1,Hc2,Hc3,Hc4 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4; 
** Hi1,Hi2,Hi3,Hi4 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (r2   = 1.4142135623730950488d0) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.4) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,2,3,4 ' 
        stop
      endif
** check on the range n1:n2 
      if ( (n1.eq.-1).and.(n2.eq.0) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) = -1  
      elseif ( (n1.eq.0).and.(n2.eq.1) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) =  1  
      elseif ( (n1.eq.-1).and.(n2.eq.1) ) then 
        infilldim =  3 
        infill(1) =  0 
        infill(2) = -1  
        infill(3) =  1  
      else 
        print*, ' illegal call of eval1dhpl with the two last ', 
     $          'arguments = (',n1,',',n2,')' 
        print*, ' the allowed values are (-1,0), (0,1), (-1,1) ' 
        stop 
      endif 
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2.d0) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2p1 ) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.le.-r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.-1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      if (n2.eq.0) return
** correct the ill-defined entries
      HY2(1,0) = - HY2(0,1) 
      Hi2(1,0) = 0d0 
      H2(1,0) = dcmplx(HY2(1,0),Hi2(1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(1,0,0) = HY3(0,0,1) 
      Hi3(1,0,0) = 0d0 
      H3(1,0,0) = dcmplx(HY3(1,0,0),Hi3(1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(1,0,0,0) = -HY4(0,0,0,1) 
      Hi4(1,0,0,0) = 0d0 
      H4(1,0,0,0) = dcmplx(HY4(1,0,0,0),Hi4(1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplat1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 1-range  (r2-1) < y <= (r2+1) 
** evaluating first the H(..,r=(1-y)/(1+y)) by calling eval1dhplat0(r)  
** and then expressing H(..,y=(1-r)/(1+r)) in terms of H(..,r) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      r = (1.d0-y)/(1.d0+y) 
*      print*,' eval1dhplat1: y = ',y,', r = ',r 
** the whole (-1,1) range is in general needed for any pair (n1,n2)
      call fillirr1dhplat0(r,nw,HR1,HR2,HR3,HR4,-1,1) 
** fillirr1dhplat1 takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                          HY1,HY2,HY3,HY4, 
     $                          Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1.d0/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                            HY1,HY2,HY3,HY4, 
     $                            Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplinm1(y,nw,H1,H2,H3,H4, 
     $                           HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=-1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplin1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      if (n1.eq.0) return
** correct the ill-defined entries
      HY2(-1,0) = - HY2(0,-1) 
      Hi2(-1,0) = Hi1(0)*HY1(-1)
      H2(-1,0) = dcmplx(HY2(-1,0),Hi2(-1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(-1,0,0) = HY1(-1)*HY2(0,0)+HY3(0,0,-1) 
      Hi3(-1,0,0) = HY1(-1)*Hi2(0,0)-HY2(0,-1)*Hi1(0)
      H3(-1,0,0) = dcmplx(HY3(-1,0,0),Hi3(-1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(-1,0,0,0) = -HY2(0,-1)*HY2(0,0)-HY4(0,0,0,-1) 
      Hi4(-1,0,0,0) = HY1(-1)*Hi3(0,0,0)+HY3(0,0,-1)*Hi1(0)
      H4(-1,0,0,0) = dcmplx(HY4(-1,0,0,0),Hi4(-1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatm1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range  -(r2+1) < y <= -(r2-1) 
** evaluating first the H(..,-y) by calling eval1dhplat1(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplat1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 


      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit double precision (a-h,o-z) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      do k1=n1,n2 
        Hi1(k1) = 0.d0 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            Hi2(k1,k2) = 0.d0 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                Hi3(k1,k2,k3) = 0.d0 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    Hi4(k1,k2,k3,k4) = 0.d0 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
* fills the reducible 1dhpl from the irreducible set
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (pinv = 0.318309886183790672d0) 
      parameter (pi   = 3.14159265358979324d0) 
** combining real and immaginary into the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        H2(k1,k2) = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            H3(k1,k2,k3) = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                H4(k1,k2,k3,k4) = 
     $               dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      iflag = 0 
      do ia =  1,infilldim 
      do ib = ia,infilldim 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,infilldim 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,infilldim 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  dble(H2(k1,k2)) 
        Hi2(k1,k2) = dimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  dble(H3(k1,k2,k3)) 
            Hi3(k1,k2,k3) = dimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  dble(H4(k1,k2,k3,k4)) 
                Hi4(k1,k2,k3,k4) = dimag(H4(k1,k2,k3,k4))*pinv 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL2(iflag,H1,H2,i1,i2,na,nb) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2 
      dimension H1(i1:i2),H2(i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with ordered indices na <= nb 
*      print*,' FILLREDHPL2, iflag =',iflag 
      if ( na.eq.nb ) then 
        H2(na,na) = 1.d0/2*( H1(na) )**2 
      else 
        H2(nb,na) = + H1(na)*H1(nb) - H2(na,nb) 
        if ( iflag.eq.1 ) then 
          call printer2(na,nb) 
        endif 
      endif 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL3(iflag,H1,H2,H3,i1,i2,ia,ib,ic) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na) 
        H3(na,na,na) = 1.d0/6*( H1(na) )**3 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL3, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ia.eq.ib ) then 
* case (na,na,nb) 
        nb = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,na,nb) 
        endif 
        H3(na,nb,na) = + H1(na)*H2(na,nb) - 2*H3(na,na,nb) 
        H3(nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb) 
     $                 - H1(na)*H2(na,nb) + H3(na,na,nb) 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL3, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb) 
        nb = ib 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nb) 
        endif 
        H3(nb,na,nb) = + H1(nb)*H2(na,nb) - 2*H3(na,nb,nb) 
        H3(nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb) 
     $                 - H1(nb)*H2(na,nb) + H3(na,nb,nb) 
* no need to protect against ic.eq.ib 
* when arriving here all indices are different 
      else 
* case (na,nb,nc)    all indices are different 
        nb = ib 
        nc = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nc) 
          call printer3(na,nc,nb) 
        endif 
        H3(nb,na,nc) = + H1(nb)*H2(na,nc) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nb,nc,na) = + H1(na)*H2(nb,nc) 
     $                 - H1(nb)*H2(na,nc) + H3(na,nc,nb) 
        H3(nc,na,nb) = + H1(nc)*H2(na,nb) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nc,nb,na) = + H1(na)*H1(nb)*H1(nc) - H1(na)*H2(nb,nc) 
     $                 - H1(nc)*H2(na,nb) + H3(na,nb,nc) 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL4(iflag,H1,H2,H3,H4,i1,i2,ia,ib,ic,id) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
      dimension H4(i1:i2,i1:i2,i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,na,na,na) 
        H4(na,na,na,na) = 1.d0/24*( H1(na) )**4 
* id cannot be anymore equal to ia 
      else if ( id.eq.ia ) then 
        print*,' FILLREDHPL4, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na,nb) 
        nb = id 
        H4(na,na,nb,na) = + H1(na)*H3(na,na,nb) - 3*H4(na,na,na,nb) 
        H4(na,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,na,nb) + 3*H4(na,na,na,nb) 
        H4(nb,na,na,na) = + 1.d0/6*H1(na)*H1(na)*H1(na)*H1(nb) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    + H1(na)*H3(na,na,nb) - H4(na,na,na,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,na,nb) 
        endif 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL4, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ic.eq.id) ) then 
* case (na,na,nb,nb) 
        nb = ic 
        H4(na,nb,na,nb) = + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    - 2*H4(na,na,nb,nb) 
        H4(na,nb,nb,na) = + H1(na)*H3(na,nb,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,na,nb) = + H1(nb)*H3(na,na,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,nb,na) = + H1(na)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,nb,nb) 
     $                    - 2*H1(nb)*H3(na,na,nb) 
     $                    + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    + 2*H4(na,na,nb,nb) 
        H4(nb,nb,na,na) = + 1.d0/4*H1(na)*H1(na)*H1(nb)*H1(nb) 
     $                    - H1(na)*H1(nb)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nb) 
     $                    + H1(nb)*H3(na,na,nb) - H4(na,na,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nb) 
        endif 
      else if ( ia.eq.ib ) then 
* case (na,na,nb,nc) 
        nb = ic 
        nc = id 
        H4(na,nb,nc,na) = + H1(na)*H3(na,nb,nc) - 2*H4(na,na,nb,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(na,nc,na,nb) = + H2(na,nb)*H2(na,nc) - 2*H4(na,na,nb,nc) 
     $                    - 2*H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(na,nc,nb,na) = + H1(na)*H3(na,nc,nb) - H2(na,nb)*H2(na,nc) 
     $                    + 2*H4(na,na,nb,nc) + H4(na,nb,na,nc) 
        H4(nb,na,na,nc) = + H1(nb)*H3(na,na,nc) - H4(na,na,nb,nc) 
     $                    - H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(nb,na,nc,na) = + H1(na)*H1(nb)*H2(na,nc) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nb)*H3(na,na,nc) + 2*H4(na,na,nb,nc) 
     $                    + 2*H4(na,na,nc,nb) + H4(na,nb,na,nc) 
        H4(nb,nc,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nb)*H2(na,nc) 
     $                    + H1(na)*H3(na,nc,nb) + H1(nb)*H3(na,na,nc) 
     $                    - H4(na,na,nc,nb) 
        H4(nc,na,na,nb) = + H1(nc)*H3(na,na,nb) - H2(na,nb)*H2(na,nc) 
     $                    + H4(na,na,nb,nc) + H4(na,na,nc,nb) 
     $                    + H4(na,nb,na,nc) 
        H4(nc,na,nb,na) = + H1(na)*H1(nc)*H2(na,nb) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nc)*H3(na,na,nb) + H2(na,nb)*H2(na,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(nc,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb)*H1(nc) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nc)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nc) + H1(nc)*H3(na,na,nb) 
     $                    - H4(na,na,nb,nc)  
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nc) 
          call printer4(na,na,nc,nb) 
          call printer4(na,nb,na,nc) 
        endif 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL4, error 3, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,nb,nb,nb) 
        nb = ib 
        H4(nb,na,nb,nb) = + H1(nb)*H3(na,nb,nb) - 3*H4(na,nb,nb,nb) 
        H4(nb,nb,na,nb) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(nb)*H3(na,nb,nb) + 3*H4(na,nb,nb,nb) 
        H4(nb,nb,nb,na) = + 1.d0/6*H1(na)*H1(nb)*H1(nb)*H1(nb) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    + H1(nb)*H3(na,nb,nb) - H4(na,nb,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nb) 
        endif 
* id cannot be anymore equal to ib 
      else if ( id.eq.ib ) then 
        print*,' FILLREDHPL4, error 4, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb,nc) 
        nb = ib 
        nc = id 
        H4(nb,na,nb,nc) = + H1(nb)*H3(na,nb,nc) 
     $                    - 2*H4(na,nb,nb,nc) - H4(na,nb,nc,nb) 
        H4(nb,na,nc,nb) = + H1(nb)*H3(na,nc,nb) - H4(na,nb,nc,nb) 
     $                    - 2*H4(na,nc,nb,nb) 
        H4(nb,nb,na,nc) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H4(na,nb,nb,nc) + H4(na,nb,nc,nb) 
     $                    + H4(na,nc,nb,nb) 
        H4(nb,nb,nc,na) = + H1(na)*H3(nb,nb,nc) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    + H1(nb)*H3(na,nc,nb) - H4(na,nc,nb,nb) 
        H4(nb,nc,na,nb) = - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H2(na,nb)*H2(nb,nc) + H4(na,nb,nc,nb) 
     $                    + 2*H4(na,nc,nb,nb) 
        H4(nb,nc,nb,na) = + H1(na)*H1(nb)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nb,nc) 
     $                    + H1(nb)*H3(na,nb,nc) 
     $                    - H2(na,nb)*H2(nb,nc) - H4(na,nb,nc,nb) 
        H4(nc,na,nb,nb) = + H1(nc)*H3(na,nb,nb) - H4(na,nb,nb,nc) 
     $                    - H4(na,nb,nc,nb) - H4(na,nc,nb,nb) 
        H4(nc,nb,na,nb) = + H1(nb)*H1(nc)*H2(na,nb) 
     $                    - 2*H1(nc)*H3(na,nb,nb) 
     $                    - H2(na,nb)*H2(nb,nc) + 2*H4(na,nb,nb,nc) 
     $                    + H4(na,nb,nc,nb) 
        H4(nc,nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb)*H1(nc) 
     $                    - H1(na)*H1(nb)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nb,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nb) + H2(na,nb)*H2(nb,nc) 
     $                    - H4(na,nb,nb,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nc) 
          call printer4(na,nb,nc,nb) 
          call printer4(na,nc,nb,nb) 
        endif 
* ic cannot be anymore equal to ib 
      else if ( ic.eq.ib ) then 
        print*,' FILLREDHPL4, error 5, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ic.eq.id ) then 
* case (na,nb,nc,nc) 
        nb = ib 
        nc = ic 
        H4(nb,na,nc,nc) = + H1(nb)*H3(na,nc,nc) - H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) - H4(na,nc,nc,nb) 
        H4(nb,nc,na,nc) = - 2*H1(nb)*H3(na,nc,nc) + H2(na,nc)*H2(nb,nc) 
     $                    + H4(na,nc,nb,nc) + 2*H4(na,nc,nc,nb) 
        H4(nb,nc,nc,na) = + H1(na)*H3(nb,nc,nc) + H1(nb)*H3(na,nc,nc) 
     $                    - H2(na,nc)*H2(nb,nc) - H4(na,nc,nc,nb) 
        H4(nc,na,nb,nc) = + H1(nc)*H3(na,nb,nc) - 2*H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,na,nc,nb) = + H1(nc)*H3(na,nc,nb) - H4(na,nc,nb,nc) 
     $                    - 2*H4(na,nc,nc,nb) 
        H4(nc,nb,na,nc) = + H1(nb)*H1(nc)*H2(na,nc) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nc) + 2*H4(na,nb,nc,nc) 
     $                    + H4(na,nc,nb,nc) 
        H4(nc,nb,nc,na) = + H1(na)*H1(nc)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nc,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nc) 
     $                    + H1(nc)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,nc,na,nb) = + 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    + H4(na,nb,nc,nc) + H4(na,nc,nb,nc) 
     $                    + H4(na,nc,nc,nb) 
        H4(nc,nc,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nc)*H1(nc) 
     $                    - H1(na)*H1(nc)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nc) 
     $                    - 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nc) - H4(na,nb,nc,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nc) 
          call printer4(na,nc,nb,nc) 
          call printer4(na,nc,nc,nb) 
        endif 
* no need to protect against id.eq.ic 
* when arriving here all indices are different 
      else 
* case (na,nb,nc,nd) all indices are different 
        nb = ib 
        nc = ic 
        nd = id 
        H4(nb,na,nc,nd) = + H1(nb)*H3(na,nc,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nc,nb,nd) - H4(na,nc,nd,nb) 
        H4(nb,na,nd,nc) = + H1(nb)*H3(na,nd,nc) - H4(na,nb,nd,nc) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nc,na,nd) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nd)*H2(nb,nc) + H4(na,nc,nb,nd) 
     $                    + H4(na,nc,nd,nb) + H4(na,nd,nc,nb) 
        H4(nb,nc,nd,na) = + H1(na)*H3(nb,nc,nd) + H1(nb)*H3(na,nd,nc) 
     $                    - H2(na,nd)*H2(nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nd,na,nc) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nc)*H2(nb,nd) + H4(na,nc,nd,nb) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nb,nd,nc,na) = + H1(na)*H3(nb,nd,nc) + H1(nb)*H3(na,nc,nd) 
     $                    - H2(na,nc)*H2(nb,nd) - H4(na,nc,nd,nb) 
        H4(nc,na,nb,nd) = + H1(nc)*H3(na,nb,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nc,nb,nd) 
        H4(nc,na,nd,nb) = + H1(nc)*H3(na,nd,nb) - H4(na,nc,nd,nb) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nc,nb,na,nd) = + H1(nb)*H1(nc)*H2(na,nd) 
     $                    - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    - H2(na,nd)*H2(nb,nc) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nd,nb,nc) 
        H4(nc,nb,nd,na) = + H1(na)*H1(nc)*H2(nb,nd) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nd) 
     $                    + H1(nc)*H3(na,nd,nb) + H2(na,nd)*H2(nb,nc) 
     $                    - H4(na,nd,nb,nc) 
        H4(nc,nd,na,nb) = - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    + H2(na,nb)*H2(nc,nd) + H4(na,nb,nd,nc) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nc,nd,nb,na) = + H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nc)*H2(nb,nd) 
     $                    + H1(na)*H3(nb,nd,nc) + H1(nc)*H3(na,nb,nd) 
     $                    - H2(na,nb)*H2(nc,nd) - H4(na,nb,nd,nc) 
        H4(nd,na,nb,nc) = + H1(nd)*H3(na,nb,nc) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nd,nb,nc) 
        H4(nd,na,nc,nb) = + H1(nd)*H3(na,nc,nb) - H4(na,nc,nb,nd) 
     $                    - H4(na,nc,nd,nb) - H4(na,nd,nc,nb) 
        H4(nd,nb,na,nc) = + H1(nb)*H1(nd)*H2(na,nc) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nc,nb,nd) 
        H4(nd,nb,nc,na) = + H1(na)*H1(nd)*H2(nb,nc) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nd)*H2(na,nc) 
     $                    + H1(nd)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nd) 
     $                    - H4(na,nc,nb,nd) 
        H4(nd,nc,na,nb) = + H1(nc)*H1(nd)*H2(na,nb) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nb)*H2(nc,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nc,nb,nd) + H4(na,nc,nd,nb) 
        H4(nd,nc,nb,na) = + H1(na)*H1(nb)*H1(nc)*H1(nd) 
     $                    - H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nd)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nd) 
     $                    - H1(nc)*H1(nd)*H2(na,nb) 
     $                    + H1(nd)*H3(na,nb,nc) 
     $                    + H2(na,nb)*H2(nc,nd) - H4(na,nb,nc,nd) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nd) 
          call printer4(na,nb,nd,nc) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nd,nb,nc) 
          call printer4(na,nd,nc,nb) 
        endif 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine printer2(na,nb) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer3(na,nb,nc) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer4(na,nb,nc,nd) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine subprint(n,na) 
      if ( na.lt.0 ) then 
        write (n,102) na 
      else 
        write (n,101) na 
      endif 
      return 
  101 format(i1,$) 
  102 format(i2,$) 
      end 

************************************************************************
** the following routines contain th set of routines evaluating 
** irreducible 1dhpl's for various values of the arguments 
************************************************************************ 
      subroutine fillh1(y,H1,HY1,Hi1,n1,n2) 
** fillh1 evaluates the 1dhpl's of weight 1 
      implicit double precision (a-h,o-z) 
      complex*16 H1 
      dimension H1(n1:n2) 
      dimension HY1(n1:n2) 
      dimension Hi1(n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
      if ( n1.eq.-1) then 
        if ( y.ge.-1.d0 ) then 
          HY1(-1) = log(1.d0+y) 
          Hi1(-1) = 0.d0 
        elseif ( y.lt.-1.d0 ) then 
          HY1(-1) = log(-1.d0-y) 
          Hi1(-1) = 1.d0 
        endif 
        H1(-1) = dcmplx(HY1(-1),pi*Hi1(-1)) 
      endif 
      if ( y.ge.0.d0 ) then 
        HY1(0) = log(y) 
*        Hi1(0) = 0.d0 
      elseif ( y.lt.0.d0 ) then 
        HY1(0) = log(-y) 
        Hi1(0) = 1.d0 
      endif 
      H1(0) = dcmplx(HY1(0),pi*Hi1(0)) 
      if ( n2.eq.1 ) then 
        if ( y.ge.1.d0 ) then 
          HY1(1) = - log(-1.d0+y) 
          Hi1(1) = 1.d0 
        elseif ( y.lt.1.d0 ) then 
          HY1(1) = - log(1.d0-y) 
          Hi1(1) = 0.d0 
        endif 
        H1(1) = dcmplx(HY1(1),pi*Hi1(1)) 
      endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluate the HPL from their power series expansions
** fillirr1dhplat0 is called by eval1dhplat0; 
** it is guaranteed that nw is in the range 1:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
** 
** for y < 0 DOES NOT evaluates the immaginary part of H(0,y) = log(y) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluating the required 1dHPL of weight 1 
      if ( n1.eq.-1) then 
** 1+y = (1+ep)/(1-ep), ep = y/(2+y) 
** log(1+y) = log((1+y)/(1-y)) = 2*ep*(1+ep^2/3+ep^4/5+.....) 
** at y= -(r2-1) = - 0.4142135624, ep = - 0.26120387496 
** ep2 = 0.068227464296, ep2^13 = 6.9 x 10^(-16) 
         ep = y/(2.d0+y) 
         e2 = ep*ep 
*         v = log(1.d0+y) 
         v = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(-1) = v 
      endif 
      if (y.ge.0d0) then 
         HY1(0) = log(y) 
      else 
         HY1(0) = log(-y) 
** the immaginary part is evaluated in the calling routine eval1dhplat0 
**       Hi1(0) = 1d0 
      endif 
      if ( n2.eq.1) then 
** 1-y = (1-ep)/(1+ep), ep = y/(2-y) 
         ep = y/(2.d0-y) 
         e2 = ep*ep 
*         u = - log(1.d0-y) 
         u = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(1) = u 
      endif 
      if ( nw.eq.1 ) return 
** from now on nw > 1 
** evaluating the Cebyshev polynomials for the expansions 
      ep = y 
      if ( n2.eq.1) then 
        tu01 = 20d0/11d0*u 
        tu02 = 2d0*tu01*tu01 - 1d0 
        tu03 = 2d0*tu01*tu02 - tu01 
        tu04 = 2d0*tu01*tu03 - tu02 
        tu05 = 2d0*tu01*tu04 - tu03 
        tu06 = 2d0*tu01*tu05 - tu04 
        tu07 = 2d0*tu01*tu06 - tu05 
        tu08 = 2d0*tu01*tu07 - tu06 
        tu09 = 2d0*tu01*tu08 - tu07 
        tu10 = 2d0*tu01*tu09 - tu08 
        tu11 = 2d0*tu01*tu10 - tu09 
        tu12 = 2d0*tu01*tu11 - tu10 
        u01 = u
        u02 = u01*u01
        u03 = u01*u02
      endif 
      if ( n1.eq.-1 ) then 
        tv01 = 20d0/11d0*v 
        tv02 = 2d0*tv01*tv01 - 1d0 
        tv03 = 2d0*tv01*tv02 - tv01 
        tv04 = 2d0*tv01*tv03 - tv02 
        tv05 = 2d0*tv01*tv04 - tv03 
        tv06 = 2d0*tv01*tv05 - tv04 
        tv07 = 2d0*tv01*tv06 - tv05 
        tv08 = 2d0*tv01*tv07 - tv06 
        tv09 = 2d0*tv01*tv08 - tv07 
        tv10 = 2d0*tv01*tv09 - tv08 
        tv11 = 2d0*tv01*tv10 - tv09 
        tv12 = 2d0*tv01*tv11 - tv10 
        v01 = v
        v02 = v01*v01
        v03 = v01*v02
      endif 
** evaluating the expansions 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = u01*(
     $  - 1.3750000000000000d-01*tu01
     $  + 4.1887406497219320d-03*tu02
     $  - 3.1529487739207918d-06*tu04
     $  + 4.0387979720925612d-09*tu06
     $  - 5.9161728033102941d-12*tu08
     $  + 9.2089553676892591d-15*tu10
     $  + 1.0041918976432192d+00)
**    it would be wrong to write 
**    if ( nw.eq.2 ) return 
**    because the (n1.eq.-1).and.(n2.eq.1) case is not yet complete 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = u01*(
     $  - 2.0733063311847688d-01*tu01
     $  + 1.1909822093713566d-02*tu02
     $  - 3.5978922715974371d-04*tu03
     $  + 1.4651506635246111d-06*tu04
     $  + 2.5264978015885375d-07*tu05
     $  - 2.9129600146853522d-09*tu06
     $  - 3.1200969105278629d-10*tu07
     $  + 5.5608760227994999d-12*tu08
     $  + 4.4683018298608872d-13*tu09
     $  - 1.0377741049678742d-14*tu10
     $  - 6.8492668900504254d-16*tu11
     $  + 1.0119083540245187d+00)
      HY3(0,1,1) = u02*(
     $  - 4.5833333333333333d-02*tu01
     $  + 1.5702519995611332d-03*tu02
     $  - 1.3132234112581665d-06*tu04
     $  + 1.7663819193128645d-09*tu06
     $  - 2.6615094984230225d-12*tu08
     $  + 4.2197085610629510d-15*tu10
     $  + 2.5157156699202004d-01)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = u01*(
     $  - 2.4309920936894915d-01*tu01
     $  + 1.7710499311295207d-02*tu02
     $  - 8.2489740705617284d-04*tu03
     $  + 2.1971556431143249d-05*tu04
     $  - 9.6292059718056785d-08*tu05
     $  - 1.3396044855444690d-08*tu06
     $  + 1.9835004778127507d-10*tu07
     $  + 1.4796431466871027d-11*tu08
     $  - 3.8471818563928138d-13*tu09
     $  - 1.8899753406004423d-14*tu10
     $  + 7.2332583052852315d-16*tu11
     $  + 2.5542041547393113d-17*tu12
     $  + 1.0176885143440038d+00)
      HY4(0,0,1,1) = u02*(
     $  - 3.8496955551090814d-02*tu01
     $  + 2.7602283087853423d-03*tu02
     $  - 1.0070794072930197d-04*tu03
     $  + 7.6339942564322393d-07*tu04
     $  + 7.7318248534321445d-08*tu05
     $  - 1.4676109918929733d-09*tu06
     $  - 9.8845338864317649d-11*tu07
     $  + 2.7576552760575231d-12*tu08
     $  + 1.4399601068889699d-13*tu09
     $  - 5.1074353206574952d-15*tu10
     $  - 2.2287715896178803d-16*tu11
     $  + 1.2775946343898593d-01)
      HY4(0,1,1,1) = u03*(
     $  - 1.1458333333333333d-02*tu01
     $  + 4.1863378899004153d-04*tu02
     $  - 3.7509447620205836d-07*tu04
     $  + 5.2322893016315673d-10*tu06
     $  - 8.0632112665315483d-13*tu08
     $  + 1.2980885707671547d-15*tu10
     $  + 5.5974564963058350d-02)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = v01*(
     $  + 1.3750000000000000d-01*tv01
     $  + 4.1887406497219320d-03*tv02
     $  - 3.1529487739207918d-06*tv04
     $  + 4.0387979720925612d-09*tv06
     $  - 5.9161728033102941d-12*tv08
     $  + 9.2089553676892591d-15*tv10
     $  + 1.0041918976432192d+00)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = v01*(
     $  + 2.0733063311847688d-01*tv01
     $  + 1.1909822093713566d-02*tv02
     $  + 3.5978922715974371d-04*tv03
     $  + 1.4651506635246111d-06*tv04
     $  - 2.5264978015885375d-07*tv05
     $  - 2.9129600146853522d-09*tv06
     $  + 3.1200969105278629d-10*tv07
     $  + 5.5608760227994999d-12*tv08
     $  - 4.4683018298608872d-13*tv09
     $  - 1.0377741049678742d-14*tv10
     $  + 6.8492668900504254d-16*tv11
     $  + 1.0119083540245187d+00)
      HY3(0,-1,-1) = v02*(
     $  + 4.5833333333333333d-02*tv01
     $  + 1.5702519995611332d-03*tv02
     $  - 1.3132234112581665d-06*tv04
     $  + 1.7663819193128645d-09*tv06
     $  - 2.6615094984230225d-12*tv08
     $  + 4.2197085610629510d-15*tv10
     $  + 2.5157156699202004d-01)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = v01*(
     $  + 2.4309920936894915d-01*tv01
     $  + 1.7710499311295207d-02*tv02
     $  + 8.2489740705617284d-04*tv03
     $  + 2.1971556431143249d-05*tv04
     $  + 9.6292059718056785d-08*tv05
     $  - 1.3396044855444690d-08*tv06
     $  - 1.9835004778127507d-10*tv07
     $  + 1.4796431466871027d-11*tv08
     $  + 3.8471818563928138d-13*tv09
     $  - 1.8899753406004423d-14*tv10
     $  - 7.2332583052852315d-16*tv11
     $  + 2.5542041547393113d-17*tv12
     $  + 1.0176885143440038d+00)
      HY4(0,0,-1,-1) = v02*(
     $  + 3.8496955551090814d-02*tv01
     $  + 2.7602283087853423d-03*tv02
     $  + 1.0070794072930197d-04*tv03
     $  + 7.6339942564322393d-07*tv04
     $  - 7.7318248534321445d-08*tv05
     $  - 1.4676109918929733d-09*tv06
     $  + 9.8845338864317649d-11*tv07
     $  + 2.7576552760575231d-12*tv08
     $  - 1.4399601068889699d-13*tv09
     $  - 5.1074353206574952d-15*tv10
     $  + 2.2287715896178803d-16*tv11
     $  + 1.2775946343898593d-01)
      HY4(0,-1,-1,-1) = v03*(
     $  + 1.1458333333333333d-02*tv01
     $  + 4.1863378899004153d-04*tv02
     $  - 3.7509447620205836d-07*tv04
     $  + 5.2322893016315673d-10*tv06
     $  - 8.0632112665315483d-13*tv08
     $  + 1.2980885707671547d-15*tv10
     $  + 5.5974564963058350d-02)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = u01*(
     $  - 1.0634379058775795d-01*tu01
     $  + 3.9945736734466125d-03*tu02
     $  - 3.7507406061620431d-05*tu03
     $  - 2.6446877835643132d-06*tu04
     $  + 6.3655681886837401d-08*tu05
     $  + 2.8057793781419732d-09*tu06
     $  - 1.1163462030090596d-10*tu07
     $  - 3.1053577733469452d-12*tu08
     $  + 1.9419310124261954d-13*tu09
     $  + 3.0929819630518251d-15*tu10
     $  - 3.3172586119612366d-16*tu11
     $  + 6.9714440173006331d-01)
     $  - 6.9314718055994530d-01*HY1(-1)
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = u01*(
     $  - 1.4936500293026433d-01*tu01
     $  + 9.1294996590286675d-03*tu02
     $  - 3.1364411079828567d-04*tu03
     $  + 3.3241884419265611d-06*tu04
     $  + 1.7149504451979047d-07*tu05
     $  - 4.9072777567462259d-09*tu06
     $  - 1.6287492110060996d-10*tu07
     $  + 7.7880352395550587d-12*tu08
     $  + 1.7144761374259933d-13*tu09
     $  - 1.2581978726037913d-14*tu10
     $  - 1.7794760437458505d-16*tu11
     $  + 2.0410050886562810d-17*tu12
     $  + 7.0227335111545365d-01)
     $  - 6.9314718055994530d-01*HY2(0,-1)
      HY3(0,1,-1) = v01*(
     $  - 1.4936500293026433d-01*tv01
     $  - 9.1294996590286675d-03*tv02
     $  - 3.1364411079828567d-04*tv03
     $  - 3.3241884419265611d-06*tv04
     $  + 1.7149504451979047d-07*tv05
     $  + 4.9072777567462259d-09*tv06
     $  - 1.6287492110060996d-10*tv07
     $  - 7.7880352395550587d-12*tv08
     $  + 1.7144761374259933d-13*tv09
     $  + 1.2581978726037913d-14*tv10
     $  - 1.7794760437458505d-16*tv11
     $  - 2.0410050886562810d-17*tv12
     $  - 7.0227335111545365d-01)
     $  + 6.9314718055994530d-01*HY2(0,1)
      HY3(-1,-1,1) = u01*(
     $  - 1.3057686804989822d-01*tu01
     $  + 8.4523830620987958d-03*tu02
     $  - 3.1975649807736323d-04*tu03
     $  + 4.6539175216815103d-06*tu04
     $  + 1.5659574963291279d-07*tu05
     $  - 7.2562180394741360d-09*tu06
     $  - 9.6461310543268369d-11*tu07
     $  + 1.1508160166123477d-11*tu08
     $  - 8.9467904347194841d-15*tu09
     $  - 1.7738601968192296d-14*tu10
     $  + 2.3555578892527793d-16*tu11
     $  + 2.6064475394479255d-17*tu12
     $  + 5.9068824834184565d-01)
     $  - 5.8224052646501250d-01*HY1(-1)
     $  - 3.4657359027997265d-01*HY1(-1)*HY1(-1)
      HY3(-1,1,1) = u01*(
     $  + 1.3339975063287861d-01*tu01
     $  - 1.1138781596983798d-02*tu02
     $  + 4.2451528270924965d-04*tu03
     $  - 3.1992340565454975d-06*tu04
     $  - 3.2466452956593668d-07*tu05
     $  + 6.5128138606318713d-09*tu06
     $  + 3.7579027495485995d-10*tu07
     $  - 1.2538868973543979d-11*tu08
     $  - 4.5112950725555316d-13*tu09
     $  + 2.3153072699638146d-14*tu10
     $  + 5.0451421432830974d-16*tu11
     $  - 4.1321853554982899d-17*tu12
     $  - 2.5136208279665204d-01)
     $  + 2.4022650695910071d-01*HY1(-1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = u01*(
     $  - 1.7139866732942541d-01*tu01
     $  + 1.2827932191861988d-02*tu02
     $  - 6.2674225522140709d-04*tu03
     $  + 1.8556075778365746d-05*tu04
     $  - 1.8020791893375308d-07*tu05
     $  - 8.5359006198704024d-09*tu06
     $  + 2.4439311365943823d-10*tu07
     $  + 7.2062003807073264d-12*tu08
     $  - 3.7410898512861779d-13*tu09
     $  - 6.6835951661525073d-15*tu10
     $  + 5.9500622408716246d-16*tu11
     $  + 7.0595654813291542d-01)
     $  - 6.9314718055994530d-01*HY3(0,0,-1)
      HY4(0,0,1,-1) = v01*(
     $  - 1.7139866732942541d-01*tv01
     $  - 1.2827932191861988d-02*tv02
     $  - 6.2674225522140709d-04*tv03
     $  - 1.8556075778365746d-05*tv04
     $  - 1.8020791893375308d-07*tv05
     $  + 8.5359006198704024d-09*tv06
     $  + 2.4439311365943823d-10*tv07
     $  - 7.2062003807073264d-12*tv08
     $  - 3.7410898512861779d-13*tv09
     $  + 6.6835951661525073d-15*tv10
     $  + 5.9500622408716246d-16*tv11
     $  - 7.0595654813291542d-01)
     $  + 6.9314718055994530d-01*HY3(0,0,1)
      HY4(0,-1,0,1) = u01*(
     $  - 2.0402091563423462d-01*tu01
     $  + 1.5326796026789534d-02*tu02
     $  - 7.5141752370299855d-04*tu03
     $  + 2.2195862858607303d-05*tu04
     $  - 1.9717216111888252d-07*tu05
     $  - 1.1758559785586910d-08*tu06
     $  + 3.3692089706460695d-10*tu07
     $  + 9.7746481899849247d-12*tu08
     $  - 5.6212130740537968d-13*tu09
     $  - 7.4090460514738114d-15*tu10
     $  + 9.1690209320756988d-16*tu11
     $  + 8.3777162181970230d-01)
     $  - 8.2246703342411321d-01*HY2(0,-1)
      HY4(0,-1,-1,1) = u01*(
     $  - 1.4659170690032865d-01*tu01
     $  + 1.1272254950417420d-02*tu02
     $  - 5.7548355241851144d-04*tu03
     $  + 1.8491926476063684d-05*tu04
     $  - 2.4426966611375586d-07*tu05
     $  - 7.0735008480283508d-09*tu06
     $  + 3.2062539799800555d-10*tu07
     $  + 3.3956059929569370d-12*tu08
     $  - 4.4853380323145144d-13*tu09
     $  + 1.2349995594214108d-15*tu10
     $  + 6.2581127147090835d-16*tu11
     $  + 5.9349428241205865d-01)
     $  - 5.8224052646501250d-01*HY2(0,-1)
     $  - 6.9314718055994530d-01*HY3(0,-1,-1)
      HY4(0,-1,1,-1) = v01*(
     $  - 1.6348235416997523d-01*tv01
     $  - 1.4075851558225201d-02*tv02
     $  - 7.0421561036652880d-04*tv03
     $  - 1.8676050687724687d-05*tv04
     $  - 6.3797094702302232d-08*tv05
     $  + 9.8853199509129313d-09*tv06
     $  + 5.9389227320253100d-11*tv07
     $  - 1.0707718572496392d-11*tv08
     $  - 7.4475662801536059d-14*tv09
     $  + 1.4186776038855675d-14*tv10
     $  + 1.0406065160682246d-16*tv11
     $  - 2.0606613212542347d-17*tv12
     $  - 4.9451017952969702d-01)
     $  + 4.8045301391820142d-01*HY2(0,-1)
     $  + 6.9314718055994530d-01*HY3(0,-1,1)
      HY4(0,1,-1,-1) = v01*(
     $  - 1.0118780392493172d-01*tv01
     $  - 1.0875902950044192d-02*tv02
     $  - 6.9860578351752513d-04*tv03
     $  - 2.5530533667092047d-05*tv04
     $  - 2.8880336859297735d-07*tv05
     $  + 1.5974880084273227d-08*tv06
     $  + 5.0308622711054577d-10*tv07
     $  - 1.5805265017025195d-11*tv08
     $  - 8.6523820415066307d-13*tv09
     $  + 1.6459560669863922d-14*tv10
     $  + 1.4691445597298165d-15*tv11
     $  - 2.5107686338477598d-01)
     $  + 2.4022650695910071d-01*HY2(0,1)
      HY4(0,-1,1,1) = u01*(
     $  + 1.0118780392493172d-01*tu01
     $  - 1.0875902950044192d-02*tu02
     $  + 6.9860578351752513d-04*tu03
     $  - 2.5530533667092047d-05*tu04
     $  + 2.8880336859297735d-07*tu05
     $  + 1.5974880084273227d-08*tu06
     $  - 5.0308622711054577d-10*tu07
     $  - 1.5805265017025195d-11*tu08
     $  + 8.6523820415066307d-13*tu09
     $  + 1.6459560669863922d-14*tu10
     $  - 1.4691445597298165d-15*tu11
     $  - 2.5107686338477598d-01)
     $  + 2.4022650695910071d-01*HY2(0,-1)
      HY4(0,1,-1,1) = u01*(
     $  + 1.6348235416997523d-01*tu01
     $  - 1.4075851558225201d-02*tu02
     $  + 7.0421561036652880d-04*tu03
     $  - 1.8676050687724687d-05*tu04
     $  + 6.3797094702302232d-08*tu05
     $  + 9.8853199509129313d-09*tu06
     $  - 5.9389227320253100d-11*tu07
     $  - 1.0707718572496392d-11*tu08
     $  + 7.4475662801536059d-14*tu09
     $  + 1.4186776038855675d-14*tu10
     $  - 1.0406065160682246d-16*tu11
     $  - 2.0606613212542347d-17*tu12
     $  - 4.9451017952969702d-01)
     $  + 4.8045301391820142d-01*HY2(0,1)
     $  - 6.9314718055994530d-01*HY3(0,1,-1)
      HY4(0,1,1,-1) = v01*(
     $  + 1.4659170690032865d-01*tv01
     $  + 1.1272254950417420d-02*tv02
     $  + 5.7548355241851144d-04*tv03
     $  + 1.8491926476063684d-05*tv04
     $  + 2.4426966611375586d-07*tv05
     $  - 7.0735008480283508d-09*tv06
     $  - 3.2062539799800555d-10*tv07
     $  + 3.3956059929569370d-12*tv08
     $  + 4.4853380323145144d-13*tv09
     $  + 1.2349995594214108d-15*tv10
     $  - 6.2581127147090835d-16*tv11
     $  + 5.9349428241205865d-01)
     $  - 5.8224052646501250d-01*HY2(0,1)
     $  + 6.9314718055994530d-01*HY3(0,1,1)
      HY4(-1,-1,-1,1) = u01*(
     $  - 1.3704186675693694d-01*tu01
     $  + 1.0738664799045256d-02*tu02
     $  - 5.6405487546831413d-04*tu03
     $  + 1.8975289545570792d-05*tu04
     $  - 2.8136699841840261d-07*tu05
     $  - 6.9483354697236112d-09*tu06
     $  + 3.9086778571652884d-10*tu07
     $  + 1.5013470800132572d-12*tu08
     $  - 5.4857187012000446d-13*tu09
     $  + 7.0089201006248443d-15*tu10
     $  + 7.1985359614000227d-16*tu11
     $  - 2.2459200848103789d-17*tu12
     $  + 5.4793287616771010d-01)
     $  - 5.3721319360804020d-01*HY1(-1)
     $  - 2.9112026323250625d-01*HY1(-1)*HY1(-1)
     $  - 1.1552453009332421d-01*HY1(-1)*HY1(-1)*HY1(-1)
      HY4(-1,-1,1,1) = u01*(
     $  + 1.0577793416858150d-01*tu01
     $  - 1.0477240350312842d-02*tu02
     $  + 6.6271618168294058d-04*tu03
     $  - 2.5384291818094501d-05*tu04
     $  + 3.7532043652092099d-07*tu05
     $  + 1.3804832669451887d-08*tu06
     $  - 6.6481886939724744d-10*tu07
     $  - 8.3615840903091915d-12*tu08
     $  + 1.1154536519042777d-12*tu09
     $  - 2.8020941621207856d-15*tu10
     $  - 1.7701739081331563d-15*tu11
     $  + 2.7769076888919123d-17*tu12
     $  - 3.1927721734213725d-01)
     $  + 3.0882537509683393d-01*HY1(-1)
     $  + 1.2011325347955035d-01*HY1(-1)*HY1(-1)
      HY4(-1,1,1,1) = u01*(
     $  - 3.2833473259986549d-02*tu01
     $  + 8.5187654858627660d-03*tu02
     $  - 7.6913672477984726d-04*tu03
     $  + 3.0885336618882649d-05*tu04
     $  - 2.3861711695799581d-07*tu05
     $  - 2.5246408915769234d-08*tu06
     $  + 5.1217849936005379d-10*tu07
     $  + 3.0326173877457694d-11*tu08
     $  - 1.0166074598017197d-12*tu09
     $  - 3.7296649908045336d-14*tu10
     $  + 1.9148798730971194d-15*tu11
     $  + 4.2483769310186512d-17*tu12
     $  + 6.3991963537293034d-02)
     $  - 5.5504108664821579d-02*HY1(-1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for r2m1 < y < r2p1
** fillirr1dhplat1 is called by eval1dhplat1 after calling 
** fillirr1dhplat0 with argument r=(1-y)/(1+y) 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      if (r.lt.0d0) then 
      Hi2(0,1) = 
     $  - HR1( -1) 
     $  - HR1(1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
     $  - 1.6449340668482264d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1( -1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 1.6449340668482264d+00*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  - HR1(1) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + HR1(1) *HR2(0,1)
     $  - HR3( -1,-1,1)
     $  + HR3( -1,1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
     $  - HR3(0,1, -1)
     $  - HR3(0,1,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  - HR1(0) *HR2(-1,1)
     $  + HR1(0) *HR2(0,-1)
     $  + HR1(0) *HR2(0,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,1)
     $  - HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
     $  - HR3(0,0,1) 
     $  - HR3(0,1, -1)
      if (r.lt.0d0) then 
      HY3(0,1,1) = HY3(0,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
     $  + 4.9348022005446793d+00*HR1(1)
      Hi3(0,0,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1) 
      Hi3(0,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      endif
      endif
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1( -1)*HR3(0,1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 8.2246703342411321d-01*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  - HR1(1) *HR3(-1,1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + HR1(1) *HR3(0,1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  + HR4( -1,1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
     $  - HR4(0,1,1,1) 
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1( -1)*HR3(0,0,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + 2.5000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR1(1)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - HR1(0) *HR3(-1,1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + HR1(0) *HR3(0,1,1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + HR1(1) *HR3(0,0,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(0,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
     $  - HR4(0,0,1, -1)
     $  - HR4(0,0,1,1) 
     $  - 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + HR1(0) *HR3(0,0,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + 5.5504108664821579d-02*HR1(1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
     $  - HR4(0,0,0,1) 
     $  - HR4(0,0,1, -1)
     $  - HR4(0,1, -1,-1)
      if (r.lt.0d0) then 
      HY4(0,0,1,1) = HY4(0,0,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 2.4674011002723396d+00*HR1(1)*HR1(1)
      HY4(0,1,1,1) = HY4(0,1,1,1) 
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(0)*HR1(1)
     $  - 3.4205442319285582d+00*HR1(1)
     $  - 4.9348022005446793d+00*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR2(0,1)
      Hi4(0,0,0,1) = 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  - 1.666666666666666d-01*HR1(1)*HR1(1)*HR1(1) 
      Hi4(0,0,1,1) = 
     $  - 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  + 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - HR1( -1)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1( -1)*HR2(0,-1) 
     $  + HR1( -1)*HR2(0,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(1)*HR1(1) 
     $  - 3.465735902799726d-01*HR1(1)*HR1(1) 
     $  - HR1(1) *HR2(-1,1) 
     $  + HR1(1) *HR2(0,-1) 
     $  + HR1(1) *HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  + HR3( -1,1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0, -1,1) 
     $  - HR3(0,1, -1) 
     $  - HR3(0,1,1) 
      Hi4(0,1,1,1) = 
     $  + 1.404707559889125d+00*HR1(-1) 
     $  + 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(0)*HR1(0) 
     $  + HR1( -1)*HR1(0)*HR1(1) 
     $  + 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(0)*HR1(1) 
     $  - HR1(0) *HR2(-1,1) 
     $  + HR1(0) *HR2(0,-1) 
     $  + HR1(0) *HR2(0,1) 
     $  + 1.404707559889125d+00*HR1(1) 
     $  - 6.931471805599453d-01*HR2(-1,1) 
     $  + 6.931471805599453d-01*HR2(0,-1) 
     $  + 6.931471805599453d-01*HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0,0, -1) 
     $  - HR3(0,0,1)  
     $  - HR3(0,1, -1) 
      endif 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
       HY2(0,-1) = 
     $  + 8.2246703342411321d-01
     $  - 6.9314718055994530d-01*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)
     $  - HR2( -1,1)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
     $  - 8.2246703342411321d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  + HR1(1) *HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3( -1,1,1)
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  - HR3( -1,-1,1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
     $  - 9.0154267736969571d-01*HR1(-1)
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 9.0154267736969571d-01*HR1(1)
     $  + 4.1123351671205660d-01*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,1,1)
     $  - HR4( -1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4( -1,1,1,1)
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
     $  - 1.5025711289494928d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  - 1.5025711289494928d-01*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
     $  - 5.5504108664821579d-02*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 5.5504108664821579d-02*HR1(1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - HR4( -1,-1,-1,1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR2(0, -1)
      if (r.lt.0d0) then 
      Hi2(-1,1) = 
     $  - HR1( -1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      if (r.lt.0d0) then 
      HY3(-1,1,1) = HY3(-1,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
      Hi3(0,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  - HR2( -1,1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HR1(-1) 
     $  - 6.9314718055994530d-01*HR1(1) 
      Hi3(-1,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
      Hi3(-1,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(0) 
     $  - HR2(0, -1) 
      endif 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
     $  - 2.4307035167006157d-01*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  - HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(-1,1,1)
     $  - 2.4307035167006157d-01*HR1(1)
     $  + 2.9112026323250625d-01*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
     $  - 5.0821521280468485d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0821521280468485d-01*HR1(1)
     $  - 5.3134677019160696d-01*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(0,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + HR4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + HR4(0,1, -1,1)
     $  + HR4(0,1,1, -1)
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
     $  - 3.8889584616810632d-01*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(-1,1,1)
     $  - 3.8889584616810632d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(-1,1,1)
     $  + 4.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - 9.4753004230127705d-02*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
     $  - 2.1407237086670622d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 2.1407237086670622d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + HR4(0, -1,1,-1)
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
     $  + 4.7533770109129867d-01*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 4.7533770109129867d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  - 5.3721319360804020d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  - HR4( -1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
     $  + 1.4780047665430420d+00*HR1(-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  - HR1(0) *HR3(0,-1,1)
     $  - HR1(0) *HR3(0,1,-1)
     $  + 1.4780047665430420d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
     $  - 5.8224052646501250d-01*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  - HR4(0, -1,0,1)
     $  + HR4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,-1)
     $  + HR4(0,1, -1,-1)
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
     $  - 1.1073038989294665d+00*HR1(-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  - 1.1073038989294665d+00*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 1.0626935403832139d+00*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  - HR2(0, -1)*HR2(0,1)
     $  + 1.0626935403832139d+00*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,1)
     $  + HR4(0, -1,0,1)
     $  + HR4(0,0, -1,1)
     $  + HR4(0,0,1, -1)
     $  + HR4(0,1, -1,-1)
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - HR4(0, -1,-1,-1)
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0,0, -1,-1)
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
      if (r.lt.0d0) then 
      HY4(0,-1,1,1) = HY4(0,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR2(-1,1)
      HY4(0,1,1,-1) = HY4(0,1,1,-1) 
     $  + 3.4205442319285582d+00*HR1(-1)
     $  + 3.4205442319285582d+00*HR1(1)
      HY4(-1,-1,1,1) = HY4(-1,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1) 
      HY4(-1,1,1,1) = HY4(-1,1,1,1)
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR2(0,-1)
      Hi4(0,0,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1(1) *HR2(-1,1) 
     $  + HR3( -1,-1,1) 
     $  - HR3( -1,1,1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1) 
      Hi4(0,-1,0,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR1(1) *HR2(-1,1) 
     $  - 2.0000000000000000d+00*HR3(-1,-1,1) 
     $  + 2.0000000000000000d+00*HR3(-1,1,1) 
      Hi4(0,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR3( -1,-1,1) 
      Hi4(0,-1,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  - 6.9314718055994530d-01*HR2(-1,1) 
      Hi4(0,1,-1,-1) =  
     $  - 2.4022650695910071d-01*HR1(-1) 
     $  - 2.4022650695910071d-01*HR1(1) 
      Hi4(0,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      Hi4(0,1,-1,1) = 
     $  - 5.8224052646501250d-01*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      Hi4(0,1,1,-1) = 
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
      Hi4(-1,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
      Hi4(-1,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      Hi4(-1,1,1,1) = 
     $  + 1.4047075598891257d+00*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      endif  
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for y > r2p1
** fillirr1dhplatinf is called by eval1dhplatinf after calling 
** fillirr1dhplat0 with argument r=1/y 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 3.2898681336964528d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      Hi2(0,1) = 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 3.2898681336964528d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,0,1) = 
     $  + 5.000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(0,1,1) = 
     $  + 1.6449340668482264d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 2.1646464674222763d+00 
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,1) = 
     $  + 2.1646464674222763d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - 2.0000000000000000d+00*HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
      HY4(0,1,1,1) = 
     $  - 5.1410353601279064d+00 
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - HX1(0) *HX3(0,1,1) 
     $  + 4.9348022005446793d+00*HX2(0,1) 
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,0,1) = 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
      Hi4(0,0,1,1) = 
     $  - 1.2020569031595942d+00 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      Hi4(0,1,1,1) = 
     $  + 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 1.6449340668482264d+00 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      HY3(0,-1,-1) = 
     $  + 1.2020569031595942d+00 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 1.8940656589944918d+00 
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0, -1) 
      HY4(0,0,-1,-1) = 
     $  - 1.8940656589944918d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX3(0,0,-1) 
     $  + HX4(0,0, -1,-1) 
     $  + 2.0000000000000000d+00*HX4(0,0,0,-1) 
      HY4(0,-1,-1,-1) = 
     $  + 1.0823232337111381d+00 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1) 
     $  + HX1(0) *HX3(0,-1,-1) 
     $  + HX1(0) *HX3(0,0,-1) 
     $  - HX4(0, -1,-1,-1) 
     $  - HX4(0,0, -1,-1) 
     $  - HX4(0,0,0, -1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 2.4674011002723396d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      Hi2(-1,1) = 
     $  - 6.9314718055994530d-01 
     $  + HX1( -1) 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 2.5190015545588625d+00 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      HY3(0,1,-1) = 
     $  + 4.3220869092982539d+00 
     $  + 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      HY3(-1,-1,1) = 
     $  - 2.7620719062289241d+00 
     $  + 2.4674011002723396d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  - HX1( -1)*HX2(0,-1) 
     $  - HX1( -1)*HX2(0,1) 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3( -1,-1,1) 
     $  + HX3(0, -1,-1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1) 
      HY3(-1,1,1) = 
     $  + 2.7620719062289241d+00 
     $  - 4.9348022005446793d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(-1,1) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3( -1,1,1) 
     $  - HX3(0, -1,1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,-1,1) = 
     $  + 8.2246703342411321d-01 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HX1(0) 
      Hi3(-1,-1,1) = 
     $  + 2.4022650695910071d-01 
     $  - 6.9314718055994530d-01*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1) 
     $  - HX1( -1)*HX1(0) 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(-1,1,1) = 
     $  + 1.8851605738073271d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 3.9234217222028759d+00
     $  + 2.5190015545588625d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX4(0,0, -1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,-1) = 
     $  - 4.1940025306306604d+00
     $  - 4.3220869092982539d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX4(0,0,0, -1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,-1,0,1) = 
     $  + 9.4703282949724591d-01
     $  + 1.8030853547393914d+00*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 3.2898681336964528d+00*HX2(0,-1)
     $  + HX4(0, -1,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,-1,1) = 
     $  + 2.5209599327464717d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,-1,1)
     $  + HX4(0, -1,0,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,1,-1) = 
     $  - 8.5266539820739622d+00
     $  - 5.5241438124578482d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,0,1)
     $  - HX4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,1,-1,-1) = 
     $  + 5.8027584430066521d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  - HX4(0,0, -1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,-1)
     $  - HX4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 6.2689427375197987d-01
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  - HX4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      HY4(0,1,-1,1) = 
     $  - 4.3326514514433017d+00
     $  - 1.3169446513992682d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - HX4(0, -1,0,1)
     $  - 3.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
     $  - HX4(0,1, -1,1)
      HY4(0,1,1,-1) = 
     $  - 1.5001934240460787d-01
     $  + 4.0790165576281924d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,1)
     $  - HX2(0, -1)*HX2(0,1)
     $  + 2.4674011002723396d+00*HX2(0,1)
     $  - 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  + HX4(0, -1,0,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - HX4(0,0,0, -1)
     $  + HX4(0,0,1, -1)
     $  - HX4(0,1,1, -1)
      HY4(-1,-1,-1,1) = 
     $  + 2.4278628067547031d+00
     $  - 2.7620719062289241d+00*HX1(-1)
     $  + 1.2337005501361698d+00*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX3(0,-1,-1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  + HX1( -1)*HX3(0,1,-1)
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX4( -1,-1,-1,1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  - HX4(0,0,1, -1)
     $  - HX4(0,1, -1,-1)
      HY4(-1,-1,1,1) = 
     $  + 2.0293560632083841d+00
     $  + 2.7620719062289241d+00*HX1(-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(-1)
     $  + 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX1(0)*HX2(0,-1)
     $  - HX1( -1)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX3(0,-1,1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  - HX1( -1)*HX3(0,1,1)
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(-1,-1,1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  + HX4( -1,-1,1,1)
     $  + HX4(0, -1,-1,1)
     $  + HX4(0, -1,1,-1)
     $  - HX4(0,0, -1,-1)
     $  + HX4(0,0, -1,1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  - HX4(0,0,1, -1)
     $  + HX4(0,0,1,1) 
     $  + HX4(0,1, -1,1)
     $  + HX4(0,1,1, -1)
      HY4(-1,1,1,1) = 
     $  - 6.4865749331714713d+00
     $  - 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(-1,1,1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  - 4.9348022005446793d+00*HX2(-1,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  + 4.9348022005446793d+00*HX2(0,1)
     $  + HX4( -1,1,1,1)
     $  - HX4(0, -1,1,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,-1,1) = 
     $  - 9.0154267736969571d-01 
     $  - 8.2246703342411321d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
      Hi4(0,-1,0,1) = 
     $  + 1.8030853547393914d+00 
     $  + 8.2246703342411321d-01*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - 2.0000000000000000d+00*HX3(0,0,-1) 
      Hi4(0,-1,-1,1) = 
     $  + 4.8170908494321862d-01 
     $  - 2.4022650695910071d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  + 6.9314718055994530d-01*HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      Hi4(0,-1,1,-1) = 
     $  + 5.7009070532142637d-01 
     $  + 4.8045301391820142d-01*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,-1) 
      Hi4(0,1,-1,-1) = 
     $  - 2.4022650695910071d-01*HX1(0) 
      Hi4(0,-1,1,1) = 
     $  - 2.7620719062289241d+00 
     $  - 1.8851605738073271d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      Hi4(0,1,-1,1) = 
     $  + 2.6736902858507163d+00 
     $  + 1.3029200473423146d+00*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  + 1.6666666666666665d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  + 6.9314718055994530d-01*HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      Hi4(0,1,1,-1) =  
     $  + 1.1401814106428527d+00 
     $  + 5.8224052646501250d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,1) 
      Hi4(-1,-1,-1,1) = 
     $  - 5.5504108664821579d-02
     $  + 2.4022650695910071d-01*HX1(-1)
     $  - 3.4657359027997265d-01*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(-1)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX1(0)
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
      Hi4(-1,-1,1,1) = 
     $  - 2.4532465311320902d+00
     $  + 1.8851605738073271d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX2(0,-1)
     $  - HX1( -1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX3( -1,-1,1)
     $  + HX3(0, -1,-1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1)
      Hi4(-1,1,1,1) = 
     $  - 5.5504108664821579d-02
     $  - 1.6449340668482264d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX2(-1,1)
     $  - HX1(0) *HX2(0,-1)
     $  - HX1(0) *HX2(0,1)
     $  + HX3( -1,1,1)
     $  - HX3(0, -1,1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluates the irreducible HPL for y =1
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
      if (nw.gt.2) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
      endif
      if (nw.gt.3) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
      endif
      endif
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 8.2246703342411321d-01
      if (nw.gt.2) then
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
      endif
      if (nw.gt.3) then
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
      endif
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
      if (nw.gt.2) then
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
      endif
      if (nw.gt.3) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
      endif
      endif
** (n1,n2) = (-1,1) -- completion endif 
      return
      end
************************************************************************
************************************************************************
*
*          FUNCTIONS for the Cross Section
*
*************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function PoleEPm2qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm2qq =  + Nc**(-1) * (
     &     - 4.D0
     &     - 8.D0*upx**(-4)
     &     + 16.D0*upx**(-3)
     &     - 8.D0*upx**(-2)
     &     - 16.D0*y*upx**(-4)
     &     + 32.D0*y*upx**(-3)
     &     - 24.D0*y*upx**(-2)
     &     + 8.D0*y*upx**(-1)
     &     - 8.D0*y**2*upx**(-4)
     &     + 16.D0*y**2*upx**(-3)
     &     - 8.D0*y**2*upx**(-2)
     &     )
      PoleEPm2qq = PoleEPm2qq + Nc * (
     &     + 8.D0
     &     + 16.D0*upx**(-4)
     &     - 32.D0*upx**(-3)
     &     + 16.D0*upx**(-2)
     &     + 32.D0*y*upx**(-4)
     &     - 64.D0*y*upx**(-3)
     &     + 48.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEPm2qq = PoleEPm2qq + Nc**3 * (
     &     - 4.D0
     &     - 8.D0*upx**(-4)
     &     + 16.D0*upx**(-3)
     &     - 8.D0*upx**(-2)
     &     - 16.D0*y*upx**(-4)
     &     + 32.D0*y*upx**(-3)
     &     - 24.D0*y*upx**(-2)
     &     + 8.D0*y*upx**(-1)
     &     - 8.D0*y**2*upx**(-4)
     &     + 16.D0*y**2*upx**(-3)
     &     - 8.D0*y**2*upx**(-2)
     &     )
      return
      end
*
*
*
      double precision function PoleEPm1qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm1qq =  + Nc**(-1) * (
     &     - 6.D0
     &     - 20.D0*upx**(-4)
     &     + 40.D0*upx**(-3)
     &     - 20.D0*upx**(-2)
     &     - 40.D0*y*upx**(-4)
     &     + 80.D0*y*upx**(-3)
     &     - 60.D0*y*upx**(-2)
     &     + 20.D0*y*upx**(-1)
     &     - 20.D0*y**2*upx**(-4)
     &     + 40.D0*y**2*upx**(-3)
     &     - 20.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + Nc * (
     &     + 12.D0
     &     + 40.D0*upx**(-4)
     &     - 80.D0*upx**(-3)
     &     + 40.D0*upx**(-2)
     &     + 80.D0*y*upx**(-4)
     &     - 160.D0*y*upx**(-3)
     &     + 120.D0*y*upx**(-2)
     &     - 40.D0*y*upx**(-1)
     &     + 40.D0*y**2*upx**(-4)
     &     - 80.D0*y**2*upx**(-3)
     &     + 40.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + Nc**3 * (
     &     - 6.D0
     &     - 20.D0*upx**(-4)
     &     + 40.D0*upx**(-3)
     &     - 20.D0*upx**(-2)
     &     - 40.D0*y*upx**(-4)
     &     + 80.D0*y*upx**(-3)
     &     - 60.D0*y*upx**(-2)
     &     + 20.D0*y*upx**(-1)
     &     - 20.D0*y**2*upx**(-4)
     &     + 40.D0*y**2*upx**(-3)
     &     - 20.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(-1)*Nc**(-1) * (
     &     + 8.D0
     &     + 16.D0*upx**(-4)
     &     - 32.D0*upx**(-3)
     &     + 16.D0*upx**(-2)
     &     + 32.D0*y*upx**(-4)
     &     - 64.D0*y*upx**(-3)
     &     + 48.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(-1)*Nc * (
     &     - 8.D0
     &     - 16.D0*upx**(-4)
     &     + 32.D0*upx**(-3)
     &     - 16.D0*upx**(-2)
     &     - 32.D0*y*upx**(-4)
     &     + 64.D0*y*upx**(-3)
     &     - 48.D0*y*upx**(-2)
     &     + 16.D0*y*upx**(-1)
     &     - 16.D0*y**2*upx**(-4)
     &     + 32.D0*y**2*upx**(-3)
     &     - 16.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(0)*Nc**(-1) * (
     &     - 8.D0*upx**(-5)
     &     + 12.D0*upx**(-4)
     &     - 2.D0*upx**(-3)
     &     - upx**(-2)
     &     - 9.D0/2.D0*upx**(-1)
     &     - 9.D0/2.D0*umx**(-1)
     &     - 16.D0*y*upx**(-5)
     &     + 24.D0*y*upx**(-4)
     &     - 12.D0*y*upx**(-3)
     &     + 2.D0*y*upx**(-2)
     &     + y*upx**(-1)
     &     + y*umx**(-1)
     &     - 8.D0*y**2*upx**(-5)
     &     + 12.D0*y**2*upx**(-4)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(0)*Nc**(-1) * (
     &     - 2.D0*y**2*upx**(-3)
     &     - y**2*upx**(-2)
     &     - 1.D0/2.D0*y**2*upx**(-1)
     &     - 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(0)*Nc * (
     &     + 8.D0*upx**(-5)
     &     - 12.D0*upx**(-4)
     &     + 2.D0*upx**(-3)
     &     + upx**(-2)
     &     + 9.D0/2.D0*upx**(-1)
     &     + 9.D0/2.D0*umx**(-1)
     &     + 16.D0*y*upx**(-5)
     &     - 24.D0*y*upx**(-4)
     &     + 12.D0*y*upx**(-3)
     &     - 2.D0*y*upx**(-2)
     &     - y*upx**(-1)
     &     - y*umx**(-1)
     &     + 8.D0*y**2*upx**(-5)
     &     - 12.D0*y**2*upx**(-4)
     &     )
      PoleEPm1qq = PoleEPm1qq + Hr1(0)*Nc * (
     &     + 2.D0*y**2*upx**(-3)
     &     + y**2*upx**(-2)
     &     + 1.D0/2.D0*y**2*upx**(-1)
     &     + 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEPm1qq = PoleEPm1qq + HAr1(-1)*Nc**(-1) * (
     &     + 16.D0
     &     + 32.D0*upx**(-4)
     &     - 64.D0*upx**(-3)
     &     + 32.D0*upx**(-2)
     &     + 64.D0*y*upx**(-4)
     &     - 128.D0*y*upx**(-3)
     &     + 96.D0*y*upx**(-2)
     &     - 32.D0*y*upx**(-1)
     &     + 32.D0*y**2*upx**(-4)
     &     - 64.D0*y**2*upx**(-3)
     &     + 32.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + HAr1(-1)*Nc * (
     &     - 24.D0
     &     - 48.D0*upx**(-4)
     &     + 96.D0*upx**(-3)
     &     - 48.D0*upx**(-2)
     &     - 96.D0*y*upx**(-4)
     &     + 192.D0*y*upx**(-3)
     &     - 144.D0*y*upx**(-2)
     &     + 48.D0*y*upx**(-1)
     &     - 48.D0*y**2*upx**(-4)
     &     + 96.D0*y**2*upx**(-3)
     &     - 48.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + HAr1(-1)*Nc**3 * (
     &     + 8.D0
     &     + 16.D0*upx**(-4)
     &     - 32.D0*upx**(-3)
     &     + 16.D0*upx**(-2)
     &     + 32.D0*y*upx**(-4)
     &     - 64.D0*y*upx**(-3)
     &     + 48.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + HBr1(-1)*Nc**(-1) * (
     &     - 16.D0
     &     - 32.D0*upx**(-4)
     &     + 64.D0*upx**(-3)
     &     - 32.D0*upx**(-2)
     &     - 64.D0*y*upx**(-4)
     &     + 128.D0*y*upx**(-3)
     &     - 96.D0*y*upx**(-2)
     &     + 32.D0*y*upx**(-1)
     &     - 32.D0*y**2*upx**(-4)
     &     + 64.D0*y**2*upx**(-3)
     &     - 32.D0*y**2*upx**(-2)
     &     )
      PoleEPm1qq = PoleEPm1qq + HBr1(-1)*Nc * (
     &     + 16.D0
     &     + 32.D0*upx**(-4)
     &     - 64.D0*upx**(-3)
     &     + 32.D0*upx**(-2)
     &     + 64.D0*y*upx**(-4)
     &     - 128.D0*y*upx**(-3)
     &     + 96.D0*y*upx**(-2)
     &     - 32.D0*y*upx**(-1)
     &     + 32.D0*y**2*upx**(-4)
     &     - 64.D0*y**2*upx**(-3)
     &     + 32.D0*y**2*upx**(-2)
     &     )
      return
      end
*
*
*
      double precision function PoleEP0qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1d0

* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEP0qq =  + Nc**(-1) * (
     &     - 14.D0
     &     - 48.D0*upx**(-4)
     &     + 96.D0*upx**(-3)
     &     - 48.D0*upx**(-2)
     &     + 32.D0*z2*upx**(-5)
     &     - 48.D0*z2*upx**(-4)
     &     - 8.D0*z2*upx**(-3)
     &     + 12.D0*z2*upx**(-2)
     &     + 10.D0*z2*upx**(-1)
     &     + 16.D0*z2*umx**(-3)
     &     - 24.D0*z2*umx**(-2)
     &     + 26.D0*z2*umx**(-1)
     &     - 96.D0*y*upx**(-4)
     &     + 192.D0*y*upx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**(-1) * (
     &     - 144.D0*y*upx**(-2)
     &     + 48.D0*y*upx**(-1)
     &     + 64.D0*y*z2*upx**(-5)
     &     - 96.D0*y*z2*upx**(-4)
     &     + 32.D0*y*z2*upx**(-3)
     &     - 4.D0*y*z2*upx**(-1)
     &     - 16.D0*y*z2*umx**(-3)
     &     + 24.D0*y*z2*umx**(-2)
     &     - 4.D0*y*z2*umx**(-1)
     &     - 48.D0*y**2*upx**(-4)
     &     + 96.D0*y**2*upx**(-3)
     &     - 48.D0*y**2*upx**(-2)
     &     + 32.D0*y**2*z2*upx**(-5)
     &     - 48.D0*y**2*z2*upx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**(-1) * (
     &     + 8.D0*y**2*z2*upx**(-3)
     &     + 4.D0*y**2*z2*upx**(-2)
     &     + 2.D0*y**2*z2*upx**(-1)
     &     + 2.D0*y**2*z2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Nl * (
     &     + 40.D0/9.D0
     &     + 80.D0/9.D0*upx**(-4)
     &     - 160.D0/9.D0*upx**(-3)
     &     + 80.D0/9.D0*upx**(-2)
     &     + 160.D0/9.D0*y*upx**(-4)
     &     - 320.D0/9.D0*y*upx**(-3)
     &     + 80.D0/3.D0*y*upx**(-2)
     &     - 80.D0/9.D0*y*upx**(-1)
     &     + 80.D0/9.D0*y**2*upx**(-4)
     &     - 160.D0/9.D0*y**2*upx**(-3)
     &     + 80.D0/9.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Nh * (
     &     + 40.D0/9.D0
     &     - 64.D0/3.D0*upx**(-6)
     &     + 64.D0*upx**(-5)
     &     - 496.D0/9.D0*upx**(-4)
     &     + 32.D0/9.D0*upx**(-3)
     &     - 16.D0/9.D0*upx**(-2)
     &     + 32.D0/3.D0*upx**(-1)
     &     - 128.D0/3.D0*y*upx**(-6)
     &     + 128.D0*y*upx**(-5)
     &     - 1184.D0/9.D0*y*upx**(-4)
     &     + 448.D0/9.D0*y*upx**(-3)
     &     + 16.D0/3.D0*y*upx**(-2)
     &     - 80.D0/9.D0*y*upx**(-1)
     &     - 64.D0/3.D0*y**2*upx**(-6)
     &     )
      PoleEP0qq = PoleEP0qq + Nh * (
     &     + 64.D0*y**2*upx**(-5)
     &     - 496.D0/9.D0*y**2*upx**(-4)
     &     + 32.D0/9.D0*y**2*upx**(-3)
     &     + 80.D0/9.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Nc * (
     &     - 88.D0/9.D0
     &     + 256.D0/9.D0*upx**(-4)
     &     - 512.D0/9.D0*upx**(-3)
     &     + 238.D0/9.D0*upx**(-2)
     &     + 2.D0*upx**(-1)
     &     - 2.D0*umx**(-2)
     &     + 2.D0*umx**(-1)
     &     - 32.D0*z2*upx**(-5)
     &     + 40.D0*z2*upx**(-4)
     &     + 22.D0*z2*upx**(-3)
     &     - 17.D0*z2*upx**(-2)
     &     - 33.D0/4.D0*z2*upx**(-1)
     &     - 12.D0*z2*umx**(-5)
     &     + 30.D0*z2*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Nc * (
     &     - 47.D0*z2*umx**(-3)
     &     + 81.D0/2.D0*z2*umx**(-2)
     &     - 145.D0/4.D0*z2*umx**(-1)
     &     + 512.D0/9.D0*y*upx**(-4)
     &     - 1024.D0/9.D0*y*upx**(-3)
     &     + 244.D0/3.D0*y*upx**(-2)
     &     - 220.D0/9.D0*y*upx**(-1)
     &     + 4.D0*y*umx**(-2)
     &     - 4.D0*y*umx**(-1)
     &     - 64.D0*y*z2*upx**(-5)
     &     + 80.D0*y*z2*upx**(-4)
     &     - 8.D0*y*z2*upx**(-3)
     &     - 8.D0*y*z2*upx**(-2)
     &     + 11.D0/2.D0*y*z2*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Nc * (
     &     + 24.D0*y*z2*umx**(-5)
     &     - 60.D0*y*z2*umx**(-4)
     &     + 70.D0*y*z2*umx**(-3)
     &     - 45.D0*y*z2*umx**(-2)
     &     + 11.D0/2.D0*y*z2*umx**(-1)
     &     + 256.D0/9.D0*y**2*upx**(-4)
     &     - 512.D0/9.D0*y**2*upx**(-3)
     &     + 238.D0/9.D0*y**2*upx**(-2)
     &     + 2.D0*y**2*upx**(-1)
     &     - 2.D0*y**2*umx**(-2)
     &     + 2.D0*y**2*umx**(-1)
     &     - 32.D0*y**2*z2*upx**(-5)
     &     + 40.D0*y**2*z2*upx**(-4)
     &     + 2.D0*y**2*z2*upx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Nc * (
     &     - 3.D0*y**2*z2*upx**(-2)
     &     - 9.D0/4.D0*y**2*z2*upx**(-1)
     &     - 12.D0*y**2*z2*umx**(-5)
     &     + 30.D0*y**2*z2*umx**(-4)
     &     - 19.D0*y**2*z2*umx**(-3)
     &     - 3.D0/2.D0*y**2*z2*umx**(-2)
     &     - 9.D0/4.D0*y**2*z2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**2*Nl * (
     &     - 40.D0/9.D0
     &     - 80.D0/9.D0*upx**(-4)
     &     + 160.D0/9.D0*upx**(-3)
     &     - 80.D0/9.D0*upx**(-2)
     &     - 160.D0/9.D0*y*upx**(-4)
     &     + 320.D0/9.D0*y*upx**(-3)
     &     - 80.D0/3.D0*y*upx**(-2)
     &     + 80.D0/9.D0*y*upx**(-1)
     &     - 80.D0/9.D0*y**2*upx**(-4)
     &     + 160.D0/9.D0*y**2*upx**(-3)
     &     - 80.D0/9.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**2*Nh * (
     &     - 40.D0/9.D0
     &     + 64.D0/3.D0*upx**(-6)
     &     - 64.D0*upx**(-5)
     &     + 496.D0/9.D0*upx**(-4)
     &     - 32.D0/9.D0*upx**(-3)
     &     + 16.D0/9.D0*upx**(-2)
     &     - 32.D0/3.D0*upx**(-1)
     &     + 128.D0/3.D0*y*upx**(-6)
     &     - 128.D0*y*upx**(-5)
     &     + 1184.D0/9.D0*y*upx**(-4)
     &     - 448.D0/9.D0*y*upx**(-3)
     &     - 16.D0/3.D0*y*upx**(-2)
     &     + 80.D0/9.D0*y*upx**(-1)
     &     + 64.D0/3.D0*y**2*upx**(-6)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**2*Nh * (
     &     - 64.D0*y**2*upx**(-5)
     &     + 496.D0/9.D0*y**2*upx**(-4)
     &     - 32.D0/9.D0*y**2*upx**(-3)
     &     - 80.D0/9.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**3 * (
     &     + 214.D0/9.D0
     &     + 176.D0/9.D0*upx**(-4)
     &     - 352.D0/9.D0*upx**(-3)
     &     + 194.D0/9.D0*upx**(-2)
     &     - 2.D0*upx**(-1)
     &     + 2.D0*umx**(-2)
     &     - 2.D0*umx**(-1)
     &     + 8.D0*z2*upx**(-4)
     &     - 14.D0*z2*upx**(-3)
     &     + 5.D0*z2*upx**(-2)
     &     - 7.D0/4.D0*z2*upx**(-1)
     &     + 12.D0*z2*umx**(-5)
     &     - 30.D0*z2*umx**(-4)
     &     + 31.D0*z2*umx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**3 * (
     &     - 33.D0/2.D0*z2*umx**(-2)
     &     + 41.D0/4.D0*z2*umx**(-1)
     &     + 352.D0/9.D0*y*upx**(-4)
     &     - 704.D0/9.D0*y*upx**(-3)
     &     + 188.D0/3.D0*y*upx**(-2)
     &     - 212.D0/9.D0*y*upx**(-1)
     &     - 4.D0*y*umx**(-2)
     &     + 4.D0*y*umx**(-1)
     &     + 16.D0*y*z2*upx**(-4)
     &     - 24.D0*y*z2*upx**(-3)
     &     + 8.D0*y*z2*upx**(-2)
     &     - 3.D0/2.D0*y*z2*upx**(-1)
     &     - 24.D0*y*z2*umx**(-5)
     &     + 60.D0*y*z2*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**3 * (
     &     - 54.D0*y*z2*umx**(-3)
     &     + 21.D0*y*z2*umx**(-2)
     &     - 3.D0/2.D0*y*z2*umx**(-1)
     &     + 176.D0/9.D0*y**2*upx**(-4)
     &     - 352.D0/9.D0*y**2*upx**(-3)
     &     + 194.D0/9.D0*y**2*upx**(-2)
     &     - 2.D0*y**2*upx**(-1)
     &     + 2.D0*y**2*umx**(-2)
     &     - 2.D0*y**2*umx**(-1)
     &     + 8.D0*y**2*z2*upx**(-4)
     &     - 10.D0*y**2*z2*upx**(-3)
     &     - y**2*z2*upx**(-2)
     &     + 1.D0/4.D0*y**2*z2*upx**(-1)
     &     + 12.D0*y**2*z2*umx**(-5)
     &     )
      PoleEP0qq = PoleEP0qq + Nc**3 * (
     &     - 30.D0*y**2*z2*umx**(-4)
     &     + 19.D0*y**2*z2*umx**(-3)
     &     + 3.D0/2.D0*y**2*z2*umx**(-2)
     &     + 1.D0/4.D0*y**2*z2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + DLog(tm**(-2)*w**2)*Nl * (
     &     + 8.D0/3.D0
     &     + 16.D0/3.D0*upx**(-4)
     &     - 32.D0/3.D0*upx**(-3)
     &     + 16.D0/3.D0*upx**(-2)
     &     + 32.D0/3.D0*y*upx**(-4)
     &     - 64.D0/3.D0*y*upx**(-3)
     &     + 16.D0*y*upx**(-2)
     &     - 16.D0/3.D0*y*upx**(-1)
     &     + 16.D0/3.D0*y**2*upx**(-4)
     &     - 32.D0/3.D0*y**2*upx**(-3)
     &     + 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + DLog(tm**(-2)*w**2)*Nc * (
     &     - 44.D0/3.D0
     &     - 88.D0/3.D0*upx**(-4)
     &     + 176.D0/3.D0*upx**(-3)
     &     - 88.D0/3.D0*upx**(-2)
     &     - 176.D0/3.D0*y*upx**(-4)
     &     + 352.D0/3.D0*y*upx**(-3)
     &     - 88.D0*y*upx**(-2)
     &     + 88.D0/3.D0*y*upx**(-1)
     &     - 88.D0/3.D0*y**2*upx**(-4)
     &     + 176.D0/3.D0*y**2*upx**(-3)
     &     - 88.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + DLog(tm**(-2)*w**2)*Nc**2*Nl * (
     &     - 8.D0/3.D0
     &     - 16.D0/3.D0*upx**(-4)
     &     + 32.D0/3.D0*upx**(-3)
     &     - 16.D0/3.D0*upx**(-2)
     &     - 32.D0/3.D0*y*upx**(-4)
     &     + 64.D0/3.D0*y*upx**(-3)
     &     - 16.D0*y*upx**(-2)
     &     + 16.D0/3.D0*y*upx**(-1)
     &     - 16.D0/3.D0*y**2*upx**(-4)
     &     + 32.D0/3.D0*y**2*upx**(-3)
     &     - 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + DLog(tm**(-2)*w**2)*Nc**3 * (
     &     + 44.D0/3.D0
     &     + 88.D0/3.D0*upx**(-4)
     &     - 176.D0/3.D0*upx**(-3)
     &     + 88.D0/3.D0*upx**(-2)
     &     + 176.D0/3.D0*y*upx**(-4)
     &     - 352.D0/3.D0*y*upx**(-3)
     &     + 88.D0*y*upx**(-2)
     &     - 88.D0/3.D0*y*upx**(-1)
     &     + 88.D0/3.D0*y**2*upx**(-4)
     &     - 176.D0/3.D0*y**2*upx**(-3)
     &     + 88.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc**(-1) * (
     &     - 12.D0
     &     + 24.D0*upx**(-4)
     &     - 48.D0*upx**(-3)
     &     + 24.D0*upx**(-2)
     &     - 32.D0*umx**(-2)
     &     + 32.D0*umx**(-1)
     &     + 48.D0*y*upx**(-4)
     &     - 96.D0*y*upx**(-3)
     &     + 72.D0*y*upx**(-2)
     &     - 24.D0*y*upx**(-1)
     &     + 32.D0*y*umx**(-2)
     &     - 32.D0*y*umx**(-1)
     &     + 24.D0*y**2*upx**(-4)
     &     - 48.D0*y**2*upx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc**(-1) * (
     &     + 24.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nl * (
     &     - 16.D0/3.D0
     &     - 32.D0/3.D0*upx**(-4)
     &     + 64.D0/3.D0*upx**(-3)
     &     - 32.D0/3.D0*upx**(-2)
     &     - 64.D0/3.D0*y*upx**(-4)
     &     + 128.D0/3.D0*y*upx**(-3)
     &     - 32.D0*y*upx**(-2)
     &     + 32.D0/3.D0*y*upx**(-1)
     &     - 32.D0/3.D0*y**2*upx**(-4)
     &     + 64.D0/3.D0*y**2*upx**(-3)
     &     - 32.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc * (
     &     + 76.D0/3.D0
     &     - 16.D0/3.D0*upx**(-4)
     &     + 32.D0/3.D0*upx**(-3)
     &     - 40.D0/3.D0*upx**(-2)
     &     + 8.D0*upx**(-1)
     &     + 24.D0*umx**(-4)
     &     - 48.D0*umx**(-3)
     &     + 72.D0*umx**(-2)
     &     - 48.D0*umx**(-1)
     &     - 32.D0/3.D0*y*upx**(-4)
     &     + 64.D0/3.D0*y*upx**(-3)
     &     - 32.D0*y*upx**(-2)
     &     + 64.D0/3.D0*y*upx**(-1)
     &     - 48.D0*y*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc * (
     &     + 96.D0*y*umx**(-3)
     &     - 96.D0*y*umx**(-2)
     &     + 48.D0*y*umx**(-1)
     &     - 16.D0/3.D0*y**2*upx**(-4)
     &     + 32.D0/3.D0*y**2*upx**(-3)
     &     - 40.D0/3.D0*y**2*upx**(-2)
     &     + 8.D0*y**2*upx**(-1)
     &     + 24.D0*y**2*umx**(-4)
     &     - 48.D0*y**2*umx**(-3)
     &     + 16.D0*y**2*umx**(-2)
     &     + 8.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc**2*Nl * (
     &     + 16.D0/3.D0
     &     + 32.D0/3.D0*upx**(-4)
     &     - 64.D0/3.D0*upx**(-3)
     &     + 32.D0/3.D0*upx**(-2)
     &     + 64.D0/3.D0*y*upx**(-4)
     &     - 128.D0/3.D0*y*upx**(-3)
     &     + 32.D0*y*upx**(-2)
     &     - 32.D0/3.D0*y*upx**(-1)
     &     + 32.D0/3.D0*y**2*upx**(-4)
     &     - 64.D0/3.D0*y**2*upx**(-3)
     &     + 32.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc**3 * (
     &     - 40.D0/3.D0
     &     - 56.D0/3.D0*upx**(-4)
     &     + 112.D0/3.D0*upx**(-3)
     &     - 32.D0/3.D0*upx**(-2)
     &     - 8.D0*upx**(-1)
     &     - 24.D0*umx**(-4)
     &     + 48.D0*umx**(-3)
     &     - 40.D0*umx**(-2)
     &     + 16.D0*umx**(-1)
     &     - 112.D0/3.D0*y*upx**(-4)
     &     + 224.D0/3.D0*y*upx**(-3)
     &     - 40.D0*y*upx**(-2)
     &     + 8.D0/3.D0*y*upx**(-1)
     &     + 48.D0*y*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*Nc**3 * (
     &     - 96.D0*y*umx**(-3)
     &     + 64.D0*y*umx**(-2)
     &     - 16.D0*y*umx**(-1)
     &     - 56.D0/3.D0*y**2*upx**(-4)
     &     + 112.D0/3.D0*y**2*upx**(-3)
     &     - 32.D0/3.D0*y**2*upx**(-2)
     &     - 8.D0*y**2*upx**(-1)
     &     - 24.D0*y**2*umx**(-4)
     &     + 48.D0*y**2*umx**(-3)
     &     - 16.D0*y**2*umx**(-2)
     &     - 8.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*HAr1(-1)*Nc**(-1) * (
     &     - 16.D0
     &     - 64.D0*upx**(-4)
     &     + 128.D0*upx**(-3)
     &     - 64.D0*upx**(-2)
     &     - 128.D0*y*upx**(-4)
     &     + 256.D0*y*upx**(-3)
     &     - 160.D0*y*upx**(-2)
     &     + 32.D0*y*upx**(-1)
     &     - 64.D0*y**2*upx**(-4)
     &     + 128.D0*y**2*upx**(-3)
     &     - 64.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*HAr1(-1)*Nc * (
     &     + 24.D0
     &     + 96.D0*upx**(-4)
     &     - 192.D0*upx**(-3)
     &     + 96.D0*upx**(-2)
     &     + 192.D0*y*upx**(-4)
     &     - 384.D0*y*upx**(-3)
     &     + 240.D0*y*upx**(-2)
     &     - 48.D0*y*upx**(-1)
     &     + 96.D0*y**2*upx**(-4)
     &     - 192.D0*y**2*upx**(-3)
     &     + 96.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*HAr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 32.D0*upx**(-4)
     &     + 64.D0*upx**(-3)
     &     - 32.D0*upx**(-2)
     &     - 64.D0*y*upx**(-4)
     &     + 128.D0*y*upx**(-3)
     &     - 80.D0*y*upx**(-2)
     &     + 16.D0*y*upx**(-1)
     &     - 32.D0*y**2*upx**(-4)
     &     + 64.D0*y**2*upx**(-3)
     &     - 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*HBr1(-1)*Nc**(-1) * (
     &     + 48.D0
     &     + 64.D0*upx**(-4)
     &     - 128.D0*upx**(-3)
     &     + 128.D0*upx**(-2)
     &     - 64.D0*upx**(-1)
     &     + 128.D0*y*upx**(-4)
     &     - 256.D0*y*upx**(-3)
     &     + 224.D0*y*upx**(-2)
     &     - 96.D0*y*upx**(-1)
     &     + 64.D0*y**2*upx**(-4)
     &     - 128.D0*y**2*upx**(-3)
     &     + 64.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(-1)*HBr1(-1)*Nc * (
     &     - 48.D0
     &     - 64.D0*upx**(-4)
     &     + 128.D0*upx**(-3)
     &     - 128.D0*upx**(-2)
     &     + 64.D0*upx**(-1)
     &     - 128.D0*y*upx**(-4)
     &     + 256.D0*y*upx**(-3)
     &     - 224.D0*y*upx**(-2)
     &     + 96.D0*y*upx**(-1)
     &     - 64.D0*y**2*upx**(-4)
     &     + 128.D0*y**2*upx**(-3)
     &     - 64.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**(-1) * (
     &     + 8.D0
     &     - 16.D0*upx**(-5)
     &     + 28.D0*upx**(-4)
     &     - 10.D0*upx**(-3)
     &     - upx**(-2)
     &     - 17.D0/2.D0*upx**(-1)
     &     + 16.D0*umx**(-2)
     &     - 25.D0/2.D0*umx**(-1)
     &     - 32.D0*y*upx**(-5)
     &     + 56.D0*y*upx**(-4)
     &     - 36.D0*y*upx**(-3)
     &     + 10.D0*y*upx**(-2)
     &     + y*upx**(-1)
     &     - 16.D0*y*umx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**(-1) * (
     &     + 17.D0*y*umx**(-1)
     &     - 16.D0*y**2*upx**(-5)
     &     + 28.D0*y**2*upx**(-4)
     &     - 10.D0*y**2*upx**(-3)
     &     - y**2*upx**(-2)
     &     - 1.D0/2.D0*y**2*upx**(-1)
     &     - 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nl * (
     &     + 8.D0/3.D0
     &     + 16.D0/3.D0*upx**(-4)
     &     - 32.D0/3.D0*upx**(-3)
     &     + 16.D0/3.D0*upx**(-2)
     &     + 32.D0/3.D0*y*upx**(-4)
     &     - 64.D0/3.D0*y*upx**(-3)
     &     + 16.D0*y*upx**(-2)
     &     - 16.D0/3.D0*y*upx**(-1)
     &     + 16.D0/3.D0*y**2*upx**(-4)
     &     - 32.D0/3.D0*y**2*upx**(-3)
     &     + 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nh * (
     &     - 8.D0/3.D0
     &     - 64.D0/3.D0*upx**(-7)
     &     + 224.D0/3.D0*upx**(-6)
     &     - 256.D0/3.D0*upx**(-5)
     &     + 80.D0/3.D0*upx**(-4)
     &     + 32.D0/3.D0*upx**(-2)
     &     - 128.D0/3.D0*y*upx**(-7)
     &     + 448.D0/3.D0*y*upx**(-6)
     &     - 192.D0*y*upx**(-5)
     &     + 320.D0/3.D0*y*upx**(-4)
     &     - 32.D0/3.D0*y*upx**(-3)
     &     - 16.D0*y*upx**(-2)
     &     + 16.D0/3.D0*y*upx**(-1)
     &     - 64.D0/3.D0*y**2*upx**(-7)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nh * (
     &     + 224.D0/3.D0*y**2*upx**(-6)
     &     - 256.D0/3.D0*y**2*upx**(-5)
     &     + 80.D0/3.D0*y**2*upx**(-4)
     &     + 32.D0/3.D0*y**2*upx**(-3)
     &     - 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc * (
     &     - 44.D0/3.D0
     &     + 16.D0*upx**(-5)
     &     - 112.D0/3.D0*upx**(-4)
     &     + 86.D0/3.D0*upx**(-3)
     &     - 13.D0/3.D0*upx**(-2)
     &     + 9.D0/2.D0*upx**(-1)
     &     - 12.D0*umx**(-4)
     &     + 24.D0*umx**(-3)
     &     - 36.D0*umx**(-2)
     &     + 41.D0/2.D0*umx**(-1)
     &     + 32.D0*y*upx**(-5)
     &     - 224.D0/3.D0*y*upx**(-4)
     &     + 220.D0/3.D0*y*upx**(-3)
     &     - 30.D0*y*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc * (
     &     + 1.D0/3.D0*y*upx**(-1)
     &     + 24.D0*y*umx**(-4)
     &     - 48.D0*y*umx**(-3)
     &     + 48.D0*y*umx**(-2)
     &     - 25.D0*y*umx**(-1)
     &     + 16.D0*y**2*upx**(-5)
     &     - 112.D0/3.D0*y**2*upx**(-4)
     &     + 86.D0/3.D0*y**2*upx**(-3)
     &     - 13.D0/3.D0*y**2*upx**(-2)
     &     - 7.D0/2.D0*y**2*upx**(-1)
     &     - 12.D0*y**2*umx**(-4)
     &     + 24.D0*y**2*umx**(-3)
     &     - 8.D0*y**2*umx**(-2)
     &     - 7.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**2*Nl * (
     &     - 8.D0/3.D0
     &     - 16.D0/3.D0*upx**(-4)
     &     + 32.D0/3.D0*upx**(-3)
     &     - 16.D0/3.D0*upx**(-2)
     &     - 32.D0/3.D0*y*upx**(-4)
     &     + 64.D0/3.D0*y*upx**(-3)
     &     - 16.D0*y*upx**(-2)
     &     + 16.D0/3.D0*y*upx**(-1)
     &     - 16.D0/3.D0*y**2*upx**(-4)
     &     + 32.D0/3.D0*y**2*upx**(-3)
     &     - 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**2*Nh * (
     &     + 8.D0/3.D0
     &     + 64.D0/3.D0*upx**(-7)
     &     - 224.D0/3.D0*upx**(-6)
     &     + 256.D0/3.D0*upx**(-5)
     &     - 80.D0/3.D0*upx**(-4)
     &     - 32.D0/3.D0*upx**(-2)
     &     + 128.D0/3.D0*y*upx**(-7)
     &     - 448.D0/3.D0*y*upx**(-6)
     &     + 192.D0*y*upx**(-5)
     &     - 320.D0/3.D0*y*upx**(-4)
     &     + 32.D0/3.D0*y*upx**(-3)
     &     + 16.D0*y*upx**(-2)
     &     - 16.D0/3.D0*y*upx**(-1)
     &     + 64.D0/3.D0*y**2*upx**(-7)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**2*Nh * (
     &     - 224.D0/3.D0*y**2*upx**(-6)
     &     + 256.D0/3.D0*y**2*upx**(-5)
     &     - 80.D0/3.D0*y**2*upx**(-4)
     &     - 32.D0/3.D0*y**2*upx**(-3)
     &     + 16.D0/3.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**3 * (
     &     + 20.D0/3.D0
     &     + 28.D0/3.D0*upx**(-4)
     &     - 56.D0/3.D0*upx**(-3)
     &     + 16.D0/3.D0*upx**(-2)
     &     + 4.D0*upx**(-1)
     &     + 12.D0*umx**(-4)
     &     - 24.D0*umx**(-3)
     &     + 20.D0*umx**(-2)
     &     - 8.D0*umx**(-1)
     &     + 56.D0/3.D0*y*upx**(-4)
     &     - 112.D0/3.D0*y*upx**(-3)
     &     + 20.D0*y*upx**(-2)
     &     - 4.D0/3.D0*y*upx**(-1)
     &     - 24.D0*y*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*Nc**3 * (
     &     + 48.D0*y*umx**(-3)
     &     - 32.D0*y*umx**(-2)
     &     + 8.D0*y*umx**(-1)
     &     + 28.D0/3.D0*y**2*upx**(-4)
     &     - 56.D0/3.D0*y**2*upx**(-3)
     &     + 16.D0/3.D0*y**2*upx**(-2)
     &     + 4.D0*y**2*upx**(-1)
     &     + 12.D0*y**2*umx**(-4)
     &     - 24.D0*y**2*umx**(-3)
     &     + 8.D0*y**2*umx**(-2)
     &     + 4.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*HAr1(-1)*Nc**(-1) * (
     &     + 8.D0
     &     + 32.D0*upx**(-4)
     &     - 64.D0*upx**(-3)
     &     + 32.D0*upx**(-2)
     &     + 64.D0*y*upx**(-4)
     &     - 128.D0*y*upx**(-3)
     &     + 80.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 32.D0*y**2*upx**(-4)
     &     - 64.D0*y**2*upx**(-3)
     &     + 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*HAr1(-1)*Nc * (
     &     - 12.D0
     &     - 48.D0*upx**(-4)
     &     + 96.D0*upx**(-3)
     &     - 48.D0*upx**(-2)
     &     - 96.D0*y*upx**(-4)
     &     + 192.D0*y*upx**(-3)
     &     - 120.D0*y*upx**(-2)
     &     + 24.D0*y*upx**(-1)
     &     - 48.D0*y**2*upx**(-4)
     &     + 96.D0*y**2*upx**(-3)
     &     - 48.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     + 4.D0
     &     + 16.D0*upx**(-4)
     &     - 32.D0*upx**(-3)
     &     + 16.D0*upx**(-2)
     &     + 32.D0*y*upx**(-4)
     &     - 64.D0*y*upx**(-3)
     &     + 40.D0*y*upx**(-2)
     &     - 8.D0*y*upx**(-1)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*HBr1(-1)*Nc**(-1) * (
     &     - 24.D0
     &     - 32.D0*upx**(-4)
     &     + 64.D0*upx**(-3)
     &     - 64.D0*upx**(-2)
     &     + 32.D0*upx**(-1)
     &     - 64.D0*y*upx**(-4)
     &     + 128.D0*y*upx**(-3)
     &     - 112.D0*y*upx**(-2)
     &     + 48.D0*y*upx**(-1)
     &     - 32.D0*y**2*upx**(-4)
     &     + 64.D0*y**2*upx**(-3)
     &     - 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr1(0)*HBr1(-1)*Nc * (
     &     + 24.D0
     &     + 32.D0*upx**(-4)
     &     - 64.D0*upx**(-3)
     &     + 64.D0*upx**(-2)
     &     - 32.D0*upx**(-1)
     &     + 64.D0*y*upx**(-4)
     &     - 128.D0*y*upx**(-3)
     &     + 112.D0*y*upx**(-2)
     &     - 48.D0*y*upx**(-1)
     &     + 32.D0*y**2*upx**(-4)
     &     - 64.D0*y**2*upx**(-3)
     &     + 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,-1)*Nc**(-1) * (
     &     - 48.D0
     &     - 32.D0*upx**(-4)
     &     + 64.D0*upx**(-3)
     &     - 96.D0*upx**(-2)
     &     + 64.D0*upx**(-1)
     &     - 64.D0*y*upx**(-4)
     &     + 128.D0*y*upx**(-3)
     &     - 160.D0*y*upx**(-2)
     &     + 96.D0*y*upx**(-1)
     &     - 32.D0*y**2*upx**(-4)
     &     + 64.D0*y**2*upx**(-3)
     &     - 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,-1)*Nc * (
     &     + 40.D0
     &     + 64.D0*upx**(-2)
     &     - 64.D0*upx**(-1)
     &     + 80.D0*y*upx**(-2)
     &     - 80.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,-1)*Nc**3 * (
     &     + 8.D0
     &     + 32.D0*upx**(-4)
     &     - 64.D0*upx**(-3)
     &     + 32.D0*upx**(-2)
     &     + 64.D0*y*upx**(-4)
     &     - 128.D0*y*upx**(-3)
     &     + 80.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 32.D0*y**2*upx**(-4)
     &     - 64.D0*y**2*upx**(-3)
     &     + 32.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,0)*Nc**(-1) * (
     &     + 24.D0
     &     + 16.D0*upx**(-4)
     &     - 32.D0*upx**(-3)
     &     + 48.D0*upx**(-2)
     &     - 32.D0*upx**(-1)
     &     + 32.D0*y*upx**(-4)
     &     - 64.D0*y*upx**(-3)
     &     + 80.D0*y*upx**(-2)
     &     - 48.D0*y*upx**(-1)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,0)*Nc * (
     &     - 20.D0
     &     - 32.D0*upx**(-2)
     &     + 32.D0*upx**(-1)
     &     - 40.D0*y*upx**(-2)
     &     + 40.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(-1,0)*Nc**3 * (
     &     - 4.D0
     &     - 16.D0*upx**(-4)
     &     + 32.D0*upx**(-3)
     &     - 16.D0*upx**(-2)
     &     - 32.D0*y*upx**(-4)
     &     + 64.D0*y*upx**(-3)
     &     - 40.D0*y*upx**(-2)
     &     + 8.D0*y*upx**(-1)
     &     - 16.D0*y**2*upx**(-4)
     &     + 32.D0*y**2*upx**(-3)
     &     - 16.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc**(-1) * (
     &     + 8.D0
     &     + 16.D0*upx**(-4)
     &     + 16.D0*upx**(-1)
     &     - 32.D0*umx**(-3)
     &     + 48.D0*umx**(-2)
     &     - 16.D0*umx**(-1)
     &     + 32.D0*y*upx**(-4)
     &     - 32.D0*y*upx**(-3)
     &     + 32.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     + 32.D0*y*umx**(-3)
     &     - 48.D0*y*umx**(-2)
     &     + 16.D0*y**2*upx**(-4)
     &     - 32.D0*y**2*upx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc**(-1) * (
     &     + 16.D0*y**2*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc * (
     &     - 8.D0
     &     - 28.D0*upx**(-3)
     &     + 10.D0*upx**(-2)
     &     - 39.D0/2.D0*upx**(-1)
     &     + 24.D0*umx**(-5)
     &     - 60.D0*umx**(-4)
     &     + 94.D0*umx**(-3)
     &     - 81.D0*umx**(-2)
     &     + 73.D0/2.D0*umx**(-1)
     &     - 16.D0*y*upx**(-3)
     &     - 16.D0*y*upx**(-2)
     &     + 13.D0*y*upx**(-1)
     &     - 48.D0*y*umx**(-5)
     &     + 120.D0*y*umx**(-4)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc * (
     &     - 140.D0*y*umx**(-3)
     &     + 90.D0*y*umx**(-2)
     &     - 3.D0*y*umx**(-1)
     &     + 12.D0*y**2*upx**(-3)
     &     - 18.D0*y**2*upx**(-2)
     &     + 1.D0/2.D0*y**2*upx**(-1)
     &     + 24.D0*y**2*umx**(-5)
     &     - 60.D0*y**2*umx**(-4)
     &     + 38.D0*y**2*umx**(-3)
     &     + 3.D0*y**2*umx**(-2)
     &     + 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc**3 * (
     &     - 16.D0*upx**(-4)
     &     + 28.D0*upx**(-3)
     &     - 10.D0*upx**(-2)
     &     + 7.D0/2.D0*upx**(-1)
     &     - 24.D0*umx**(-5)
     &     + 60.D0*umx**(-4)
     &     - 62.D0*umx**(-3)
     &     + 33.D0*umx**(-2)
     &     - 41.D0/2.D0*umx**(-1)
     &     - 32.D0*y*upx**(-4)
     &     + 48.D0*y*upx**(-3)
     &     - 16.D0*y*upx**(-2)
     &     + 3.D0*y*upx**(-1)
     &     + 48.D0*y*umx**(-5)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,-1)*Nc**3 * (
     &     - 120.D0*y*umx**(-4)
     &     + 108.D0*y*umx**(-3)
     &     - 42.D0*y*umx**(-2)
     &     + 3.D0*y*umx**(-1)
     &     - 16.D0*y**2*upx**(-4)
     &     + 20.D0*y**2*upx**(-3)
     &     + 2.D0*y**2*upx**(-2)
     &     - 1.D0/2.D0*y**2*upx**(-1)
     &     - 24.D0*y**2*umx**(-5)
     &     + 60.D0*y**2*umx**(-4)
     &     - 38.D0*y**2*umx**(-3)
     &     - 3.D0*y**2*umx**(-2)
     &     - 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc**(-1) * (
     &     - 8.D0*upx**(-5)
     &     + 12.D0*upx**(-4)
     &     - 18.D0*upx**(-3)
     &     + 7.D0*upx**(-2)
     &     - 25.D0/2.D0*upx**(-1)
     &     + 16.D0*umx**(-3)
     &     - 24.D0*umx**(-2)
     &     + 7.D0/2.D0*umx**(-1)
     &     - 16.D0*y*upx**(-5)
     &     + 24.D0*y*upx**(-4)
     &     - 28.D0*y*upx**(-3)
     &     + 10.D0*y*upx**(-2)
     &     + y*upx**(-1)
     &     - 16.D0*y*umx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc**(-1) * (
     &     + 24.D0*y*umx**(-2)
     &     + y*umx**(-1)
     &     - 8.D0*y**2*upx**(-5)
     &     + 12.D0*y**2*upx**(-4)
     &     - 2.D0*y**2*upx**(-3)
     &     - y**2*upx**(-2)
     &     - 1.D0/2.D0*y**2*upx**(-1)
     &     - 1.D0/2.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc * (
     &     + 8.D0*upx**(-5)
     &     - 20.D0*upx**(-4)
     &     + 32.D0*upx**(-3)
     &     - 12.D0*upx**(-2)
     &     + 57.D0/4.D0*upx**(-1)
     &     - 12.D0*umx**(-5)
     &     + 30.D0*umx**(-4)
     &     - 47.D0*umx**(-3)
     &     + 81.D0/2.D0*umx**(-2)
     &     - 55.D0/4.D0*umx**(-1)
     &     + 16.D0*y*upx**(-5)
     &     - 40.D0*y*upx**(-4)
     &     + 52.D0*y*upx**(-3)
     &     - 18.D0*y*upx**(-2)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc * (
     &     + 1.D0/2.D0*y*upx**(-1)
     &     + 24.D0*y*umx**(-5)
     &     - 60.D0*y*umx**(-4)
     &     + 70.D0*y*umx**(-3)
     &     - 45.D0*y*umx**(-2)
     &     + 1.D0/2.D0*y*umx**(-1)
     &     + 8.D0*y**2*upx**(-5)
     &     - 20.D0*y**2*upx**(-4)
     &     + 12.D0*y**2*upx**(-3)
     &     + 2.D0*y**2*upx**(-2)
     &     + 1.D0/4.D0*y**2*upx**(-1)
     &     - 12.D0*y**2*umx**(-5)
     &     + 30.D0*y**2*umx**(-4)
     &     - 19.D0*y**2*umx**(-3)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc * (
     &     - 3.D0/2.D0*y**2*umx**(-2)
     &     + 1.D0/4.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc**3 * (
     &     + 8.D0*upx**(-4)
     &     - 14.D0*upx**(-3)
     &     + 5.D0*upx**(-2)
     &     - 7.D0/4.D0*upx**(-1)
     &     + 12.D0*umx**(-5)
     &     - 30.D0*umx**(-4)
     &     + 31.D0*umx**(-3)
     &     - 33.D0/2.D0*umx**(-2)
     &     + 41.D0/4.D0*umx**(-1)
     &     + 16.D0*y*upx**(-4)
     &     - 24.D0*y*upx**(-3)
     &     + 8.D0*y*upx**(-2)
     &     - 3.D0/2.D0*y*upx**(-1)
     &     - 24.D0*y*umx**(-5)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(0,0)*Nc**3 * (
     &     + 60.D0*y*umx**(-4)
     &     - 54.D0*y*umx**(-3)
     &     + 21.D0*y*umx**(-2)
     &     - 3.D0/2.D0*y*umx**(-1)
     &     + 8.D0*y**2*upx**(-4)
     &     - 10.D0*y**2*upx**(-3)
     &     - y**2*upx**(-2)
     &     + 1.D0/4.D0*y**2*upx**(-1)
     &     + 12.D0*y**2*umx**(-5)
     &     - 30.D0*y**2*umx**(-4)
     &     + 19.D0*y**2*umx**(-3)
     &     + 3.D0/2.D0*y**2*umx**(-2)
     &     + 1.D0/4.D0*y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(1,0)*Nc**(-1) * (
     &     + 8.D0
     &     - 16.D0*upx**(-5)
     &     + 40.D0*upx**(-4)
     &     - 36.D0*upx**(-3)
     &     + 14.D0*upx**(-2)
     &     - 9.D0*upx**(-1)
     &     - 9.D0*umx**(-1)
     &     - 32.D0*y*upx**(-5)
     &     + 80.D0*y*upx**(-4)
     &     - 88.D0*y*upx**(-3)
     &     + 52.D0*y*upx**(-2)
     &     - 14.D0*y*upx**(-1)
     &     + 2.D0*y*umx**(-1)
     &     - 16.D0*y**2*upx**(-5)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(1,0)*Nc**(-1) * (
     &     + 40.D0*y**2*upx**(-4)
     &     - 36.D0*y**2*upx**(-3)
     &     + 14.D0*y**2*upx**(-2)
     &     - y**2*upx**(-1)
     &     - y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(1,0)*Nc * (
     &     - 8.D0
     &     + 16.D0*upx**(-5)
     &     - 40.D0*upx**(-4)
     &     + 36.D0*upx**(-3)
     &     - 14.D0*upx**(-2)
     &     + 9.D0*upx**(-1)
     &     + 9.D0*umx**(-1)
     &     + 32.D0*y*upx**(-5)
     &     - 80.D0*y*upx**(-4)
     &     + 88.D0*y*upx**(-3)
     &     - 52.D0*y*upx**(-2)
     &     + 14.D0*y*upx**(-1)
     &     - 2.D0*y*umx**(-1)
     &     + 16.D0*y**2*upx**(-5)
     &     )
      PoleEP0qq = PoleEP0qq + Hr2(1,0)*Nc * (
     &     - 40.D0*y**2*upx**(-4)
     &     + 36.D0*y**2*upx**(-3)
     &     - 14.D0*y**2*upx**(-2)
     &     + y**2*upx**(-1)
     &     + y**2*umx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr1(-1)*Nc**(-1) * (
     &     - 8.D0
     &     + 8.D0*y**(-1)*upx**(-2)
     &     - 8.D0*y**(-1)*upx**(-1)
     &     - 8.D0*y**(-1)
     &     + 16.D0*upx**(-2)
     &     - 16.D0*upx**(-1)
     &     + 8.D0*y*upx**(-2)
     &     - 8.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr1(-1)*Nc * (
     &     + 12.D0
     &     - 12.D0*y**(-1)*upx**(-2)
     &     + 12.D0*y**(-1)*upx**(-1)
     &     + 12.D0*y**(-1)
     &     - 24.D0*upx**(-2)
     &     + 24.D0*upx**(-1)
     &     - 12.D0*y*upx**(-2)
     &     + 12.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr1(-1)*Nc**3 * (
     &     - 4.D0
     &     + 4.D0*y**(-1)*upx**(-2)
     &     - 4.D0*y**(-1)*upx**(-1)
     &     - 4.D0*y**(-1)
     &     + 8.D0*upx**(-2)
     &     - 8.D0*upx**(-1)
     &     + 4.D0*y*upx**(-2)
     &     - 4.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(-1,-1)*Nc**(-1) * (
     &     - 16.D0
     &     - 32.D0*y*upx**(-2)
     &     + 32.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(-1,-1)*Nc * (
     &     + 24.D0
     &     + 48.D0*y*upx**(-2)
     &     - 48.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(-1,-1)*Nc**3 * (
     &     - 8.D0
     &     - 16.D0*y*upx**(-2)
     &     + 16.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(0,-1)*Nc**(-1) * (
     &     + 8.D0
     &     + 16.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(0,-1)*Nc * (
     &     - 12.D0
     &     - 24.D0*y*upx**(-2)
     &     + 24.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HAr2(0,-1)*Nc**3 * (
     &     + 4.D0
     &     + 8.D0*y*upx**(-2)
     &     - 8.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr1(-1)*Nc**(-1) * (
     &     + 16.D0
     &     - 8.D0*z**(-1)*upx**(-2)
     &     + 8.D0*z**(-1)*upx**(-1)
     &     + 8.D0*z**(-1)
     &     + 8.D0*y*upx**(-2)
     &     - 8.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr1(-1)*Nc * (
     &     - 16.D0
     &     + 8.D0*z**(-1)*upx**(-2)
     &     - 8.D0*z**(-1)*upx**(-1)
     &     - 8.D0*z**(-1)
     &     - 8.D0*y*upx**(-2)
     &     + 8.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr2(-1,-1)*Nc**(-1) * (
     &     - 16.D0
     &     - 64.D0*upx**(-2)
     &     + 64.D0*upx**(-1)
     &     - 32.D0*y*upx**(-2)
     &     + 32.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr2(-1,-1)*Nc * (
     &     + 16.D0
     &     + 64.D0*upx**(-2)
     &     - 64.D0*upx**(-1)
     &     + 32.D0*y*upx**(-2)
     &     - 32.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr2(0,-1)*Nc**(-1) * (
     &     + 8.D0
     &     + 32.D0*upx**(-2)
     &     - 32.D0*upx**(-1)
     &     + 16.D0*y*upx**(-2)
     &     - 16.D0*y*upx**(-1)
     &     )
      PoleEP0qq = PoleEP0qq + HBr2(0,-1)*Nc * (
     &     - 8.D0
     &     - 32.D0*upx**(-2)
     &     + 32.D0*upx**(-1)
     &     - 16.D0*y*upx**(-2)
     &     + 16.D0*y*upx**(-1)
     &     )
      return
      end

*
*******************************************************************************

************************************************************************
************************************************************************
*
*          FUNCTIONS for the Cross Section
*
*************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function PoleEPm2gg(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm2gg =  + Nc**2 * (
     &     - 16.D0
     &     - 64.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 80.D0*upy**(-1)
     &     - 64.D0*upz**(-2)
     &     + 80.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEPm2gg = PoleEPm2gg + Nc**4 * (
     &     + 16.D0
     &     + 32.D0*upy**(-2)
     &     - 40.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm2gg = PoleEPm2gg + 32.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
      return
      end
*
*
*
      double precision function PoleEPm1gg(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm1gg =  + Nc**(-2) * (
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Nc**(-1)*Nl * (
     &     - 32.D0/3.D0*upy**(-2)
     &     - 64.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 40.D0/3.D0*upy**(-1)
     &     - 32.D0/3.D0*upz**(-2)
     &     + 40.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Nc*Nl * (
     &     + 16.D0/3.D0
     &     + 64.D0/3.D0*upy**(-2)
     &     + 64.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 80.D0/3.D0*upy**(-1)
     &     + 64.D0/3.D0*upz**(-2)
     &     - 80.D0/3.D0*upz**(-1)
     &     + 32.D0/3.D0*ypzp2**(-2)
     &     + 32.D0/3.D0*ypzp2**(-1)
     &     - 16.D0/3.D0*z*upy**(-1)
     &     + 64.D0/3.D0*z*ypzp2**(-2)
     &     - 32.D0/3.D0*z*ypzp2**(-1)
     &     + 32.D0/3.D0*z**2*ypzp2**(-2)
     &     - 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Nc**2 * (
     &     - 136.D0/3.D0
     &     - 496.D0/3.D0*upy**(-2)
     &     - 448.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 524.D0/3.D0*upy**(-1)
     &     - 496.D0/3.D0*upz**(-2)
     &     + 524.D0/3.D0*upz**(-1)
     &     - 176.D0/3.D0*ypzp2**(-2)
     &     - 368.D0/3.D0*ypzp2**(-1)
     &     + 28.D0/3.D0*z*upy**(-1)
     &     - 352.D0/3.D0*z*ypzp2**(-2)
     &     + 176.D0/3.D0*z*ypzp2**(-1)
     &     - 176.D0/3.D0*z**2*ypzp2**(-2)
     &     + 28.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Nc**3*Nl * (
     &     - 16.D0/3.D0
     &     - 32.D0/3.D0*upy**(-2)
     &     + 40.D0/3.D0*upy**(-1)
     &     - 32.D0/3.D0*upz**(-2)
     &     + 40.D0/3.D0*upz**(-1)
     &     - 32.D0/3.D0*ypzp2**(-2)
     &     - 32.D0/3.D0*ypzp2**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     - 64.D0/3.D0*z*ypzp2**(-2)
     &     + 32.D0/3.D0*z*ypzp2**(-1)
     &     - 32.D0/3.D0*z**2*ypzp2**(-2)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Nc**4 * (
     &     + 64.D0/3.D0
     &     + 224.D0/3.D0*upy**(-2)
     &     - 232.D0/3.D0*upy**(-1)
     &     + 224.D0/3.D0*upz**(-2)
     &     - 232.D0/3.D0*upz**(-1)
     &     + 128.D0/3.D0*ypzp2**(-2)
     &     + 320.D0/3.D0*ypzp2**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     + 256.D0/3.D0*z*ypzp2**(-2)
     &     - 128.D0/3.D0*z*ypzp2**(-1)
     &     + 128.D0/3.D0*z**2*ypzp2**(-2)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(-1)*Nc**2 * (
     &     + 16.D0
     &     + 32.D0*upy**(-2)
     &     - 40.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(-1)*Nc**4 * (
     &     - 16.D0
     &     - 32.D0*upy**(-2)
     &     + 40.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0)*Nc**(-2) * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0)*Nc**2 * (
     &     - 8.D0
     &     - 16.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 32.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0)*Nc**2 * (
     &     + 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     - 64.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 16.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0)*Nc**4 * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + Hr1(0) * (
     &     + 16.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + HAr1(-1)*Nc**2 * (
     &     + 64.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 80.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + HAr1(-1)*Nc**4 * (
     &     - 32.D0*upy**(-2)
     &     + 40.D0*upy**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1gg = PoleEPm1gg + HAr1(-1) * (
     &     - 32.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + HBr1(-1)*Nc**2 * (
     &     + 16.D0
     &     + 32.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 64.D0*upz**(-2)
     &     - 80.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 8.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + HBr1(-1)*Nc**4 * (
     &     - 16.D0
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + HBr1(-1) * (
     &     - 32.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1gg = PoleEPm1gg + 24.D0
     &     + 320.D0/3.D0*upy**(-2)
     &     + 544.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 352.D0/3.D0*upy**(-1)
     &     + 320.D0/3.D0*upz**(-2)
     &     - 352.D0/3.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 32.D0/3.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 32.D0/3.D0*y*upz**(-1)
      return
      end
*
*
*
      double precision function PoleEP0gg(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,zp2,zm3,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	zp2 = 2.d0+z
	zm3 = z-3.d0
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEP0gg =  + Nc**(-2) * (
     &     - 8.D0
     &     - 4.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)
     &     - 24.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 42.D0*upy**(-1)
     &     - 24.D0*upz**(-2)
     &     + 42.D0*upz**(-1)
     &     - 8.D0*z2*upy**(-3)
     &     - 48.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*upy**(-2)
     &     - 192.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**(-2) * (
     &     + 186.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z2*upy**(-1)
     &     - 8.D0*z2*upz**(-3)
     &     - 48.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*upz**(-2)
     &     + 186.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z2*upz**(-1)
     &     - 124.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2
     &     + 8.D0*z*upy**(-2)
     &     - 2.D0*z*upy**(-1)
     &     - 8.D0*z*z2*upy**(-3)
     &     + 48.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**(-2) * (
     &     + 6.D0*z*z2*upy**(-1)
     &     - 2.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)
     &     - 2.D0*y*upz**(-1)
     &     - 8.D0*y*z2*upz**(-3)
     &     + 48.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y*z2*upz**(-1)
     &     - 2.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**(-1)*Nl * (
     &     - 16.D0/3.D0
     &     - 16.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-1)
     &     - 16.D0/3.D0*z*upy**(-1)
     &     - 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc*Nl * (
     &     + 28.D0/3.D0*upy**(-1)
     &     + 28.D0/3.D0*upz**(-1)
     &     - 32.D0/3.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     + 32.D0/3.D0*z*upy**(-1)
     &     - 64.D0/3.D0*z*ypzp2**(-2)
     &     + 32.D0/3.D0*z*ypzp2**(-1)
     &     - 32.D0/3.D0*z**2*ypzp2**(-2)
     &     + 32.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc*Nh * (
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 4.D0/3.D0*upy**(-1)
     &     - 4.D0/3.D0*upz**(-1)
     &     + 128.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*ypzp2**(-1)
     &     + 24.D0*z2*upy**(-1)*upz**(-1)
     &     - 96.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     - 4.D0*y**(-1)*upz**(-1)
     &     - 4.D0*y**(-1)*zp2**(-1)
     &     - 6.D0*y**(-1)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     - 4.D0*z**(-1)*ypzp2**(-1)
     &     - 6.D0*z**(-1)
     &     - 88.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)*zm3**(-1)
     &     + 170.D0/3.D0*upy**(-1)
     &     - 88.D0*upz**(-2)
     &     - 16.D0*upz**(-1)*ypzm2**(-1)
     &     + 170.D0/3.D0*upz**(-1)
     &     + 4.D0*zp2**(-1)*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     + 16.D0*zm3**(-1)*ypzm2**(-1)
     &     + 20.D0*ypzm2**(-1)
     &     + 80.D0/3.D0*ypzp2**(-2)
     &     - 172.D0*ypzp2**(-1)
     &     - 8.D0*z2*upy**(-3)
     &     + 200.D0*z2*upy**(-2)
     &     + 192.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 120.D0*z2*upy**(-1)*upz**(-1)
     &     - 192.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 96.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 78.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 170.D0*z2*upy**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     - 8.D0*z2*upz**(-3)
     &     + 200.D0*z2*upz**(-2)
     &     - 192.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 96.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 78.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 170.D0*z2*upz**(-1)
     &     + 192.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 192.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 96.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     + 240.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 88.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z2*ypzp2**(-2)
     &     - 72.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z2*ypzp2**(-1)
     &     + 52.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 38.D0*z2
     &     + 8.D0*z*upy**(-2)
     &     - 146.D0/3.D0*z*upy**(-1)
     &     - 8.D0*z*ypzm2**(-1)
     &     + 160.D0/3.D0*z*ypzp2**(-2)
     &     - 104.D0/3.D0*z*ypzp2**(-1)
     &     - 8.D0*z*z2*upy**(-3)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     - 16.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 46.D0*z*z2*upy**(-1)
     &     - 96.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 256.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 160.D0*z*z2*ypzp2**(-2)
     &     + 304.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 80.D0*z*z2*ypzp2**(-1)
     &     - 48.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzm2**(-1)
     &     + 80.D0/3.D0*z**2*ypzp2**(-2)
     &     - 4.D0*z**2*ypzp2**(-1)
     &     - 6.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**2 * (
     &     - 8.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z**2*z2*ypzp2**(-2)
     &     + 88.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)
     &     - 146.D0/3.D0*y*upz**(-1)
     &     - 8.D0*y*z2*upz**(-3)
     &     - 16.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 46.D0*y*z2*upz**(-1)
     &     + 32.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**3*Nl * (
     &     + 16.D0/3.D0
     &     - 4.D0*upy**(-1)
     &     - 4.D0*upz**(-1)
     &     + 32.D0/3.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 16.D0/3.D0*z*upy**(-1)
     &     + 64.D0/3.D0*z*ypzp2**(-2)
     &     - 32.D0/3.D0*z*ypzp2**(-1)
     &     + 32.D0/3.D0*z**2*ypzp2**(-2)
     &     - 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**3*Nh * (
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 4.D0/3.D0*upy**(-1)
     &     + 4.D0/3.D0*upz**(-1)
     &     - 128.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*ypzp2**(-1)
     &     - 24.D0*z2*upy**(-1)*upz**(-1)
     &     + 96.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**4 * (
     &     - 82.D0/3.D0
     &     + 2.D0*y**(-1)*zp2**(-1)
     &     + 2.D0*y**(-1)
     &     + 2.D0*z**(-1)*ypzp2**(-1)
     &     + 2.D0*z**(-1)
     &     + 32.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*zm3**(-1)
     &     - 10.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     + 16.D0*upz**(-1)*ypzm2**(-1)
     &     - 10.D0*upz**(-1)
     &     - 2.D0*zp2**(-1)*ypzp2**(-1)
     &     - 16.D0*zm3**(-1)*ypzm2**(-1)
     &     - 20.D0*ypzm2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**4 * (
     &     - 128.D0/3.D0*ypzp2**(-2)
     &     + 120.D0*ypzp2**(-1)
     &     - 80.D0*z2*upy**(-2)
     &     + 192.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 96.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 90.D0*z2*upy**(-1)
     &     - 80.D0*z2*upz**(-2)
     &     + 192.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 96.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 10.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**4 * (
     &     + 90.D0*z2*upz**(-1)
     &     - 192.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 192.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 96.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 240.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 88.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 80.D0*z2*ypzp2**(-2)
     &     - 24.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z2*ypzp2**(-1)
     &     - 8.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z2
     &     )
      PoleEP0gg = PoleEP0gg + Nc**4 * (
     &     + 70.D0/3.D0*z*upy**(-1)
     &     + 8.D0*z*ypzm2**(-1)
     &     - 256.D0/3.D0*z*ypzp2**(-2)
     &     + 152.D0/3.D0*z*ypzp2**(-1)
     &     + 8.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*z*z2*upy**(-1)
     &     + 96.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 160.D0*z*z2*ypzp2**(-2)
     &     - 48.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z*z2*ypzp2**(-1)
     &     + 10.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzm2**(-1)
     &     - 128.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Nc**4 * (
     &     + 4.D0*z**2*ypzp2**(-1)
     &     + 2.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 80.D0*z**2*z2*ypzp2**(-2)
     &     - 24.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 70.D0/3.D0*y*upz**(-1)
     &     + 8.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*y*z2*upz**(-1)
     &     - 6.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2)*Nc**(-1)*Nl * (
     &     + 32.D0/3.D0*upy**(-2)
     &     + 64.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 40.D0/3.D0*upy**(-1)
     &     + 32.D0/3.D0*upz**(-2)
     &     - 40.D0/3.D0*upz**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2)*Nc*Nl * (
     &     - 16.D0/3.D0
     &     - 64.D0/3.D0*upy**(-2)
     &     - 64.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 80.D0/3.D0*upy**(-1)
     &     - 64.D0/3.D0*upz**(-2)
     &     + 80.D0/3.D0*upz**(-1)
     &     - 32.D0/3.D0*ypzp2**(-2)
     &     - 32.D0/3.D0*ypzp2**(-1)
     &     + 16.D0/3.D0*z*upy**(-1)
     &     - 64.D0/3.D0*z*ypzp2**(-2)
     &     + 32.D0/3.D0*z*ypzp2**(-1)
     &     - 32.D0/3.D0*z**2*ypzp2**(-2)
     &     + 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2)*Nc**2 * (
     &     + 88.D0/3.D0
     &     + 352.D0/3.D0*upy**(-2)
     &     + 352.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 440.D0/3.D0*upy**(-1)
     &     + 352.D0/3.D0*upz**(-2)
     &     - 440.D0/3.D0*upz**(-1)
     &     + 176.D0/3.D0*ypzp2**(-2)
     &     + 176.D0/3.D0*ypzp2**(-1)
     &     - 88.D0/3.D0*z*upy**(-1)
     &     + 352.D0/3.D0*z*ypzp2**(-2)
     &     - 176.D0/3.D0*z*ypzp2**(-1)
     &     + 176.D0/3.D0*z**2*ypzp2**(-2)
     &     - 88.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2)*Nc**3*Nl * (
     &     + 16.D0/3.D0
     &     + 32.D0/3.D0*upy**(-2)
     &     - 40.D0/3.D0*upy**(-1)
     &     + 32.D0/3.D0*upz**(-2)
     &     - 40.D0/3.D0*upz**(-1)
     &     + 32.D0/3.D0*ypzp2**(-2)
     &     + 32.D0/3.D0*ypzp2**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     + 64.D0/3.D0*z*ypzp2**(-2)
     &     - 32.D0/3.D0*z*ypzp2**(-1)
     &     + 32.D0/3.D0*z**2*ypzp2**(-2)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2)*Nc**4 * (
     &     - 88.D0/3.D0
     &     - 176.D0/3.D0*upy**(-2)
     &     + 220.D0/3.D0*upy**(-1)
     &     - 176.D0/3.D0*upz**(-2)
     &     + 220.D0/3.D0*upz**(-1)
     &     - 176.D0/3.D0*ypzp2**(-2)
     &     - 176.D0/3.D0*ypzp2**(-1)
     &     + 44.D0/3.D0*z*upy**(-1)
     &     - 352.D0/3.D0*z*ypzp2**(-2)
     &     + 176.D0/3.D0*z*ypzp2**(-1)
     &     - 176.D0/3.D0*z**2*ypzp2**(-2)
     &     + 44.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Dlog(tm**(-2)*w**2) * (
     &     - 176.D0/3.D0*upy**(-2)
     &     - 352.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 220.D0/3.D0*upy**(-1)
     &     - 176.D0/3.D0*upz**(-2)
     &     + 220.D0/3.D0*upz**(-1)
     &     + 44.D0/3.D0*z*upy**(-1)
     &     + 44.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*Nc**2 * (
     &     - 12.D0
     &     + 192.D0*upy**(-1)*zm3**(-2)
     &     + 64.D0*upy**(-1)*zm3**(-1)
     &     + 12.D0*upy**(-1)
     &     + 192.D0*upz**(-1)*ypzm2**(-2)
     &     + 64.D0*upz**(-1)*ypzm2**(-1)
     &     + 12.D0*upz**(-1)
     &     - 192.D0*zm3**(-2)*ypzm2**(-1)
     &     - 192.D0*zm3**(-1)*ypzm2**(-2)
     &     - 64.D0*zm3**(-1)*ypzm2**(-1)
     &     - 240.D0*ypzm2**(-2)
     &     - 48.D0*ypzm2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*Nc**2 * (
     &     + 96.D0*z*ypzm2**(-2)
     &     + 16.D0*z*ypzm2**(-1)
     &     - 32.D0*z*ypzp2**(-1)
     &     - 48.D0*z**2*ypzm2**(-2)
     &     + 16.D0*z**2*ypzm2**(-1)
     &     - 16.D0*z**2*ypzp2**(-1)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*Nc**4 * (
     &     + 12.D0
     &     - 192.D0*upy**(-1)*zm3**(-2)
     &     - 64.D0*upy**(-1)*zm3**(-1)
     &     - 12.D0*upy**(-1)
     &     - 192.D0*upz**(-1)*ypzm2**(-2)
     &     - 64.D0*upz**(-1)*ypzm2**(-1)
     &     - 12.D0*upz**(-1)
     &     + 192.D0*zm3**(-2)*ypzm2**(-1)
     &     + 192.D0*zm3**(-1)*ypzm2**(-2)
     &     + 64.D0*zm3**(-1)*ypzm2**(-1)
     &     + 240.D0*ypzm2**(-2)
     &     + 48.D0*ypzm2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*Nc**4 * (
     &     - 96.D0*z*ypzm2**(-2)
     &     - 16.D0*z*ypzm2**(-1)
     &     + 32.D0*z*ypzp2**(-1)
     &     + 48.D0*z**2*ypzm2**(-2)
     &     - 16.D0*z**2*ypzm2**(-1)
     &     + 16.D0*z**2*ypzp2**(-1)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*HAr1(-1)*Nc**2 * (
     &     - 16.D0
     &     - 64.D0*upy**(-2)
     &     + 40.D0*upy**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     + 8.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*HAr1(-1)*Nc**4 * (
     &     + 16.D0
     &     + 64.D0*upy**(-2)
     &     - 40.D0*upy**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     - 8.D0*z*upy**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*HBr1(-1)*Nc**2 * (
     &     - 16.D0
     &     - 64.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(-1)*HBr1(-1)*Nc**4 * (
     &     + 16.D0
     &     + 64.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**(-2) * (
     &     + 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc*Nh * (
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 256.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**2 * (
     &     + 6.D0
     &     - 96.D0*upy**(-1)*zm3**(-2)
     &     - 32.D0*upy**(-1)*zm3**(-1)
     &     - 14.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*upy**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-2)
     &     - 32.D0*upz**(-1)*ypzm2**(-1)
     &     - 14.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*upz**(-1)
     &     + 96.D0*zm3**(-2)*ypzm2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-2)
     &     + 32.D0*zm3**(-1)*ypzm2**(-1)
     &     + 120.D0*ypzm2**(-2)
     &     + 24.D0*ypzm2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**2 * (
     &     + 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)
     &     + 2.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 48.D0*z*ypzm2**(-2)
     &     - 8.D0*z*ypzm2**(-1)
     &     + 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 22.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*ypzm2**(-2)
     &     - 8.D0*z**2*ypzm2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**2 * (
     &     - 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)
     &     + 2.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     + 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**3*Nh * (
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 256.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**4 * (
     &     - 6.D0
     &     + 96.D0*upy**(-1)*zm3**(-2)
     &     + 32.D0*upy**(-1)*zm3**(-1)
     &     + 6.D0*upy**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-2)
     &     + 32.D0*upz**(-1)*ypzm2**(-1)
     &     + 6.D0*upz**(-1)
     &     - 96.D0*zm3**(-2)*ypzm2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-2)
     &     - 32.D0*zm3**(-1)*ypzm2**(-1)
     &     - 120.D0*ypzm2**(-2)
     &     - 24.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*Nc**4 * (
     &     + 48.D0*z*ypzm2**(-2)
     &     + 8.D0*z*ypzm2**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 24.D0*z**2*ypzm2**(-2)
     &     + 8.D0*z**2*ypzm2**(-1)
     &     - 8.D0*z**2*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HAr1(-1)*Nc**(-2) * (
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 84.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 88.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HAr1(-1)*Nc**2 * (
     &     + 8.D0
     &     + 32.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)
     &     - 88.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HAr1(-1)*Nc**2 * (
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HAr1(-1)*Nc**4 * (
     &     - 8.D0
     &     - 32.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HAr1(-1) * (
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 52.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HBr1(-1)*Nc**(-2) * (
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 88.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 84.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HBr1(-1)*Nc**2 * (
     &     + 8.D0
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 88.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)
     &     - 32.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HBr1(-1)*Nc**2 * (
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HBr1(-1)*Nc**4 * (
     &     - 8.D0
     &     - 32.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0)*HBr1(-1) * (
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 52.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0) * (
     &     - 2.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr1(0) * (
     &     + 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,-1)*Nc**2 * (
     &     + 40.D0*upy**(-1)
     &     + 40.D0*upz**(-1)
     &     - 64.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,-1)*Nc**4 * (
     &     - 40.D0*upy**(-1)
     &     - 40.D0*upz**(-1)
     &     + 64.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,0)*Nc**(-2) * (
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 172.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 172.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 72.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,0)*Nc**2 * (
     &     - 128.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 120.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)
     &     + 120.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     - 64.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)
     &     - 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,0)*Nc**2 * (
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,0)*Nc**4 * (
     &     + 20.D0*upy**(-1)
     &     + 20.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(-1,0) * (
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 52.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 52.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**2 * (
     &     + 384.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)
     &     + 384.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     - 384.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 384.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 480.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 176.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**2 * (
     &     + 32.D0*ypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 192.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**2 * (
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**4 * (
     &     - 384.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upy**(-1)
     &     - 384.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 384.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 384.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 480.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 176.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**4 * (
     &     - 32.D0*ypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 192.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,-1)*Nc**4 * (
     &     + 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**(-2) * (
     &     - 8.D0
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     + 18.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*upy**(-1)
     &     + 18.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*upz**(-1)
     &     + 4.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z*upy**(-1)
     &     - 10.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y*upz**(-1)
     &     - 10.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc*Nh * (
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**2 * (
     &     + 2.D0
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     - 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 38.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upy**(-1)
     &     - 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 38.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upz**(-1)
     &     + 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 240.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**2 * (
     &     + 88.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 72.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 12.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z*upy**(-1)
     &     - 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**2 * (
     &     + 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y*upz**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**3*Nh * (
     &     + 8.D0*upy**(-1)*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**4 * (
     &     + 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)
     &     + 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upz**(-1)
     &     - 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 240.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 88.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**4 * (
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     + 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0)*Nc**4 * (
     &     - 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0) * (
     &     + 6.D0
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     + 10.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*upy**(-1)
     &     + 10.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(0,0) * (
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     + 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(1,0)*Nc**(-2) * (
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 136.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 136.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(1,0)*Nc**2 * (
     &     - 128.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(1,0) * (
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 72.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 72.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + Hr2(1,0) * (
     &     + 32.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*Nc**(-2) * (
     &     - 4.D0
     &     + 4.D0*y**(-2)*upz**(-1)
     &     + 2.D0*y**(-2)
     &     + 4.D0*y**(-1)*upz**(-1)
     &     + 4.D0*y**(-1)
     &     + 2.D0*y**(-1)*z
     &     + 24.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     - 8.D0*upz**(-1)
     &     - 8.D0*z*upy**(-2)
     &     - 8.D0*z*upy**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*Nc**2 * (
     &     + 20.D0
     &     + 4.D0*y**(-2)*upz**(-1)
     &     + 4.D0*y**(-2)*zp2**(-1)
     &     + 6.D0*y**(-2)
     &     - 12.D0*y**(-1)*upz**(-1)
     &     - 4.D0*y**(-1)*zp2**(-2)
     &     - 4.D0*y**(-1)*zp2**(-1)
     &     - 8.D0*y**(-1)
     &     + 6.D0*y**(-1)*z
     &     + 88.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 48.D0*upy**(-1)
     &     + 4.D0*zp2**(-2)*ypzp2**(-1)
     &     + 4.D0*zp2**(-1)*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*Nc**2 * (
     &     - 16.D0*ypzp2**(-2)
     &     + 36.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-2)
     &     + 8.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 12.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*Nc**4 * (
     &     - 2.D0*y**(-2)*zp2**(-1)
     &     - 2.D0*y**(-2)
     &     + 2.D0*y**(-1)*zp2**(-2)
     &     + 6.D0*y**(-1)*zp2**(-1)
     &     + 6.D0*y**(-1)
     &     - 2.D0*y**(-1)*z
     &     - 32.D0*upy**(-2)
     &     + 8.D0*upy**(-1)
     &     - 2.D0*zp2**(-2)*ypzp2**(-1)
     &     - 6.D0*zp2**(-1)*ypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 18.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*Nc**4 * (
     &     - 10.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*HBr1(-1)*Nc**2 * (
     &     - 64.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     - 64.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1)*HBr1(-1) * (
     &     + 64.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 64.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1) * (
     &     - 16.D0
     &     - 8.D0*y**(-2)*upz**(-1)
     &     - 2.D0*y**(-2)*zp2**(-1)
     &     - 6.D0*y**(-2)
     &     + 8.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)*zp2**(-2)
     &     - 2.D0*y**(-1)*zp2**(-1)
     &     - 2.D0*y**(-1)
     &     - 6.D0*y**(-1)*z
     &     - 80.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 80.D0*upy**(-1)
     &     + 8.D0*upz**(-1)
     &     - 2.D0*zp2**(-2)*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr1(-1) * (
     &     + 2.D0*zp2**(-1)*ypzp2**(-1)
     &     - 18.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-2)
     &     + 8.D0*z*upy**(-1)
     &     - 2.D0*z*ypzp2**(-1)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(-1,-1)*Nc**2 * (
     &     + 16.D0
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 80.D0*upy**(-1)
     &     + 40.D0*upz**(-1)
     &     - 64.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(-1,-1)*Nc**4 * (
     &     - 16.D0
     &     - 40.D0*upy**(-1)
     &     + 64.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(-1,-1) * (
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     - 40.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(0,-1)*Nc**(-2) * (
     &     - 8.D0
     &     - 8.D0*upy**(-3)
     &     - 16.D0*upy**(-2)
     &     + 24.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upy**(-1)
     &     - 16.D0*upz**(-1)
     &     - 8.D0*z*upy**(-3)
     &     - 4.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(0,-1)*Nc**2 * (
     &     - 8.D0
     &     - 8.D0*upy**(-3)
     &     - 8.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 36.D0*upy**(-1)
     &     - 32.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-3)
     &     - 12.D0*z*upy**(-1)
     &     + 20.D0*z*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(0,-1)*Nc**4 * (
     &     + 8.D0
     &     + 20.D0*upy**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HAr2(0,-1) * (
     &     + 8.D0
     &     + 16.D0*upy**(-3)
     &     + 24.D0*upy**(-2)
     &     - 56.D0*upy**(-1)*upz**(-1)
     &     + 12.D0*upy**(-1)
     &     + 32.D0*upz**(-1)*ypzp2**(-1)
     &     + 16.D0*upz**(-1)
     &     + 12.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-3)
     &     + 12.D0*z*upy**(-1)
     &     - 4.D0*z*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1)*Nc**(-2) * (
     &     - 4.D0
     &     + 4.D0*z**(-2)*upy**(-1)
     &     + 2.D0*z**(-2)
     &     + 4.D0*z**(-1)*upy**(-1)
     &     + 4.D0*z**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 8.D0*upy**(-1)
     &     + 24.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 2.D0*y*z**(-1)
     &     - 8.D0*y*upz**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1)*Nc**2 * (
     &     + 16.D0
     &     + 4.D0*z**(-2)*upy**(-1)
     &     + 4.D0*z**(-2)*ypzp2**(-1)
     &     + 6.D0*z**(-2)
     &     - 12.D0*z**(-1)*upy**(-1)
     &     - 4.D0*z**(-1)*ypzp2**(-1)
     &     - 8.D0*z**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 88.D0*upz**(-2)
     &     - 48.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 44.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1)*Nc**2 * (
     &     + 20.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 6.D0*y*z**(-1)
     &     - 8.D0*y*upz**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1)*Nc**4 * (
     &     + 6.D0
     &     - 2.D0*z**(-2)*ypzp2**(-1)
     &     - 2.D0*z**(-2)
     &     + 6.D0*z**(-1)*ypzp2**(-1)
     &     + 6.D0*z**(-1)
     &     - 32.D0*upz**(-2)
     &     + 8.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 30.D0*ypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 22.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 2.D0*y*z**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1) * (
     &     - 18.D0
     &     - 8.D0*z**(-2)*upy**(-1)
     &     - 2.D0*z**(-2)*ypzp2**(-1)
     &     - 6.D0*z**(-2)
     &     + 8.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)*ypzp2**(-1)
     &     - 2.D0*z**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)
     &     - 80.D0*upz**(-2)
     &     + 80.D0*upz**(-1)
     &     - 14.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     + 2.D0*z*ypzp2**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr1(-1) * (
     &     - 6.D0*y*z**(-1)
     &     + 16.D0*y*upz**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(-1,-1)*Nc**2 * (
     &     - 16.D0
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     + 80.D0*upz**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-1)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(-1,-1)*Nc**4 * (
     &     + 16.D0
     &     - 40.D0*upz**(-1)
     &     - 32.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(-1,-1) * (
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     - 40.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(0,-1)*Nc**(-2) * (
     &     - 8.D0
     &     + 24.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)
     &     - 8.D0*upz**(-3)
     &     - 16.D0*upz**(-2)
     &     + 4.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-3)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(0,-1)*Nc**2 * (
     &     + 12.D0
     &     - 8.D0*upz**(-3)
     &     - 8.D0*upz**(-2)
     &     + 32.D0*upz**(-1)*ypzp2**(-1)
     &     - 36.D0*upz**(-1)
     &     - 20.D0*ypzp2**(-1)
     &     - 20.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-3)
     &     - 12.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(0,-1)*Nc**4 * (
     &     - 8.D0
     &     + 20.D0*upz**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + HBr2(0,-1) * (
     &     + 4.D0
     &     - 24.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*upy**(-1)
     &     + 16.D0*upz**(-3)
     &     + 24.D0*upz**(-2)
     &     - 32.D0*upz**(-1)*ypzp2**(-1)
     &     + 12.D0*upz**(-1)
     &     + 20.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 4.D0*z*ypzp2**(-1)
     &     + 16.D0*y*upz**(-3)
     &     + 12.D0*y*upz**(-1)
     &     )
      PoleEP0gg = PoleEP0gg + 106.D0/3.D0
     &     + 8.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)*zp2**(-1)
     &     + 6.D0*y**(-1)
     &     + 8.D0*z**(-1)*upy**(-1)
     &     + 2.D0*z**(-1)*ypzp2**(-1)
     &     + 6.D0*z**(-1)
     &     + 80.D0*upy**(-2)
     &     + 128.D0*upy**(-1)*upz**(-1)
     &     - 266.D0/3.D0*upy**(-1)
     &     + 80.D0*upz**(-2)
     &     - 266.D0/3.D0*upz**(-1)
     &     - 2.D0*zp2**(-1)*ypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 52.D0*ypzp2**(-1)
     &
      PoleEP0gg = PoleEP0gg + 16.D0*z2*upy**(-3)
     &     + 48.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 104.D0*z2*upy**(-2)
     &     - 120.D0*z2*upy**(-1)*upz**(-1)
     &     - 118.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 74.D0*z2*upy**(-1)
     &     + 16.D0*z2*upz**(-3)
     &     + 48.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 104.D0*z2*upz**(-2)
     &     - 118.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 74.D0*z2*upz**(-1)
     &     + 128.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*ypzp2**(-1)
     &     + 80.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &
      PoleEP0gg = PoleEP0gg - 6.D0*z2
     &     - 16.D0*z*upy**(-2)
     &     + 82.D0/3.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z*z2*upy**(-3)
     &     - 48.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 56.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 22.D0*z*z2*upy**(-1)
     &     + 256.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 256.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     + 14.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &
      PoleEP0gg = PoleEP0gg - 64.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*
     & sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)
     &     + 82.D0/3.D0*y*upz**(-1)
     &     + 16.D0*y*z2*upz**(-3)
     &     - 48.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 56.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 22.D0*y*z2*upz**(-1)
     &     - 24.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
      return
      end

*
*******************************************************************************
************************************************************************
************************************************************************
*
*          FUNCTIONS for the Cross Section
*
*************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function T34PoleEPm2qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T34PoleEPm2qq =  + Nc**(-2) * (
     &     - 8.D0
     &     - 16.D0*upx**(-4)
     &     - 32.D0*upx**(-4)*y
     &     - 16.D0*upx**(-4)*y**2
     &     + 32.D0*upx**(-3)
     &     + 64.D0*upx**(-3)*y
     &     + 32.D0*upx**(-3)*y**2
     &     - 16.D0*upx**(-2)
     &     - 48.D0*upx**(-2)*y
     &     - 16.D0*upx**(-2)*y**2
     &     + 16.D0*upx**(-1)*y
     &     )/4.D0
      T34PoleEPm2qq = T34PoleEPm2qq + (
     &     + 16.D0
     &     + 32.D0*upx**(-4)
     &     + 64.D0*upx**(-4)*y
     &     + 32.D0*upx**(-4)*y**2
     &     - 64.D0*upx**(-3)
     &     - 128.D0*upx**(-3)*y
     &     - 64.D0*upx**(-3)*y**2
     &     + 32.D0*upx**(-2)
     &     + 96.D0*upx**(-2)*y
     &     + 32.D0*upx**(-2)*y**2
     &     - 32.D0*upx**(-1)*y
     &     )/4.D0
      T34PoleEPm2qq = T34PoleEPm2qq + Nc**2 * (
     &     - 8.D0
     &     - 16.D0*upx**(-4)
     &     - 32.D0*upx**(-4)*y
     &     - 16.D0*upx**(-4)*y**2
     &     + 32.D0*upx**(-3)
     &     + 64.D0*upx**(-3)*y
     &     + 32.D0*upx**(-3)*y**2
     &     - 16.D0*upx**(-2)
     &     - 48.D0*upx**(-2)*y
     &     - 16.D0*upx**(-2)*y**2
     &     + 16.D0*upx**(-1)*y
     &     )/4.D0
      return
      end
*
*
*
      double precision function T34PoleEPm1qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T34PoleEPm1qq =  + Nc**(-2) * (
     &     + 6.D0
     &     + 20.D0*upx**(-4)
     &     + 40.D0*upx**(-4)*y
     &     + 20.D0*upx**(-4)*y**2
     &     - 40.D0*upx**(-3)
     &     - 80.D0*upx**(-3)*y
     &     - 40.D0*upx**(-3)*y**2
     &     + 20.D0*upx**(-2)
     &     + 60.D0*upx**(-2)*y
     &     + 20.D0*upx**(-2)*y**2
     &     - 20.D0*upx**(-1)*y
     &     - 8.D0*Hr1(-1)
     &     - 16.D0*Hr1(-1)*upx**(-4)
     &     - 32.D0*Hr1(-1)*upx**(-4)*y
     &     - 16.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 32.D0*Hr1(-1)*upx**(-3)
     &     + 64.D0*Hr1(-1)*upx**(-3)*y
     &     + 32.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 16.D0*Hr1(-1)*upx**(-2)
     &     - 48.D0*Hr1(-1)*upx**(-2)*y
     &     - 16.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 16.D0*Hr1(-1)*upx**(-1)*y
     &     - 16.D0*HAr1(-1)
     &     - 32.D0*HAr1(-1)*upx**(-4)
     &     - 64.D0*HAr1(-1)*upx**(-4)*y
     &     - 32.D0*HAr1(-1)*upx**(-4)*y**2
     &     + 64.D0*HAr1(-1)*upx**(-3)
     &     + 128.D0*HAr1(-1)*upx**(-3)*y
     &     + 64.D0*HAr1(-1)*upx**(-3)*y**2
     &     - 32.D0*HAr1(-1)*upx**(-2)
     &     - 96.D0*HAr1(-1)*upx**(-2)*y
     &     - 32.D0*HAr1(-1)*upx**(-2)*y**2
     &     + 32.D0*HAr1(-1)*upx**(-1)*y
     &     + 16.D0*HBr1(-1)
     &     + 32.D0*HBr1(-1)*upx**(-4)
     &     + 64.D0*HBr1(-1)*upx**(-4)*y
     &     + 32.D0*HBr1(-1)*upx**(-4)*y**2
     &     - 64.D0*HBr1(-1)*upx**(-3)
     &     - 128.D0*HBr1(-1)*upx**(-3)*y
     &     - 64.D0*HBr1(-1)*upx**(-3)*y**2
     &     + 32.D0*HBr1(-1)*upx**(-2)
     &     + 96.D0*HBr1(-1)*upx**(-2)*y
     &     + 32.D0*HBr1(-1)*upx**(-2)*y**2
     &     - 32.D0*HBr1(-1)*upx**(-1)*y
     &     + 9/2.D0*Hr1(0)*umx**(-1)
     &     - Hr1(0)*umx**(-1)*y
     &     + 1/2.D0*Hr1(0)*umx**(-1)*y**2
     &     + 8.D0*Hr1(0)*upx**(-5)
     &     + 16.D0*Hr1(0)*upx**(-5)*y
     &     + 8.D0*Hr1(0)*upx**(-5)*y**2
     &     - 12.D0*Hr1(0)*upx**(-4)
     &     - 24.D0*Hr1(0)*upx**(-4)*y
     &     - 12.D0*Hr1(0)*upx**(-4)*y**2
     &     + 2.D0*Hr1(0)*upx**(-3)
     &     + 12.D0*Hr1(0)*upx**(-3)*y
     &     + 2.D0*Hr1(0)*upx**(-3)*y**2
     &     + Hr1(0)*upx**(-2)
     &     - 2.D0*Hr1(0)*upx**(-2)*y
     &     + Hr1(0)*upx**(-2)*y**2
     &     + 9/2.D0*Hr1(0)*upx**(-1)
     &     - Hr1(0)*upx**(-1)*y
     &     + 1/2.D0*Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T34PoleEPm1qq = T34PoleEPm1qq + (
     &     - 12.D0
     &     - 40.D0*upx**(-4)
     &     - 80.D0*upx**(-4)*y
     &     - 40.D0*upx**(-4)*y**2
     &     + 80.D0*upx**(-3)
     &     + 160.D0*upx**(-3)*y
     &     + 80.D0*upx**(-3)*y**2
     &     - 40.D0*upx**(-2)
     &     - 120.D0*upx**(-2)*y
     &     - 40.D0*upx**(-2)*y**2
     &     + 40.D0*upx**(-1)*y
     &     + 8.D0*Hr1(-1)
     &     + 16.D0*Hr1(-1)*upx**(-4)
     &     + 32.D0*Hr1(-1)*upx**(-4)*y
     &     + 16.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 32.D0*Hr1(-1)*upx**(-3)
     &     - 64.D0*Hr1(-1)*upx**(-3)*y
     &     - 32.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 16.D0*Hr1(-1)*upx**(-2)
     &     + 48.D0*Hr1(-1)*upx**(-2)*y
     &     + 16.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 16.D0*Hr1(-1)*upx**(-1)*y
     &     + 24.D0*HAr1(-1)
     &     + 48.D0*HAr1(-1)*upx**(-4)
     &     + 96.D0*HAr1(-1)*upx**(-4)*y
     &     + 48.D0*HAr1(-1)*upx**(-4)*y**2
     &     - 96.D0*HAr1(-1)*upx**(-3)
     &     - 192.D0*HAr1(-1)*upx**(-3)*y
     &     - 96.D0*HAr1(-1)*upx**(-3)*y**2
     &     + 48.D0*HAr1(-1)*upx**(-2)
     &     + 144.D0*HAr1(-1)*upx**(-2)*y
     &     + 48.D0*HAr1(-1)*upx**(-2)*y**2
     &     - 48.D0*HAr1(-1)*upx**(-1)*y
     &     - 16.D0*HBr1(-1)
     &     - 32.D0*HBr1(-1)*upx**(-4)
     &     - 64.D0*HBr1(-1)*upx**(-4)*y
     &     - 32.D0*HBr1(-1)*upx**(-4)*y**2
     &     + 64.D0*HBr1(-1)*upx**(-3)
     &     + 128.D0*HBr1(-1)*upx**(-3)*y
     &     + 64.D0*HBr1(-1)*upx**(-3)*y**2
     &     - 32.D0*HBr1(-1)*upx**(-2)
     &     - 96.D0*HBr1(-1)*upx**(-2)*y
     &     - 32.D0*HBr1(-1)*upx**(-2)*y**2
     &     + 32.D0*HBr1(-1)*upx**(-1)*y
     &     - 9/2.D0*Hr1(0)*umx**(-1)
     &     + Hr1(0)*umx**(-1)*y
     &     - 1/2.D0*Hr1(0)*umx**(-1)*y**2
     &     - 8.D0*Hr1(0)*upx**(-5)
     &     - 16.D0*Hr1(0)*upx**(-5)*y
     &     - 8.D0*Hr1(0)*upx**(-5)*y**2
     &     + 12.D0*Hr1(0)*upx**(-4)
     &     + 24.D0*Hr1(0)*upx**(-4)*y
     &     + 12.D0*Hr1(0)*upx**(-4)*y**2
     &     - 2.D0*Hr1(0)*upx**(-3)
     &     - 12.D0*Hr1(0)*upx**(-3)*y
     &     - 2.D0*Hr1(0)*upx**(-3)*y**2
     &     - Hr1(0)*upx**(-2)
     &     + 2.D0*Hr1(0)*upx**(-2)*y
     &     - Hr1(0)*upx**(-2)*y**2
     &     - 9/2.D0*Hr1(0)*upx**(-1)
     &     + Hr1(0)*upx**(-1)*y
     &     - 1/2.D0*Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T34PoleEPm1qq = T34PoleEPm1qq + Nc**2 * (
     &     + 6.D0
     &     + 20.D0*upx**(-4)
     &     + 40.D0*upx**(-4)*y
     &     + 20.D0*upx**(-4)*y**2
     &     - 40.D0*upx**(-3)
     &     - 80.D0*upx**(-3)*y
     &     - 40.D0*upx**(-3)*y**2
     &     + 20.D0*upx**(-2)
     &     + 60.D0*upx**(-2)*y
     &     + 20.D0*upx**(-2)*y**2
     &     - 20.D0*upx**(-1)*y
     &     - 8.D0*HAr1(-1)
     &     - 16.D0*HAr1(-1)*upx**(-4)
     &     - 32.D0*HAr1(-1)*upx**(-4)*y
     &     - 16.D0*HAr1(-1)*upx**(-4)*y**2
     &     + 32.D0*HAr1(-1)*upx**(-3)
     &     + 64.D0*HAr1(-1)*upx**(-3)*y
     &     + 32.D0*HAr1(-1)*upx**(-3)*y**2
     &     - 16.D0*HAr1(-1)*upx**(-2)
     &     - 48.D0*HAr1(-1)*upx**(-2)*y
     &     - 16.D0*HAr1(-1)*upx**(-2)*y**2
     &     + 16.D0*HAr1(-1)*upx**(-1)*y
     &     )/(-2.D0)
      
      return
      end
*
*
*
      double precision function T34PoleEP0qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1d0

* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T34PoleEP0qq =  + Nc**(-2) * (
     &     - 7.D0
     &     + 8.D0*umx**(-3)*z2
     &     - 8.D0*umx**(-3)*y*z2
     &     - 12.D0*umx**(-2)*z2
     &     + 12.D0*umx**(-2)*y*z2
     &     + 13.D0*umx**(-1)*z2
     &     - 2.D0*umx**(-1)*y*z2
     &     + umx**(-1)*y**2*z2
     &     + 16.D0*upx**(-5)*z2
     &     + 32.D0*upx**(-5)*y*z2
     &     + 16.D0*upx**(-5)*y**2*z2
     &     - 24.D0*upx**(-4)
     &     - 24.D0*upx**(-4)*z2
     &     - 48.D0*upx**(-4)*y
     &     - 48.D0*upx**(-4)*y*z2
     &     - 24.D0*upx**(-4)*y**2
     &     - 24.D0*upx**(-4)*y**2*z2
     &     + 48.D0*upx**(-3)
     &     - 4.D0*upx**(-3)*z2
     &     + 96.D0*upx**(-3)*y
     &     + 16.D0*upx**(-3)*y*z2
     &     + 48.D0*upx**(-3)*y**2
     &     + 4.D0*upx**(-3)*y**2*z2
     &     - 24.D0*upx**(-2)
     &     + 6.D0*upx**(-2)*z2
     &     - 72.D0*upx**(-2)*y
     &     - 24.D0*upx**(-2)*y**2
     &     + 2.D0*upx**(-2)*y**2*z2
     &     + 5.D0*upx**(-1)*z2
     &     + 24.D0*upx**(-1)*y
     &     - 2.D0*upx**(-1)*y*z2
     &     + upx**(-1)*y**2*z2
     &     - 6.D0*Hr1(-1)
     &     - 16.D0*Hr1(-1)*umx**(-2)
     &     + 16.D0*Hr1(-1)*umx**(-2)*y
     &     + 16.D0*Hr1(-1)*umx**(-1)
     &     - 16.D0*Hr1(-1)*umx**(-1)*y
     &     + 12.D0*Hr1(-1)*upx**(-4)
     &     + 24.D0*Hr1(-1)*upx**(-4)*y
     &     + 12.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 24.D0*Hr1(-1)*upx**(-3)
     &     - 48.D0*Hr1(-1)*upx**(-3)*y
     &     - 24.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 12.D0*Hr1(-1)*upx**(-2)
     &     + 36.D0*Hr1(-1)*upx**(-2)*y
     &     + 12.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 12.D0*Hr1(-1)*upx**(-1)*y
     &     - 8.D0*Hr1(-1)*HAr1(-1)
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 64.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 128.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 64.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 80.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 16.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 24.D0*Hr1(-1)*HBr1(-1)
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 64.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 128.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 64.D0*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 112.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 48.D0*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 4.D0*HAr1(-1)
     &     + 4.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 8.D0*HAr1(-1)*upx**(-2)
     &     + 4.D0*HAr1(-1)*upx**(-2)*y
     &     - 4.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 8.D0*HAr1(-1)*upx**(-1)
     &     - 4.D0*HAr1(-1)*upx**(-1)*y
     &     - 4.D0*HAr1(-1)*y**(-1)
     &     + 4.D0*HAr1(-1)*Hr1(0)
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 64.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 32.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 40.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 8.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 8.D0*HBr1(-1)
     &     - 4.D0*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 4.D0*HBr1(-1)*upx**(-2)*y
     &     + 4.D0*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 4.D0*HBr1(-1)*upx**(-1)*y
     &     + 4.D0*HBr1(-1)*z**(-1)
     &     - 12.D0*HBr1(-1)*Hr1(0)
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 32.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 64.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 32.D0*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 56.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 24.D0*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 24.D0*Hr2(-1,-1)
     &     - 16.D0*Hr2(-1,-1)*upx**(-4)
     &     - 32.D0*Hr2(-1,-1)*upx**(-4)*y
     &     - 16.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 32.D0*Hr2(-1,-1)*upx**(-3)
     &     + 64.D0*Hr2(-1,-1)*upx**(-3)*y
     &     + 32.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 48.D0*Hr2(-1,-1)*upx**(-2)
     &     - 80.D0*Hr2(-1,-1)*upx**(-2)*y
     &     - 16.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 32.D0*Hr2(-1,-1)*upx**(-1)
     &     + 48.D0*Hr2(-1,-1)*upx**(-1)*y
     &     - 8.D0*HAr2(-1,-1)
     &     - 16.D0*HAr2(-1,-1)*upx**(-2)*y
     &     + 16.D0*HAr2(-1,-1)*upx**(-1)*y
     &     - 8.D0*HBr2(-1,-1)
     &     - 32.D0*HBr2(-1,-1)*upx**(-2)
     &     - 16.D0*HBr2(-1,-1)*upx**(-2)*y
     &     + 32.D0*HBr2(-1,-1)*upx**(-1)
     &     + 16.D0*HBr2(-1,-1)*upx**(-1)*y
     &     + 12.D0*Hr2(-1,0)
     &     + 8.D0*Hr2(-1,0)*upx**(-4)
     &     + 16.D0*Hr2(-1,0)*upx**(-4)*y
     &     + 8.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     - 16.D0*Hr2(-1,0)*upx**(-3)
     &     - 32.D0*Hr2(-1,0)*upx**(-3)*y
     &     - 16.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     + 24.D0*Hr2(-1,0)*upx**(-2)
     &     + 40.D0*Hr2(-1,0)*upx**(-2)*y
     &     + 8.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     - 16.D0*Hr2(-1,0)*upx**(-1)
     &     - 24.D0*Hr2(-1,0)*upx**(-1)*y
     &     + 4.D0*Hr1(0)
     &     + 8.D0*Hr1(0)*umx**(-2)
     &     - 8.D0*Hr1(0)*umx**(-2)*y
     &     - 25/4.D0*Hr1(0)*umx**(-1)
     &     + 17/2.D0*Hr1(0)*umx**(-1)*y
     &     - 1/4.D0*Hr1(0)*umx**(-1)*y**2
     &     - 8.D0*Hr1(0)*upx**(-5)
     &     - 16.D0*Hr1(0)*upx**(-5)*y
     &     - 8.D0*Hr1(0)*upx**(-5)*y**2
     &     + 14.D0*Hr1(0)*upx**(-4)
     &     + 28.D0*Hr1(0)*upx**(-4)*y
     &     + 14.D0*Hr1(0)*upx**(-4)*y**2
     &     - 5.D0*Hr1(0)*upx**(-3)
     &     - 18.D0*Hr1(0)*upx**(-3)*y
     &     - 5.D0*Hr1(0)*upx**(-3)*y**2
     &     - 1/2.D0*Hr1(0)*upx**(-2)
     &     + 5.D0*Hr1(0)*upx**(-2)*y
     &     - 1/2.D0*Hr1(0)*upx**(-2)*y**2
     &     - 17/4.D0*Hr1(0)*upx**(-1)
     &     + 1/2.D0*Hr1(0)*upx**(-1)*y
     &     - 1/4.D0*Hr1(0)*upx**(-1)*y**2
     &     + 4.D0*Hr2(0,-1)
     &     - 16.D0*Hr2(0,-1)*umx**(-3)
     &     + 16.D0*Hr2(0,-1)*umx**(-3)*y
     &     + 24.D0*Hr2(0,-1)*umx**(-2)
     &     - 24.D0*Hr2(0,-1)*umx**(-2)*y
     &     - 8.D0*Hr2(0,-1)*umx**(-1)
     &     + 8.D0*Hr2(0,-1)*upx**(-4)
     &     + 16.D0*Hr2(0,-1)*upx**(-4)*y
     &     + 8.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     - 16.D0*Hr2(0,-1)*upx**(-3)*y
     &     - 16.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     + 16.D0*Hr2(0,-1)*upx**(-2)*y
     &     + 8.D0*Hr2(0,-1)*upx**(-2)*y**2
     &     + 8.D0*Hr2(0,-1)*upx**(-1)
     &     - 8.D0*Hr2(0,-1)*upx**(-1)*y
     &     + 4.D0*HAr2(0,-1)
     &     + 8.D0*HAr2(0,-1)*upx**(-2)*y
     &     - 8.D0*HAr2(0,-1)*upx**(-1)*y
     &     + 4.D0*HBr2(0,-1)
     &     + 16.D0*HBr2(0,-1)*upx**(-2)
     &     + 8.D0*HBr2(0,-1)*upx**(-2)*y
     &     - 16.D0*HBr2(0,-1)*upx**(-1)
     &     - 8.D0*HBr2(0,-1)*upx**(-1)*y
     &     + 8.D0*Hr2(0,0)*umx**(-3)
     &     - 8.D0*Hr2(0,0)*umx**(-3)*y
     &     - 12.D0*Hr2(0,0)*umx**(-2)
     &     + 12.D0*Hr2(0,0)*umx**(-2)*y
     &     + 7/4.D0*Hr2(0,0)*umx**(-1)
     &     + 1/2.D0*Hr2(0,0)*umx**(-1)*y
     &     - 1/4.D0*Hr2(0,0)*umx**(-1)*y**2
     &     - 4.D0*Hr2(0,0)*upx**(-5)
     &     - 8.D0*Hr2(0,0)*upx**(-5)*y
     &     - 4.D0*Hr2(0,0)*upx**(-5)*y**2
     &     + 6.D0*Hr2(0,0)*upx**(-4)
     &     + 12.D0*Hr2(0,0)*upx**(-4)*y
     &     + 6.D0*Hr2(0,0)*upx**(-4)*y**2
     &     - 9.D0*Hr2(0,0)*upx**(-3)
     &     - 14.D0*Hr2(0,0)*upx**(-3)*y
     &     - Hr2(0,0)*upx**(-3)*y**2
     &     + 7/2.D0*Hr2(0,0)*upx**(-2)
     &     + 5.D0*Hr2(0,0)*upx**(-2)*y
     &     - 1/2.D0*Hr2(0,0)*upx**(-2)*y**2
     &     - 25/4.D0*Hr2(0,0)*upx**(-1)
     &     + 1/2.D0*Hr2(0,0)*upx**(-1)*y
     &     - 1/4.D0*Hr2(0,0)*upx**(-1)*y**2
     &     + 4.D0*Hr2(1,0)
     &     - 9/2.D0*Hr2(1,0)*umx**(-1)
     &     + Hr2(1,0)*umx**(-1)*y
     &     - 1/2.D0*Hr2(1,0)*umx**(-1)*y**2
     &     - 8.D0*Hr2(1,0)*upx**(-5)
     &     - 16.D0*Hr2(1,0)*upx**(-5)*y
     &     - 8.D0*Hr2(1,0)*upx**(-5)*y**2
     &     + 20.D0*Hr2(1,0)*upx**(-4)
     &     + 40.D0*Hr2(1,0)*upx**(-4)*y
     &     + 20.D0*Hr2(1,0)*upx**(-4)*y**2
     &     - 18.D0*Hr2(1,0)*upx**(-3)
     &     - 44.D0*Hr2(1,0)*upx**(-3)*y
     &     - 18.D0*Hr2(1,0)*upx**(-3)*y**2
     &     + 7.D0*Hr2(1,0)*upx**(-2)
     &     + 26.D0*Hr2(1,0)*upx**(-2)*y
     &     + 7.D0*Hr2(1,0)*upx**(-2)*y**2
     &     - 9/2.D0*Hr2(1,0)*upx**(-1)
     &     - 7.D0*Hr2(1,0)*upx**(-1)*y
     &     - 1/2.D0*Hr2(1,0)*upx**(-1)*y**2
     &     )
      T34PoleEP0qq = T34PoleEP0qq + Nc**(-1)*Nh * (
     &     + 20.D0/9.D0
     &     - 32/3.D0*upx**(-6)
     &     - 64/3.D0*upx**(-6)*y
     &     - 32/3.D0*upx**(-6)*y**2
     &     + 32.D0*upx**(-5)
     &     + 64.D0*upx**(-5)*y
     &     + 32.D0*upx**(-5)*y**2
     &     - 248/9.D0*upx**(-4)
     &     - 592/9.D0*upx**(-4)*y
     &     - 248/9.D0*upx**(-4)*y**2
     &     + 16/9.D0*upx**(-3)
     &     + 224/9.D0*upx**(-3)*y
     &     + 16/9.D0*upx**(-3)*y**2
     &     - 8/9.D0*upx**(-2)
     &     + 8/3.D0*upx**(-2)*y
     &     + 40/9.D0*upx**(-2)*y**2
     &     + 16/3.D0*upx**(-1)
     &     - 40/9.D0*upx**(-1)*y
     &     - 4/3.D0*Hr1(0)
     &     - 32/3.D0*Hr1(0)*upx**(-7)
     &     - 64/3.D0*Hr1(0)*upx**(-7)*y
     &     - 32/3.D0*Hr1(0)*upx**(-7)*y**2
     &     + 112/3.D0*Hr1(0)*upx**(-6)
     &     + 224/3.D0*Hr1(0)*upx**(-6)*y
     &     + 112/3.D0*Hr1(0)*upx**(-6)*y**2
     &     - 128/3.D0*Hr1(0)*upx**(-5)
     &     - 96.D0*Hr1(0)*upx**(-5)*y
     &     - 128/3.D0*Hr1(0)*upx**(-5)*y**2
     &     + 40/3.D0*Hr1(0)*upx**(-4)
     &     + 160/3.D0*Hr1(0)*upx**(-4)*y
     &     + 40/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-3)*y
     &     + 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-2)
     &     - 8.D0*Hr1(0)*upx**(-2)*y
     &     - 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-1)*y
     &     )
      T34PoleEP0qq = T34PoleEP0qq + Nc**(-1)*Nl * (
     &     + 20.D0/9.D0
     &     + 40/9.D0*upx**(-4)
     &     + 80/9.D0*upx**(-4)*y
     &     + 40/9.D0*upx**(-4)*y**2
     &     - 80/9.D0*upx**(-3)
     &     - 160/9.D0*upx**(-3)*y
     &     - 80/9.D0*upx**(-3)*y**2
     &     + 40/9.D0*upx**(-2)
     &     + 40/3.D0*upx**(-2)*y
     &     + 40/9.D0*upx**(-2)*y**2
     &     - 40/9.D0*upx**(-1)*y
     &     - 8/3.D0*Hr1(-1)
     &     - 16/3.D0*Hr1(-1)*upx**(-4)
     &     - 32/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 16/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 32/3.D0*Hr1(-1)*upx**(-3)
     &     + 64/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 32/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 16/3.D0*Hr1(-1)*upx**(-2)
     &     - 16.D0*Hr1(-1)*upx**(-2)*y
     &     - 16/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 16/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 4/3.D0*Hr1(0)
     &     + 8/3.D0*Hr1(0)*upx**(-4)
     &     + 16/3.D0*Hr1(0)*upx**(-4)*y
     &     + 8/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-3)
     &     - 32/3.D0*Hr1(0)*upx**(-3)*y
     &     - 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-2)
     &     + 8.D0*Hr1(0)*upx**(-2)*y
     &     + 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 8/3.D0*Hr1(0)*upx**(-1)*y
     &     + 4/3.D0*DLog(tm**(-2)*w**2)
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 8.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T34PoleEP0qq = T34PoleEP0qq + (
     &     - 44.D0/9.D0
     &     - 6.D0*umx**(-5)*z2
     &     + 12.D0*umx**(-5)*y*z2
     &     - 6.D0*umx**(-5)*y**2*z2
     &     + 15.D0*umx**(-4)*z2
     &     - 30.D0*umx**(-4)*y*z2
     &     + 15.D0*umx**(-4)*y**2*z2
     &     - 47/2.D0*umx**(-3)*z2
     &     + 35.D0*umx**(-3)*y*z2
     &     - 19/2.D0*umx**(-3)*y**2*z2
     &     - umx**(-2)
     &     + 81/4.D0*umx**(-2)*z2
     &     + 2.D0*umx**(-2)*y
     &     - 45/2.D0*umx**(-2)*y*z2
     &     - umx**(-2)*y**2
     &     - 3/4.D0*umx**(-2)*y**2*z2
     &     + umx**(-1)
     &     - 145/8.D0*umx**(-1)*z2
     &     - 2.D0*umx**(-1)*y
     &     + 11/4.D0*umx**(-1)*y*z2
     &     + umx**(-1)*y**2
     &     - 9/8.D0*umx**(-1)*y**2*z2
     &     - 16.D0*upx**(-5)*z2
     &     - 32.D0*upx**(-5)*y*z2
     &     - 16.D0*upx**(-5)*y**2*z2
     &     + 128/9.D0*upx**(-4)
     &     + 20.D0*upx**(-4)*z2
     &     + 256/9.D0*upx**(-4)*y
     &     + 40.D0*upx**(-4)*y*z2
     &     + 128/9.D0*upx**(-4)*y**2
     &     + 20.D0*upx**(-4)*y**2*z2
     &     - 256/9.D0*upx**(-3)
     &     + 11.D0*upx**(-3)*z2
     &     - 512/9.D0*upx**(-3)*y
     &     - 4.D0*upx**(-3)*y*z2
     &     - 256/9.D0*upx**(-3)*y**2
     &     + upx**(-3)*y**2*z2
     &     + 119/9.D0*upx**(-2)
     &     - 17/2.D0*upx**(-2)*z2
     &     + 122/3.D0*upx**(-2)*y
     &     - 4.D0*upx**(-2)*y*z2
     &     + 119/9.D0*upx**(-2)*y**2
     &     - 3/2.D0*upx**(-2)*y**2*z2
     &     + upx**(-1)
     &     - 33/8.D0*upx**(-1)*z2
     &     - 110/9.D0*upx**(-1)*y
     &     + 11/4.D0*upx**(-1)*y*z2
     &     + upx**(-1)*y**2
     &     - 9/8.D0*upx**(-1)*y**2*z2
     &     + 38/3.D0*Hr1(-1)
     &     + 12.D0*Hr1(-1)*umx**(-4)
     &     - 24.D0*Hr1(-1)*umx**(-4)*y
     &     + 12.D0*Hr1(-1)*umx**(-4)*y**2
     &     - 24.D0*Hr1(-1)*umx**(-3)
     &     + 48.D0*Hr1(-1)*umx**(-3)*y
     &     - 24.D0*Hr1(-1)*umx**(-3)*y**2
     &     + 36.D0*Hr1(-1)*umx**(-2)
     &     - 48.D0*Hr1(-1)*umx**(-2)*y
     &     + 8.D0*Hr1(-1)*umx**(-2)*y**2
     &     - 24.D0*Hr1(-1)*umx**(-1)
     &     + 24.D0*Hr1(-1)*umx**(-1)*y
     &     + 4.D0*Hr1(-1)*umx**(-1)*y**2
     &     - 8/3.D0*Hr1(-1)*upx**(-4)
     &     - 16/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 8/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 16/3.D0*Hr1(-1)*upx**(-3)
     &     + 32/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 16/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 20/3.D0*Hr1(-1)*upx**(-2)
     &     - 16.D0*Hr1(-1)*upx**(-2)*y
     &     - 20/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 4.D0*Hr1(-1)*upx**(-1)
     &     + 32/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 4.D0*Hr1(-1)*upx**(-1)*y**2
     &     + 12.D0*Hr1(-1)*HAr1(-1)
     &     + 48.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 48.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 96.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 192.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 96.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 48.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 120.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 48.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 24.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     - 24.D0*Hr1(-1)*HBr1(-1)
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     + 64.D0*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     + 128.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     + 64.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     - 112.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     + 48.D0*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     + 6.D0*HAr1(-1)
     &     - 6.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 12.D0*HAr1(-1)*upx**(-2)
     &     - 6.D0*HAr1(-1)*upx**(-2)*y
     &     + 6.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 12.D0*HAr1(-1)*upx**(-1)
     &     + 6.D0*HAr1(-1)*upx**(-1)*y
     &     + 6.D0*HAr1(-1)*y**(-1)
     &     - 6.D0*HAr1(-1)*Hr1(0)
     &     - 24.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 24.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 48.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 96.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 48.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 24.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 60.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 24.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 12.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 8.D0*HBr1(-1)
     &     + 4.D0*HBr1(-1)*upx**(-2)*z**(-1)
     &     - 4.D0*HBr1(-1)*upx**(-2)*y
     &     - 4.D0*HBr1(-1)*upx**(-1)*z**(-1)
     &     + 4.D0*HBr1(-1)*upx**(-1)*y
     &     - 4.D0*HBr1(-1)*z**(-1)
     &     + 12.D0*HBr1(-1)*Hr1(0)
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-4)
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0*HBr1(-1)*Hr1(0)*upx**(-3)
     &     - 64.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 32.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-2)
     &     + 56.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-1)
     &     - 24.D0*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 20.D0*Hr2(-1,-1)
     &     + 32.D0*Hr2(-1,-1)*upx**(-2)
     &     + 40.D0*Hr2(-1,-1)*upx**(-2)*y
     &     - 32.D0*Hr2(-1,-1)*upx**(-1)
     &     - 40.D0*Hr2(-1,-1)*upx**(-1)*y
     &     + 12.D0*HAr2(-1,-1)
     &     + 24.D0*HAr2(-1,-1)*upx**(-2)*y
     &     - 24.D0*HAr2(-1,-1)*upx**(-1)*y
     &     + 8.D0*HBr2(-1,-1)
     &     + 32.D0*HBr2(-1,-1)*upx**(-2)
     &     + 16.D0*HBr2(-1,-1)*upx**(-2)*y
     &     - 32.D0*HBr2(-1,-1)*upx**(-1)
     &     - 16.D0*HBr2(-1,-1)*upx**(-1)*y
     &     - 10.D0*Hr2(-1,0)
     &     - 16.D0*Hr2(-1,0)*upx**(-2)
     &     - 20.D0*Hr2(-1,0)*upx**(-2)*y
     &     + 16.D0*Hr2(-1,0)*upx**(-1)
     &     + 20.D0*Hr2(-1,0)*upx**(-1)*y
     &     - 22/3.D0*Hr1(0)
     &     - 6.D0*Hr1(0)*umx**(-4)
     &     + 12.D0*Hr1(0)*umx**(-4)*y
     &     - 6.D0*Hr1(0)*umx**(-4)*y**2
     &     + 12.D0*Hr1(0)*umx**(-3)
     &     - 24.D0*Hr1(0)*umx**(-3)*y
     &     + 12.D0*Hr1(0)*umx**(-3)*y**2
     &     - 18.D0*Hr1(0)*umx**(-2)
     &     + 24.D0*Hr1(0)*umx**(-2)*y
     &     - 4.D0*Hr1(0)*umx**(-2)*y**2
     &     + 41/4.D0*Hr1(0)*umx**(-1)
     &     - 25/2.D0*Hr1(0)*umx**(-1)*y
     &     - 7/4.D0*Hr1(0)*umx**(-1)*y**2
     &     + 8.D0*Hr1(0)*upx**(-5)
     &     + 16.D0*Hr1(0)*upx**(-5)*y
     &     + 8.D0*Hr1(0)*upx**(-5)*y**2
     &     - 56/3.D0*Hr1(0)*upx**(-4)
     &     - 112/3.D0*Hr1(0)*upx**(-4)*y
     &     - 56/3.D0*Hr1(0)*upx**(-4)*y**2
     &     + 43/3.D0*Hr1(0)*upx**(-3)
     &     + 110/3.D0*Hr1(0)*upx**(-3)*y
     &     + 43/3.D0*Hr1(0)*upx**(-3)*y**2
     &     - 13/6.D0*Hr1(0)*upx**(-2)
     &     - 15.D0*Hr1(0)*upx**(-2)*y
     &     - 13/6.D0*Hr1(0)*upx**(-2)*y**2
     &     + 9/4.D0*Hr1(0)*upx**(-1)
     &     + 1/6.D0*Hr1(0)*upx**(-1)*y
     &     - 7/4.D0*Hr1(0)*upx**(-1)*y**2
     &     - 4.D0*Hr2(0,-1)
     &     + 12.D0*Hr2(0,-1)*umx**(-5)
     &     - 24.D0*Hr2(0,-1)*umx**(-5)*y
     &     + 12.D0*Hr2(0,-1)*umx**(-5)*y**2
     &     - 30.D0*Hr2(0,-1)*umx**(-4)
     &     + 60.D0*Hr2(0,-1)*umx**(-4)*y
     &     - 30.D0*Hr2(0,-1)*umx**(-4)*y**2
     &     + 47.D0*Hr2(0,-1)*umx**(-3)
     &     - 70.D0*Hr2(0,-1)*umx**(-3)*y
     &     + 19.D0*Hr2(0,-1)*umx**(-3)*y**2
     &     - 81/2.D0*Hr2(0,-1)*umx**(-2)
     &     + 45.D0*Hr2(0,-1)*umx**(-2)*y
     &     + 3/2.D0*Hr2(0,-1)*umx**(-2)*y**2
     &     + 73/4.D0*Hr2(0,-1)*umx**(-1)
     &     - 3/2.D0*Hr2(0,-1)*umx**(-1)*y
     &     + 1/4.D0*Hr2(0,-1)*umx**(-1)*y**2
     &     - 14.D0*Hr2(0,-1)*upx**(-3)
     &     - 8.D0*Hr2(0,-1)*upx**(-3)*y
     &     + 6.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     + 5.D0*Hr2(0,-1)*upx**(-2)
     &     - 8.D0*Hr2(0,-1)*upx**(-2)*y
     &     - 9.D0*Hr2(0,-1)*upx**(-2)*y**2
     &     - 39/4.D0*Hr2(0,-1)*upx**(-1)
     &     + 13/2.D0*Hr2(0,-1)*upx**(-1)*y
     &     + 1/4.D0*Hr2(0,-1)*upx**(-1)*y**2
     &     - 6.D0*HAr2(0,-1)
     &     - 12.D0*HAr2(0,-1)*upx**(-2)*y
     &     + 12.D0*HAr2(0,-1)*upx**(-1)*y
     &     - 4.D0*HBr2(0,-1)
     &     - 16.D0*HBr2(0,-1)*upx**(-2)
     &     - 8.D0*HBr2(0,-1)*upx**(-2)*y
     &     + 16.D0*HBr2(0,-1)*upx**(-1)
     &     + 8.D0*HBr2(0,-1)*upx**(-1)*y
     &     - 6.D0*Hr2(0,0)*umx**(-5)
     &     + 12.D0*Hr2(0,0)*umx**(-5)*y
     &     - 6.D0*Hr2(0,0)*umx**(-5)*y**2
     &     + 15.D0*Hr2(0,0)*umx**(-4)
     &     - 30.D0*Hr2(0,0)*umx**(-4)*y
     &     + 15.D0*Hr2(0,0)*umx**(-4)*y**2
     &     - 47/2.D0*Hr2(0,0)*umx**(-3)
     &     + 35.D0*Hr2(0,0)*umx**(-3)*y
     &     - 19/2.D0*Hr2(0,0)*umx**(-3)*y**2
     &     + 81/4.D0*Hr2(0,0)*umx**(-2)
     &     - 45/2.D0*Hr2(0,0)*umx**(-2)*y
     &     - 3/4.D0*Hr2(0,0)*umx**(-2)*y**2
     &     - 55/8.D0*Hr2(0,0)*umx**(-1)
     &     + 1/4.D0*Hr2(0,0)*umx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*umx**(-1)*y**2
     &     + 4.D0*Hr2(0,0)*upx**(-5)
     &     + 8.D0*Hr2(0,0)*upx**(-5)*y
     &     + 4.D0*Hr2(0,0)*upx**(-5)*y**2
     &     - 10.D0*Hr2(0,0)*upx**(-4)
     &     - 20.D0*Hr2(0,0)*upx**(-4)*y
     &     - 10.D0*Hr2(0,0)*upx**(-4)*y**2
     &     + 16.D0*Hr2(0,0)*upx**(-3)
     &     + 26.D0*Hr2(0,0)*upx**(-3)*y
     &     + 6.D0*Hr2(0,0)*upx**(-3)*y**2
     &     - 6.D0*Hr2(0,0)*upx**(-2)
     &     - 9.D0*Hr2(0,0)*upx**(-2)*y
     &     + Hr2(0,0)*upx**(-2)*y**2
     &     + 57/8.D0*Hr2(0,0)*upx**(-1)
     &     + 1/4.D0*Hr2(0,0)*upx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*upx**(-1)*y**2
     &     - 4.D0*Hr2(1,0)
     &     + 9/2.D0*Hr2(1,0)*umx**(-1)
     &     - Hr2(1,0)*umx**(-1)*y
     &     + 1/2.D0*Hr2(1,0)*umx**(-1)*y**2
     &     + 8.D0*Hr2(1,0)*upx**(-5)
     &     + 16.D0*Hr2(1,0)*upx**(-5)*y
     &     + 8.D0*Hr2(1,0)*upx**(-5)*y**2
     &     - 20.D0*Hr2(1,0)*upx**(-4)
     &     - 40.D0*Hr2(1,0)*upx**(-4)*y
     &     - 20.D0*Hr2(1,0)*upx**(-4)*y**2
     &     + 18.D0*Hr2(1,0)*upx**(-3)
     &     + 44.D0*Hr2(1,0)*upx**(-3)*y
     &     + 18.D0*Hr2(1,0)*upx**(-3)*y**2
     &     - 7.D0*Hr2(1,0)*upx**(-2)
     &     - 26.D0*Hr2(1,0)*upx**(-2)*y
     &     - 7.D0*Hr2(1,0)*upx**(-2)*y**2
     &     + 9/2.D0*Hr2(1,0)*upx**(-1)
     &     + 7.D0*Hr2(1,0)*upx**(-1)*y
     &     + 1/2.D0*Hr2(1,0)*upx**(-1)*y**2
     &     - 22/3.D0*DLog(tm**(-2)*w**2)
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 44.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T34PoleEP0qq = T34PoleEP0qq + Nc*Nh * (
     &     - 20.D0/9.D0
     &     + 32/3.D0*upx**(-6)
     &     + 64/3.D0*upx**(-6)*y
     &     + 32/3.D0*upx**(-6)*y**2
     &     - 32.D0*upx**(-5)
     &     - 64.D0*upx**(-5)*y
     &     - 32.D0*upx**(-5)*y**2
     &     + 248/9.D0*upx**(-4)
     &     + 592/9.D0*upx**(-4)*y
     &     + 248/9.D0*upx**(-4)*y**2
     &     - 16/9.D0*upx**(-3)
     &     - 224/9.D0*upx**(-3)*y
     &     - 16/9.D0*upx**(-3)*y**2
     &     + 8/9.D0*upx**(-2)
     &     - 8/3.D0*upx**(-2)*y
     &     - 40/9.D0*upx**(-2)*y**2
     &     - 16/3.D0*upx**(-1)
     &     + 40/9.D0*upx**(-1)*y
     &     + 4/3.D0*Hr1(0)
     &     + 32/3.D0*Hr1(0)*upx**(-7)
     &     + 64/3.D0*Hr1(0)*upx**(-7)*y
     &     + 32/3.D0*Hr1(0)*upx**(-7)*y**2
     &     - 112/3.D0*Hr1(0)*upx**(-6)
     &     - 224/3.D0*Hr1(0)*upx**(-6)*y
     &     - 112/3.D0*Hr1(0)*upx**(-6)*y**2
     &     + 128/3.D0*Hr1(0)*upx**(-5)
     &     + 96.D0*Hr1(0)*upx**(-5)*y
     &     + 128/3.D0*Hr1(0)*upx**(-5)*y**2
     &     - 40/3.D0*Hr1(0)*upx**(-4)
     &     - 160/3.D0*Hr1(0)*upx**(-4)*y
     &     - 40/3.D0*Hr1(0)*upx**(-4)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-3)*y
     &     - 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-2)
     &     + 8.D0*Hr1(0)*upx**(-2)*y
     &     + 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 8/3.D0*Hr1(0)*upx**(-1)*y
     &     )
      T34PoleEP0qq = T34PoleEP0qq + Nc*Nl * (
     &     - 20.D0/9.D0
     &     - 40/9.D0*upx**(-4)
     &     - 80/9.D0*upx**(-4)*y
     &     - 40/9.D0*upx**(-4)*y**2
     &     + 80/9.D0*upx**(-3)
     &     + 160/9.D0*upx**(-3)*y
     &     + 80/9.D0*upx**(-3)*y**2
     &     - 40/9.D0*upx**(-2)
     &     - 40/3.D0*upx**(-2)*y
     &     - 40/9.D0*upx**(-2)*y**2
     &     + 40/9.D0*upx**(-1)*y
     &     + 8/3.D0*Hr1(-1)
     &     + 16/3.D0*Hr1(-1)*upx**(-4)
     &     + 32/3.D0*Hr1(-1)*upx**(-4)*y
     &     + 16/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 32/3.D0*Hr1(-1)*upx**(-3)
     &     - 64/3.D0*Hr1(-1)*upx**(-3)*y
     &     - 32/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 16/3.D0*Hr1(-1)*upx**(-2)
     &     + 16.D0*Hr1(-1)*upx**(-2)*y
     &     + 16/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 16/3.D0*Hr1(-1)*upx**(-1)*y
     &     - 4/3.D0*Hr1(0)
     &     - 8/3.D0*Hr1(0)*upx**(-4)
     &     - 16/3.D0*Hr1(0)*upx**(-4)*y
     &     - 8/3.D0*Hr1(0)*upx**(-4)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-3)
     &     + 32/3.D0*Hr1(0)*upx**(-3)*y
     &     + 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     - 8/3.D0*Hr1(0)*upx**(-2)
     &     - 8.D0*Hr1(0)*upx**(-2)*y
     &     - 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-1)*y
     &     - 4/3.D0*DLog(tm**(-2)*w**2)
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 8.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T34PoleEP0qq = T34PoleEP0qq + Nc**2 * (
     &     + 107.D0/9.D0
     &     + 6.D0*umx**(-5)*z2
     &     - 12.D0*umx**(-5)*y*z2
     &     + 6.D0*umx**(-5)*y**2*z2
     &     - 15.D0*umx**(-4)*z2
     &     + 30.D0*umx**(-4)*y*z2
     &     - 15.D0*umx**(-4)*y**2*z2
     &     + 31/2.D0*umx**(-3)*z2
     &     - 27.D0*umx**(-3)*y*z2
     &     + 19/2.D0*umx**(-3)*y**2*z2
     &     + umx**(-2)
     &     - 33/4.D0*umx**(-2)*z2
     &     - 2.D0*umx**(-2)*y
     &     + 21/2.D0*umx**(-2)*y*z2
     &     + umx**(-2)*y**2
     &     + 3/4.D0*umx**(-2)*y**2*z2
     &     - umx**(-1)
     &     + 41/8.D0*umx**(-1)*z2
     &     + 2.D0*umx**(-1)*y
     &     - 3/4.D0*umx**(-1)*y*z2
     &     - umx**(-1)*y**2
     &     + 1/8.D0*umx**(-1)*y**2*z2
     &     + 88/9.D0*upx**(-4)
     &     + 4.D0*upx**(-4)*z2
     &     + 176/9.D0*upx**(-4)*y
     &     + 8.D0*upx**(-4)*y*z2
     &     + 88/9.D0*upx**(-4)*y**2
     &     + 4.D0*upx**(-4)*y**2*z2
     &     - 176/9.D0*upx**(-3)
     &     - 7.D0*upx**(-3)*z2
     &     - 352/9.D0*upx**(-3)*y
     &     - 12.D0*upx**(-3)*y*z2
     &     - 176/9.D0*upx**(-3)*y**2
     &     - 5.D0*upx**(-3)*y**2*z2
     &     + 97/9.D0*upx**(-2)
     &     + 5/2.D0*upx**(-2)*z2
     &     + 94/3.D0*upx**(-2)*y
     &     + 4.D0*upx**(-2)*y*z2
     &     + 97/9.D0*upx**(-2)*y**2
     &     - 1/2.D0*upx**(-2)*y**2*z2
     &     - upx**(-1)
     &     - 7/8.D0*upx**(-1)*z2
     &     - 106/9.D0*upx**(-1)*y
     &     - 3/4.D0*upx**(-1)*y*z2
     &     - upx**(-1)*y**2
     &     + 1/8.D0*upx**(-1)*y**2*z2
     &     - 20/3.D0*Hr1(-1)
     &     - 12.D0*Hr1(-1)*umx**(-4)
     &     + 24.D0*Hr1(-1)*umx**(-4)*y
     &     - 12.D0*Hr1(-1)*umx**(-4)*y**2
     &     + 24.D0*Hr1(-1)*umx**(-3)
     &     - 48.D0*Hr1(-1)*umx**(-3)*y
     &     + 24.D0*Hr1(-1)*umx**(-3)*y**2
     &     - 20.D0*Hr1(-1)*umx**(-2)
     &     + 32.D0*Hr1(-1)*umx**(-2)*y
     &     - 8.D0*Hr1(-1)*umx**(-2)*y**2
     &     + 8.D0*Hr1(-1)*umx**(-1)
     &     - 8.D0*Hr1(-1)*umx**(-1)*y
     &     - 4.D0*Hr1(-1)*umx**(-1)*y**2
     &     - 28/3.D0*Hr1(-1)*upx**(-4)
     &     - 56/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 28/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 56/3.D0*Hr1(-1)*upx**(-3)
     &     + 112/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 56/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 16/3.D0*Hr1(-1)*upx**(-2)
     &     - 20.D0*Hr1(-1)*upx**(-2)*y
     &     - 16/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 4.D0*Hr1(-1)*upx**(-1)
     &     + 4/3.D0*Hr1(-1)*upx**(-1)*y
     &     - 4.D0*Hr1(-1)*upx**(-1)*y**2
     &     - 4.D0*Hr1(-1)*HAr1(-1)
     &     - 16.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 16.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 32.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 64.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 32.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 16.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 40.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 16.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 8.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     - 2.D0*HAr1(-1)
     &     + 2.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 4.D0*HAr1(-1)*upx**(-2)
     &     + 2.D0*HAr1(-1)*upx**(-2)*y
     &     - 2.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 4.D0*HAr1(-1)*upx**(-1)
     &     - 2.D0*HAr1(-1)*upx**(-1)*y
     &     - 2.D0*HAr1(-1)*y**(-1)
     &     + 2.D0*HAr1(-1)*Hr1(0)
     &     + 8.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 8.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 16.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 32.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 16.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 8.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 20.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 8.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 4.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 4.D0*Hr2(-1,-1)
     &     + 16.D0*Hr2(-1,-1)*upx**(-4)
     &     + 32.D0*Hr2(-1,-1)*upx**(-4)*y
     &     + 16.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 32.D0*Hr2(-1,-1)*upx**(-3)
     &     - 64.D0*Hr2(-1,-1)*upx**(-3)*y
     &     - 32.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 16.D0*Hr2(-1,-1)*upx**(-2)
     &     + 40.D0*Hr2(-1,-1)*upx**(-2)*y
     &     + 16.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     - 8.D0*Hr2(-1,-1)*upx**(-1)*y
     &     - 4.D0*HAr2(-1,-1)
     &     - 8.D0*HAr2(-1,-1)*upx**(-2)*y
     &     + 8.D0*HAr2(-1,-1)*upx**(-1)*y
     &     - 2.D0*Hr2(-1,0)
     &     - 8.D0*Hr2(-1,0)*upx**(-4)
     &     - 16.D0*Hr2(-1,0)*upx**(-4)*y
     &     - 8.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     + 16.D0*Hr2(-1,0)*upx**(-3)
     &     + 32.D0*Hr2(-1,0)*upx**(-3)*y
     &     + 16.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     - 8.D0*Hr2(-1,0)*upx**(-2)
     &     - 20.D0*Hr2(-1,0)*upx**(-2)*y
     &     - 8.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     + 4.D0*Hr2(-1,0)*upx**(-1)*y
     &     + 10/3.D0*Hr1(0)
     &     + 6.D0*Hr1(0)*umx**(-4)
     &     - 12.D0*Hr1(0)*umx**(-4)*y
     &     + 6.D0*Hr1(0)*umx**(-4)*y**2
     &     - 12.D0*Hr1(0)*umx**(-3)
     &     + 24.D0*Hr1(0)*umx**(-3)*y
     &     - 12.D0*Hr1(0)*umx**(-3)*y**2
     &     + 10.D0*Hr1(0)*umx**(-2)
     &     - 16.D0*Hr1(0)*umx**(-2)*y
     &     + 4.D0*Hr1(0)*umx**(-2)*y**2
     &     - 4.D0*Hr1(0)*umx**(-1)
     &     + 4.D0*Hr1(0)*umx**(-1)*y
     &     + 2.D0*Hr1(0)*umx**(-1)*y**2
     &     + 14/3.D0*Hr1(0)*upx**(-4)
     &     + 28/3.D0*Hr1(0)*upx**(-4)*y
     &     + 14/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 28/3.D0*Hr1(0)*upx**(-3)
     &     - 56/3.D0*Hr1(0)*upx**(-3)*y
     &     - 28/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-2)
     &     + 10.D0*Hr1(0)*upx**(-2)*y
     &     + 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     + 2.D0*Hr1(0)*upx**(-1)
     &     - 2/3.D0*Hr1(0)*upx**(-1)*y
     &     + 2.D0*Hr1(0)*upx**(-1)*y**2
     &     - 12.D0*Hr2(0,-1)*umx**(-5)
     &     + 24.D0*Hr2(0,-1)*umx**(-5)*y
     &     - 12.D0*Hr2(0,-1)*umx**(-5)*y**2
     &     + 30.D0*Hr2(0,-1)*umx**(-4)
     &     - 60.D0*Hr2(0,-1)*umx**(-4)*y
     &     + 30.D0*Hr2(0,-1)*umx**(-4)*y**2
     &     - 31.D0*Hr2(0,-1)*umx**(-3)
     &     + 54.D0*Hr2(0,-1)*umx**(-3)*y
     &     - 19.D0*Hr2(0,-1)*umx**(-3)*y**2
     &     + 33/2.D0*Hr2(0,-1)*umx**(-2)
     &     - 21.D0*Hr2(0,-1)*umx**(-2)*y
     &     - 3/2.D0*Hr2(0,-1)*umx**(-2)*y**2
     &     - 41/4.D0*Hr2(0,-1)*umx**(-1)
     &     + 3/2.D0*Hr2(0,-1)*umx**(-1)*y
     &     - 1/4.D0*Hr2(0,-1)*umx**(-1)*y**2
     &     - 8.D0*Hr2(0,-1)*upx**(-4)
     &     - 16.D0*Hr2(0,-1)*upx**(-4)*y
     &     - 8.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     + 14.D0*Hr2(0,-1)*upx**(-3)
     &     + 24.D0*Hr2(0,-1)*upx**(-3)*y
     &     + 10.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     - 5.D0*Hr2(0,-1)*upx**(-2)
     &     - 8.D0*Hr2(0,-1)*upx**(-2)*y
     &     + Hr2(0,-1)*upx**(-2)*y**2
     &     + 7/4.D0*Hr2(0,-1)*upx**(-1)
     &     + 3/2.D0*Hr2(0,-1)*upx**(-1)*y
     &     - 1/4.D0*Hr2(0,-1)*upx**(-1)*y**2
     &     + 2.D0*HAr2(0,-1)
     &     + 4.D0*HAr2(0,-1)*upx**(-2)*y
     &     - 4.D0*HAr2(0,-1)*upx**(-1)*y
     &     + 6.D0*Hr2(0,0)*umx**(-5)
     &     - 12.D0*Hr2(0,0)*umx**(-5)*y
     &     + 6.D0*Hr2(0,0)*umx**(-5)*y**2
     &     - 15.D0*Hr2(0,0)*umx**(-4)
     &     + 30.D0*Hr2(0,0)*umx**(-4)*y
     &     - 15.D0*Hr2(0,0)*umx**(-4)*y**2
     &     + 31/2.D0*Hr2(0,0)*umx**(-3)
     &     - 27.D0*Hr2(0,0)*umx**(-3)*y
     &     + 19/2.D0*Hr2(0,0)*umx**(-3)*y**2
     &     - 33/4.D0*Hr2(0,0)*umx**(-2)
     &     + 21/2.D0*Hr2(0,0)*umx**(-2)*y
     &     + 3/4.D0*Hr2(0,0)*umx**(-2)*y**2
     &     + 41/8.D0*Hr2(0,0)*umx**(-1)
     &     - 3/4.D0*Hr2(0,0)*umx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*umx**(-1)*y**2
     &     + 4.D0*Hr2(0,0)*upx**(-4)
     &     + 8.D0*Hr2(0,0)*upx**(-4)*y
     &     + 4.D0*Hr2(0,0)*upx**(-4)*y**2
     &     - 7.D0*Hr2(0,0)*upx**(-3)
     &     - 12.D0*Hr2(0,0)*upx**(-3)*y
     &     - 5.D0*Hr2(0,0)*upx**(-3)*y**2
     &     + 5/2.D0*Hr2(0,0)*upx**(-2)
     &     + 4.D0*Hr2(0,0)*upx**(-2)*y
     &     - 1/2.D0*Hr2(0,0)*upx**(-2)*y**2
     &     - 7/8.D0*Hr2(0,0)*upx**(-1)
     &     - 3/4.D0*Hr2(0,0)*upx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*upx**(-1)*y**2
     &     + 22/3.D0*DLog(tm**(-2)*w**2)
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 44.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      return
      end

*
*******************************************************************************
************************************************************************
************************************************************************
*
*          FUNCTIONS for the Cross Section
*
*************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function T13PoleEPm2qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T13PoleEPm2qq =  + Nc**(-2) * (
     &     - 16.D0
     &     - 32.D0*upx**(-4)
     &     - 64.D0*upx**(-4)*y
     &     - 32.D0*upx**(-4)*y**2
     &     + 64.D0*upx**(-3)
     &     + 128.D0*upx**(-3)*y
     &     + 64.D0*upx**(-3)*y**2
     &     - 32.D0*upx**(-2)
     &     - 96.D0*upx**(-2)*y
     &     - 32.D0*upx**(-2)*y**2
     &     + 32.D0*upx**(-1)*y
     &     )/4.D0
      T13PoleEPm2qq = T13PoleEPm2qq + (
     &     + 40.D0
     &     + 80.D0*upx**(-4)
     &     + 160.D0*upx**(-4)*y
     &     + 80.D0*upx**(-4)*y**2
     &     - 160.D0*upx**(-3)
     &     - 320.D0*upx**(-3)*y
     &     - 160.D0*upx**(-3)*y**2
     &     + 80.D0*upx**(-2)
     &     + 240.D0*upx**(-2)*y
     &     + 80.D0*upx**(-2)*y**2
     &     - 80.D0*upx**(-1)*y
     &     )/4.D0
      T13PoleEPm2qq = T13PoleEPm2qq + Nc**2 * (
     &     - 32.D0
     &     - 64.D0*upx**(-4)
     &     - 128.D0*upx**(-4)*y
     &     - 64.D0*upx**(-4)*y**2
     &     + 128.D0*upx**(-3)
     &     + 256.D0*upx**(-3)*y
     &     + 128.D0*upx**(-3)*y**2
     &     - 64.D0*upx**(-2)
     &     - 192.D0*upx**(-2)*y
     &     - 64.D0*upx**(-2)*y**2
     &     + 64.D0*upx**(-1)*y
     &     )/4.D0
      T13PoleEPm2qq = T13PoleEPm2qq + Nc**4 * (
     &     + 8.D0
     &     + 16.D0*upx**(-4)
     &     + 32.D0*upx**(-4)*y
     &     + 16.D0*upx**(-4)*y**2
     &     - 32.D0*upx**(-3)
     &     - 64.D0*upx**(-3)*y
     &     - 32.D0*upx**(-3)*y**2
     &     + 16.D0*upx**(-2)
     &     + 48.D0*upx**(-2)*y
     &     + 16.D0*upx**(-2)*y**2
     &     - 16.D0*upx**(-1)*y
     &     )/4.D0
      return
      end
*
*
*
      double precision function T13PoleEPm1qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T13PoleEPm1qq =  + Nc**(-2) * (
     &     + 12.D0
     &     + 40.D0*upx**(-4)
     &     + 80.D0*upx**(-4)*y
     &     + 40.D0*upx**(-4)*y**2
     &     - 80.D0*upx**(-3)
     &     - 160.D0*upx**(-3)*y
     &     - 80.D0*upx**(-3)*y**2
     &     + 40.D0*upx**(-2)
     &     + 120.D0*upx**(-2)*y
     &     + 40.D0*upx**(-2)*y**2
     &     - 40.D0*upx**(-1)*y
     &     - 16.D0*Hr1(-1)
     &     - 32.D0*Hr1(-1)*upx**(-4)
     &     - 64.D0*Hr1(-1)*upx**(-4)*y
     &     - 32.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 64.D0*Hr1(-1)*upx**(-3)
     &     + 128.D0*Hr1(-1)*upx**(-3)*y
     &     + 64.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 32.D0*Hr1(-1)*upx**(-2)
     &     - 96.D0*Hr1(-1)*upx**(-2)*y
     &     - 32.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 32.D0*Hr1(-1)*upx**(-1)*y
     &     - 24.D0*HAr1(-1)
     &     - 48.D0*HAr1(-1)*upx**(-4)
     &     - 96.D0*HAr1(-1)*upx**(-4)*y
     &     - 48.D0*HAr1(-1)*upx**(-4)*y**2
     &     + 96.D0*HAr1(-1)*upx**(-3)
     &     + 192.D0*HAr1(-1)*upx**(-3)*y
     &     + 96.D0*HAr1(-1)*upx**(-3)*y**2
     &     - 48.D0*HAr1(-1)*upx**(-2)
     &     - 144.D0*HAr1(-1)*upx**(-2)*y
     &     - 48.D0*HAr1(-1)*upx**(-2)*y**2
     &     + 48.D0*HAr1(-1)*upx**(-1)*y
     &     + 24.D0*HBr1(-1)
     &     + 48.D0*HBr1(-1)*upx**(-4)
     &     + 96.D0*HBr1(-1)*upx**(-4)*y
     &     + 48.D0*HBr1(-1)*upx**(-4)*y**2
     &     - 96.D0*HBr1(-1)*upx**(-3)
     &     - 192.D0*HBr1(-1)*upx**(-3)*y
     &     - 96.D0*HBr1(-1)*upx**(-3)*y**2
     &     + 48.D0*HBr1(-1)*upx**(-2)
     &     + 144.D0*HBr1(-1)*upx**(-2)*y
     &     + 48.D0*HBr1(-1)*upx**(-2)*y**2
     &     - 48.D0*HBr1(-1)*upx**(-1)*y
     &     + 9.D0*Hr1(0)*umx**(-1)
     &     - 2.D0*Hr1(0)*umx**(-1)*y
     &     + Hr1(0)*umx**(-1)*y**2
     &     + 16.D0*Hr1(0)*upx**(-5)
     &     + 32.D0*Hr1(0)*upx**(-5)*y
     &     + 16.D0*Hr1(0)*upx**(-5)*y**2
     &     - 24.D0*Hr1(0)*upx**(-4)
     &     - 48.D0*Hr1(0)*upx**(-4)*y
     &     - 24.D0*Hr1(0)*upx**(-4)*y**2
     &     + 4.D0*Hr1(0)*upx**(-3)
     &     + 24.D0*Hr1(0)*upx**(-3)*y
     &     + 4.D0*Hr1(0)*upx**(-3)*y**2
     &     + 2.D0*Hr1(0)*upx**(-2)
     &     - 4.D0*Hr1(0)*upx**(-2)*y
     &     + 2.D0*Hr1(0)*upx**(-2)*y**2
     &     + 9.D0*Hr1(0)*upx**(-1)
     &     - 2.D0*Hr1(0)*upx**(-1)*y
     &     + Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T13PoleEPm1qq = T13PoleEPm1qq + (
     &     - 30.D0
     &     - 100.D0*upx**(-4)
     &     - 200.D0*upx**(-4)*y
     &     - 100.D0*upx**(-4)*y**2
     &     + 200.D0*upx**(-3)
     &     + 400.D0*upx**(-3)*y
     &     + 200.D0*upx**(-3)*y**2
     &     - 100.D0*upx**(-2)
     &     - 300.D0*upx**(-2)*y
     &     - 100.D0*upx**(-2)*y**2
     &     + 100.D0*upx**(-1)*y
     &     + 24.D0*Hr1(-1)
     &     + 48.D0*Hr1(-1)*upx**(-4)
     &     + 96.D0*Hr1(-1)*upx**(-4)*y
     &     + 48.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 96.D0*Hr1(-1)*upx**(-3)
     &     - 192.D0*Hr1(-1)*upx**(-3)*y
     &     - 96.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 48.D0*Hr1(-1)*upx**(-2)
     &     + 144.D0*Hr1(-1)*upx**(-2)*y
     &     + 48.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 48.D0*Hr1(-1)*upx**(-1)*y
     &     + 48.D0*HAr1(-1)
     &     + 96.D0*HAr1(-1)*upx**(-4)
     &     + 192.D0*HAr1(-1)*upx**(-4)*y
     &     + 96.D0*HAr1(-1)*upx**(-4)*y**2
     &     - 192.D0*HAr1(-1)*upx**(-3)
     &     - 384.D0*HAr1(-1)*upx**(-3)*y
     &     - 192.D0*HAr1(-1)*upx**(-3)*y**2
     &     + 96.D0*HAr1(-1)*upx**(-2)
     &     + 288.D0*HAr1(-1)*upx**(-2)*y
     &     + 96.D0*HAr1(-1)*upx**(-2)*y**2
     &     - 96.D0*HAr1(-1)*upx**(-1)*y
     &     - 32.D0*HBr1(-1)
     &     - 64.D0*HBr1(-1)*upx**(-4)
     &     - 128.D0*HBr1(-1)*upx**(-4)*y
     &     - 64.D0*HBr1(-1)*upx**(-4)*y**2
     &     + 128.D0*HBr1(-1)*upx**(-3)
     &     + 256.D0*HBr1(-1)*upx**(-3)*y
     &     + 128.D0*HBr1(-1)*upx**(-3)*y**2
     &     - 64.D0*HBr1(-1)*upx**(-2)
     &     - 192.D0*HBr1(-1)*upx**(-2)*y
     &     - 64.D0*HBr1(-1)*upx**(-2)*y**2
     &     + 64.D0*HBr1(-1)*upx**(-1)*y
     &     - 27/2.D0*Hr1(0)*umx**(-1)
     &     + 3.D0*Hr1(0)*umx**(-1)*y
     &     - 3/2.D0*Hr1(0)*umx**(-1)*y**2
     &     - 24.D0*Hr1(0)*upx**(-5)
     &     - 48.D0*Hr1(0)*upx**(-5)*y
     &     - 24.D0*Hr1(0)*upx**(-5)*y**2
     &     + 36.D0*Hr1(0)*upx**(-4)
     &     + 72.D0*Hr1(0)*upx**(-4)*y
     &     + 36.D0*Hr1(0)*upx**(-4)*y**2
     &     - 6.D0*Hr1(0)*upx**(-3)
     &     - 36.D0*Hr1(0)*upx**(-3)*y
     &     - 6.D0*Hr1(0)*upx**(-3)*y**2
     &     - 3.D0*Hr1(0)*upx**(-2)
     &     + 6.D0*Hr1(0)*upx**(-2)*y
     &     - 3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 27/2.D0*Hr1(0)*upx**(-1)
     &     + 3.D0*Hr1(0)*upx**(-1)*y
     &     - 3/2.D0*Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T13PoleEPm1qq = T13PoleEPm1qq + Nc**2 * (
     &     + 24.D0
     &     + 80.D0*upx**(-4)
     &     + 160.D0*upx**(-4)*y
     &     + 80.D0*upx**(-4)*y**2
     &     - 160.D0*upx**(-3)
     &     - 320.D0*upx**(-3)*y
     &     - 160.D0*upx**(-3)*y**2
     &     + 80.D0*upx**(-2)
     &     + 240.D0*upx**(-2)*y
     &     + 80.D0*upx**(-2)*y**2
     &     - 80.D0*upx**(-1)*y
     &     - 8.D0*Hr1(-1)
     &     - 16.D0*Hr1(-1)*upx**(-4)
     &     - 32.D0*Hr1(-1)*upx**(-4)*y
     &     - 16.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 32.D0*Hr1(-1)*upx**(-3)
     &     + 64.D0*Hr1(-1)*upx**(-3)*y
     &     + 32.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 16.D0*Hr1(-1)*upx**(-2)
     &     - 48.D0*Hr1(-1)*upx**(-2)*y
     &     - 16.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 16.D0*Hr1(-1)*upx**(-1)*y
     &     - 32.D0*HAr1(-1)
     &     - 64.D0*HAr1(-1)*upx**(-4)
     &     - 128.D0*HAr1(-1)*upx**(-4)*y
     &     - 64.D0*HAr1(-1)*upx**(-4)*y**2
     &     + 128.D0*HAr1(-1)*upx**(-3)
     &     + 256.D0*HAr1(-1)*upx**(-3)*y
     &     + 128.D0*HAr1(-1)*upx**(-3)*y**2
     &     - 64.D0*HAr1(-1)*upx**(-2)
     &     - 192.D0*HAr1(-1)*upx**(-2)*y
     &     - 64.D0*HAr1(-1)*upx**(-2)*y**2
     &     + 64.D0*HAr1(-1)*upx**(-1)*y
     &     + 8.D0*HBr1(-1)
     &     + 16.D0*HBr1(-1)*upx**(-4)
     &     + 32.D0*HBr1(-1)*upx**(-4)*y
     &     + 16.D0*HBr1(-1)*upx**(-4)*y**2
     &     - 32.D0*HBr1(-1)*upx**(-3)
     &     - 64.D0*HBr1(-1)*upx**(-3)*y
     &     - 32.D0*HBr1(-1)*upx**(-3)*y**2
     &     + 16.D0*HBr1(-1)*upx**(-2)
     &     + 48.D0*HBr1(-1)*upx**(-2)*y
     &     + 16.D0*HBr1(-1)*upx**(-2)*y**2
     &     - 16.D0*HBr1(-1)*upx**(-1)*y
     &     + 9/2.D0*Hr1(0)*umx**(-1)
     &     - Hr1(0)*umx**(-1)*y
     &     + 1/2.D0*Hr1(0)*umx**(-1)*y**2
     &     + 8.D0*Hr1(0)*upx**(-5)
     &     + 16.D0*Hr1(0)*upx**(-5)*y
     &     + 8.D0*Hr1(0)*upx**(-5)*y**2
     &     - 12.D0*Hr1(0)*upx**(-4)
     &     - 24.D0*Hr1(0)*upx**(-4)*y
     &     - 12.D0*Hr1(0)*upx**(-4)*y**2
     &     + 2.D0*Hr1(0)*upx**(-3)
     &     + 12.D0*Hr1(0)*upx**(-3)*y
     &     + 2.D0*Hr1(0)*upx**(-3)*y**2
     &     + Hr1(0)*upx**(-2)
     &     - 2.D0*Hr1(0)*upx**(-2)*y
     &     + Hr1(0)*upx**(-2)*y**2
     &     + 9/2.D0*Hr1(0)*upx**(-1)
     &     - Hr1(0)*upx**(-1)*y
     &     + 1/2.D0*Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T13PoleEPm1qq = T13PoleEPm1qq + Nc**4 * (
     &     - 6.D0
     &     - 20.D0*upx**(-4)
     &     - 40.D0*upx**(-4)*y
     &     - 20.D0*upx**(-4)*y**2
     &     + 40.D0*upx**(-3)
     &     + 80.D0*upx**(-3)*y
     &     + 40.D0*upx**(-3)*y**2
     &     - 20.D0*upx**(-2)
     &     - 60.D0*upx**(-2)*y
     &     - 20.D0*upx**(-2)*y**2
     &     + 20.D0*upx**(-1)*y
     &     + 8.D0*HAr1(-1)
     &     + 16.D0*HAr1(-1)*upx**(-4)
     &     + 32.D0*HAr1(-1)*upx**(-4)*y
     &     + 16.D0*HAr1(-1)*upx**(-4)*y**2
     &     - 32.D0*HAr1(-1)*upx**(-3)
     &     - 64.D0*HAr1(-1)*upx**(-3)*y
     &     - 32.D0*HAr1(-1)*upx**(-3)*y**2
     &     + 16.D0*HAr1(-1)*upx**(-2)
     &     + 48.D0*HAr1(-1)*upx**(-2)*y
     &     + 16.D0*HAr1(-1)*upx**(-2)*y**2
     &     - 16.D0*HAr1(-1)*upx**(-1)*y
     &     )/(-2.D0)
      
      return
      end
*
*
*
      double precision function T13PoleEP0qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1d0

* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T13PoleEP0qq =  + Nc**(-2) * (
     &     - 14.D0
     &     + 12.D0*umx**(-3)*z2
     &     - 12.D0*umx**(-3)*y*z2
     &     - 18.D0*umx**(-2)*z2
     &     + 18.D0*umx**(-2)*y*z2
     &     + 24.D0*umx**(-1)*z2
     &     - 4.D0*umx**(-1)*y*z2
     &     + 2.D0*umx**(-1)*y**2*z2
     &     + 32.D0*upx**(-5)*z2
     &     + 64.D0*upx**(-5)*y*z2
     &     + 32.D0*upx**(-5)*y**2*z2
     &     - 48.D0*upx**(-4)
     &     - 48.D0*upx**(-4)*z2
     &     - 96.D0*upx**(-4)*y
     &     - 96.D0*upx**(-4)*y*z2
     &     - 48.D0*upx**(-4)*y**2
     &     - 48.D0*upx**(-4)*y**2*z2
     &     + 96.D0*upx**(-3)
     &     - 4.D0*upx**(-3)*z2
     &     + 192.D0*upx**(-3)*y
     &     + 36.D0*upx**(-3)*y*z2
     &     + 96.D0*upx**(-3)*y**2
     &     + 8.D0*upx**(-3)*y**2*z2
     &     - 48.D0*upx**(-2)
     &     + 10.D0*upx**(-2)*z2
     &     - 144.D0*upx**(-2)*y
     &     - 2.D0*upx**(-2)*y*z2
     &     - 48.D0*upx**(-2)*y**2
     &     + 4.D0*upx**(-2)*y**2*z2
     &     + 12.D0*upx**(-1)*z2
     &     + 48.D0*upx**(-1)*y
     &     - 4.D0*upx**(-1)*y*z2
     &     + 2.D0*upx**(-1)*y**2*z2
     &     - 8.D0*Hr1(-1)
     &     - 24.D0*Hr1(-1)*umx**(-2)
     &     + 24.D0*Hr1(-1)*umx**(-2)*y
     &     + 24.D0*Hr1(-1)*umx**(-1)
     &     - 24.D0*Hr1(-1)*umx**(-1)*y
     &     + 24.D0*Hr1(-1)*upx**(-4)
     &     + 48.D0*Hr1(-1)*upx**(-4)*y
     &     + 24.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 48.D0*Hr1(-1)*upx**(-3)
     &     - 96.D0*Hr1(-1)*upx**(-3)*y
     &     - 48.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 24.D0*Hr1(-1)*upx**(-2)
     &     + 72.D0*Hr1(-1)*upx**(-2)*y
     &     + 24.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 24.D0*Hr1(-1)*upx**(-1)*y
     &     - 12.D0*Hr1(-1)*HAr1(-1)
     &     - 48.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 96.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 48.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 192.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 48.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 120.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 48.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 24.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 36.D0*Hr1(-1)*HBr1(-1)
     &     + 48.D0*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 96.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 48.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 96.D0*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 192.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 96.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 96.D0*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 168.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 48.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 48.D0*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 72.D0*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 6.D0*HAr1(-1)
     &     + 6.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 12.D0*HAr1(-1)*upx**(-2)
     &     + 6.D0*HAr1(-1)*upx**(-2)*y
     &     - 6.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 12.D0*HAr1(-1)*upx**(-1)
     &     - 6.D0*HAr1(-1)*upx**(-1)*y
     &     - 6.D0*HAr1(-1)*y**(-1)
     &     + 6.D0*HAr1(-1)*Hr1(0)
     &     + 24.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 48.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 24.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 96.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 24.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 60.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 24.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 12.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 12.D0*HBr1(-1)
     &     - 6.D0*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 6.D0*HBr1(-1)*upx**(-2)*y
     &     + 6.D0*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 6.D0*HBr1(-1)*upx**(-1)*y
     &     + 6.D0*HBr1(-1)*z**(-1)
     &     - 18.D0*HBr1(-1)*Hr1(0)
     &     - 24.D0*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 48.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 24.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 48.D0*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 96.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 48.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 48.D0*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 84.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 24.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 24.D0*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 36.D0*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 40.D0*Hr2(-1,-1)
     &     - 32.D0*Hr2(-1,-1)*upx**(-4)
     &     - 64.D0*Hr2(-1,-1)*upx**(-4)*y
     &     - 32.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 64.D0*Hr2(-1,-1)*upx**(-3)
     &     + 128.D0*Hr2(-1,-1)*upx**(-3)*y
     &     + 64.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 80.D0*Hr2(-1,-1)*upx**(-2)
     &     - 144.D0*Hr2(-1,-1)*upx**(-2)*y
     &     - 32.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 48.D0*Hr2(-1,-1)*upx**(-1)
     &     + 80.D0*Hr2(-1,-1)*upx**(-1)*y
     &     - 12.D0*HAr2(-1,-1)
     &     - 24.D0*HAr2(-1,-1)*upx**(-2)*y
     &     + 24.D0*HAr2(-1,-1)*upx**(-1)*y
     &     - 12.D0*HBr2(-1,-1)
     &     - 48.D0*HBr2(-1,-1)*upx**(-2)
     &     - 24.D0*HBr2(-1,-1)*upx**(-2)*y
     &     + 48.D0*HBr2(-1,-1)*upx**(-1)
     &     + 24.D0*HBr2(-1,-1)*upx**(-1)*y
     &     + 20.D0*Hr2(-1,0)
     &     + 16.D0*Hr2(-1,0)*upx**(-4)
     &     + 32.D0*Hr2(-1,0)*upx**(-4)*y
     &     + 16.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     - 32.D0*Hr2(-1,0)*upx**(-3)
     &     - 64.D0*Hr2(-1,0)*upx**(-3)*y
     &     - 32.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     + 40.D0*Hr2(-1,0)*upx**(-2)
     &     + 72.D0*Hr2(-1,0)*upx**(-2)*y
     &     + 16.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     - 24.D0*Hr2(-1,0)*upx**(-1)
     &     - 40.D0*Hr2(-1,0)*upx**(-1)*y
     &     + 6.D0*Hr1(0)
     &     + 12.D0*Hr1(0)*umx**(-2)
     &     - 12.D0*Hr1(0)*umx**(-2)*y
     &     - 17/2.D0*Hr1(0)*umx**(-1)
     &     + 13.D0*Hr1(0)*umx**(-1)*y
     &     - 1/2.D0*Hr1(0)*umx**(-1)*y**2
     &     - 16.D0*Hr1(0)*upx**(-5)
     &     - 32.D0*Hr1(0)*upx**(-5)*y
     &     - 16.D0*Hr1(0)*upx**(-5)*y**2
     &     + 28.D0*Hr1(0)*upx**(-4)
     &     + 56.D0*Hr1(0)*upx**(-4)*y
     &     + 28.D0*Hr1(0)*upx**(-4)*y**2
     &     - 10.D0*Hr1(0)*upx**(-3)
     &     - 36.D0*Hr1(0)*upx**(-3)*y
     &     - 10.D0*Hr1(0)*upx**(-3)*y**2
     &     - Hr1(0)*upx**(-2)
     &     + 10.D0*Hr1(0)*upx**(-2)*y
     &     - Hr1(0)*upx**(-2)*y**2
     &     - 17/2.D0*Hr1(0)*upx**(-1)
     &     + Hr1(0)*upx**(-1)*y
     &     - 1/2.D0*Hr1(0)*upx**(-1)*y**2
     &     + 8.D0*Hr2(0,-1)
     &     - 24.D0*Hr2(0,-1)*umx**(-3)
     &     + 24.D0*Hr2(0,-1)*umx**(-3)*y
     &     + 36.D0*Hr2(0,-1)*umx**(-2)
     &     - 36.D0*Hr2(0,-1)*umx**(-2)*y
     &     - 12.D0*Hr2(0,-1)*umx**(-1)
     &     + 16.D0*Hr2(0,-1)*upx**(-4)
     &     + 32.D0*Hr2(0,-1)*upx**(-4)*y
     &     + 16.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     - 8.D0*Hr2(0,-1)*upx**(-3)
     &     - 40.D0*Hr2(0,-1)*upx**(-3)*y
     &     - 32.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     + 4.D0*Hr2(0,-1)*upx**(-2)
     &     + 36.D0*Hr2(0,-1)*upx**(-2)*y
     &     + 16.D0*Hr2(0,-1)*upx**(-2)*y**2
     &     + 12.D0*Hr2(0,-1)*upx**(-1)
     &     - 16.D0*Hr2(0,-1)*upx**(-1)*y
     &     + 6.D0*HAr2(0,-1)
     &     + 12.D0*HAr2(0,-1)*upx**(-2)*y
     &     - 12.D0*HAr2(0,-1)*upx**(-1)*y
     &     + 6.D0*HBr2(0,-1)
     &     + 24.D0*HBr2(0,-1)*upx**(-2)
     &     + 12.D0*HBr2(0,-1)*upx**(-2)*y
     &     - 24.D0*HBr2(0,-1)*upx**(-1)
     &     - 12.D0*HBr2(0,-1)*upx**(-1)*y
     &     + 12.D0*Hr2(0,0)*umx**(-3)
     &     - 12.D0*Hr2(0,0)*umx**(-3)*y
     &     - 18.D0*Hr2(0,0)*umx**(-2)
     &     + 18.D0*Hr2(0,0)*umx**(-2)*y
     &     + 3/2.D0*Hr2(0,0)*umx**(-1)
     &     + Hr2(0,0)*umx**(-1)*y
     &     - 1/2.D0*Hr2(0,0)*umx**(-1)*y**2
     &     - 8.D0*Hr2(0,0)*upx**(-5)
     &     - 16.D0*Hr2(0,0)*upx**(-5)*y
     &     - 8.D0*Hr2(0,0)*upx**(-5)*y**2
     &     + 12.D0*Hr2(0,0)*upx**(-4)
     &     + 24.D0*Hr2(0,0)*upx**(-4)*y
     &     + 12.D0*Hr2(0,0)*upx**(-4)*y**2
     &     - 14.D0*Hr2(0,0)*upx**(-3)
     &     - 24.D0*Hr2(0,0)*upx**(-3)*y
     &     - 2.D0*Hr2(0,0)*upx**(-3)*y**2
     &     + 5.D0*Hr2(0,0)*upx**(-2)
     &     + 8.D0*Hr2(0,0)*upx**(-2)*y
     &     - Hr2(0,0)*upx**(-2)*y**2
     &     - 21/2.D0*Hr2(0,0)*upx**(-1)
     &     + Hr2(0,0)*upx**(-1)*y
     &     - 1/2.D0*Hr2(0,0)*upx**(-1)*y**2
     &     + 8.D0*Hr2(1,0)
     &     - 9.D0*Hr2(1,0)*umx**(-1)
     &     + 2.D0*Hr2(1,0)*umx**(-1)*y
     &     - Hr2(1,0)*umx**(-1)*y**2
     &     - 16.D0*Hr2(1,0)*upx**(-5)
     &     - 32.D0*Hr2(1,0)*upx**(-5)*y
     &     - 16.D0*Hr2(1,0)*upx**(-5)*y**2
     &     + 40.D0*Hr2(1,0)*upx**(-4)
     &     + 80.D0*Hr2(1,0)*upx**(-4)*y
     &     + 40.D0*Hr2(1,0)*upx**(-4)*y**2
     &     - 36.D0*Hr2(1,0)*upx**(-3)
     &     - 88.D0*Hr2(1,0)*upx**(-3)*y
     &     - 36.D0*Hr2(1,0)*upx**(-3)*y**2
     &     + 14.D0*Hr2(1,0)*upx**(-2)
     &     + 52.D0*Hr2(1,0)*upx**(-2)*y
     &     + 14.D0*Hr2(1,0)*upx**(-2)*y**2
     &     - 9.D0*Hr2(1,0)*upx**(-1)
     &     - 14.D0*Hr2(1,0)*upx**(-1)*y
     &     - Hr2(1,0)*upx**(-1)*y**2
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**(-1)*Nh * (
     &     + 40.D0/9.D0
     &     - 64/3.D0*upx**(-6)
     &     - 128/3.D0*upx**(-6)*y
     &     - 64/3.D0*upx**(-6)*y**2
     &     + 64.D0*upx**(-5)
     &     + 128.D0*upx**(-5)*y
     &     + 64.D0*upx**(-5)*y**2
     &     - 496/9.D0*upx**(-4)
     &     - 1184/9.D0*upx**(-4)*y
     &     - 496/9.D0*upx**(-4)*y**2
     &     + 32/9.D0*upx**(-3)
     &     + 448/9.D0*upx**(-3)*y
     &     + 32/9.D0*upx**(-3)*y**2
     &     - 16/9.D0*upx**(-2)
     &     + 16/3.D0*upx**(-2)*y
     &     + 80/9.D0*upx**(-2)*y**2
     &     + 32/3.D0*upx**(-1)
     &     - 80/9.D0*upx**(-1)*y
     &     - 8/3.D0*Hr1(0)
     &     - 64/3.D0*Hr1(0)*upx**(-7)
     &     - 128/3.D0*Hr1(0)*upx**(-7)*y
     &     - 64/3.D0*Hr1(0)*upx**(-7)*y**2
     &     + 224/3.D0*Hr1(0)*upx**(-6)
     &     + 448/3.D0*Hr1(0)*upx**(-6)*y
     &     + 224/3.D0*Hr1(0)*upx**(-6)*y**2
     &     - 256/3.D0*Hr1(0)*upx**(-5)
     &     - 192.D0*Hr1(0)*upx**(-5)*y
     &     - 256/3.D0*Hr1(0)*upx**(-5)*y**2
     &     + 80/3.D0*Hr1(0)*upx**(-4)
     &     + 320/3.D0*Hr1(0)*upx**(-4)*y
     &     + 80/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 32/3.D0*Hr1(0)*upx**(-3)*y
     &     + 32/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 32/3.D0*Hr1(0)*upx**(-2)
     &     - 16.D0*Hr1(0)*upx**(-2)*y
     &     - 16/3.D0*Hr1(0)*upx**(-2)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**(-1)*Nl * (
     &     + 40.D0/9.D0
     &     + 80/9.D0*upx**(-4)
     &     + 160/9.D0*upx**(-4)*y
     &     + 80/9.D0*upx**(-4)*y**2
     &     - 160/9.D0*upx**(-3)
     &     - 320/9.D0*upx**(-3)*y
     &     - 160/9.D0*upx**(-3)*y**2
     &     + 80/9.D0*upx**(-2)
     &     + 80/3.D0*upx**(-2)*y
     &     + 80/9.D0*upx**(-2)*y**2
     &     - 80/9.D0*upx**(-1)*y
     &     - 16/3.D0*Hr1(-1)
     &     - 32/3.D0*Hr1(-1)*upx**(-4)
     &     - 64/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 32/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 64/3.D0*Hr1(-1)*upx**(-3)
     &     + 128/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 64/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 32/3.D0*Hr1(-1)*upx**(-2)
     &     - 32.D0*Hr1(-1)*upx**(-2)*y
     &     - 32/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 32/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 8/3.D0*Hr1(0)
     &     + 16/3.D0*Hr1(0)*upx**(-4)
     &     + 32/3.D0*Hr1(0)*upx**(-4)*y
     &     + 16/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 32/3.D0*Hr1(0)*upx**(-3)
     &     - 64/3.D0*Hr1(0)*upx**(-3)*y
     &     - 32/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-2)
     &     + 16.D0*Hr1(0)*upx**(-2)*y
     &     + 16/3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-1)*y
     &     + 8/3.D0*DLog(tm**(-2)*w**2)
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 64/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 16.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + (
     &     - 25.D0/9.D0
     &     - 12.D0*umx**(-5)*z2
     &     + 24.D0*umx**(-5)*y*z2
     &     - 12.D0*umx**(-5)*y**2*z2
     &     + 30.D0*umx**(-4)*z2
     &     - 60.D0*umx**(-4)*y*z2
     &     + 30.D0*umx**(-4)*y**2*z2
     &     - 47.D0*umx**(-3)*z2
     &     + 70.D0*umx**(-3)*y*z2
     &     - 19.D0*umx**(-3)*y**2*z2
     &     - 2.D0*umx**(-2)
     &     + 81/2.D0*umx**(-2)*z2
     &     + 4.D0*umx**(-2)*y
     &     - 45.D0*umx**(-2)*y*z2
     &     - 2.D0*umx**(-2)*y**2
     &     - 3/2.D0*umx**(-2)*y**2*z2
     &     + 2.D0*umx**(-1)
     &     - 181/4.D0*umx**(-1)*z2
     &     - 4.D0*umx**(-1)*y
     &     + 15/2.D0*umx**(-1)*y*z2
     &     + 2.D0*umx**(-1)*y**2
     &     - 13/4.D0*umx**(-1)*y**2*z2
     &     - 48.D0*upx**(-5)*z2
     &     - 96.D0*upx**(-5)*y*z2
     &     - 48.D0*upx**(-5)*y**2*z2
     &     + 472/9.D0*upx**(-4)
     &     + 64.D0*upx**(-4)*z2
     &     + 944/9.D0*upx**(-4)*y
     &     + 128.D0*upx**(-4)*y*z2
     &     + 472/9.D0*upx**(-4)*y**2
     &     + 64.D0*upx**(-4)*y**2*z2
     &     - 944/9.D0*upx**(-3)
     &     + 18.D0*upx**(-3)*z2
     &     - 1888/9.D0*upx**(-3)*y
     &     - 32.D0*upx**(-3)*y*z2
     &     - 944/9.D0*upx**(-3)*y**2
     &     - 2.D0*upx**(-3)*y**2*z2
     &     + 454/9.D0*upx**(-2)
     &     - 19.D0*upx**(-2)*z2
     &     + 460/3.D0*upx**(-2)*y
     &     - 4.D0*upx**(-2)*y*z2
     &     + 454/9.D0*upx**(-2)*y**2
     &     - 5.D0*upx**(-2)*y**2*z2
     &     + 2.D0*upx**(-1)
     &     - 69/4.D0*upx**(-1)*z2
     &     - 436/9.D0*upx**(-1)*y
     &     + 15/2.D0*upx**(-1)*y*z2
     &     + 2.D0*upx**(-1)*y**2
     &     - 13/4.D0*upx**(-1)*y**2*z2
     &     + 70/3.D0*Hr1(-1)
     &     + 24.D0*Hr1(-1)*umx**(-4)
     &     - 48.D0*Hr1(-1)*umx**(-4)*y
     &     + 24.D0*Hr1(-1)*umx**(-4)*y**2
     &     - 48.D0*Hr1(-1)*umx**(-3)
     &     + 96.D0*Hr1(-1)*umx**(-3)*y
     &     - 48.D0*Hr1(-1)*umx**(-3)*y**2
     &     + 72.D0*Hr1(-1)*umx**(-2)
     &     - 96.D0*Hr1(-1)*umx**(-2)*y
     &     + 16.D0*Hr1(-1)*umx**(-2)*y**2
     &     - 48.D0*Hr1(-1)*umx**(-1)
     &     + 48.D0*Hr1(-1)*umx**(-1)*y
     &     + 8.D0*Hr1(-1)*umx**(-1)*y**2
     &     - 52/3.D0*Hr1(-1)*upx**(-4)
     &     - 104/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 52/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 104/3.D0*Hr1(-1)*upx**(-3)
     &     + 208/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 104/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 76/3.D0*Hr1(-1)*upx**(-2)
     &     - 68.D0*Hr1(-1)*upx**(-2)*y
     &     - 76/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 8.D0*Hr1(-1)*upx**(-1)
     &     + 100/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 8.D0*Hr1(-1)*upx**(-1)*y**2
     &     + 24.D0*Hr1(-1)*HAr1(-1)
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 192.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 192.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 384.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 192.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 240.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 96.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 48.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     - 48.D0*Hr1(-1)*HBr1(-1)
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     - 128.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     + 128.D0*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     + 256.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     + 128.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     - 128.D0*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     - 224.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     + 64.D0*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     + 96.D0*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     + 12.D0*HAr1(-1)
     &     - 12.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 24.D0*HAr1(-1)*upx**(-2)
     &     - 12.D0*HAr1(-1)*upx**(-2)*y
     &     + 12.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 24.D0*HAr1(-1)*upx**(-1)
     &     + 12.D0*HAr1(-1)*upx**(-1)*y
     &     + 12.D0*HAr1(-1)*y**(-1)
     &     - 12.D0*HAr1(-1)*Hr1(0)
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 96.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 96.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 192.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 96.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 120.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 48.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 24.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 16.D0*HBr1(-1)
     &     + 8.D0*HBr1(-1)*upx**(-2)*z**(-1)
     &     - 8.D0*HBr1(-1)*upx**(-2)*y
     &     - 8.D0*HBr1(-1)*upx**(-1)*z**(-1)
     &     + 8.D0*HBr1(-1)*upx**(-1)*y
     &     - 8.D0*HBr1(-1)*z**(-1)
     &     + 24.D0*HBr1(-1)*Hr1(0)
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-4)
     &     + 64.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 64.D0*HBr1(-1)*Hr1(0)*upx**(-3)
     &     - 128.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 64.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 64.D0*HBr1(-1)*Hr1(0)*upx**(-2)
     &     + 112.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 32.D0*HBr1(-1)*Hr1(0)*upx**(-1)
     &     - 48.D0*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 48.D0*Hr2(-1,-1)
     &     + 16.D0*Hr2(-1,-1)*upx**(-4)
     &     + 32.D0*Hr2(-1,-1)*upx**(-4)*y
     &     + 16.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 32.D0*Hr2(-1,-1)*upx**(-3)
     &     - 64.D0*Hr2(-1,-1)*upx**(-3)*y
     &     - 32.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 80.D0*Hr2(-1,-1)*upx**(-2)
     &     + 128.D0*Hr2(-1,-1)*upx**(-2)*y
     &     + 16.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     - 64.D0*Hr2(-1,-1)*upx**(-1)
     &     - 96.D0*Hr2(-1,-1)*upx**(-1)*y
     &     + 24.D0*HAr2(-1,-1)
     &     + 48.D0*HAr2(-1,-1)*upx**(-2)*y
     &     - 48.D0*HAr2(-1,-1)*upx**(-1)*y
     &     + 16.D0*HBr2(-1,-1)
     &     + 64.D0*HBr2(-1,-1)*upx**(-2)
     &     + 32.D0*HBr2(-1,-1)*upx**(-2)*y
     &     - 64.D0*HBr2(-1,-1)*upx**(-1)
     &     - 32.D0*HBr2(-1,-1)*upx**(-1)*y
     &     - 24.D0*Hr2(-1,0)
     &     - 8.D0*Hr2(-1,0)*upx**(-4)
     &     - 16.D0*Hr2(-1,0)*upx**(-4)*y
     &     - 8.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     + 16.D0*Hr2(-1,0)*upx**(-3)
     &     + 32.D0*Hr2(-1,0)*upx**(-3)*y
     &     + 16.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     - 40.D0*Hr2(-1,0)*upx**(-2)
     &     - 64.D0*Hr2(-1,0)*upx**(-2)*y
     &     - 8.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     + 32.D0*Hr2(-1,0)*upx**(-1)
     &     + 48.D0*Hr2(-1,0)*upx**(-1)*y
     &     - 44/3.D0*Hr1(0)
     &     - 12.D0*Hr1(0)*umx**(-4)
     &     + 24.D0*Hr1(0)*umx**(-4)*y
     &     - 12.D0*Hr1(0)*umx**(-4)*y**2
     &     + 24.D0*Hr1(0)*umx**(-3)
     &     - 48.D0*Hr1(0)*umx**(-3)*y
     &     + 24.D0*Hr1(0)*umx**(-3)*y**2
     &     - 36.D0*Hr1(0)*umx**(-2)
     &     + 48.D0*Hr1(0)*umx**(-2)*y
     &     - 8.D0*Hr1(0)*umx**(-2)*y**2
     &     + 75/4.D0*Hr1(0)*umx**(-1)
     &     - 51/2.D0*Hr1(0)*umx**(-1)*y
     &     - 13/4.D0*Hr1(0)*umx**(-1)*y**2
     &     + 24.D0*Hr1(0)*upx**(-5)
     &     + 48.D0*Hr1(0)*upx**(-5)*y
     &     + 24.D0*Hr1(0)*upx**(-5)*y**2
     &     - 154/3.D0*Hr1(0)*upx**(-4)
     &     - 308/3.D0*Hr1(0)*upx**(-4)*y
     &     - 154/3.D0*Hr1(0)*upx**(-4)*y**2
     &     + 101/3.D0*Hr1(0)*upx**(-3)
     &     + 274/3.D0*Hr1(0)*upx**(-3)*y
     &     + 101/3.D0*Hr1(0)*upx**(-3)*y**2
     &     - 23/6.D0*Hr1(0)*upx**(-2)
     &     - 35.D0*Hr1(0)*upx**(-2)*y
     &     - 23/6.D0*Hr1(0)*upx**(-2)*y**2
     &     + 35/4.D0*Hr1(0)*upx**(-1)
     &     - 1/6.D0*Hr1(0)*upx**(-1)*y
     &     - 13/4.D0*Hr1(0)*upx**(-1)*y**2
     &     - 12.D0*Hr2(0,-1)
     &     + 24.D0*Hr2(0,-1)*umx**(-5)
     &     - 48.D0*Hr2(0,-1)*umx**(-5)*y
     &     + 24.D0*Hr2(0,-1)*umx**(-5)*y**2
     &     - 60.D0*Hr2(0,-1)*umx**(-4)
     &     + 120.D0*Hr2(0,-1)*umx**(-4)*y
     &     - 60.D0*Hr2(0,-1)*umx**(-4)*y**2
     &     + 94.D0*Hr2(0,-1)*umx**(-3)
     &     - 140.D0*Hr2(0,-1)*umx**(-3)*y
     &     + 38.D0*Hr2(0,-1)*umx**(-3)*y**2
     &     - 81.D0*Hr2(0,-1)*umx**(-2)
     &     + 90.D0*Hr2(0,-1)*umx**(-2)*y
     &     + 3.D0*Hr2(0,-1)*umx**(-2)*y**2
     &     + 73/2.D0*Hr2(0,-1)*umx**(-1)
     &     - 3.D0*Hr2(0,-1)*umx**(-1)*y
     &     + 1/2.D0*Hr2(0,-1)*umx**(-1)*y**2
     &     - 8.D0*Hr2(0,-1)*upx**(-4)
     &     - 16.D0*Hr2(0,-1)*upx**(-4)*y
     &     - 8.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     - 12.D0*Hr2(0,-1)*upx**(-3)
     &     + 16.D0*Hr2(0,-1)*upx**(-3)*y
     &     + 28.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     + 2.D0*Hr2(0,-1)*upx**(-2)
     &     - 40.D0*Hr2(0,-1)*upx**(-2)*y
     &     - 26.D0*Hr2(0,-1)*upx**(-2)*y**2
     &     - 39/2.D0*Hr2(0,-1)*upx**(-1)
     &     + 21.D0*Hr2(0,-1)*upx**(-1)*y
     &     + 1/2.D0*Hr2(0,-1)*upx**(-1)*y**2
     &     - 12.D0*HAr2(0,-1)
     &     - 24.D0*HAr2(0,-1)*upx**(-2)*y
     &     + 24.D0*HAr2(0,-1)*upx**(-1)*y
     &     - 8.D0*HBr2(0,-1)
     &     - 32.D0*HBr2(0,-1)*upx**(-2)
     &     - 16.D0*HBr2(0,-1)*upx**(-2)*y
     &     + 32.D0*HBr2(0,-1)*upx**(-1)
     &     + 16.D0*HBr2(0,-1)*upx**(-1)*y
     &     - 12.D0*Hr2(0,0)*umx**(-5)
     &     + 24.D0*Hr2(0,0)*umx**(-5)*y
     &     - 12.D0*Hr2(0,0)*umx**(-5)*y**2
     &     + 30.D0*Hr2(0,0)*umx**(-4)
     &     - 60.D0*Hr2(0,0)*umx**(-4)*y
     &     + 30.D0*Hr2(0,0)*umx**(-4)*y**2
     &     - 47.D0*Hr2(0,0)*umx**(-3)
     &     + 70.D0*Hr2(0,0)*umx**(-3)*y
     &     - 19.D0*Hr2(0,0)*umx**(-3)*y**2
     &     + 81/2.D0*Hr2(0,0)*umx**(-2)
     &     - 45.D0*Hr2(0,0)*umx**(-2)*y
     &     - 3/2.D0*Hr2(0,0)*umx**(-2)*y**2
     &     - 23/2.D0*Hr2(0,0)*umx**(-1)
     &     + 1/2.D0*Hr2(0,0)*umx**(-1)*y**2
     &     + 12.D0*Hr2(0,0)*upx**(-5)
     &     + 24.D0*Hr2(0,0)*upx**(-5)*y
     &     + 12.D0*Hr2(0,0)*upx**(-5)*y**2
     &     - 26.D0*Hr2(0,0)*upx**(-4)
     &     - 52.D0*Hr2(0,0)*upx**(-4)*y
     &     - 26.D0*Hr2(0,0)*upx**(-4)*y**2
     &     + 33.D0*Hr2(0,0)*upx**(-3)
     &     + 58.D0*Hr2(0,0)*upx**(-3)*y
     &     + 13.D0*Hr2(0,0)*upx**(-3)*y**2
     &     - 23/2.D0*Hr2(0,0)*upx**(-2)
     &     - 19.D0*Hr2(0,0)*upx**(-2)*y
     &     + 5/2.D0*Hr2(0,0)*upx**(-2)*y**2
     &     + 33/2.D0*Hr2(0,0)*upx**(-1)
     &     + 1/2.D0*Hr2(0,0)*upx**(-1)*y**2
     &     - 12.D0*Hr2(1,0)
     &     + 27/2.D0*Hr2(1,0)*umx**(-1)
     &     - 3.D0*Hr2(1,0)*umx**(-1)*y
     &     + 3/2.D0*Hr2(1,0)*umx**(-1)*y**2
     &     + 24.D0*Hr2(1,0)*upx**(-5)
     &     + 48.D0*Hr2(1,0)*upx**(-5)*y
     &     + 24.D0*Hr2(1,0)*upx**(-5)*y**2
     &     - 60.D0*Hr2(1,0)*upx**(-4)
     &     - 120.D0*Hr2(1,0)*upx**(-4)*y
     &     - 60.D0*Hr2(1,0)*upx**(-4)*y**2
     &     + 54.D0*Hr2(1,0)*upx**(-3)
     &     + 132.D0*Hr2(1,0)*upx**(-3)*y
     &     + 54.D0*Hr2(1,0)*upx**(-3)*y**2
     &     - 21.D0*Hr2(1,0)*upx**(-2)
     &     - 78.D0*Hr2(1,0)*upx**(-2)*y
     &     - 21.D0*Hr2(1,0)*upx**(-2)*y**2
     &     + 27/2.D0*Hr2(1,0)*upx**(-1)
     &     + 21.D0*Hr2(1,0)*upx**(-1)*y
     &     + 3/2.D0*Hr2(1,0)*upx**(-1)*y**2
     &     - 44/3.D0*DLog(tm**(-2)*w**2)
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 352/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 88.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc*Nh * (
     &     - 20.D0/3.D0
     &     + 32.D0*upx**(-6)
     &     + 64.D0*upx**(-6)*y
     &     + 32.D0*upx**(-6)*y**2
     &     - 96.D0*upx**(-5)
     &     - 192.D0*upx**(-5)*y
     &     - 96.D0*upx**(-5)*y**2
     &     + 248/3.D0*upx**(-4)
     &     + 592/3.D0*upx**(-4)*y
     &     + 248/3.D0*upx**(-4)*y**2
     &     - 16/3.D0*upx**(-3)
     &     - 224/3.D0*upx**(-3)*y
     &     - 16/3.D0*upx**(-3)*y**2
     &     + 8/3.D0*upx**(-2)
     &     - 8.D0*upx**(-2)*y
     &     - 40/3.D0*upx**(-2)*y**2
     &     - 16.D0*upx**(-1)
     &     + 40/3.D0*upx**(-1)*y
     &     + 4.D0*Hr1(0)
     &     + 32.D0*Hr1(0)*upx**(-7)
     &     + 64.D0*Hr1(0)*upx**(-7)*y
     &     + 32.D0*Hr1(0)*upx**(-7)*y**2
     &     - 112.D0*Hr1(0)*upx**(-6)
     &     - 224.D0*Hr1(0)*upx**(-6)*y
     &     - 112.D0*Hr1(0)*upx**(-6)*y**2
     &     + 128.D0*Hr1(0)*upx**(-5)
     &     + 288.D0*Hr1(0)*upx**(-5)*y
     &     + 128.D0*Hr1(0)*upx**(-5)*y**2
     &     - 40.D0*Hr1(0)*upx**(-4)
     &     - 160.D0*Hr1(0)*upx**(-4)*y
     &     - 40.D0*Hr1(0)*upx**(-4)*y**2
     &     + 16.D0*Hr1(0)*upx**(-3)*y
     &     - 16.D0*Hr1(0)*upx**(-3)*y**2
     &     - 16.D0*Hr1(0)*upx**(-2)
     &     + 24.D0*Hr1(0)*upx**(-2)*y
     &     + 8.D0*Hr1(0)*upx**(-2)*y**2
     &     - 8.D0*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc*Nl * (
     &     - 20.D0/3.D0
     &     - 40/3.D0*upx**(-4)
     &     - 80/3.D0*upx**(-4)*y
     &     - 40/3.D0*upx**(-4)*y**2
     &     + 80/3.D0*upx**(-3)
     &     + 160/3.D0*upx**(-3)*y
     &     + 80/3.D0*upx**(-3)*y**2
     &     - 40/3.D0*upx**(-2)
     &     - 40.D0*upx**(-2)*y
     &     - 40/3.D0*upx**(-2)*y**2
     &     + 40/3.D0*upx**(-1)*y
     &     + 8.D0*Hr1(-1)
     &     + 16.D0*Hr1(-1)*upx**(-4)
     &     + 32.D0*Hr1(-1)*upx**(-4)*y
     &     + 16.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 32.D0*Hr1(-1)*upx**(-3)
     &     - 64.D0*Hr1(-1)*upx**(-3)*y
     &     - 32.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 16.D0*Hr1(-1)*upx**(-2)
     &     + 48.D0*Hr1(-1)*upx**(-2)*y
     &     + 16.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 16.D0*Hr1(-1)*upx**(-1)*y
     &     - 4.D0*Hr1(0)
     &     - 8.D0*Hr1(0)*upx**(-4)
     &     - 16.D0*Hr1(0)*upx**(-4)*y
     &     - 8.D0*Hr1(0)*upx**(-4)*y**2
     &     + 16.D0*Hr1(0)*upx**(-3)
     &     + 32.D0*Hr1(0)*upx**(-3)*y
     &     + 16.D0*Hr1(0)*upx**(-3)*y**2
     &     - 8.D0*Hr1(0)*upx**(-2)
     &     - 24.D0*Hr1(0)*upx**(-2)*y
     &     - 8.D0*Hr1(0)*upx**(-2)*y**2
     &     + 8.D0*Hr1(0)*upx**(-1)*y
     &     - 4.D0*DLog(tm**(-2)*w**2)
     &     - 8.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 16.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 8.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 16.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 32.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 16.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 8.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 24.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 8.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 8.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**2 * (
     &     + 86.D0/3.D0
     &     + 18.D0*umx**(-5)*z2
     &     - 36.D0*umx**(-5)*y*z2
     &     + 18.D0*umx**(-5)*y**2*z2
     &     - 45.D0*umx**(-4)*z2
     &     + 90.D0*umx**(-4)*y*z2
     &     - 45.D0*umx**(-4)*y**2*z2
     &     + 101/2.D0*umx**(-3)*z2
     &     - 85.D0*umx**(-3)*y*z2
     &     + 57/2.D0*umx**(-3)*y**2*z2
     &     + 3.D0*umx**(-2)
     &     - 123/4.D0*umx**(-2)*z2
     &     - 6.D0*umx**(-2)*y
     &     + 75/2.D0*umx**(-2)*y*z2
     &     + 3.D0*umx**(-2)*y**2
     &     + 9/4.D0*umx**(-2)*y**2*z2
     &     - 3.D0*umx**(-1)
     &     + 211/8.D0*umx**(-1)*z2
     &     + 6.D0*umx**(-1)*y
     &     - 17/4.D0*umx**(-1)*y*z2
     &     - 3.D0*umx**(-1)*y**2
     &     + 11/8.D0*umx**(-1)*y**2*z2
     &     + 16.D0*upx**(-5)*z2
     &     + 32.D0*upx**(-5)*y*z2
     &     + 16.D0*upx**(-5)*y**2*z2
     &     + 16/3.D0*upx**(-4)
     &     - 12.D0*upx**(-4)*z2
     &     + 32/3.D0*upx**(-4)*y
     &     - 24.D0*upx**(-4)*y*z2
     &     + 16/3.D0*upx**(-4)*y**2
     &     - 12.D0*upx**(-4)*y**2*z2
     &     - 32/3.D0*upx**(-3)
     &     - 21.D0*upx**(-3)*z2
     &     - 64/3.D0*upx**(-3)*y
     &     - 16.D0*upx**(-3)*y*z2
     &     - 32/3.D0*upx**(-3)*y**2
     &     - 11.D0*upx**(-3)*y**2*z2
     &     + 25/3.D0*upx**(-2)
     &     + 23/2.D0*upx**(-2)*z2
     &     + 22.D0*upx**(-2)*y
     &     + 10.D0*upx**(-2)*y*z2
     &     + 25/3.D0*upx**(-2)*y**2
     &     + 1/2.D0*upx**(-2)*y**2*z2
     &     - 3.D0*upx**(-1)
     &     + 35/8.D0*upx**(-1)*z2
     &     - 34/3.D0*upx**(-1)*y
     &     - 17/4.D0*upx**(-1)*y*z2
     &     - 3.D0*upx**(-1)*y**2
     &     + 11/8.D0*upx**(-1)*y**2*z2
     &     - 22.D0*Hr1(-1)
     &     - 36.D0*Hr1(-1)*umx**(-4)
     &     + 72.D0*Hr1(-1)*umx**(-4)*y
     &     - 36.D0*Hr1(-1)*umx**(-4)*y**2
     &     + 72.D0*Hr1(-1)*umx**(-3)
     &     - 144.D0*Hr1(-1)*umx**(-3)*y
     &     + 72.D0*Hr1(-1)*umx**(-3)*y**2
     &     - 68.D0*Hr1(-1)*umx**(-2)
     &     + 104.D0*Hr1(-1)*umx**(-2)*y
     &     - 24.D0*Hr1(-1)*umx**(-2)*y**2
     &     + 32.D0*Hr1(-1)*umx**(-1)
     &     - 32.D0*Hr1(-1)*umx**(-1)*y
     &     - 12.D0*Hr1(-1)*umx**(-1)*y**2
     &     - 16.D0*Hr1(-1)*upx**(-4)
     &     - 32.D0*Hr1(-1)*upx**(-4)*y
     &     - 16.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 32.D0*Hr1(-1)*upx**(-3)
     &     + 64.D0*Hr1(-1)*upx**(-3)*y
     &     + 32.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 4.D0*Hr1(-1)*upx**(-2)
     &     - 24.D0*Hr1(-1)*upx**(-2)*y
     &     - 4.D0*Hr1(-1)*upx**(-2)*y**2
     &     - 12.D0*Hr1(-1)*upx**(-1)
     &     - 8.D0*Hr1(-1)*upx**(-1)*y
     &     - 12.D0*Hr1(-1)*upx**(-1)*y**2
     &     - 16.D0*Hr1(-1)*HAr1(-1)
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 128.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 128.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 256.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 128.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 160.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 32.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 12.D0*Hr1(-1)*HBr1(-1)
     &     + 16.D0*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 16.D0*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 64.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 32.D0*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 32.D0*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 56.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 16.D0*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 16.D0*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 24.D0*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 8.D0*HAr1(-1)
     &     + 8.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 16.D0*HAr1(-1)*upx**(-2)
     &     + 8.D0*HAr1(-1)*upx**(-2)*y
     &     - 8.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 16.D0*HAr1(-1)*upx**(-1)
     &     - 8.D0*HAr1(-1)*upx**(-1)*y
     &     - 8.D0*HAr1(-1)*y**(-1)
     &     + 8.D0*HAr1(-1)*Hr1(0)
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 64.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 64.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 128.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 64.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 80.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 16.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 4.D0*HBr1(-1)
     &     - 2.D0*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 2.D0*HBr1(-1)*upx**(-2)*y
     &     + 2.D0*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 2.D0*HBr1(-1)*upx**(-1)*y
     &     + 2.D0*HBr1(-1)*z**(-1)
     &     - 6.D0*HBr1(-1)*Hr1(0)
     &     - 8.D0*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8.D0*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 32.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16.D0*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 16.D0*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 28.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8.D0*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 8.D0*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 12.D0*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 4.D0*Hr2(-1,-1)
     &     + 32.D0*Hr2(-1,-1)*upx**(-4)
     &     + 64.D0*Hr2(-1,-1)*upx**(-4)*y
     &     + 32.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 64.D0*Hr2(-1,-1)*upx**(-3)
     &     - 128.D0*Hr2(-1,-1)*upx**(-3)*y
     &     - 64.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 16.D0*Hr2(-1,-1)*upx**(-2)
     &     + 56.D0*Hr2(-1,-1)*upx**(-2)*y
     &     + 32.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 16.D0*Hr2(-1,-1)*upx**(-1)
     &     + 8.D0*Hr2(-1,-1)*upx**(-1)*y
     &     - 16.D0*HAr2(-1,-1)
     &     - 32.D0*HAr2(-1,-1)*upx**(-2)*y
     &     + 32.D0*HAr2(-1,-1)*upx**(-1)*y
     &     - 4.D0*HBr2(-1,-1)
     &     - 16.D0*HBr2(-1,-1)*upx**(-2)
     &     - 8.D0*HBr2(-1,-1)*upx**(-2)*y
     &     + 16.D0*HBr2(-1,-1)*upx**(-1)
     &     + 8.D0*HBr2(-1,-1)*upx**(-1)*y
     &     + 2.D0*Hr2(-1,0)
     &     - 16.D0*Hr2(-1,0)*upx**(-4)
     &     - 32.D0*Hr2(-1,0)*upx**(-4)*y
     &     - 16.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     + 32.D0*Hr2(-1,0)*upx**(-3)
     &     + 64.D0*Hr2(-1,0)*upx**(-3)*y
     &     + 32.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     - 8.D0*Hr2(-1,0)*upx**(-2)
     &     - 28.D0*Hr2(-1,0)*upx**(-2)*y
     &     - 16.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     - 8.D0*Hr2(-1,0)*upx**(-1)
     &     - 4.D0*Hr2(-1,0)*upx**(-1)*y
     &     + 12.D0*Hr1(0)
     &     + 18.D0*Hr1(0)*umx**(-4)
     &     - 36.D0*Hr1(0)*umx**(-4)*y
     &     + 18.D0*Hr1(0)*umx**(-4)*y**2
     &     - 36.D0*Hr1(0)*umx**(-3)
     &     + 72.D0*Hr1(0)*umx**(-3)*y
     &     - 36.D0*Hr1(0)*umx**(-3)*y**2
     &     + 34.D0*Hr1(0)*umx**(-2)
     &     - 52.D0*Hr1(0)*umx**(-2)*y
     &     + 12.D0*Hr1(0)*umx**(-2)*y**2
     &     - 57/4.D0*Hr1(0)*umx**(-1)
     &     + 33/2.D0*Hr1(0)*umx**(-1)*y
     &     + 23/4.D0*Hr1(0)*umx**(-1)*y**2
     &     - 8.D0*Hr1(0)*upx**(-5)
     &     - 16.D0*Hr1(0)*upx**(-5)*y
     &     - 8.D0*Hr1(0)*upx**(-5)*y**2
     &     + 28.D0*Hr1(0)*upx**(-4)
     &     + 56.D0*Hr1(0)*upx**(-4)*y
     &     + 28.D0*Hr1(0)*upx**(-4)*y**2
     &     - 33.D0*Hr1(0)*upx**(-3)
     &     - 74.D0*Hr1(0)*upx**(-3)*y
     &     - 33.D0*Hr1(0)*upx**(-3)*y**2
     &     + 15/2.D0*Hr1(0)*upx**(-2)
     &     + 35.D0*Hr1(0)*upx**(-2)*y
     &     + 15/2.D0*Hr1(0)*upx**(-2)*y**2
     &     + 7/4.D0*Hr1(0)*upx**(-1)
     &     - 3/2.D0*Hr1(0)*upx**(-1)*y
     &     + 23/4.D0*Hr1(0)*upx**(-1)*y**2
     &     + 4.D0*Hr2(0,-1)
     &     - 36.D0*Hr2(0,-1)*umx**(-5)
     &     + 72.D0*Hr2(0,-1)*umx**(-5)*y
     &     - 36.D0*Hr2(0,-1)*umx**(-5)*y**2
     &     + 90.D0*Hr2(0,-1)*umx**(-4)
     &     - 180.D0*Hr2(0,-1)*umx**(-4)*y
     &     + 90.D0*Hr2(0,-1)*umx**(-4)*y**2
     &     - 101.D0*Hr2(0,-1)*umx**(-3)
     &     + 170.D0*Hr2(0,-1)*umx**(-3)*y
     &     - 57.D0*Hr2(0,-1)*umx**(-3)*y**2
     &     + 123/2.D0*Hr2(0,-1)*umx**(-2)
     &     - 75.D0*Hr2(0,-1)*umx**(-2)*y
     &     - 9/2.D0*Hr2(0,-1)*umx**(-2)*y**2
     &     - 139/4.D0*Hr2(0,-1)*umx**(-1)
     &     + 9/2.D0*Hr2(0,-1)*umx**(-1)*y
     &     - 3/4.D0*Hr2(0,-1)*umx**(-1)*y**2
     &     - 16.D0*Hr2(0,-1)*upx**(-4)
     &     - 32.D0*Hr2(0,-1)*upx**(-4)*y
     &     - 16.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     + 34.D0*Hr2(0,-1)*upx**(-3)
     &     + 48.D0*Hr2(0,-1)*upx**(-3)*y
     &     + 14.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     - 11.D0*Hr2(0,-1)*upx**(-2)
     &     - 4.D0*Hr2(0,-1)*upx**(-2)*y
     &     + 11.D0*Hr2(0,-1)*upx**(-2)*y**2
     &     + 37/4.D0*Hr2(0,-1)*upx**(-1)
     &     - 7/2.D0*Hr2(0,-1)*upx**(-1)*y
     &     - 3/4.D0*Hr2(0,-1)*upx**(-1)*y**2
     &     + 8.D0*HAr2(0,-1)
     &     + 16.D0*HAr2(0,-1)*upx**(-2)*y
     &     - 16.D0*HAr2(0,-1)*upx**(-1)*y
     &     + 2.D0*HBr2(0,-1)
     &     + 8.D0*HBr2(0,-1)*upx**(-2)
     &     + 4.D0*HBr2(0,-1)*upx**(-2)*y
     &     - 8.D0*HBr2(0,-1)*upx**(-1)
     &     - 4.D0*HBr2(0,-1)*upx**(-1)*y
     &     + 18.D0*Hr2(0,0)*umx**(-5)
     &     - 36.D0*Hr2(0,0)*umx**(-5)*y
     &     + 18.D0*Hr2(0,0)*umx**(-5)*y**2
     &     - 45.D0*Hr2(0,0)*umx**(-4)
     &     + 90.D0*Hr2(0,0)*umx**(-4)*y
     &     - 45.D0*Hr2(0,0)*umx**(-4)*y**2
     &     + 101/2.D0*Hr2(0,0)*umx**(-3)
     &     - 85.D0*Hr2(0,0)*umx**(-3)*y
     &     + 57/2.D0*Hr2(0,0)*umx**(-3)*y**2
     &     - 123/4.D0*Hr2(0,0)*umx**(-2)
     &     + 75/2.D0*Hr2(0,0)*umx**(-2)*y
     &     + 9/4.D0*Hr2(0,0)*umx**(-2)*y**2
     &     + 121/8.D0*Hr2(0,0)*umx**(-1)
     &     - 7/4.D0*Hr2(0,0)*umx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*umx**(-1)*y**2
     &     - 4.D0*Hr2(0,0)*upx**(-5)
     &     - 8.D0*Hr2(0,0)*upx**(-5)*y
     &     - 4.D0*Hr2(0,0)*upx**(-5)*y**2
     &     + 18.D0*Hr2(0,0)*upx**(-4)
     &     + 36.D0*Hr2(0,0)*upx**(-4)*y
     &     + 18.D0*Hr2(0,0)*upx**(-4)*y**2
     &     - 26.D0*Hr2(0,0)*upx**(-3)
     &     - 46.D0*Hr2(0,0)*upx**(-3)*y
     &     - 16.D0*Hr2(0,0)*upx**(-3)*y**2
     &     + 9.D0*Hr2(0,0)*upx**(-2)
     &     + 15.D0*Hr2(0,0)*upx**(-2)*y
     &     - 2.D0*Hr2(0,0)*upx**(-2)*y**2
     &     - 55/8.D0*Hr2(0,0)*upx**(-1)
     &     - 7/4.D0*Hr2(0,0)*upx**(-1)*y
     &     + 1/8.D0*Hr2(0,0)*upx**(-1)*y**2
     &     + 4.D0*Hr2(1,0)
     &     - 9/2.D0*Hr2(1,0)*umx**(-1)
     &     + Hr2(1,0)*umx**(-1)*y
     &     - 1/2.D0*Hr2(1,0)*umx**(-1)*y**2
     &     - 8.D0*Hr2(1,0)*upx**(-5)
     &     - 16.D0*Hr2(1,0)*upx**(-5)*y
     &     - 8.D0*Hr2(1,0)*upx**(-5)*y**2
     &     + 20.D0*Hr2(1,0)*upx**(-4)
     &     + 40.D0*Hr2(1,0)*upx**(-4)*y
     &     + 20.D0*Hr2(1,0)*upx**(-4)*y**2
     &     - 18.D0*Hr2(1,0)*upx**(-3)
     &     - 44.D0*Hr2(1,0)*upx**(-3)*y
     &     - 18.D0*Hr2(1,0)*upx**(-3)*y**2
     &     + 7.D0*Hr2(1,0)*upx**(-2)
     &     + 26.D0*Hr2(1,0)*upx**(-2)*y
     &     + 7.D0*Hr2(1,0)*upx**(-2)*y**2
     &     - 9/2.D0*Hr2(1,0)*upx**(-1)
     &     - 7.D0*Hr2(1,0)*upx**(-1)*y
     &     - 1/2.D0*Hr2(1,0)*upx**(-1)*y**2
     &     + 22.D0*DLog(tm**(-2)*w**2)
     &     + 44.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 88.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 44.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 88.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 176.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 88.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 44.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 132.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 44.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 44.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**3*Nh * (
     &     + 20.D0/9.D0
     &     - 32/3.D0*upx**(-6)
     &     - 64/3.D0*upx**(-6)*y
     &     - 32/3.D0*upx**(-6)*y**2
     &     + 32.D0*upx**(-5)
     &     + 64.D0*upx**(-5)*y
     &     + 32.D0*upx**(-5)*y**2
     &     - 248/9.D0*upx**(-4)
     &     - 592/9.D0*upx**(-4)*y
     &     - 248/9.D0*upx**(-4)*y**2
     &     + 16/9.D0*upx**(-3)
     &     + 224/9.D0*upx**(-3)*y
     &     + 16/9.D0*upx**(-3)*y**2
     &     - 8/9.D0*upx**(-2)
     &     + 8/3.D0*upx**(-2)*y
     &     + 40/9.D0*upx**(-2)*y**2
     &     + 16/3.D0*upx**(-1)
     &     - 40/9.D0*upx**(-1)*y
     &     - 4/3.D0*Hr1(0)
     &     - 32/3.D0*Hr1(0)*upx**(-7)
     &     - 64/3.D0*Hr1(0)*upx**(-7)*y
     &     - 32/3.D0*Hr1(0)*upx**(-7)*y**2
     &     + 112/3.D0*Hr1(0)*upx**(-6)
     &     + 224/3.D0*Hr1(0)*upx**(-6)*y
     &     + 112/3.D0*Hr1(0)*upx**(-6)*y**2
     &     - 128/3.D0*Hr1(0)*upx**(-5)
     &     - 96.D0*Hr1(0)*upx**(-5)*y
     &     - 128/3.D0*Hr1(0)*upx**(-5)*y**2
     &     + 40/3.D0*Hr1(0)*upx**(-4)
     &     + 160/3.D0*Hr1(0)*upx**(-4)*y
     &     + 40/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-3)*y
     &     + 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 16/3.D0*Hr1(0)*upx**(-2)
     &     - 8.D0*Hr1(0)*upx**(-2)*y
     &     - 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**3*Nl * (
     &     + 20.D0/9.D0
     &     + 40/9.D0*upx**(-4)
     &     + 80/9.D0*upx**(-4)*y
     &     + 40/9.D0*upx**(-4)*y**2
     &     - 80/9.D0*upx**(-3)
     &     - 160/9.D0*upx**(-3)*y
     &     - 80/9.D0*upx**(-3)*y**2
     &     + 40/9.D0*upx**(-2)
     &     + 40/3.D0*upx**(-2)*y
     &     + 40/9.D0*upx**(-2)*y**2
     &     - 40/9.D0*upx**(-1)*y
     &     - 8/3.D0*Hr1(-1)
     &     - 16/3.D0*Hr1(-1)*upx**(-4)
     &     - 32/3.D0*Hr1(-1)*upx**(-4)*y
     &     - 16/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     + 32/3.D0*Hr1(-1)*upx**(-3)
     &     + 64/3.D0*Hr1(-1)*upx**(-3)*y
     &     + 32/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     - 16/3.D0*Hr1(-1)*upx**(-2)
     &     - 16.D0*Hr1(-1)*upx**(-2)*y
     &     - 16/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 16/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 4/3.D0*Hr1(0)
     &     + 8/3.D0*Hr1(0)*upx**(-4)
     &     + 16/3.D0*Hr1(0)*upx**(-4)*y
     &     + 8/3.D0*Hr1(0)*upx**(-4)*y**2
     &     - 16/3.D0*Hr1(0)*upx**(-3)
     &     - 32/3.D0*Hr1(0)*upx**(-3)*y
     &     - 16/3.D0*Hr1(0)*upx**(-3)*y**2
     &     + 8/3.D0*Hr1(0)*upx**(-2)
     &     + 8.D0*Hr1(0)*upx**(-2)*y
     &     + 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 8/3.D0*Hr1(0)*upx**(-1)*y
     &     + 4/3.D0*DLog(tm**(-2)*w**2)
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 32/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 16/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 8.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 8/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qq = T13PoleEP0qq + Nc**4 * (
     &     - 107.D0/9.D0
     &     - 6.D0*umx**(-5)*z2
     &     + 12.D0*umx**(-5)*y*z2
     &     - 6.D0*umx**(-5)*y**2*z2
     &     + 15.D0*umx**(-4)*z2
     &     - 30.D0*umx**(-4)*y*z2
     &     + 15.D0*umx**(-4)*y**2*z2
     &     - 31/2.D0*umx**(-3)*z2
     &     + 27.D0*umx**(-3)*y*z2
     &     - 19/2.D0*umx**(-3)*y**2*z2
     &     - umx**(-2)
     &     + 33/4.D0*umx**(-2)*z2
     &     + 2.D0*umx**(-2)*y
     &     - 21/2.D0*umx**(-2)*y*z2
     &     - umx**(-2)*y**2
     &     - 3/4.D0*umx**(-2)*y**2*z2
     &     + umx**(-1)
     &     - 41/8.D0*umx**(-1)*z2
     &     - 2.D0*umx**(-1)*y
     &     + 3/4.D0*umx**(-1)*y*z2
     &     + umx**(-1)*y**2
     &     - 1/8.D0*umx**(-1)*y**2*z2
     &     - 88/9.D0*upx**(-4)
     &     - 4.D0*upx**(-4)*z2
     &     - 176/9.D0*upx**(-4)*y
     &     - 8.D0*upx**(-4)*y*z2
     &     - 88/9.D0*upx**(-4)*y**2
     &     - 4.D0*upx**(-4)*y**2*z2
     &     + 176/9.D0*upx**(-3)
     &     + 7.D0*upx**(-3)*z2
     &     + 352/9.D0*upx**(-3)*y
     &     + 12.D0*upx**(-3)*y*z2
     &     + 176/9.D0*upx**(-3)*y**2
     &     + 5.D0*upx**(-3)*y**2*z2
     &     - 97/9.D0*upx**(-2)
     &     - 5/2.D0*upx**(-2)*z2
     &     - 94/3.D0*upx**(-2)*y
     &     - 4.D0*upx**(-2)*y*z2
     &     - 97/9.D0*upx**(-2)*y**2
     &     + 1/2.D0*upx**(-2)*y**2*z2
     &     + upx**(-1)
     &     + 7/8.D0*upx**(-1)*z2
     &     + 106/9.D0*upx**(-1)*y
     &     + 3/4.D0*upx**(-1)*y*z2
     &     + upx**(-1)*y**2
     &     - 1/8.D0*upx**(-1)*y**2*z2
     &     + 20/3.D0*Hr1(-1)
     &     + 12.D0*Hr1(-1)*umx**(-4)
     &     - 24.D0*Hr1(-1)*umx**(-4)*y
     &     + 12.D0*Hr1(-1)*umx**(-4)*y**2
     &     - 24.D0*Hr1(-1)*umx**(-3)
     &     + 48.D0*Hr1(-1)*umx**(-3)*y
     &     - 24.D0*Hr1(-1)*umx**(-3)*y**2
     &     + 20.D0*Hr1(-1)*umx**(-2)
     &     - 32.D0*Hr1(-1)*umx**(-2)*y
     &     + 8.D0*Hr1(-1)*umx**(-2)*y**2
     &     - 8.D0*Hr1(-1)*umx**(-1)
     &     + 8.D0*Hr1(-1)*umx**(-1)*y
     &     + 4.D0*Hr1(-1)*umx**(-1)*y**2
     &     + 28/3.D0*Hr1(-1)*upx**(-4)
     &     + 56/3.D0*Hr1(-1)*upx**(-4)*y
     &     + 28/3.D0*Hr1(-1)*upx**(-4)*y**2
     &     - 56/3.D0*Hr1(-1)*upx**(-3)
     &     - 112/3.D0*Hr1(-1)*upx**(-3)*y
     &     - 56/3.D0*Hr1(-1)*upx**(-3)*y**2
     &     + 16/3.D0*Hr1(-1)*upx**(-2)
     &     + 20.D0*Hr1(-1)*upx**(-2)*y
     &     + 16/3.D0*Hr1(-1)*upx**(-2)*y**2
     &     + 4.D0*Hr1(-1)*upx**(-1)
     &     - 4/3.D0*Hr1(-1)*upx**(-1)*y
     &     + 4.D0*Hr1(-1)*upx**(-1)*y**2
     &     + 4.D0*Hr1(-1)*HAr1(-1)
     &     + 16.D0*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 32.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 16.D0*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 64.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 32.D0*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 16.D0*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 40.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 16.D0*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 8.D0*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 2.D0*HAr1(-1)
     &     - 2.D0*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 4.D0*HAr1(-1)*upx**(-2)
     &     - 2.D0*HAr1(-1)*upx**(-2)*y
     &     + 2.D0*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 4.D0*HAr1(-1)*upx**(-1)
     &     + 2.D0*HAr1(-1)*upx**(-1)*y
     &     + 2.D0*HAr1(-1)*y**(-1)
     &     - 2.D0*HAr1(-1)*Hr1(0)
     &     - 8.D0*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 16.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8.D0*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 32.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16.D0*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 8.D0*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 20.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8.D0*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 4.D0*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 4.D0*Hr2(-1,-1)
     &     - 16.D0*Hr2(-1,-1)*upx**(-4)
     &     - 32.D0*Hr2(-1,-1)*upx**(-4)*y
     &     - 16.D0*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 32.D0*Hr2(-1,-1)*upx**(-3)
     &     + 64.D0*Hr2(-1,-1)*upx**(-3)*y
     &     + 32.D0*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 16.D0*Hr2(-1,-1)*upx**(-2)
     &     - 40.D0*Hr2(-1,-1)*upx**(-2)*y
     &     - 16.D0*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 8.D0*Hr2(-1,-1)*upx**(-1)*y
     &     + 4.D0*HAr2(-1,-1)
     &     + 8.D0*HAr2(-1,-1)*upx**(-2)*y
     &     - 8.D0*HAr2(-1,-1)*upx**(-1)*y
     &     + 2.D0*Hr2(-1,0)
     &     + 8.D0*Hr2(-1,0)*upx**(-4)
     &     + 16.D0*Hr2(-1,0)*upx**(-4)*y
     &     + 8.D0*Hr2(-1,0)*upx**(-4)*y**2
     &     - 16.D0*Hr2(-1,0)*upx**(-3)
     &     - 32.D0*Hr2(-1,0)*upx**(-3)*y
     &     - 16.D0*Hr2(-1,0)*upx**(-3)*y**2
     &     + 8.D0*Hr2(-1,0)*upx**(-2)
     &     + 20.D0*Hr2(-1,0)*upx**(-2)*y
     &     + 8.D0*Hr2(-1,0)*upx**(-2)*y**2
     &     - 4.D0*Hr2(-1,0)*upx**(-1)*y
     &     - 10/3.D0*Hr1(0)
     &     - 6.D0*Hr1(0)*umx**(-4)
     &     + 12.D0*Hr1(0)*umx**(-4)*y
     &     - 6.D0*Hr1(0)*umx**(-4)*y**2
     &     + 12.D0*Hr1(0)*umx**(-3)
     &     - 24.D0*Hr1(0)*umx**(-3)*y
     &     + 12.D0*Hr1(0)*umx**(-3)*y**2
     &     - 10.D0*Hr1(0)*umx**(-2)
     &     + 16.D0*Hr1(0)*umx**(-2)*y
     &     - 4.D0*Hr1(0)*umx**(-2)*y**2
     &     + 4.D0*Hr1(0)*umx**(-1)
     &     - 4.D0*Hr1(0)*umx**(-1)*y
     &     - 2.D0*Hr1(0)*umx**(-1)*y**2
     &     - 14/3.D0*Hr1(0)*upx**(-4)
     &     - 28/3.D0*Hr1(0)*upx**(-4)*y
     &     - 14/3.D0*Hr1(0)*upx**(-4)*y**2
     &     + 28/3.D0*Hr1(0)*upx**(-3)
     &     + 56/3.D0*Hr1(0)*upx**(-3)*y
     &     + 28/3.D0*Hr1(0)*upx**(-3)*y**2
     &     - 8/3.D0*Hr1(0)*upx**(-2)
     &     - 10.D0*Hr1(0)*upx**(-2)*y
     &     - 8/3.D0*Hr1(0)*upx**(-2)*y**2
     &     - 2.D0*Hr1(0)*upx**(-1)
     &     + 2/3.D0*Hr1(0)*upx**(-1)*y
     &     - 2.D0*Hr1(0)*upx**(-1)*y**2
     &     + 12.D0*Hr2(0,-1)*umx**(-5)
     &     - 24.D0*Hr2(0,-1)*umx**(-5)*y
     &     + 12.D0*Hr2(0,-1)*umx**(-5)*y**2
     &     - 30.D0*Hr2(0,-1)*umx**(-4)
     &     + 60.D0*Hr2(0,-1)*umx**(-4)*y
     &     - 30.D0*Hr2(0,-1)*umx**(-4)*y**2
     &     + 31.D0*Hr2(0,-1)*umx**(-3)
     &     - 54.D0*Hr2(0,-1)*umx**(-3)*y
     &     + 19.D0*Hr2(0,-1)*umx**(-3)*y**2
     &     - 33/2.D0*Hr2(0,-1)*umx**(-2)
     &     + 21.D0*Hr2(0,-1)*umx**(-2)*y
     &     + 3/2.D0*Hr2(0,-1)*umx**(-2)*y**2
     &     + 41/4.D0*Hr2(0,-1)*umx**(-1)
     &     - 3/2.D0*Hr2(0,-1)*umx**(-1)*y
     &     + 1/4.D0*Hr2(0,-1)*umx**(-1)*y**2
     &     + 8.D0*Hr2(0,-1)*upx**(-4)
     &     + 16.D0*Hr2(0,-1)*upx**(-4)*y
     &     + 8.D0*Hr2(0,-1)*upx**(-4)*y**2
     &     - 14.D0*Hr2(0,-1)*upx**(-3)
     &     - 24.D0*Hr2(0,-1)*upx**(-3)*y
     &     - 10.D0*Hr2(0,-1)*upx**(-3)*y**2
     &     + 5.D0*Hr2(0,-1)*upx**(-2)
     &     + 8.D0*Hr2(0,-1)*upx**(-2)*y
     &     - Hr2(0,-1)*upx**(-2)*y**2
     &     - 7/4.D0*Hr2(0,-1)*upx**(-1)
     &     - 3/2.D0*Hr2(0,-1)*upx**(-1)*y
     &     + 1/4.D0*Hr2(0,-1)*upx**(-1)*y**2
     &     - 2.D0*HAr2(0,-1)
     &     - 4.D0*HAr2(0,-1)*upx**(-2)*y
     &     + 4.D0*HAr2(0,-1)*upx**(-1)*y
     &     - 6.D0*Hr2(0,0)*umx**(-5)
     &     + 12.D0*Hr2(0,0)*umx**(-5)*y
     &     - 6.D0*Hr2(0,0)*umx**(-5)*y**2
     &     + 15.D0*Hr2(0,0)*umx**(-4)
     &     - 30.D0*Hr2(0,0)*umx**(-4)*y
     &     + 15.D0*Hr2(0,0)*umx**(-4)*y**2
     &     - 31/2.D0*Hr2(0,0)*umx**(-3)
     &     + 27.D0*Hr2(0,0)*umx**(-3)*y
     &     - 19/2.D0*Hr2(0,0)*umx**(-3)*y**2
     &     + 33/4.D0*Hr2(0,0)*umx**(-2)
     &     - 21/2.D0*Hr2(0,0)*umx**(-2)*y
     &     - 3/4.D0*Hr2(0,0)*umx**(-2)*y**2
     &     - 41/8.D0*Hr2(0,0)*umx**(-1)
     &     + 3/4.D0*Hr2(0,0)*umx**(-1)*y
     &     - 1/8.D0*Hr2(0,0)*umx**(-1)*y**2
     &     - 4.D0*Hr2(0,0)*upx**(-4)
     &     - 8.D0*Hr2(0,0)*upx**(-4)*y
     &     - 4.D0*Hr2(0,0)*upx**(-4)*y**2
     &     + 7.D0*Hr2(0,0)*upx**(-3)
     &     + 12.D0*Hr2(0,0)*upx**(-3)*y
     &     + 5.D0*Hr2(0,0)*upx**(-3)*y**2
     &     - 5/2.D0*Hr2(0,0)*upx**(-2)
     &     - 4.D0*Hr2(0,0)*upx**(-2)*y
     &     + 1/2.D0*Hr2(0,0)*upx**(-2)*y**2
     &     + 7/8.D0*Hr2(0,0)*upx**(-1)
     &     + 3/4.D0*Hr2(0,0)*upx**(-1)*y
     &     - 1/8.D0*Hr2(0,0)*upx**(-1)*y**2
     &     - 22/3.D0*DLog(tm**(-2)*w**2)
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 176/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 88/3.D0*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 44.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 44/3.D0*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      return
      end

      double precision function T13PoleEP0qqc(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1d0

* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T13PoleEP0qqc =  + Nc**(-2) * (
     &     - 14
     &     + 12*umx**(-3)*z2
     &     - 12*umx**(-3)*y*z2
     &     - 18*umx**(-2)*z2
     &     + 18*umx**(-2)*y*z2
     &     + 24*umx**(-1)*z2
     &     - 4*umx**(-1)*y*z2
     &     + 2*umx**(-1)*y**2*z2
     &     + 32*upx**(-5)*z2
     &     + 64*upx**(-5)*y*z2
     &     + 32*upx**(-5)*y**2*z2
     &     - 48*upx**(-4)
     &     - 48*upx**(-4)*z2
     &     - 96*upx**(-4)*y
     &     - 96*upx**(-4)*y*z2
     &     - 48*upx**(-4)*y**2
     &     - 48*upx**(-4)*y**2*z2
     &     + 96*upx**(-3)
     &     - 4*upx**(-3)*z2
     &     + 192*upx**(-3)*y
     &     + 36*upx**(-3)*y*z2
     &     + 96*upx**(-3)*y**2
     &     + 8*upx**(-3)*y**2*z2
     &     - 48*upx**(-2)
     &     + 10*upx**(-2)*z2
     &     - 144*upx**(-2)*y
     &     - 2*upx**(-2)*y*z2
     &     - 48*upx**(-2)*y**2
     &     + 4*upx**(-2)*y**2*z2
     &     + 12*upx**(-1)*z2
     &     + 48*upx**(-1)*y
     &     - 4*upx**(-1)*y*z2
     &     + 2*upx**(-1)*y**2*z2
     &     - 8*Hr1(-1)
     &     - 24*Hr1(-1)*umx**(-2)
     &     + 24*Hr1(-1)*umx**(-2)*y
     &     + 24*Hr1(-1)*umx**(-1)
     &     - 24*Hr1(-1)*umx**(-1)*y
     &     + 24*Hr1(-1)*upx**(-4)
     &     + 48*Hr1(-1)*upx**(-4)*y
     &     + 24*Hr1(-1)*upx**(-4)*y**2
     &     - 48*Hr1(-1)*upx**(-3)
     &     - 96*Hr1(-1)*upx**(-3)*y
     &     - 48*Hr1(-1)*upx**(-3)*y**2
     &     + 24*Hr1(-1)*upx**(-2)
     &     + 72*Hr1(-1)*upx**(-2)*y
     &     + 24*Hr1(-1)*upx**(-2)*y**2
     &     - 24*Hr1(-1)*upx**(-1)*y
     &     - 12*Hr1(-1)*HAr1(-1)
     &     - 48*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 96*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 48*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 192*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 48*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 120*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 48*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 24*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 36*Hr1(-1)*HBr1(-1)
     &     + 48*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 96*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 48*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 96*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 192*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 96*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 96*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 168*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 48*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 48*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 72*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 6*HAr1(-1)
     &     + 6*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 12*HAr1(-1)*upx**(-2)
     &     + 6*HAr1(-1)*upx**(-2)*y
     &     - 6*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 12*HAr1(-1)*upx**(-1)
     &     - 6*HAr1(-1)*upx**(-1)*y
     &     - 6*HAr1(-1)*y**(-1)
     &     + 6*HAr1(-1)*Hr1(0)
     &     + 24*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 48*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 24*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 96*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 24*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 60*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 24*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 12*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 12*HBr1(-1)
     &     - 6*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 6*HBr1(-1)*upx**(-2)*y
     &     + 6*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 6*HBr1(-1)*upx**(-1)*y
     &     + 6*HBr1(-1)*z**(-1)
     &     - 18*HBr1(-1)*Hr1(0)
     &     - 24*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 48*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 24*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 48*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 96*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 48*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 48*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 84*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 24*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 24*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 36*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 40*Hr2(-1,-1)
     &     - 32*Hr2(-1,-1)*upx**(-4)
     &     - 64*Hr2(-1,-1)*upx**(-4)*y
     &     - 32*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 64*Hr2(-1,-1)*upx**(-3)
     &     + 128*Hr2(-1,-1)*upx**(-3)*y
     &     + 64*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 80*Hr2(-1,-1)*upx**(-2)
     &     - 144*Hr2(-1,-1)*upx**(-2)*y
     &     - 32*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 48*Hr2(-1,-1)*upx**(-1)
     &     + 80*Hr2(-1,-1)*upx**(-1)*y
     &     - 12*HAr2(-1,-1)
     &     - 24*HAr2(-1,-1)*upx**(-2)*y
     &     + 24*HAr2(-1,-1)*upx**(-1)*y
     &     - 12*HBr2(-1,-1)
     &     - 48*HBr2(-1,-1)*upx**(-2)
     &     - 24*HBr2(-1,-1)*upx**(-2)*y
     &     + 48*HBr2(-1,-1)*upx**(-1)
     &     + 24*HBr2(-1,-1)*upx**(-1)*y
     &     + 20*Hr2(-1,0)
     &     + 16*Hr2(-1,0)*upx**(-4)
     &     + 32*Hr2(-1,0)*upx**(-4)*y
     &     + 16*Hr2(-1,0)*upx**(-4)*y**2
     &     - 32*Hr2(-1,0)*upx**(-3)
     &     - 64*Hr2(-1,0)*upx**(-3)*y
     &     - 32*Hr2(-1,0)*upx**(-3)*y**2
     &     + 40*Hr2(-1,0)*upx**(-2)
     &     + 72*Hr2(-1,0)*upx**(-2)*y
     &     + 16*Hr2(-1,0)*upx**(-2)*y**2
     &     - 24*Hr2(-1,0)*upx**(-1)
     &     - 40*Hr2(-1,0)*upx**(-1)*y
     &     + 6*Hr1(0)
     &     + 12*Hr1(0)*umx**(-2)
     &     - 12*Hr1(0)*umx**(-2)*y
     &     - 17.D0/2*Hr1(0)*umx**(-1)
     &     + 13*Hr1(0)*umx**(-1)*y
     &     - 1.D0/2*Hr1(0)*umx**(-1)*y**2
     &     - 16*Hr1(0)*upx**(-5)
     &     - 32*Hr1(0)*upx**(-5)*y
     &     - 16*Hr1(0)*upx**(-5)*y**2
     &     + 28*Hr1(0)*upx**(-4)
     &     + 56*Hr1(0)*upx**(-4)*y
     &     + 28*Hr1(0)*upx**(-4)*y**2
     &     - 10*Hr1(0)*upx**(-3)
     &     - 36*Hr1(0)*upx**(-3)*y
     &     - 10*Hr1(0)*upx**(-3)*y**2
     &     - Hr1(0)*upx**(-2)
     &     + 10*Hr1(0)*upx**(-2)*y
     &     - Hr1(0)*upx**(-2)*y**2
     &     - 17.D0/2*Hr1(0)*upx**(-1)
     &     + Hr1(0)*upx**(-1)*y
     &     - 1.D0/2*Hr1(0)*upx**(-1)*y**2
     &     + 8*Hr2(0,-1)
     &     - 24*Hr2(0,-1)*umx**(-3)
     &     + 24*Hr2(0,-1)*umx**(-3)*y
     &     + 36*Hr2(0,-1)*umx**(-2)
     &     - 36*Hr2(0,-1)*umx**(-2)*y
     &     - 12*Hr2(0,-1)*umx**(-1)
     &     + 16*Hr2(0,-1)*upx**(-4)
     &     + 32*Hr2(0,-1)*upx**(-4)*y
     &     + 16*Hr2(0,-1)*upx**(-4)*y**2
     &     - 8*Hr2(0,-1)*upx**(-3)
     &     - 40*Hr2(0,-1)*upx**(-3)*y
     &     - 32*Hr2(0,-1)*upx**(-3)*y**2
     &     + 4*Hr2(0,-1)*upx**(-2)
     &     + 36*Hr2(0,-1)*upx**(-2)*y
     &     + 16*Hr2(0,-1)*upx**(-2)*y**2
     &     + 12*Hr2(0,-1)*upx**(-1)
     &     - 16*Hr2(0,-1)*upx**(-1)*y
     &     + 6*HAr2(0,-1)
     &     + 12*HAr2(0,-1)*upx**(-2)*y
     &     - 12*HAr2(0,-1)*upx**(-1)*y
     &     + 6*HBr2(0,-1)
     &     + 24*HBr2(0,-1)*upx**(-2)
     &     + 12*HBr2(0,-1)*upx**(-2)*y
     &     - 24*HBr2(0,-1)*upx**(-1)
     &     - 12*HBr2(0,-1)*upx**(-1)*y
     &     + 12*Hr2(0,0)*umx**(-3)
     &     - 12*Hr2(0,0)*umx**(-3)*y
     &     - 18*Hr2(0,0)*umx**(-2)
     &     + 18*Hr2(0,0)*umx**(-2)*y
     &     + 3.D0/2*Hr2(0,0)*umx**(-1)
     &     + Hr2(0,0)*umx**(-1)*y
     &     - 1.D0/2*Hr2(0,0)*umx**(-1)*y**2
     &     - 8*Hr2(0,0)*upx**(-5)
     &     - 16*Hr2(0,0)*upx**(-5)*y
     &     - 8*Hr2(0,0)*upx**(-5)*y**2
     &     + 12*Hr2(0,0)*upx**(-4)
     &     + 24*Hr2(0,0)*upx**(-4)*y
     &     + 12*Hr2(0,0)*upx**(-4)*y**2
     &     - 14*Hr2(0,0)*upx**(-3)
     &     - 24*Hr2(0,0)*upx**(-3)*y
     &     - 2*Hr2(0,0)*upx**(-3)*y**2
     &     + 5*Hr2(0,0)*upx**(-2)
     &     + 8*Hr2(0,0)*upx**(-2)*y
     &     - Hr2(0,0)*upx**(-2)*y**2
     &     - 21.D0/2*Hr2(0,0)*upx**(-1)
     &     + Hr2(0,0)*upx**(-1)*y
     &     - 1.D0/2*Hr2(0,0)*upx**(-1)*y**2
     &     + 8*Hr2(1,0)
     &     - 9*Hr2(1,0)*umx**(-1)
     &     + 2*Hr2(1,0)*umx**(-1)*y
     &     - Hr2(1,0)*umx**(-1)*y**2
     &     - 16*Hr2(1,0)*upx**(-5)
     &     - 32*Hr2(1,0)*upx**(-5)*y
     &     - 16*Hr2(1,0)*upx**(-5)*y**2
     &     + 40*Hr2(1,0)*upx**(-4)
     &     + 80*Hr2(1,0)*upx**(-4)*y
     &     + 40*Hr2(1,0)*upx**(-4)*y**2
     &     - 36*Hr2(1,0)*upx**(-3)
     &     - 88*Hr2(1,0)*upx**(-3)*y
     &     - 36*Hr2(1,0)*upx**(-3)*y**2
     &     + 14*Hr2(1,0)*upx**(-2)
     &     + 52*Hr2(1,0)*upx**(-2)*y
     &     + 14*Hr2(1,0)*upx**(-2)*y**2
     &     - 9*Hr2(1,0)*upx**(-1)
     &     - 14*Hr2(1,0)*upx**(-1)*y
     &     - Hr2(1,0)*upx**(-1)*y**2
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**(-1)*Nh * (
     &     + 40.D0/9
     &     - 64.D0/3*upx**(-6)
     &     - 128.D0/3*upx**(-6)*y
     &     - 64.D0/3*upx**(-6)*y**2
     &     + 64*upx**(-5)
     &     + 128*upx**(-5)*y
     &     + 64*upx**(-5)*y**2
     &     - 496.D0/9*upx**(-4)
     &     - 1184.D0/9*upx**(-4)*y
     &     - 496.D0/9*upx**(-4)*y**2
     &     + 32.D0/9*upx**(-3)
     &     + 448.D0/9*upx**(-3)*y
     &     + 32.D0/9*upx**(-3)*y**2
     &     - 16.D0/9*upx**(-2)
     &     + 16.D0/3*upx**(-2)*y
     &     + 80.D0/9*upx**(-2)*y**2
     &     + 32.D0/3*upx**(-1)
     &     - 80.D0/9*upx**(-1)*y
     &     - 8.D0/3*Hr1(0)
     &     - 64.D0/3*Hr1(0)*upx**(-7)
     &     - 128.D0/3*Hr1(0)*upx**(-7)*y
     &     - 64.D0/3*Hr1(0)*upx**(-7)*y**2
     &     + 224.D0/3*Hr1(0)*upx**(-6)
     &     + 448.D0/3*Hr1(0)*upx**(-6)*y
     &     + 224.D0/3*Hr1(0)*upx**(-6)*y**2
     &     - 256.D0/3*Hr1(0)*upx**(-5)
     &     - 192*Hr1(0)*upx**(-5)*y
     &     - 256.D0/3*Hr1(0)*upx**(-5)*y**2
     &     + 80.D0/3*Hr1(0)*upx**(-4)
     &     + 320.D0/3*Hr1(0)*upx**(-4)*y
     &     + 80.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y
     &     + 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 32.D0/3*Hr1(0)*upx**(-2)
     &     - 16*Hr1(0)*upx**(-2)*y
     &     - 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**(-1)*Nl * (
     &     + 40.D0/9
     &     + 80.D0/9*upx**(-4)
     &     + 160.D0/9*upx**(-4)*y
     &     + 80.D0/9*upx**(-4)*y**2
     &     - 160.D0/9*upx**(-3)
     &     - 320.D0/9*upx**(-3)*y
     &     - 160.D0/9*upx**(-3)*y**2
     &     + 80.D0/9*upx**(-2)
     &     + 80.D0/3*upx**(-2)*y
     &     + 80.D0/9*upx**(-2)*y**2
     &     - 80.D0/9*upx**(-1)*y
     &     - 16.D0/3*Hr1(-1)
     &     - 32.D0/3*Hr1(-1)*upx**(-4)
     &     - 64.D0/3*Hr1(-1)*upx**(-4)*y
     &     - 32.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     + 64.D0/3*Hr1(-1)*upx**(-3)
     &     + 128.D0/3*Hr1(-1)*upx**(-3)*y
     &     + 64.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     - 32.D0/3*Hr1(-1)*upx**(-2)
     &     - 32*Hr1(-1)*upx**(-2)*y
     &     - 32.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 32.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 8.D0/3*Hr1(0)
     &     + 16.D0/3*Hr1(0)*upx**(-4)
     &     + 32.D0/3*Hr1(0)*upx**(-4)*y
     &     + 16.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0/3*Hr1(0)*upx**(-3)
     &     - 64.D0/3*Hr1(0)*upx**(-3)*y
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-2)
     &     + 16*Hr1(0)*upx**(-2)*y
     &     + 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-1)*y
     &     + 8.D0/3*DLog(tm**(-2)*w**2)
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 64.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 16*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + (
     &     - 25.D0/9
     &     - 12*umx**(-5)*z2
     &     + 24*umx**(-5)*y*z2
     &     - 12*umx**(-5)*y**2*z2
     &     + 30*umx**(-4)*z2
     &     - 60*umx**(-4)*y*z2
     &     + 30*umx**(-4)*y**2*z2
     &     - 47*umx**(-3)*z2
     &     + 70*umx**(-3)*y*z2
     &     - 19*umx**(-3)*y**2*z2
     &     - 2*umx**(-2)
     &     + 81.D0/2*umx**(-2)*z2
     &     + 4*umx**(-2)*y
     &     - 45*umx**(-2)*y*z2
     &     - 2*umx**(-2)*y**2
     &     - 3.D0/2*umx**(-2)*y**2*z2
     &     + 2*umx**(-1)
     &     - 181.D0/4*umx**(-1)*z2
     &     - 4*umx**(-1)*y
     &     + 15.D0/2*umx**(-1)*y*z2
     &     + 2*umx**(-1)*y**2
     &     - 13.D0/4*umx**(-1)*y**2*z2
     &     - 48*upx**(-5)*z2
     &     - 96*upx**(-5)*y*z2
     &     - 48*upx**(-5)*y**2*z2
     &     + 472.D0/9*upx**(-4)
     &     + 64*upx**(-4)*z2
     &     + 944.D0/9*upx**(-4)*y
     &     + 128*upx**(-4)*y*z2
     &     + 472.D0/9*upx**(-4)*y**2
     &     + 64*upx**(-4)*y**2*z2
     &     - 944.D0/9*upx**(-3)
     &     + 18*upx**(-3)*z2
     &     - 1888.D0/9*upx**(-3)*y
     &     - 32*upx**(-3)*y*z2
     &     - 944.D0/9*upx**(-3)*y**2
     &     - 2*upx**(-3)*y**2*z2
     &     + 454.D0/9*upx**(-2)
     &     - 19*upx**(-2)*z2
     &     + 460.D0/3*upx**(-2)*y
     &     - 4*upx**(-2)*y*z2
     &     + 454.D0/9*upx**(-2)*y**2
     &     - 5*upx**(-2)*y**2*z2
     &     + 2*upx**(-1)
     &     - 69.D0/4*upx**(-1)*z2
     &     - 436.D0/9*upx**(-1)*y
     &     + 15.D0/2*upx**(-1)*y*z2
     &     + 2*upx**(-1)*y**2
     &     - 13.D0/4*upx**(-1)*y**2*z2
     &     + 70.D0/3*Hr1(-1)
     &     + 24*Hr1(-1)*umx**(-4)
     &     - 48*Hr1(-1)*umx**(-4)*y
     &     + 24*Hr1(-1)*umx**(-4)*y**2
     &     - 48*Hr1(-1)*umx**(-3)
     &     + 96*Hr1(-1)*umx**(-3)*y
     &     - 48*Hr1(-1)*umx**(-3)*y**2
     &     + 72*Hr1(-1)*umx**(-2)
     &     - 96*Hr1(-1)*umx**(-2)*y
     &     + 16*Hr1(-1)*umx**(-2)*y**2
     &     - 48*Hr1(-1)*umx**(-1)
     &     + 48*Hr1(-1)*umx**(-1)*y
     &     + 8*Hr1(-1)*umx**(-1)*y**2
     &     - 52.D0/3*Hr1(-1)*upx**(-4)
     &     - 104.D0/3*Hr1(-1)*upx**(-4)*y
     &     - 52.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     + 104.D0/3*Hr1(-1)*upx**(-3)
     &     + 208.D0/3*Hr1(-1)*upx**(-3)*y
     &     + 104.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     - 76.D0/3*Hr1(-1)*upx**(-2)
     &     - 68*Hr1(-1)*upx**(-2)*y
     &     - 76.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 8*Hr1(-1)*upx**(-1)
     &     + 100.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 8*Hr1(-1)*upx**(-1)*y**2
     &     + 24*Hr1(-1)*HAr1(-1)
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 192*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 192*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 384*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 192*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 240*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 48*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     - 48*Hr1(-1)*HBr1(-1)
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     - 128*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     + 128*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     + 256*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     + 128*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     - 128*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     - 224*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     + 64*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     + 96*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     + 12*HAr1(-1)
     &     - 12*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 24*HAr1(-1)*upx**(-2)
     &     - 12*HAr1(-1)*upx**(-2)*y
     &     + 12*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 24*HAr1(-1)*upx**(-1)
     &     + 12*HAr1(-1)*upx**(-1)*y
     &     + 12*HAr1(-1)*y**(-1)
     &     - 12*HAr1(-1)*Hr1(0)
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 96*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 96*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 192*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 96*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 120*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 24*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 16*HBr1(-1)
     &     + 8*HBr1(-1)*upx**(-2)*z**(-1)
     &     - 8*HBr1(-1)*upx**(-2)*y
     &     - 8*HBr1(-1)*upx**(-1)*z**(-1)
     &     + 8*HBr1(-1)*upx**(-1)*y
     &     - 8*HBr1(-1)*z**(-1)
     &     + 24*HBr1(-1)*Hr1(0)
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-4)
     &     + 64*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 64*HBr1(-1)*Hr1(0)*upx**(-3)
     &     - 128*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 64*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 64*HBr1(-1)*Hr1(0)*upx**(-2)
     &     + 112*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 32*HBr1(-1)*Hr1(0)*upx**(-1)
     &     - 48*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 48*Hr2(-1,-1)
     &     + 16*Hr2(-1,-1)*upx**(-4)
     &     + 32*Hr2(-1,-1)*upx**(-4)*y
     &     + 16*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 32*Hr2(-1,-1)*upx**(-3)
     &     - 64*Hr2(-1,-1)*upx**(-3)*y
     &     - 32*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 80*Hr2(-1,-1)*upx**(-2)
     &     + 128*Hr2(-1,-1)*upx**(-2)*y
     &     + 16*Hr2(-1,-1)*upx**(-2)*y**2
     &     - 64*Hr2(-1,-1)*upx**(-1)
     &     - 96*Hr2(-1,-1)*upx**(-1)*y
     &     + 24*HAr2(-1,-1)
     &     + 48*HAr2(-1,-1)*upx**(-2)*y
     &     - 48*HAr2(-1,-1)*upx**(-1)*y
     &     + 16*HBr2(-1,-1)
     &     + 64*HBr2(-1,-1)*upx**(-2)
     &     + 32*HBr2(-1,-1)*upx**(-2)*y
     &     - 64*HBr2(-1,-1)*upx**(-1)
     &     - 32*HBr2(-1,-1)*upx**(-1)*y
     &     - 24*Hr2(-1,0)
     &     - 8*Hr2(-1,0)*upx**(-4)
     &     - 16*Hr2(-1,0)*upx**(-4)*y
     &     - 8*Hr2(-1,0)*upx**(-4)*y**2
     &     + 16*Hr2(-1,0)*upx**(-3)
     &     + 32*Hr2(-1,0)*upx**(-3)*y
     &     + 16*Hr2(-1,0)*upx**(-3)*y**2
     &     - 40*Hr2(-1,0)*upx**(-2)
     &     - 64*Hr2(-1,0)*upx**(-2)*y
     &     - 8*Hr2(-1,0)*upx**(-2)*y**2
     &     + 32*Hr2(-1,0)*upx**(-1)
     &     + 48*Hr2(-1,0)*upx**(-1)*y
     &     - 44.D0/3*Hr1(0)
     &     - 12*Hr1(0)*umx**(-4)
     &     + 24*Hr1(0)*umx**(-4)*y
     &     - 12*Hr1(0)*umx**(-4)*y**2
     &     + 24*Hr1(0)*umx**(-3)
     &     - 48*Hr1(0)*umx**(-3)*y
     &     + 24*Hr1(0)*umx**(-3)*y**2
     &     - 36*Hr1(0)*umx**(-2)
     &     + 48*Hr1(0)*umx**(-2)*y
     &     - 8*Hr1(0)*umx**(-2)*y**2
     &     + 75.D0/4*Hr1(0)*umx**(-1)
     &     - 51.D0/2*Hr1(0)*umx**(-1)*y
     &     - 13.D0/4*Hr1(0)*umx**(-1)*y**2
     &     + 24*Hr1(0)*upx**(-5)
     &     + 48*Hr1(0)*upx**(-5)*y
     &     + 24*Hr1(0)*upx**(-5)*y**2
     &     - 154.D0/3*Hr1(0)*upx**(-4)
     &     - 308.D0/3*Hr1(0)*upx**(-4)*y
     &     - 154.D0/3*Hr1(0)*upx**(-4)*y**2
     &     + 101.D0/3*Hr1(0)*upx**(-3)
     &     + 274.D0/3*Hr1(0)*upx**(-3)*y
     &     + 101.D0/3*Hr1(0)*upx**(-3)*y**2
     &     - 23.D0/6*Hr1(0)*upx**(-2)
     &     - 35*Hr1(0)*upx**(-2)*y
     &     - 23.D0/6*Hr1(0)*upx**(-2)*y**2
     &     + 35.D0/4*Hr1(0)*upx**(-1)
     &     - 1.D0/6*Hr1(0)*upx**(-1)*y
     &     - 13.D0/4*Hr1(0)*upx**(-1)*y**2
     &     - 12*Hr2(0,-1)
     &     + 24*Hr2(0,-1)*umx**(-5)
     &     - 48*Hr2(0,-1)*umx**(-5)*y
     &     + 24*Hr2(0,-1)*umx**(-5)*y**2
     &     - 60*Hr2(0,-1)*umx**(-4)
     &     + 120*Hr2(0,-1)*umx**(-4)*y
     &     - 60*Hr2(0,-1)*umx**(-4)*y**2
     &     + 94*Hr2(0,-1)*umx**(-3)
     &     - 140*Hr2(0,-1)*umx**(-3)*y
     &     + 38*Hr2(0,-1)*umx**(-3)*y**2
     &     - 81*Hr2(0,-1)*umx**(-2)
     &     + 90*Hr2(0,-1)*umx**(-2)*y
     &     + 3*Hr2(0,-1)*umx**(-2)*y**2
     &     + 73.D0/2*Hr2(0,-1)*umx**(-1)
     &     - 3*Hr2(0,-1)*umx**(-1)*y
     &     + 1.D0/2*Hr2(0,-1)*umx**(-1)*y**2
     &     - 8*Hr2(0,-1)*upx**(-4)
     &     - 16*Hr2(0,-1)*upx**(-4)*y
     &     - 8*Hr2(0,-1)*upx**(-4)*y**2
     &     - 12*Hr2(0,-1)*upx**(-3)
     &     + 16*Hr2(0,-1)*upx**(-3)*y
     &     + 28*Hr2(0,-1)*upx**(-3)*y**2
     &     + 2*Hr2(0,-1)*upx**(-2)
     &     - 40*Hr2(0,-1)*upx**(-2)*y
     &     - 26*Hr2(0,-1)*upx**(-2)*y**2
     &     - 39.D0/2*Hr2(0,-1)*upx**(-1)
     &     + 21*Hr2(0,-1)*upx**(-1)*y
     &     + 1.D0/2*Hr2(0,-1)*upx**(-1)*y**2
     &     - 12*HAr2(0,-1)
     &     - 24*HAr2(0,-1)*upx**(-2)*y
     &     + 24*HAr2(0,-1)*upx**(-1)*y
     &     - 8*HBr2(0,-1)
     &     - 32*HBr2(0,-1)*upx**(-2)
     &     - 16*HBr2(0,-1)*upx**(-2)*y
     &     + 32*HBr2(0,-1)*upx**(-1)
     &     + 16*HBr2(0,-1)*upx**(-1)*y
     &     - 12*Hr2(0,0)*umx**(-5)
     &     + 24*Hr2(0,0)*umx**(-5)*y
     &     - 12*Hr2(0,0)*umx**(-5)*y**2
     &     + 30*Hr2(0,0)*umx**(-4)
     &     - 60*Hr2(0,0)*umx**(-4)*y
     &     + 30*Hr2(0,0)*umx**(-4)*y**2
     &     - 47*Hr2(0,0)*umx**(-3)
     &     + 70*Hr2(0,0)*umx**(-3)*y
     &     - 19*Hr2(0,0)*umx**(-3)*y**2
     &     + 81.D0/2*Hr2(0,0)*umx**(-2)
     &     - 45*Hr2(0,0)*umx**(-2)*y
     &     - 3.D0/2*Hr2(0,0)*umx**(-2)*y**2
     &     - 23.D0/2*Hr2(0,0)*umx**(-1)
     &     + 1.D0/2*Hr2(0,0)*umx**(-1)*y**2
     &     + 12*Hr2(0,0)*upx**(-5)
     &     + 24*Hr2(0,0)*upx**(-5)*y
     &     + 12*Hr2(0,0)*upx**(-5)*y**2
     &     - 26*Hr2(0,0)*upx**(-4)
     &     - 52*Hr2(0,0)*upx**(-4)*y
     &     - 26*Hr2(0,0)*upx**(-4)*y**2
     &     + 33*Hr2(0,0)*upx**(-3)
     &     + 58*Hr2(0,0)*upx**(-3)*y
     &     + 13*Hr2(0,0)*upx**(-3)*y**2
     &     - 23.D0/2*Hr2(0,0)*upx**(-2)
     &     - 19*Hr2(0,0)*upx**(-2)*y
     &     + 5.D0/2*Hr2(0,0)*upx**(-2)*y**2
     &     + 33.D0/2*Hr2(0,0)*upx**(-1)
     &     + 1.D0/2*Hr2(0,0)*upx**(-1)*y**2
     &     - 12*Hr2(1,0)
     &     + 27.D0/2*Hr2(1,0)*umx**(-1)
     &     - 3*Hr2(1,0)*umx**(-1)*y
     &     + 3.D0/2*Hr2(1,0)*umx**(-1)*y**2
     &     + 24*Hr2(1,0)*upx**(-5)
     &     + 48*Hr2(1,0)*upx**(-5)*y
     &     + 24*Hr2(1,0)*upx**(-5)*y**2
     &     - 60*Hr2(1,0)*upx**(-4)
     &     - 120*Hr2(1,0)*upx**(-4)*y
     &     - 60*Hr2(1,0)*upx**(-4)*y**2
     &     + 54*Hr2(1,0)*upx**(-3)
     &     + 132*Hr2(1,0)*upx**(-3)*y
     &     + 54*Hr2(1,0)*upx**(-3)*y**2
     &     - 21*Hr2(1,0)*upx**(-2)
     &     - 78*Hr2(1,0)*upx**(-2)*y
     &     - 21*Hr2(1,0)*upx**(-2)*y**2
     &     + 27.D0/2*Hr2(1,0)*upx**(-1)
     &     + 21*Hr2(1,0)*upx**(-1)*y
     &     + 3.D0/2*Hr2(1,0)*upx**(-1)*y**2
     &     - 44.D0/3*DLog(tm**(-2)*w**2)
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 352.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 88*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc*Nh * (
     &     - 20.D0/3
     &     + 32*upx**(-6)
     &     + 64*upx**(-6)*y
     &     + 32*upx**(-6)*y**2
     &     - 96*upx**(-5)
     &     - 192*upx**(-5)*y
     &     - 96*upx**(-5)*y**2
     &     + 248.D0/3*upx**(-4)
     &     + 592.D0/3*upx**(-4)*y
     &     + 248.D0/3*upx**(-4)*y**2
     &     - 16.D0/3*upx**(-3)
     &     - 224.D0/3*upx**(-3)*y
     &     - 16.D0/3*upx**(-3)*y**2
     &     + 8.D0/3*upx**(-2)
     &     - 8*upx**(-2)*y
     &     - 40.D0/3*upx**(-2)*y**2
     &     - 16*upx**(-1)
     &     + 40.D0/3*upx**(-1)*y
     &     + 4*Hr1(0)
     &     + 32*Hr1(0)*upx**(-7)
     &     + 64*Hr1(0)*upx**(-7)*y
     &     + 32*Hr1(0)*upx**(-7)*y**2
     &     - 112*Hr1(0)*upx**(-6)
     &     - 224*Hr1(0)*upx**(-6)*y
     &     - 112*Hr1(0)*upx**(-6)*y**2
     &     + 128*Hr1(0)*upx**(-5)
     &     + 288*Hr1(0)*upx**(-5)*y
     &     + 128*Hr1(0)*upx**(-5)*y**2
     &     - 40*Hr1(0)*upx**(-4)
     &     - 160*Hr1(0)*upx**(-4)*y
     &     - 40*Hr1(0)*upx**(-4)*y**2
     &     + 16*Hr1(0)*upx**(-3)*y
     &     - 16*Hr1(0)*upx**(-3)*y**2
     &     - 16*Hr1(0)*upx**(-2)
     &     + 24*Hr1(0)*upx**(-2)*y
     &     + 8*Hr1(0)*upx**(-2)*y**2
     &     - 8*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc*Nl * (
     &     - 20.D0/3
     &     - 40.D0/3*upx**(-4)
     &     - 80.D0/3*upx**(-4)*y
     &     - 40.D0/3*upx**(-4)*y**2
     &     + 80.D0/3*upx**(-3)
     &     + 160.D0/3*upx**(-3)*y
     &     + 80.D0/3*upx**(-3)*y**2
     &     - 40.D0/3*upx**(-2)
     &     - 40*upx**(-2)*y
     &     - 40.D0/3*upx**(-2)*y**2
     &     + 40.D0/3*upx**(-1)*y
     &     + 8*Hr1(-1)
     &     + 16*Hr1(-1)*upx**(-4)
     &     + 32*Hr1(-1)*upx**(-4)*y
     &     + 16*Hr1(-1)*upx**(-4)*y**2
     &     - 32*Hr1(-1)*upx**(-3)
     &     - 64*Hr1(-1)*upx**(-3)*y
     &     - 32*Hr1(-1)*upx**(-3)*y**2
     &     + 16*Hr1(-1)*upx**(-2)
     &     + 48*Hr1(-1)*upx**(-2)*y
     &     + 16*Hr1(-1)*upx**(-2)*y**2
     &     - 16*Hr1(-1)*upx**(-1)*y
     &     - 4*Hr1(0)
     &     - 8*Hr1(0)*upx**(-4)
     &     - 16*Hr1(0)*upx**(-4)*y
     &     - 8*Hr1(0)*upx**(-4)*y**2
     &     + 16*Hr1(0)*upx**(-3)
     &     + 32*Hr1(0)*upx**(-3)*y
     &     + 16*Hr1(0)*upx**(-3)*y**2
     &     - 8*Hr1(0)*upx**(-2)
     &     - 24*Hr1(0)*upx**(-2)*y
     &     - 8*Hr1(0)*upx**(-2)*y**2
     &     + 8*Hr1(0)*upx**(-1)*y
     &     - 4*DLog(tm**(-2)*w**2)
     &     - 8*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 16*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 8*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 16*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 32*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 16*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 8*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 24*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 8*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 8*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**2 * (
     &     + 86.D0/3
     &     + 18*umx**(-5)*z2
     &     - 36*umx**(-5)*y*z2
     &     + 18*umx**(-5)*y**2*z2
     &     - 45*umx**(-4)*z2
     &     + 90*umx**(-4)*y*z2
     &     - 45*umx**(-4)*y**2*z2
     &     + 101.D0/2*umx**(-3)*z2
     &     - 85*umx**(-3)*y*z2
     &     + 57.D0/2*umx**(-3)*y**2*z2
     &     + 3*umx**(-2)
     &     - 123.D0/4*umx**(-2)*z2
     &     - 6*umx**(-2)*y
     &     + 75.D0/2*umx**(-2)*y*z2
     &     + 3*umx**(-2)*y**2
     &     + 9.D0/4*umx**(-2)*y**2*z2
     &     - 3*umx**(-1)
     &     + 211.D0/8*umx**(-1)*z2
     &     + 6*umx**(-1)*y
     &     - 17.D0/4*umx**(-1)*y*z2
     &     - 3*umx**(-1)*y**2
     &     + 11.D0/8*umx**(-1)*y**2*z2
     &     + 16*upx**(-5)*z2
     &     + 32*upx**(-5)*y*z2
     &     + 16*upx**(-5)*y**2*z2
     &     + 16.D0/3*upx**(-4)
     &     - 12*upx**(-4)*z2
     &     + 32.D0/3*upx**(-4)*y
     &     - 24*upx**(-4)*y*z2
     &     + 16.D0/3*upx**(-4)*y**2
     &     - 12*upx**(-4)*y**2*z2
     &     - 32.D0/3*upx**(-3)
     &     - 21*upx**(-3)*z2
     &     - 64.D0/3*upx**(-3)*y
     &     - 16*upx**(-3)*y*z2
     &     - 32.D0/3*upx**(-3)*y**2
     &     - 11*upx**(-3)*y**2*z2
     &     + 25.D0/3*upx**(-2)
     &     + 23.D0/2*upx**(-2)*z2
     &     + 22*upx**(-2)*y
     &     + 10*upx**(-2)*y*z2
     &     + 25.D0/3*upx**(-2)*y**2
     &     + 1.D0/2*upx**(-2)*y**2*z2
     &     - 3*upx**(-1)
     &     + 35.D0/8*upx**(-1)*z2
     &     - 34.D0/3*upx**(-1)*y
     &     - 17.D0/4*upx**(-1)*y*z2
     &     - 3*upx**(-1)*y**2
     &     + 11.D0/8*upx**(-1)*y**2*z2
     &     - 22*Hr1(-1)
     &     - 36*Hr1(-1)*umx**(-4)
     &     + 72*Hr1(-1)*umx**(-4)*y
     &     - 36*Hr1(-1)*umx**(-4)*y**2
     &     + 72*Hr1(-1)*umx**(-3)
     &     - 144*Hr1(-1)*umx**(-3)*y
     &     + 72*Hr1(-1)*umx**(-3)*y**2
     &     - 68*Hr1(-1)*umx**(-2)
     &     + 104*Hr1(-1)*umx**(-2)*y
     &     - 24*Hr1(-1)*umx**(-2)*y**2
     &     + 32*Hr1(-1)*umx**(-1)
     &     - 32*Hr1(-1)*umx**(-1)*y
     &     - 12*Hr1(-1)*umx**(-1)*y**2
     &     - 16*Hr1(-1)*upx**(-4)
     &     - 32*Hr1(-1)*upx**(-4)*y
     &     - 16*Hr1(-1)*upx**(-4)*y**2
     &     + 32*Hr1(-1)*upx**(-3)
     &     + 64*Hr1(-1)*upx**(-3)*y
     &     + 32*Hr1(-1)*upx**(-3)*y**2
     &     - 4*Hr1(-1)*upx**(-2)
     &     - 24*Hr1(-1)*upx**(-2)*y
     &     - 4*Hr1(-1)*upx**(-2)*y**2
     &     - 12*Hr1(-1)*upx**(-1)
     &     - 8*Hr1(-1)*upx**(-1)*y
     &     - 12*Hr1(-1)*upx**(-1)*y**2
     &     - 16*Hr1(-1)*HAr1(-1)
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 128*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 128*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 256*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 128*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 160*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 32*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 12*Hr1(-1)*HBr1(-1)
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 32*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 32*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 56*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 16*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 24*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 8*HAr1(-1)
     &     + 8*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 16*HAr1(-1)*upx**(-2)
     &     + 8*HAr1(-1)*upx**(-2)*y
     &     - 8*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 16*HAr1(-1)*upx**(-1)
     &     - 8*HAr1(-1)*upx**(-1)*y
     &     - 8*HAr1(-1)*y**(-1)
     &     + 8*HAr1(-1)*Hr1(0)
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 64*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 64*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 128*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 64*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 80*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 16*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 4*HBr1(-1)
     &     - 2*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 2*HBr1(-1)*upx**(-2)*y
     &     + 2*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 2*HBr1(-1)*upx**(-1)*y
     &     + 2*HBr1(-1)*z**(-1)
     &     - 6*HBr1(-1)*Hr1(0)
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 28*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 8*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 12*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 4*Hr2(-1,-1)
     &     + 32*Hr2(-1,-1)*upx**(-4)
     &     + 64*Hr2(-1,-1)*upx**(-4)*y
     &     + 32*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 64*Hr2(-1,-1)*upx**(-3)
     &     - 128*Hr2(-1,-1)*upx**(-3)*y
     &     - 64*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 16*Hr2(-1,-1)*upx**(-2)
     &     + 56*Hr2(-1,-1)*upx**(-2)*y
     &     + 32*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 16*Hr2(-1,-1)*upx**(-1)
     &     + 8*Hr2(-1,-1)*upx**(-1)*y
     &     - 16*HAr2(-1,-1)
     &     - 32*HAr2(-1,-1)*upx**(-2)*y
     &     + 32*HAr2(-1,-1)*upx**(-1)*y
     &     - 4*HBr2(-1,-1)
     &     - 16*HBr2(-1,-1)*upx**(-2)
     &     - 8*HBr2(-1,-1)*upx**(-2)*y
     &     + 16*HBr2(-1,-1)*upx**(-1)
     &     + 8*HBr2(-1,-1)*upx**(-1)*y
     &     + 2*Hr2(-1,0)
     &     - 16*Hr2(-1,0)*upx**(-4)
     &     - 32*Hr2(-1,0)*upx**(-4)*y
     &     - 16*Hr2(-1,0)*upx**(-4)*y**2
     &     + 32*Hr2(-1,0)*upx**(-3)
     &     + 64*Hr2(-1,0)*upx**(-3)*y
     &     + 32*Hr2(-1,0)*upx**(-3)*y**2
     &     - 8*Hr2(-1,0)*upx**(-2)
     &     - 28*Hr2(-1,0)*upx**(-2)*y
     &     - 16*Hr2(-1,0)*upx**(-2)*y**2
     &     - 8*Hr2(-1,0)*upx**(-1)
     &     - 4*Hr2(-1,0)*upx**(-1)*y
     &     + 12*Hr1(0)
     &     + 18*Hr1(0)*umx**(-4)
     &     - 36*Hr1(0)*umx**(-4)*y
     &     + 18*Hr1(0)*umx**(-4)*y**2
     &     - 36*Hr1(0)*umx**(-3)
     &     + 72*Hr1(0)*umx**(-3)*y
     &     - 36*Hr1(0)*umx**(-3)*y**2
     &     + 34*Hr1(0)*umx**(-2)
     &     - 52*Hr1(0)*umx**(-2)*y
     &     + 12*Hr1(0)*umx**(-2)*y**2
     &     - 57.D0/4*Hr1(0)*umx**(-1)
     &     + 33.D0/2*Hr1(0)*umx**(-1)*y
     &     + 23.D0/4*Hr1(0)*umx**(-1)*y**2
     &     - 8*Hr1(0)*upx**(-5)
     &     - 16*Hr1(0)*upx**(-5)*y
     &     - 8*Hr1(0)*upx**(-5)*y**2
     &     + 28*Hr1(0)*upx**(-4)
     &     + 56*Hr1(0)*upx**(-4)*y
     &     + 28*Hr1(0)*upx**(-4)*y**2
     &     - 33*Hr1(0)*upx**(-3)
     &     - 74*Hr1(0)*upx**(-3)*y
     &     - 33*Hr1(0)*upx**(-3)*y**2
     &     + 15.D0/2*Hr1(0)*upx**(-2)
     &     + 35*Hr1(0)*upx**(-2)*y
     &     + 15.D0/2*Hr1(0)*upx**(-2)*y**2
     &     + 7.D0/4*Hr1(0)*upx**(-1)
     &     - 3.D0/2*Hr1(0)*upx**(-1)*y
     &     + 23.D0/4*Hr1(0)*upx**(-1)*y**2
     &     + 4*Hr2(0,-1)
     &     - 36*Hr2(0,-1)*umx**(-5)
     &     + 72*Hr2(0,-1)*umx**(-5)*y
     &     - 36*Hr2(0,-1)*umx**(-5)*y**2
     &     + 90*Hr2(0,-1)*umx**(-4)
     &     - 180*Hr2(0,-1)*umx**(-4)*y
     &     + 90*Hr2(0,-1)*umx**(-4)*y**2
     &     - 101*Hr2(0,-1)*umx**(-3)
     &     + 170*Hr2(0,-1)*umx**(-3)*y
     &     - 57*Hr2(0,-1)*umx**(-3)*y**2
     &     + 123.D0/2*Hr2(0,-1)*umx**(-2)
     &     - 75*Hr2(0,-1)*umx**(-2)*y
     &     - 9.D0/2*Hr2(0,-1)*umx**(-2)*y**2
     &     - 139.D0/4*Hr2(0,-1)*umx**(-1)
     &     + 9.D0/2*Hr2(0,-1)*umx**(-1)*y
     &     - 3.D0/4*Hr2(0,-1)*umx**(-1)*y**2
     &     - 16*Hr2(0,-1)*upx**(-4)
     &     - 32*Hr2(0,-1)*upx**(-4)*y
     &     - 16*Hr2(0,-1)*upx**(-4)*y**2
     &     + 34*Hr2(0,-1)*upx**(-3)
     &     + 48*Hr2(0,-1)*upx**(-3)*y
     &     + 14*Hr2(0,-1)*upx**(-3)*y**2
     &     - 11*Hr2(0,-1)*upx**(-2)
     &     - 4*Hr2(0,-1)*upx**(-2)*y
     &     + 11*Hr2(0,-1)*upx**(-2)*y**2
     &     + 37.D0/4*Hr2(0,-1)*upx**(-1)
     &     - 7.D0/2*Hr2(0,-1)*upx**(-1)*y
     &     - 3.D0/4*Hr2(0,-1)*upx**(-1)*y**2
     &     + 8*HAr2(0,-1)
     &     + 16*HAr2(0,-1)*upx**(-2)*y
     &     - 16*HAr2(0,-1)*upx**(-1)*y
     &     + 2*HBr2(0,-1)
     &     + 8*HBr2(0,-1)*upx**(-2)
     &     + 4*HBr2(0,-1)*upx**(-2)*y
     &     - 8*HBr2(0,-1)*upx**(-1)
     &     - 4*HBr2(0,-1)*upx**(-1)*y
     &     + 18*Hr2(0,0)*umx**(-5)
     &     - 36*Hr2(0,0)*umx**(-5)*y
     &     + 18*Hr2(0,0)*umx**(-5)*y**2
     &     - 45*Hr2(0,0)*umx**(-4)
     &     + 90*Hr2(0,0)*umx**(-4)*y
     &     - 45*Hr2(0,0)*umx**(-4)*y**2
     &     + 101.D0/2*Hr2(0,0)*umx**(-3)
     &     - 85*Hr2(0,0)*umx**(-3)*y
     &     + 57.D0/2*Hr2(0,0)*umx**(-3)*y**2
     &     - 123.D0/4*Hr2(0,0)*umx**(-2)
     &     + 75.D0/2*Hr2(0,0)*umx**(-2)*y
     &     + 9.D0/4*Hr2(0,0)*umx**(-2)*y**2
     &     + 121.D0/8*Hr2(0,0)*umx**(-1)
     &     - 7.D0/4*Hr2(0,0)*umx**(-1)*y
     &     + 1.D0/8*Hr2(0,0)*umx**(-1)*y**2
     &     - 4*Hr2(0,0)*upx**(-5)
     &     - 8*Hr2(0,0)*upx**(-5)*y
     &     - 4*Hr2(0,0)*upx**(-5)*y**2
     &     + 18*Hr2(0,0)*upx**(-4)
     &     + 36*Hr2(0,0)*upx**(-4)*y
     &     + 18*Hr2(0,0)*upx**(-4)*y**2
     &     - 26*Hr2(0,0)*upx**(-3)
     &     - 46*Hr2(0,0)*upx**(-3)*y
     &     - 16*Hr2(0,0)*upx**(-3)*y**2
     &     + 9*Hr2(0,0)*upx**(-2)
     &     + 15*Hr2(0,0)*upx**(-2)*y
     &     - 2*Hr2(0,0)*upx**(-2)*y**2
     &     - 55.D0/8*Hr2(0,0)*upx**(-1)
     &     - 7.D0/4*Hr2(0,0)*upx**(-1)*y
     &     + 1.D0/8*Hr2(0,0)*upx**(-1)*y**2
     &     + 4*Hr2(1,0)
     &     - 9.D0/2*Hr2(1,0)*umx**(-1)
     &     + Hr2(1,0)*umx**(-1)*y
     &     - 1.D0/2*Hr2(1,0)*umx**(-1)*y**2
     &     - 8*Hr2(1,0)*upx**(-5)
     &     - 16*Hr2(1,0)*upx**(-5)*y
     &     - 8*Hr2(1,0)*upx**(-5)*y**2
     &     + 20*Hr2(1,0)*upx**(-4)
     &     + 40*Hr2(1,0)*upx**(-4)*y
     &     + 20*Hr2(1,0)*upx**(-4)*y**2
     &     - 18*Hr2(1,0)*upx**(-3)
     &     - 44*Hr2(1,0)*upx**(-3)*y
     &     - 18*Hr2(1,0)*upx**(-3)*y**2
     &     + 7*Hr2(1,0)*upx**(-2)
     &     + 26*Hr2(1,0)*upx**(-2)*y
     &     + 7*Hr2(1,0)*upx**(-2)*y**2
     &     - 9.D0/2*Hr2(1,0)*upx**(-1)
     &     - 7*Hr2(1,0)*upx**(-1)*y
     &     - 1.D0/2*Hr2(1,0)*upx**(-1)*y**2
     &     + 22*DLog(tm**(-2)*w**2)
     &     + 44*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 88*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 44*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 88*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 176*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 88*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 44*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 132*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 44*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 44*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**3*Nh * (
     &     + 20.D0/9
     &     - 32.D0/3*upx**(-6)
     &     - 64.D0/3*upx**(-6)*y
     &     - 32.D0/3*upx**(-6)*y**2
     &     + 32*upx**(-5)
     &     + 64*upx**(-5)*y
     &     + 32*upx**(-5)*y**2
     &     - 248.D0/9*upx**(-4)
     &     - 592.D0/9*upx**(-4)*y
     &     - 248.D0/9*upx**(-4)*y**2
     &     + 16.D0/9*upx**(-3)
     &     + 224.D0/9*upx**(-3)*y
     &     + 16.D0/9*upx**(-3)*y**2
     &     - 8.D0/9*upx**(-2)
     &     + 8.D0/3*upx**(-2)*y
     &     + 40.D0/9*upx**(-2)*y**2
     &     + 16.D0/3*upx**(-1)
     &     - 40.D0/9*upx**(-1)*y
     &     - 4.D0/3*Hr1(0)
     &     - 32.D0/3*Hr1(0)*upx**(-7)
     &     - 64.D0/3*Hr1(0)*upx**(-7)*y
     &     - 32.D0/3*Hr1(0)*upx**(-7)*y**2
     &     + 112.D0/3*Hr1(0)*upx**(-6)
     &     + 224.D0/3*Hr1(0)*upx**(-6)*y
     &     + 112.D0/3*Hr1(0)*upx**(-6)*y**2
     &     - 128.D0/3*Hr1(0)*upx**(-5)
     &     - 96*Hr1(0)*upx**(-5)*y
     &     - 128.D0/3*Hr1(0)*upx**(-5)*y**2
     &     + 40.D0/3*Hr1(0)*upx**(-4)
     &     + 160.D0/3*Hr1(0)*upx**(-4)*y
     &     + 40.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-3)*y
     &     + 16.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-2)
     &     - 8*Hr1(0)*upx**(-2)*y
     &     - 8.D0/3*Hr1(0)*upx**(-2)*y**2
     &     + 8.D0/3*Hr1(0)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**3*Nl * (
     &     + 20.D0/9
     &     + 40.D0/9*upx**(-4)
     &     + 80.D0/9*upx**(-4)*y
     &     + 40.D0/9*upx**(-4)*y**2
     &     - 80.D0/9*upx**(-3)
     &     - 160.D0/9*upx**(-3)*y
     &     - 80.D0/9*upx**(-3)*y**2
     &     + 40.D0/9*upx**(-2)
     &     + 40.D0/3*upx**(-2)*y
     &     + 40.D0/9*upx**(-2)*y**2
     &     - 40.D0/9*upx**(-1)*y
     &     - 8.D0/3*Hr1(-1)
     &     - 16.D0/3*Hr1(-1)*upx**(-4)
     &     - 32.D0/3*Hr1(-1)*upx**(-4)*y
     &     - 16.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     + 32.D0/3*Hr1(-1)*upx**(-3)
     &     + 64.D0/3*Hr1(-1)*upx**(-3)*y
     &     + 32.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     - 16.D0/3*Hr1(-1)*upx**(-2)
     &     - 16*Hr1(-1)*upx**(-2)*y
     &     - 16.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 16.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 4.D0/3*Hr1(0)
     &     + 8.D0/3*Hr1(0)*upx**(-4)
     &     + 16.D0/3*Hr1(0)*upx**(-4)*y
     &     + 8.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-3)
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y
     &     - 16.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 8.D0/3*Hr1(0)*upx**(-2)
     &     + 8*Hr1(0)*upx**(-2)*y
     &     + 8.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 8.D0/3*Hr1(0)*upx**(-1)*y
     &     + 4.D0/3*DLog(tm**(-2)*w**2)
     &     + 8.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 8.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 8.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 8*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 8.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 8.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T13PoleEP0qqc = T13PoleEP0qqc + Nc**4 * (
     &     - 107.D0/9
     &     - 6*umx**(-5)*z2
     &     + 12*umx**(-5)*y*z2
     &     - 6*umx**(-5)*y**2*z2
     &     + 15*umx**(-4)*z2
     &     - 30*umx**(-4)*y*z2
     &     + 15*umx**(-4)*y**2*z2
     &     - 31.D0/2*umx**(-3)*z2
     &     + 27*umx**(-3)*y*z2
     &     - 19.D0/2*umx**(-3)*y**2*z2
     &     - umx**(-2)
     &     + 33.D0/4*umx**(-2)*z2
     &     + 2*umx**(-2)*y
     &     - 21.D0/2*umx**(-2)*y*z2
     &     - umx**(-2)*y**2
     &     - 3.D0/4*umx**(-2)*y**2*z2
     &     + umx**(-1)
     &     - 41.D0/8*umx**(-1)*z2
     &     - 2*umx**(-1)*y
     &     + 3.D0/4*umx**(-1)*y*z2
     &     + umx**(-1)*y**2
     &     - 1.D0/8*umx**(-1)*y**2*z2
     &     - 88.D0/9*upx**(-4)
     &     - 4*upx**(-4)*z2
     &     - 176.D0/9*upx**(-4)*y
     &     - 8*upx**(-4)*y*z2
     &     - 88.D0/9*upx**(-4)*y**2
     &     - 4*upx**(-4)*y**2*z2
     &     + 176.D0/9*upx**(-3)
     &     + 7*upx**(-3)*z2
     &     + 352.D0/9*upx**(-3)*y
     &     + 12*upx**(-3)*y*z2
     &     + 176.D0/9*upx**(-3)*y**2
     &     + 5*upx**(-3)*y**2*z2
     &     - 97.D0/9*upx**(-2)
     &     - 5.D0/2*upx**(-2)*z2
     &     - 94.D0/3*upx**(-2)*y
     &     - 4*upx**(-2)*y*z2
     &     - 97.D0/9*upx**(-2)*y**2
     &     + 1.D0/2*upx**(-2)*y**2*z2
     &     + upx**(-1)
     &     + 7.D0/8*upx**(-1)*z2
     &     + 106.D0/9*upx**(-1)*y
     &     + 3.D0/4*upx**(-1)*y*z2
     &     + upx**(-1)*y**2
     &     - 1.D0/8*upx**(-1)*y**2*z2
     &     + 20.D0/3*Hr1(-1)
     &     + 12*Hr1(-1)*umx**(-4)
     &     - 24*Hr1(-1)*umx**(-4)*y
     &     + 12*Hr1(-1)*umx**(-4)*y**2
     &     - 24*Hr1(-1)*umx**(-3)
     &     + 48*Hr1(-1)*umx**(-3)*y
     &     - 24*Hr1(-1)*umx**(-3)*y**2
     &     + 20*Hr1(-1)*umx**(-2)
     &     - 32*Hr1(-1)*umx**(-2)*y
     &     + 8*Hr1(-1)*umx**(-2)*y**2
     &     - 8*Hr1(-1)*umx**(-1)
     &     + 8*Hr1(-1)*umx**(-1)*y
     &     + 4*Hr1(-1)*umx**(-1)*y**2
     &     + 28.D0/3*Hr1(-1)*upx**(-4)
     &     + 56.D0/3*Hr1(-1)*upx**(-4)*y
     &     + 28.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     - 56.D0/3*Hr1(-1)*upx**(-3)
     &     - 112.D0/3*Hr1(-1)*upx**(-3)*y
     &     - 56.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     + 16.D0/3*Hr1(-1)*upx**(-2)
     &     + 20*Hr1(-1)*upx**(-2)*y
     &     + 16.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 4*Hr1(-1)*upx**(-1)
     &     - 4.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 4*Hr1(-1)*upx**(-1)*y**2
     &     + 4*Hr1(-1)*HAr1(-1)
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 32*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 32*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 32*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 40*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 8*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 2*HAr1(-1)
     &     - 2*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 4*HAr1(-1)*upx**(-2)
     &     - 2*HAr1(-1)*upx**(-2)*y
     &     + 2*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 4*HAr1(-1)*upx**(-1)
     &     + 2*HAr1(-1)*upx**(-1)*y
     &     + 2*HAr1(-1)*y**(-1)
     &     - 2*HAr1(-1)*Hr1(0)
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 16*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 20*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 4*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 4*Hr2(-1,-1)
     &     - 16*Hr2(-1,-1)*upx**(-4)
     &     - 32*Hr2(-1,-1)*upx**(-4)*y
     &     - 16*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 32*Hr2(-1,-1)*upx**(-3)
     &     + 64*Hr2(-1,-1)*upx**(-3)*y
     &     + 32*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 16*Hr2(-1,-1)*upx**(-2)
     &     - 40*Hr2(-1,-1)*upx**(-2)*y
     &     - 16*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 8*Hr2(-1,-1)*upx**(-1)*y
     &     + 4*HAr2(-1,-1)
     &     + 8*HAr2(-1,-1)*upx**(-2)*y
     &     - 8*HAr2(-1,-1)*upx**(-1)*y
     &     + 2*Hr2(-1,0)
     &     + 8*Hr2(-1,0)*upx**(-4)
     &     + 16*Hr2(-1,0)*upx**(-4)*y
     &     + 8*Hr2(-1,0)*upx**(-4)*y**2
     &     - 16*Hr2(-1,0)*upx**(-3)
     &     - 32*Hr2(-1,0)*upx**(-3)*y
     &     - 16*Hr2(-1,0)*upx**(-3)*y**2
     &     + 8*Hr2(-1,0)*upx**(-2)
     &     + 20*Hr2(-1,0)*upx**(-2)*y
     &     + 8*Hr2(-1,0)*upx**(-2)*y**2
     &     - 4*Hr2(-1,0)*upx**(-1)*y
     &     - 10.D0/3*Hr1(0)
     &     - 6*Hr1(0)*umx**(-4)
     &     + 12*Hr1(0)*umx**(-4)*y
     &     - 6*Hr1(0)*umx**(-4)*y**2
     &     + 12*Hr1(0)*umx**(-3)
     &     - 24*Hr1(0)*umx**(-3)*y
     &     + 12*Hr1(0)*umx**(-3)*y**2
     &     - 10*Hr1(0)*umx**(-2)
     &     + 16*Hr1(0)*umx**(-2)*y
     &     - 4*Hr1(0)*umx**(-2)*y**2
     &     + 4*Hr1(0)*umx**(-1)
     &     - 4*Hr1(0)*umx**(-1)*y
     &     - 2*Hr1(0)*umx**(-1)*y**2
     &     - 14.D0/3*Hr1(0)*upx**(-4)
     &     - 28.D0/3*Hr1(0)*upx**(-4)*y
     &     - 14.D0/3*Hr1(0)*upx**(-4)*y**2
     &     + 28.D0/3*Hr1(0)*upx**(-3)
     &     + 56.D0/3*Hr1(0)*upx**(-3)*y
     &     + 28.D0/3*Hr1(0)*upx**(-3)*y**2
     &     - 8.D0/3*Hr1(0)*upx**(-2)
     &     - 10*Hr1(0)*upx**(-2)*y
     &     - 8.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 2*Hr1(0)*upx**(-1)
     &     + 2.D0/3*Hr1(0)*upx**(-1)*y
     &     - 2*Hr1(0)*upx**(-1)*y**2
     &     + 12*Hr2(0,-1)*umx**(-5)
     &     - 24*Hr2(0,-1)*umx**(-5)*y
     &     + 12*Hr2(0,-1)*umx**(-5)*y**2
     &     - 30*Hr2(0,-1)*umx**(-4)
     &     + 60*Hr2(0,-1)*umx**(-4)*y
     &     - 30*Hr2(0,-1)*umx**(-4)*y**2
     &     + 31*Hr2(0,-1)*umx**(-3)
     &     - 54*Hr2(0,-1)*umx**(-3)*y
     &     + 19*Hr2(0,-1)*umx**(-3)*y**2
     &     - 33.D0/2*Hr2(0,-1)*umx**(-2)
     &     + 21*Hr2(0,-1)*umx**(-2)*y
     &     + 3.D0/2*Hr2(0,-1)*umx**(-2)*y**2
     &     + 41.D0/4*Hr2(0,-1)*umx**(-1)
     &     - 3.D0/2*Hr2(0,-1)*umx**(-1)*y
     &     + 1.D0/4*Hr2(0,-1)*umx**(-1)*y**2
     &     + 8*Hr2(0,-1)*upx**(-4)
     &     + 16*Hr2(0,-1)*upx**(-4)*y
     &     + 8*Hr2(0,-1)*upx**(-4)*y**2
     &     - 14*Hr2(0,-1)*upx**(-3)
     &     - 24*Hr2(0,-1)*upx**(-3)*y
     &     - 10*Hr2(0,-1)*upx**(-3)*y**2
     &     + 5*Hr2(0,-1)*upx**(-2)
     &     + 8*Hr2(0,-1)*upx**(-2)*y
     &     - Hr2(0,-1)*upx**(-2)*y**2
     &     - 7.D0/4*Hr2(0,-1)*upx**(-1)
     &     - 3.D0/2*Hr2(0,-1)*upx**(-1)*y
     &     + 1.D0/4*Hr2(0,-1)*upx**(-1)*y**2
     &     - 2*HAr2(0,-1)
     &     - 4*HAr2(0,-1)*upx**(-2)*y
     &     + 4*HAr2(0,-1)*upx**(-1)*y
     &     - 6*Hr2(0,0)*umx**(-5)
     &     + 12*Hr2(0,0)*umx**(-5)*y
     &     - 6*Hr2(0,0)*umx**(-5)*y**2
     &     + 15*Hr2(0,0)*umx**(-4)
     &     - 30*Hr2(0,0)*umx**(-4)*y
     &     + 15*Hr2(0,0)*umx**(-4)*y**2
     &     - 31.D0/2*Hr2(0,0)*umx**(-3)
     &     + 27*Hr2(0,0)*umx**(-3)*y
     &     - 19.D0/2*Hr2(0,0)*umx**(-3)*y**2
     &     + 33.D0/4*Hr2(0,0)*umx**(-2)
     &     - 21.D0/2*Hr2(0,0)*umx**(-2)*y
     &     - 3.D0/4*Hr2(0,0)*umx**(-2)*y**2
     &     - 41.D0/8*Hr2(0,0)*umx**(-1)
     &     + 3.D0/4*Hr2(0,0)*umx**(-1)*y
     &     - 1.D0/8*Hr2(0,0)*umx**(-1)*y**2
     &     - 4*Hr2(0,0)*upx**(-4)
     &     - 8*Hr2(0,0)*upx**(-4)*y
     &     - 4*Hr2(0,0)*upx**(-4)*y**2
     &     + 7*Hr2(0,0)*upx**(-3)
     &     + 12*Hr2(0,0)*upx**(-3)*y
     &     + 5*Hr2(0,0)*upx**(-3)*y**2
     &     - 5.D0/2*Hr2(0,0)*upx**(-2)
     &     - 4*Hr2(0,0)*upx**(-2)*y
     &     + 1.D0/2*Hr2(0,0)*upx**(-2)*y**2
     &     + 7.D0/8*Hr2(0,0)*upx**(-1)
     &     + 3.D0/4*Hr2(0,0)*upx**(-1)*y
     &     - 1.D0/8*Hr2(0,0)*upx**(-1)*y**2
     &     - 22.D0/3*DLog(tm**(-2)*w**2)
     &     - 44.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 44.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 44.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 44*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 44.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 44.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      return
      end

*
*******************************************************************************
************************************************************************
************************************************************************
*
*          FUNCTIONS for the Cross Section
*
*************************************************************************
**************************** 72 characters *****************************
*
*
      double precision function T23PoleEPm2qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T23PoleEPm2qq =  + Nc**(-2) * (
     &     + 16
     &     + 32*upx**(-4)
     &     + 64*upx**(-4)*y
     &     + 32*upx**(-4)*y**2
     &     - 64*upx**(-3)
     &     - 128*upx**(-3)*y
     &     - 64*upx**(-3)*y**2
     &     + 32*upx**(-2)
     &     + 96*upx**(-2)*y
     &     + 32*upx**(-2)*y**2
     &     - 32*upx**(-1)*y
     &     )/4.D0
      T23PoleEPm2qq = T23PoleEPm2qq + (
     &     - 32
     &     - 64*upx**(-4)
     &     - 128*upx**(-4)*y
     &     - 64*upx**(-4)*y**2
     &     + 128*upx**(-3)
     &     + 256*upx**(-3)*y
     &     + 128*upx**(-3)*y**2
     &     - 64*upx**(-2)
     &     - 192*upx**(-2)*y
     &     - 64*upx**(-2)*y**2
     &     + 64*upx**(-1)*y
     &     )/4.D0
      T23PoleEPm2qq = T23PoleEPm2qq + Nc**2 * (
     &     + 16
     &     + 32*upx**(-4)
     &     + 64*upx**(-4)*y
     &     + 32*upx**(-4)*y**2
     &     - 64*upx**(-3)
     &     - 128*upx**(-3)*y
     &     - 64*upx**(-3)*y**2
     &     + 32*upx**(-2)
     &     + 96*upx**(-2)*y
     &     + 32*upx**(-2)*y**2
     &     - 32*upx**(-1)*y
     &     )/4.D0
      return
      end
*
*
*
      double precision function T23PoleEPm1qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T23PoleEPm1qq =  + Nc**(-2) * (
     &     - 12
     &     - 40*upx**(-4)
     &     - 80*upx**(-4)*y
     &     - 40*upx**(-4)*y**2
     &     + 80*upx**(-3)
     &     + 160*upx**(-3)*y
     &     + 80*upx**(-3)*y**2
     &     - 40*upx**(-2)
     &     - 120*upx**(-2)*y
     &     - 40*upx**(-2)*y**2
     &     + 40*upx**(-1)*y
     &     + 16*Hr1(-1)
     &     + 32*Hr1(-1)*upx**(-4)
     &     + 64*Hr1(-1)*upx**(-4)*y
     &     + 32*Hr1(-1)*upx**(-4)*y**2
     &     - 64*Hr1(-1)*upx**(-3)
     &     - 128*Hr1(-1)*upx**(-3)*y
     &     - 64*Hr1(-1)*upx**(-3)*y**2
     &     + 32*Hr1(-1)*upx**(-2)
     &     + 96*Hr1(-1)*upx**(-2)*y
     &     + 32*Hr1(-1)*upx**(-2)*y**2
     &     - 32*Hr1(-1)*upx**(-1)*y
     &     + 24*HAr1(-1)
     &     + 48*HAr1(-1)*upx**(-4)
     &     + 96*HAr1(-1)*upx**(-4)*y
     &     + 48*HAr1(-1)*upx**(-4)*y**2
     &     - 96*HAr1(-1)*upx**(-3)
     &     - 192*HAr1(-1)*upx**(-3)*y
     &     - 96*HAr1(-1)*upx**(-3)*y**2
     &     + 48*HAr1(-1)*upx**(-2)
     &     + 144*HAr1(-1)*upx**(-2)*y
     &     + 48*HAr1(-1)*upx**(-2)*y**2
     &     - 48*HAr1(-1)*upx**(-1)*y
     &     - 24*HBr1(-1)
     &     - 48*HBr1(-1)*upx**(-4)
     &     - 96*HBr1(-1)*upx**(-4)*y
     &     - 48*HBr1(-1)*upx**(-4)*y**2
     &     + 96*HBr1(-1)*upx**(-3)
     &     + 192*HBr1(-1)*upx**(-3)*y
     &     + 96*HBr1(-1)*upx**(-3)*y**2
     &     - 48*HBr1(-1)*upx**(-2)
     &     - 144*HBr1(-1)*upx**(-2)*y
     &     - 48*HBr1(-1)*upx**(-2)*y**2
     &     + 48*HBr1(-1)*upx**(-1)*y
     &     - 9*Hr1(0)*umx**(-1)
     &     + 2*Hr1(0)*umx**(-1)*y
     &     - Hr1(0)*umx**(-1)*y**2
     &     - 16*Hr1(0)*upx**(-5)
     &     - 32*Hr1(0)*upx**(-5)*y
     &     - 16*Hr1(0)*upx**(-5)*y**2
     &     + 24*Hr1(0)*upx**(-4)
     &     + 48*Hr1(0)*upx**(-4)*y
     &     + 24*Hr1(0)*upx**(-4)*y**2
     &     - 4*Hr1(0)*upx**(-3)
     &     - 24*Hr1(0)*upx**(-3)*y
     &     - 4*Hr1(0)*upx**(-3)*y**2
     &     - 2*Hr1(0)*upx**(-2)
     &     + 4*Hr1(0)*upx**(-2)*y
     &     - 2*Hr1(0)*upx**(-2)*y**2
     &     - 9*Hr1(0)*upx**(-1)
     &     + 2*Hr1(0)*upx**(-1)*y
     &     - Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T23PoleEPm1qq = T23PoleEPm1qq + (
     &     + 24
     &     + 80*upx**(-4)
     &     + 160*upx**(-4)*y
     &     + 80*upx**(-4)*y**2
     &     - 160*upx**(-3)
     &     - 320*upx**(-3)*y
     &     - 160*upx**(-3)*y**2
     &     + 80*upx**(-2)
     &     + 240*upx**(-2)*y
     &     + 80*upx**(-2)*y**2
     &     - 80*upx**(-1)*y
     &     - 16*Hr1(-1)
     &     - 32*Hr1(-1)*upx**(-4)
     &     - 64*Hr1(-1)*upx**(-4)*y
     &     - 32*Hr1(-1)*upx**(-4)*y**2
     &     + 64*Hr1(-1)*upx**(-3)
     &     + 128*Hr1(-1)*upx**(-3)*y
     &     + 64*Hr1(-1)*upx**(-3)*y**2
     &     - 32*Hr1(-1)*upx**(-2)
     &     - 96*Hr1(-1)*upx**(-2)*y
     &     - 32*Hr1(-1)*upx**(-2)*y**2
     &     + 32*Hr1(-1)*upx**(-1)*y
     &     - 32*HAr1(-1)
     &     - 64*HAr1(-1)*upx**(-4)
     &     - 128*HAr1(-1)*upx**(-4)*y
     &     - 64*HAr1(-1)*upx**(-4)*y**2
     &     + 128*HAr1(-1)*upx**(-3)
     &     + 256*HAr1(-1)*upx**(-3)*y
     &     + 128*HAr1(-1)*upx**(-3)*y**2
     &     - 64*HAr1(-1)*upx**(-2)
     &     - 192*HAr1(-1)*upx**(-2)*y
     &     - 64*HAr1(-1)*upx**(-2)*y**2
     &     + 64*HAr1(-1)*upx**(-1)*y
     &     + 16*HBr1(-1)
     &     + 32*HBr1(-1)*upx**(-4)
     &     + 64*HBr1(-1)*upx**(-4)*y
     &     + 32*HBr1(-1)*upx**(-4)*y**2
     &     - 64*HBr1(-1)*upx**(-3)
     &     - 128*HBr1(-1)*upx**(-3)*y
     &     - 64*HBr1(-1)*upx**(-3)*y**2
     &     + 32*HBr1(-1)*upx**(-2)
     &     + 96*HBr1(-1)*upx**(-2)*y
     &     + 32*HBr1(-1)*upx**(-2)*y**2
     &     - 32*HBr1(-1)*upx**(-1)*y
     &     + 9*Hr1(0)*umx**(-1)
     &     - 2*Hr1(0)*umx**(-1)*y
     &     + Hr1(0)*umx**(-1)*y**2
     &     + 16*Hr1(0)*upx**(-5)
     &     + 32*Hr1(0)*upx**(-5)*y
     &     + 16*Hr1(0)*upx**(-5)*y**2
     &     - 24*Hr1(0)*upx**(-4)
     &     - 48*Hr1(0)*upx**(-4)*y
     &     - 24*Hr1(0)*upx**(-4)*y**2
     &     + 4*Hr1(0)*upx**(-3)
     &     + 24*Hr1(0)*upx**(-3)*y
     &     + 4*Hr1(0)*upx**(-3)*y**2
     &     + 2*Hr1(0)*upx**(-2)
     &     - 4*Hr1(0)*upx**(-2)*y
     &     + 2*Hr1(0)*upx**(-2)*y**2
     &     + 9*Hr1(0)*upx**(-1)
     &     - 2*Hr1(0)*upx**(-1)*y
     &     + Hr1(0)*upx**(-1)*y**2
     &     )/(-2.D0)
      T23PoleEPm1qq = T23PoleEPm1qq + Nc**2 * (
     &     - 12
     &     - 40*upx**(-4)
     &     - 80*upx**(-4)*y
     &     - 40*upx**(-4)*y**2
     &     + 80*upx**(-3)
     &     + 160*upx**(-3)*y
     &     + 80*upx**(-3)*y**2
     &     - 40*upx**(-2)
     &     - 120*upx**(-2)*y
     &     - 40*upx**(-2)*y**2
     &     + 40*upx**(-1)*y
     &     + 8*HAr1(-1)
     &     + 16*HAr1(-1)*upx**(-4)
     &     + 32*HAr1(-1)*upx**(-4)*y
     &     + 16*HAr1(-1)*upx**(-4)*y**2
     &     - 32*HAr1(-1)*upx**(-3)
     &     - 64*HAr1(-1)*upx**(-3)*y
     &     - 32*HAr1(-1)*upx**(-3)*y**2
     &     + 16*HAr1(-1)*upx**(-2)
     &     + 48*HAr1(-1)*upx**(-2)*y
     &     + 16*HAr1(-1)*upx**(-2)*y**2
     &     - 16*HAr1(-1)*upx**(-1)*y
     &     + 8*HBr1(-1)
     &     + 16*HBr1(-1)*upx**(-4)
     &     + 32*HBr1(-1)*upx**(-4)*y
     &     + 16*HBr1(-1)*upx**(-4)*y**2
     &     - 32*HBr1(-1)*upx**(-3)
     &     - 64*HBr1(-1)*upx**(-3)*y
     &     - 32*HBr1(-1)*upx**(-3)*y**2
     &     + 16*HBr1(-1)*upx**(-2)
     &     + 48*HBr1(-1)*upx**(-2)*y
     &     + 16*HBr1(-1)*upx**(-2)*y**2
     &     - 16*HBr1(-1)*upx**(-1)*y
     &     )/(-2.D0)
      return
      end
*
*
*
      double precision function T23PoleEP0qq(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1d0

* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      T23PoleEP0qq =  + Nc**(-2) * (
     &     + 14
     &     - 12*umx**(-3)*z2
     &     + 12*umx**(-3)*y*z2
     &     + 18*umx**(-2)*z2
     &     - 18*umx**(-2)*y*z2
     &     - 24*umx**(-1)*z2
     &     + 4*umx**(-1)*y*z2
     &     - 2*umx**(-1)*y**2*z2
     &     - 32*upx**(-5)*z2
     &     - 64*upx**(-5)*y*z2
     &     - 32*upx**(-5)*y**2*z2
     &     + 48*upx**(-4)
     &     + 48*upx**(-4)*z2
     &     + 96*upx**(-4)*y
     &     + 96*upx**(-4)*y*z2
     &     + 48*upx**(-4)*y**2
     &     + 48*upx**(-4)*y**2*z2
     &     - 96*upx**(-3)
     &     + 4*upx**(-3)*z2
     &     - 192*upx**(-3)*y
     &     - 36*upx**(-3)*y*z2
     &     - 96*upx**(-3)*y**2
     &     - 8*upx**(-3)*y**2*z2
     &     + 48*upx**(-2)
     &     - 10*upx**(-2)*z2
     &     + 144*upx**(-2)*y
     &     + 2*upx**(-2)*y*z2
     &     + 48*upx**(-2)*y**2
     &     - 4*upx**(-2)*y**2*z2
     &     - 12*upx**(-1)*z2
     &     - 48*upx**(-1)*y
     &     + 4*upx**(-1)*y*z2
     &     - 2*upx**(-1)*y**2*z2
     &     + 8*Hr1(-1)
     &     + 24*Hr1(-1)*umx**(-2)
     &     - 24*Hr1(-1)*umx**(-2)*y
     &     - 24*Hr1(-1)*umx**(-1)
     &     + 24*Hr1(-1)*umx**(-1)*y
     &     - 24*Hr1(-1)*upx**(-4)
     &     - 48*Hr1(-1)*upx**(-4)*y
     &     - 24*Hr1(-1)*upx**(-4)*y**2
     &     + 48*Hr1(-1)*upx**(-3)
     &     + 96*Hr1(-1)*upx**(-3)*y
     &     + 48*Hr1(-1)*upx**(-3)*y**2
     &     - 24*Hr1(-1)*upx**(-2)
     &     - 72*Hr1(-1)*upx**(-2)*y
     &     - 24*Hr1(-1)*upx**(-2)*y**2
     &     + 24*Hr1(-1)*upx**(-1)*y
     &     + 12*Hr1(-1)*HAr1(-1)
     &     + 48*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 96*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 48*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 96*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 192*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 96*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 48*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 120*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 48*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 24*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     - 36*Hr1(-1)*HBr1(-1)
     &     - 48*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     - 96*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     - 48*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     + 96*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     + 192*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     + 96*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     - 96*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     - 168*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     - 48*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     + 48*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     + 72*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     + 6*HAr1(-1)
     &     - 6*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 12*HAr1(-1)*upx**(-2)
     &     - 6*HAr1(-1)*upx**(-2)*y
     &     + 6*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 12*HAr1(-1)*upx**(-1)
     &     + 6*HAr1(-1)*upx**(-1)*y
     &     + 6*HAr1(-1)*y**(-1)
     &     - 6*HAr1(-1)*Hr1(0)
     &     - 24*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 48*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 24*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 48*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 96*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 48*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 24*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 60*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 24*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 12*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 12*HBr1(-1)
     &     + 6*HBr1(-1)*upx**(-2)*z**(-1)
     &     - 6*HBr1(-1)*upx**(-2)*y
     &     - 6*HBr1(-1)*upx**(-1)*z**(-1)
     &     + 6*HBr1(-1)*upx**(-1)*y
     &     - 6*HBr1(-1)*z**(-1)
     &     + 18*HBr1(-1)*Hr1(0)
     &     + 24*HBr1(-1)*Hr1(0)*upx**(-4)
     &     + 48*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 24*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 48*HBr1(-1)*Hr1(0)*upx**(-3)
     &     - 96*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 48*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 48*HBr1(-1)*Hr1(0)*upx**(-2)
     &     + 84*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 24*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 24*HBr1(-1)*Hr1(0)*upx**(-1)
     &     - 36*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 40*Hr2(-1,-1)
     &     + 32*Hr2(-1,-1)*upx**(-4)
     &     + 64*Hr2(-1,-1)*upx**(-4)*y
     &     + 32*Hr2(-1,-1)*upx**(-4)*y**2
     &     - 64*Hr2(-1,-1)*upx**(-3)
     &     - 128*Hr2(-1,-1)*upx**(-3)*y
     &     - 64*Hr2(-1,-1)*upx**(-3)*y**2
     &     + 80*Hr2(-1,-1)*upx**(-2)
     &     + 144*Hr2(-1,-1)*upx**(-2)*y
     &     + 32*Hr2(-1,-1)*upx**(-2)*y**2
     &     - 48*Hr2(-1,-1)*upx**(-1)
     &     - 80*Hr2(-1,-1)*upx**(-1)*y
     &     + 12*HAr2(-1,-1)
     &     + 24*HAr2(-1,-1)*upx**(-2)*y
     &     - 24*HAr2(-1,-1)*upx**(-1)*y
     &     + 12*HBr2(-1,-1)
     &     + 48*HBr2(-1,-1)*upx**(-2)
     &     + 24*HBr2(-1,-1)*upx**(-2)*y
     &     - 48*HBr2(-1,-1)*upx**(-1)
     &     - 24*HBr2(-1,-1)*upx**(-1)*y
     &     - 20*Hr2(-1,0)
     &     - 16*Hr2(-1,0)*upx**(-4)
     &     - 32*Hr2(-1,0)*upx**(-4)*y
     &     - 16*Hr2(-1,0)*upx**(-4)*y**2
     &     + 32*Hr2(-1,0)*upx**(-3)
     &     + 64*Hr2(-1,0)*upx**(-3)*y
     &     + 32*Hr2(-1,0)*upx**(-3)*y**2
     &     - 40*Hr2(-1,0)*upx**(-2)
     &     - 72*Hr2(-1,0)*upx**(-2)*y
     &     - 16*Hr2(-1,0)*upx**(-2)*y**2
     &     + 24*Hr2(-1,0)*upx**(-1)
     &     + 40*Hr2(-1,0)*upx**(-1)*y
     &     - 6*Hr1(0)
     &     - 12*Hr1(0)*umx**(-2)
     &     + 12*Hr1(0)*umx**(-2)*y
     &     + 17.D0/2*Hr1(0)*umx**(-1)
     &     - 13*Hr1(0)*umx**(-1)*y
     &     + 1.D0/2*Hr1(0)*umx**(-1)*y**2
     &     + 16*Hr1(0)*upx**(-5)
     &     + 32*Hr1(0)*upx**(-5)*y
     &     + 16*Hr1(0)*upx**(-5)*y**2
     &     - 28*Hr1(0)*upx**(-4)
     &     - 56*Hr1(0)*upx**(-4)*y
     &     - 28*Hr1(0)*upx**(-4)*y**2
     &     + 10*Hr1(0)*upx**(-3)
     &     + 36*Hr1(0)*upx**(-3)*y
     &     + 10*Hr1(0)*upx**(-3)*y**2
     &     + Hr1(0)*upx**(-2)
     &     - 10*Hr1(0)*upx**(-2)*y
     &     + Hr1(0)*upx**(-2)*y**2
     &     + 17.D0/2*Hr1(0)*upx**(-1)
     &     - Hr1(0)*upx**(-1)*y
     &     + 1.D0/2*Hr1(0)*upx**(-1)*y**2
     &     - 8*Hr2(0,-1)
     &     + 24*Hr2(0,-1)*umx**(-3)
     &     - 24*Hr2(0,-1)*umx**(-3)*y
     &     - 36*Hr2(0,-1)*umx**(-2)
     &     + 36*Hr2(0,-1)*umx**(-2)*y
     &     + 12*Hr2(0,-1)*umx**(-1)
     &     - 16*Hr2(0,-1)*upx**(-4)
     &     - 32*Hr2(0,-1)*upx**(-4)*y
     &     - 16*Hr2(0,-1)*upx**(-4)*y**2
     &     + 8*Hr2(0,-1)*upx**(-3)
     &     + 40*Hr2(0,-1)*upx**(-3)*y
     &     + 32*Hr2(0,-1)*upx**(-3)*y**2
     &     - 4*Hr2(0,-1)*upx**(-2)
     &     - 36*Hr2(0,-1)*upx**(-2)*y
     &     - 16*Hr2(0,-1)*upx**(-2)*y**2
     &     - 12*Hr2(0,-1)*upx**(-1)
     &     + 16*Hr2(0,-1)*upx**(-1)*y
     &     - 6*HAr2(0,-1)
     &     - 12*HAr2(0,-1)*upx**(-2)*y
     &     + 12*HAr2(0,-1)*upx**(-1)*y
     &     - 6*HBr2(0,-1)
     &     - 24*HBr2(0,-1)*upx**(-2)
     &     - 12*HBr2(0,-1)*upx**(-2)*y
     &     + 24*HBr2(0,-1)*upx**(-1)
     &     + 12*HBr2(0,-1)*upx**(-1)*y
     &     - 12*Hr2(0,0)*umx**(-3)
     &     + 12*Hr2(0,0)*umx**(-3)*y
     &     + 18*Hr2(0,0)*umx**(-2)
     &     - 18*Hr2(0,0)*umx**(-2)*y
     &     - 3.D0/2*Hr2(0,0)*umx**(-1)
     &     - Hr2(0,0)*umx**(-1)*y
     &     + 1.D0/2*Hr2(0,0)*umx**(-1)*y**2
     &     + 8*Hr2(0,0)*upx**(-5)
     &     + 16*Hr2(0,0)*upx**(-5)*y
     &     + 8*Hr2(0,0)*upx**(-5)*y**2
     &     - 12*Hr2(0,0)*upx**(-4)
     &     - 24*Hr2(0,0)*upx**(-4)*y
     &     - 12*Hr2(0,0)*upx**(-4)*y**2
     &     + 14*Hr2(0,0)*upx**(-3)
     &     + 24*Hr2(0,0)*upx**(-3)*y
     &     + 2*Hr2(0,0)*upx**(-3)*y**2
     &     - 5*Hr2(0,0)*upx**(-2)
     &     - 8*Hr2(0,0)*upx**(-2)*y
     &     + Hr2(0,0)*upx**(-2)*y**2
     &     + 21.D0/2*Hr2(0,0)*upx**(-1)
     &     - Hr2(0,0)*upx**(-1)*y
     &     + 1.D0/2*Hr2(0,0)*upx**(-1)*y**2
     &     - 8*Hr2(1,0)
     &     + 9*Hr2(1,0)*umx**(-1)
     &     - 2*Hr2(1,0)*umx**(-1)*y
     &     + Hr2(1,0)*umx**(-1)*y**2
     &     + 16*Hr2(1,0)*upx**(-5)
     &     + 32*Hr2(1,0)*upx**(-5)*y
     &     + 16*Hr2(1,0)*upx**(-5)*y**2
     &     - 40*Hr2(1,0)*upx**(-4)
     &     - 80*Hr2(1,0)*upx**(-4)*y
     &     - 40*Hr2(1,0)*upx**(-4)*y**2
     &     + 36*Hr2(1,0)*upx**(-3)
     &     + 88*Hr2(1,0)*upx**(-3)*y
     &     + 36*Hr2(1,0)*upx**(-3)*y**2
     &     - 14*Hr2(1,0)*upx**(-2)
     &     - 52*Hr2(1,0)*upx**(-2)*y
     &     - 14*Hr2(1,0)*upx**(-2)*y**2
     &     + 9*Hr2(1,0)*upx**(-1)
     &     + 14*Hr2(1,0)*upx**(-1)*y
     &     + Hr2(1,0)*upx**(-1)*y**2
     &     )
      T23PoleEP0qq = T23PoleEP0qq + Nc**(-1)*Nh * (
     &     - 40.D0/9
     &     + 64.D0/3*upx**(-6)
     &     + 128.D0/3*upx**(-6)*y
     &     + 64.D0/3*upx**(-6)*y**2
     &     - 64*upx**(-5)
     &     - 128*upx**(-5)*y
     &     - 64*upx**(-5)*y**2
     &     + 496.D0/9*upx**(-4)
     &     + 1184.D0/9*upx**(-4)*y
     &     + 496.D0/9*upx**(-4)*y**2
     &     - 32.D0/9*upx**(-3)
     &     - 448.D0/9*upx**(-3)*y
     &     - 32.D0/9*upx**(-3)*y**2
     &     + 16.D0/9*upx**(-2)
     &     - 16.D0/3*upx**(-2)*y
     &     - 80.D0/9*upx**(-2)*y**2
     &     - 32.D0/3*upx**(-1)
     &     + 80.D0/9*upx**(-1)*y
     &     + 8.D0/3*Hr1(0)
     &     + 64.D0/3*Hr1(0)*upx**(-7)
     &     + 128.D0/3*Hr1(0)*upx**(-7)*y
     &     + 64.D0/3*Hr1(0)*upx**(-7)*y**2
     &     - 224.D0/3*Hr1(0)*upx**(-6)
     &     - 448.D0/3*Hr1(0)*upx**(-6)*y
     &     - 224.D0/3*Hr1(0)*upx**(-6)*y**2
     &     + 256.D0/3*Hr1(0)*upx**(-5)
     &     + 192*Hr1(0)*upx**(-5)*y
     &     + 256.D0/3*Hr1(0)*upx**(-5)*y**2
     &     - 80.D0/3*Hr1(0)*upx**(-4)
     &     - 320.D0/3*Hr1(0)*upx**(-4)*y
     &     - 80.D0/3*Hr1(0)*upx**(-4)*y**2
     &     + 32.D0/3*Hr1(0)*upx**(-3)*y
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     - 32.D0/3*Hr1(0)*upx**(-2)
     &     + 16*Hr1(0)*upx**(-2)*y
     &     + 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-1)*y
     &     )
      T23PoleEP0qq = T23PoleEP0qq + Nc**(-1)*Nl * (
     &     - 40.D0/9
     &     - 80.D0/9*upx**(-4)
     &     - 160.D0/9*upx**(-4)*y
     &     - 80.D0/9*upx**(-4)*y**2
     &     + 160.D0/9*upx**(-3)
     &     + 320.D0/9*upx**(-3)*y
     &     + 160.D0/9*upx**(-3)*y**2
     &     - 80.D0/9*upx**(-2)
     &     - 80.D0/3*upx**(-2)*y
     &     - 80.D0/9*upx**(-2)*y**2
     &     + 80.D0/9*upx**(-1)*y
     &     + 16.D0/3*Hr1(-1)
     &     + 32.D0/3*Hr1(-1)*upx**(-4)
     &     + 64.D0/3*Hr1(-1)*upx**(-4)*y
     &     + 32.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     - 64.D0/3*Hr1(-1)*upx**(-3)
     &     - 128.D0/3*Hr1(-1)*upx**(-3)*y
     &     - 64.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     + 32.D0/3*Hr1(-1)*upx**(-2)
     &     + 32*Hr1(-1)*upx**(-2)*y
     &     + 32.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     - 32.D0/3*Hr1(-1)*upx**(-1)*y
     &     - 8.D0/3*Hr1(0)
     &     - 16.D0/3*Hr1(0)*upx**(-4)
     &     - 32.D0/3*Hr1(0)*upx**(-4)*y
     &     - 16.D0/3*Hr1(0)*upx**(-4)*y**2
     &     + 32.D0/3*Hr1(0)*upx**(-3)
     &     + 64.D0/3*Hr1(0)*upx**(-3)*y
     &     + 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-2)
     &     - 16*Hr1(0)*upx**(-2)*y
     &     - 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-1)*y
     &     - 8.D0/3*DLog(tm**(-2)*w**2)
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 64.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 16*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T23PoleEP0qq = T23PoleEP0qq + (
     &     + 88.D0/9
     &     + 12*umx**(-5)*z2
     &     - 24*umx**(-5)*y*z2
     &     + 12*umx**(-5)*y**2*z2
     &     - 30*umx**(-4)*z2
     &     + 60*umx**(-4)*y*z2
     &     - 30*umx**(-4)*y**2*z2
     &     + 39*umx**(-3)*z2
     &     - 62*umx**(-3)*y*z2
     &     + 19*umx**(-3)*y**2*z2
     &     + 2*umx**(-2)
     &     - 57.D0/2*umx**(-2)*z2
     &     - 4*umx**(-2)*y
     &     + 33*umx**(-2)*y*z2
     &     + 2*umx**(-2)*y**2
     &     + 3.D0/2*umx**(-2)*y**2*z2
     &     - 2*umx**(-1)
     &     + 129.D0/4*umx**(-1)*z2
     &     + 4*umx**(-1)*y
     &     - 11.D0/2*umx**(-1)*y*z2
     &     - 2*umx**(-1)*y**2
     &     + 9.D0/4*umx**(-1)*y**2*z2
     &     + 32*upx**(-5)*z2
     &     + 64*upx**(-5)*y*z2
     &     + 32*upx**(-5)*y**2*z2
     &     - 256.D0/9*upx**(-4)
     &     - 40*upx**(-4)*z2
     &     - 512.D0/9*upx**(-4)*y
     &     - 80*upx**(-4)*y*z2
     &     - 256.D0/9*upx**(-4)*y**2
     &     - 40*upx**(-4)*y**2*z2
     &     + 512.D0/9*upx**(-3)
     &     - 14*upx**(-3)*z2
     &     + 1024.D0/9*upx**(-3)*y
     &     + 16*upx**(-3)*y*z2
     &     + 512.D0/9*upx**(-3)*y**2
     &     - 2*upx**(-3)*y**2*z2
     &     - 238.D0/9*upx**(-2)
     &     + 13*upx**(-2)*z2
     &     - 244.D0/3*upx**(-2)*y
     &     + 4*upx**(-2)*y*z2
     &     - 238.D0/9*upx**(-2)*y**2
     &     + 3*upx**(-2)*y**2*z2
     &     - 2*upx**(-1)
     &     + 49.D0/4*upx**(-1)*z2
     &     + 220.D0/9*upx**(-1)*y
     &     - 11.D0/2*upx**(-1)*y*z2
     &     - 2*upx**(-1)*y**2
     &     + 9.D0/4*upx**(-1)*y**2*z2
     &     - 52.D0/3*Hr1(-1)
     &     - 24*Hr1(-1)*umx**(-4)
     &     + 48*Hr1(-1)*umx**(-4)*y
     &     - 24*Hr1(-1)*umx**(-4)*y**2
     &     + 48*Hr1(-1)*umx**(-3)
     &     - 96*Hr1(-1)*umx**(-3)*y
     &     + 48*Hr1(-1)*umx**(-3)*y**2
     &     - 56*Hr1(-1)*umx**(-2)
     &     + 80*Hr1(-1)*umx**(-2)*y
     &     - 16*Hr1(-1)*umx**(-2)*y**2
     &     + 32*Hr1(-1)*umx**(-1)
     &     - 32*Hr1(-1)*umx**(-1)*y
     &     - 8*Hr1(-1)*umx**(-1)*y**2
     &     + 16.D0/3*Hr1(-1)*upx**(-4)
     &     + 32.D0/3*Hr1(-1)*upx**(-4)*y
     &     + 16.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     - 32.D0/3*Hr1(-1)*upx**(-3)
     &     - 64.D0/3*Hr1(-1)*upx**(-3)*y
     &     - 32.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     + 40.D0/3*Hr1(-1)*upx**(-2)
     &     + 32*Hr1(-1)*upx**(-2)*y
     &     + 40.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     - 8*Hr1(-1)*upx**(-1)
     &     - 64.D0/3*Hr1(-1)*upx**(-1)*y
     &     - 8*Hr1(-1)*upx**(-1)*y**2
     &     - 16*Hr1(-1)*HAr1(-1)
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     - 128*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     + 128*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     + 256*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     + 128*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     - 160*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     + 32*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 24*Hr1(-1)*HBr1(-1)
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 64*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 128*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 64*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 112*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 32*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 48*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     - 8*HAr1(-1)
     &     + 8*HAr1(-1)*upx**(-2)*y**(-1)
     &     + 16*HAr1(-1)*upx**(-2)
     &     + 8*HAr1(-1)*upx**(-2)*y
     &     - 8*HAr1(-1)*upx**(-1)*y**(-1)
     &     - 16*HAr1(-1)*upx**(-1)
     &     - 8*HAr1(-1)*upx**(-1)*y
     &     - 8*HAr1(-1)*y**(-1)
     &     + 8*HAr1(-1)*Hr1(0)
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-4)
     &     + 64*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     - 64*HAr1(-1)*Hr1(0)*upx**(-3)
     &     - 128*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     - 64*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-2)
     &     + 80*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     - 16*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 8*HBr1(-1)
     &     - 4*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 4*HBr1(-1)*upx**(-2)*y
     &     + 4*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 4*HBr1(-1)*upx**(-1)*y
     &     + 4*HBr1(-1)*z**(-1)
     &     - 12*HBr1(-1)*Hr1(0)
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 32*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 64*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 32*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 56*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 16*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 24*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 24*Hr2(-1,-1)
     &     - 32*Hr2(-1,-1)*upx**(-2)
     &     - 48*Hr2(-1,-1)*upx**(-2)*y
     &     + 32*Hr2(-1,-1)*upx**(-1)
     &     + 48*Hr2(-1,-1)*upx**(-1)*y
     &     - 16*HAr2(-1,-1)
     &     - 32*HAr2(-1,-1)*upx**(-2)*y
     &     + 32*HAr2(-1,-1)*upx**(-1)*y
     &     - 8*HBr2(-1,-1)
     &     - 32*HBr2(-1,-1)*upx**(-2)
     &     - 16*HBr2(-1,-1)*upx**(-2)*y
     &     + 32*HBr2(-1,-1)*upx**(-1)
     &     + 16*HBr2(-1,-1)*upx**(-1)*y
     &     + 12*Hr2(-1,0)
     &     + 16*Hr2(-1,0)*upx**(-2)
     &     + 24*Hr2(-1,0)*upx**(-2)*y
     &     - 16*Hr2(-1,0)*upx**(-1)
     &     - 24*Hr2(-1,0)*upx**(-1)*y
     &     + 32.D0/3*Hr1(0)
     &     + 12*Hr1(0)*umx**(-4)
     &     - 24*Hr1(0)*umx**(-4)*y
     &     + 12*Hr1(0)*umx**(-4)*y**2
     &     - 24*Hr1(0)*umx**(-3)
     &     + 48*Hr1(0)*umx**(-3)*y
     &     - 24*Hr1(0)*umx**(-3)*y**2
     &     + 28*Hr1(0)*umx**(-2)
     &     - 40*Hr1(0)*umx**(-2)*y
     &     + 8*Hr1(0)*umx**(-2)*y**2
     &     - 25.D0/2*Hr1(0)*umx**(-1)
     &     + 17*Hr1(0)*umx**(-1)*y
     &     + 7.D0/2*Hr1(0)*umx**(-1)*y**2
     &     - 16*Hr1(0)*upx**(-5)
     &     - 32*Hr1(0)*upx**(-5)*y
     &     - 16*Hr1(0)*upx**(-5)*y**2
     &     + 112.D0/3*Hr1(0)*upx**(-4)
     &     + 224.D0/3*Hr1(0)*upx**(-4)*y
     &     + 112.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 86.D0/3*Hr1(0)*upx**(-3)
     &     - 220.D0/3*Hr1(0)*upx**(-3)*y
     &     - 86.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 13.D0/3*Hr1(0)*upx**(-2)
     &     + 30*Hr1(0)*upx**(-2)*y
     &     + 13.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 9.D0/2*Hr1(0)*upx**(-1)
     &     - 1.D0/3*Hr1(0)*upx**(-1)*y
     &     + 7.D0/2*Hr1(0)*upx**(-1)*y**2
     &     + 8*Hr2(0,-1)
     &     - 24*Hr2(0,-1)*umx**(-5)
     &     + 48*Hr2(0,-1)*umx**(-5)*y
     &     - 24*Hr2(0,-1)*umx**(-5)*y**2
     &     + 60*Hr2(0,-1)*umx**(-4)
     &     - 120*Hr2(0,-1)*umx**(-4)*y
     &     + 60*Hr2(0,-1)*umx**(-4)*y**2
     &     - 78*Hr2(0,-1)*umx**(-3)
     &     + 124*Hr2(0,-1)*umx**(-3)*y
     &     - 38*Hr2(0,-1)*umx**(-3)*y**2
     &     + 57*Hr2(0,-1)*umx**(-2)
     &     - 66*Hr2(0,-1)*umx**(-2)*y
     &     - 3*Hr2(0,-1)*umx**(-2)*y**2
     &     - 57.D0/2*Hr2(0,-1)*umx**(-1)
     &     + 3*Hr2(0,-1)*umx**(-1)*y
     &     - 1.D0/2*Hr2(0,-1)*umx**(-1)*y**2
     &     + 12*Hr2(0,-1)*upx**(-3)
     &     - 12*Hr2(0,-1)*upx**(-3)*y**2
     &     - 2*Hr2(0,-1)*upx**(-2)
     &     + 24*Hr2(0,-1)*upx**(-2)*y
     &     + 18*Hr2(0,-1)*upx**(-2)*y**2
     &     + 23.D0/2*Hr2(0,-1)*upx**(-1)
     &     - 13*Hr2(0,-1)*upx**(-1)*y
     &     - 1.D0/2*Hr2(0,-1)*upx**(-1)*y**2
     &     + 8*HAr2(0,-1)
     &     + 16*HAr2(0,-1)*upx**(-2)*y
     &     - 16*HAr2(0,-1)*upx**(-1)*y
     &     + 4*HBr2(0,-1)
     &     + 16*HBr2(0,-1)*upx**(-2)
     &     + 8*HBr2(0,-1)*upx**(-2)*y
     &     - 16*HBr2(0,-1)*upx**(-1)
     &     - 8*HBr2(0,-1)*upx**(-1)*y
     &     + 12*Hr2(0,0)*umx**(-5)
     &     - 24*Hr2(0,0)*umx**(-5)*y
     &     + 12*Hr2(0,0)*umx**(-5)*y**2
     &     - 30*Hr2(0,0)*umx**(-4)
     &     + 60*Hr2(0,0)*umx**(-4)*y
     &     - 30*Hr2(0,0)*umx**(-4)*y**2
     &     + 39*Hr2(0,0)*umx**(-3)
     &     - 62*Hr2(0,0)*umx**(-3)*y
     &     + 19*Hr2(0,0)*umx**(-3)*y**2
     &     - 57.D0/2*Hr2(0,0)*umx**(-2)
     &     + 33*Hr2(0,0)*umx**(-2)*y
     &     + 3.D0/2*Hr2(0,0)*umx**(-2)*y**2
     &     + 39.D0/4*Hr2(0,0)*umx**(-1)
     &     - 1.D0/2*Hr2(0,0)*umx**(-1)*y
     &     - 1.D0/4*Hr2(0,0)*umx**(-1)*y**2
     &     - 8*Hr2(0,0)*upx**(-5)
     &     - 16*Hr2(0,0)*upx**(-5)*y
     &     - 8*Hr2(0,0)*upx**(-5)*y**2
     &     + 20*Hr2(0,0)*upx**(-4)
     &     + 40*Hr2(0,0)*upx**(-4)*y
     &     + 20*Hr2(0,0)*upx**(-4)*y**2
     &     - 24*Hr2(0,0)*upx**(-3)
     &     - 44*Hr2(0,0)*upx**(-3)*y
     &     - 12*Hr2(0,0)*upx**(-3)*y**2
     &     + 8*Hr2(0,0)*upx**(-2)
     &     + 14*Hr2(0,0)*upx**(-2)*y
     &     - 2*Hr2(0,0)*upx**(-2)*y**2
     &     - 41.D0/4*Hr2(0,0)*upx**(-1)
     &     - 1.D0/2*Hr2(0,0)*upx**(-1)*y
     &     - 1.D0/4*Hr2(0,0)*upx**(-1)*y**2
     &     + 8*Hr2(1,0)
     &     - 9*Hr2(1,0)*umx**(-1)
     &     + 2*Hr2(1,0)*umx**(-1)*y
     &     - Hr2(1,0)*umx**(-1)*y**2
     &     - 16*Hr2(1,0)*upx**(-5)
     &     - 32*Hr2(1,0)*upx**(-5)*y
     &     - 16*Hr2(1,0)*upx**(-5)*y**2
     &     + 40*Hr2(1,0)*upx**(-4)
     &     + 80*Hr2(1,0)*upx**(-4)*y
     &     + 40*Hr2(1,0)*upx**(-4)*y**2
     &     - 36*Hr2(1,0)*upx**(-3)
     &     - 88*Hr2(1,0)*upx**(-3)*y
     &     - 36*Hr2(1,0)*upx**(-3)*y**2
     &     + 14*Hr2(1,0)*upx**(-2)
     &     + 52*Hr2(1,0)*upx**(-2)*y
     &     + 14*Hr2(1,0)*upx**(-2)*y**2
     &     - 9*Hr2(1,0)*upx**(-1)
     &     - 14*Hr2(1,0)*upx**(-1)*y
     &     - Hr2(1,0)*upx**(-1)*y**2
     &     + 44.D0/3*DLog(tm**(-2)*w**2)
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 352.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 88*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T23PoleEP0qq = T23PoleEP0qq + Nc*Nh * (
     &     + 40.D0/9
     &     - 64.D0/3*upx**(-6)
     &     - 128.D0/3*upx**(-6)*y
     &     - 64.D0/3*upx**(-6)*y**2
     &     + 64*upx**(-5)
     &     + 128*upx**(-5)*y
     &     + 64*upx**(-5)*y**2
     &     - 496.D0/9*upx**(-4)
     &     - 1184.D0/9*upx**(-4)*y
     &     - 496.D0/9*upx**(-4)*y**2
     &     + 32.D0/9*upx**(-3)
     &     + 448.D0/9*upx**(-3)*y
     &     + 32.D0/9*upx**(-3)*y**2
     &     - 16.D0/9*upx**(-2)
     &     + 16.D0/3*upx**(-2)*y
     &     + 80.D0/9*upx**(-2)*y**2
     &     + 32.D0/3*upx**(-1)
     &     - 80.D0/9*upx**(-1)*y
     &     - 8.D0/3*Hr1(0)
     &     - 64.D0/3*Hr1(0)*upx**(-7)
     &     - 128.D0/3*Hr1(0)*upx**(-7)*y
     &     - 64.D0/3*Hr1(0)*upx**(-7)*y**2
     &     + 224.D0/3*Hr1(0)*upx**(-6)
     &     + 448.D0/3*Hr1(0)*upx**(-6)*y
     &     + 224.D0/3*Hr1(0)*upx**(-6)*y**2
     &     - 256.D0/3*Hr1(0)*upx**(-5)
     &     - 192*Hr1(0)*upx**(-5)*y
     &     - 256.D0/3*Hr1(0)*upx**(-5)*y**2
     &     + 80.D0/3*Hr1(0)*upx**(-4)
     &     + 320.D0/3*Hr1(0)*upx**(-4)*y
     &     + 80.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y
     &     + 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 32.D0/3*Hr1(0)*upx**(-2)
     &     - 16*Hr1(0)*upx**(-2)*y
     &     - 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-1)*y
     &     )
      T23PoleEP0qq = T23PoleEP0qq + Nc*Nl * (
     &     + 40.D0/9
     &     + 80.D0/9*upx**(-4)
     &     + 160.D0/9*upx**(-4)*y
     &     + 80.D0/9*upx**(-4)*y**2
     &     - 160.D0/9*upx**(-3)
     &     - 320.D0/9*upx**(-3)*y
     &     - 160.D0/9*upx**(-3)*y**2
     &     + 80.D0/9*upx**(-2)
     &     + 80.D0/3*upx**(-2)*y
     &     + 80.D0/9*upx**(-2)*y**2
     &     - 80.D0/9*upx**(-1)*y
     &     - 16.D0/3*Hr1(-1)
     &     - 32.D0/3*Hr1(-1)*upx**(-4)
     &     - 64.D0/3*Hr1(-1)*upx**(-4)*y
     &     - 32.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     + 64.D0/3*Hr1(-1)*upx**(-3)
     &     + 128.D0/3*Hr1(-1)*upx**(-3)*y
     &     + 64.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     - 32.D0/3*Hr1(-1)*upx**(-2)
     &     - 32*Hr1(-1)*upx**(-2)*y
     &     - 32.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 32.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 8.D0/3*Hr1(0)
     &     + 16.D0/3*Hr1(0)*upx**(-4)
     &     + 32.D0/3*Hr1(0)*upx**(-4)*y
     &     + 16.D0/3*Hr1(0)*upx**(-4)*y**2
     &     - 32.D0/3*Hr1(0)*upx**(-3)
     &     - 64.D0/3*Hr1(0)*upx**(-3)*y
     &     - 32.D0/3*Hr1(0)*upx**(-3)*y**2
     &     + 16.D0/3*Hr1(0)*upx**(-2)
     &     + 16*Hr1(0)*upx**(-2)*y
     &     + 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-1)*y
     &     + 8.D0/3*DLog(tm**(-2)*w**2)
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     + 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     - 64.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     - 32.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     + 16*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     + 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     - 16.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      T23PoleEP0qq = T23PoleEP0qq + Nc**2 * (
     &     - 214.D0/9
     &     - 12*umx**(-5)*z2
     &     + 24*umx**(-5)*y*z2
     &     - 12*umx**(-5)*y**2*z2
     &     + 30*umx**(-4)*z2
     &     - 60*umx**(-4)*y*z2
     &     + 30*umx**(-4)*y**2*z2
     &     - 27*umx**(-3)*z2
     &     + 50*umx**(-3)*y*z2
     &     - 19*umx**(-3)*y**2*z2
     &     - 2*umx**(-2)
     &     + 21.D0/2*umx**(-2)*z2
     &     + 4*umx**(-2)*y
     &     - 15*umx**(-2)*y*z2
     &     - 2*umx**(-2)*y**2
     &     - 3.D0/2*umx**(-2)*y**2*z2
     &     + 2*umx**(-1)
     &     - 33.D0/4*umx**(-1)*z2
     &     - 4*umx**(-1)*y
     &     + 3.D0/2*umx**(-1)*y*z2
     &     + 2*umx**(-1)*y**2
     &     - 1.D0/4*umx**(-1)*y**2*z2
     &     - 176.D0/9*upx**(-4)
     &     - 8*upx**(-4)*z2
     &     - 352.D0/9*upx**(-4)*y
     &     - 16*upx**(-4)*y*z2
     &     - 176.D0/9*upx**(-4)*y**2
     &     - 8*upx**(-4)*y**2*z2
     &     + 352.D0/9*upx**(-3)
     &     + 10*upx**(-3)*z2
     &     + 704.D0/9*upx**(-3)*y
     &     + 20*upx**(-3)*y*z2
     &     + 352.D0/9*upx**(-3)*y**2
     &     + 10*upx**(-3)*y**2*z2
     &     - 194.D0/9*upx**(-2)
     &     - 3*upx**(-2)*z2
     &     - 188.D0/3*upx**(-2)*y
     &     - 6*upx**(-2)*y*z2
     &     - 194.D0/9*upx**(-2)*y**2
     &     + upx**(-2)*y**2*z2
     &     + 2*upx**(-1)
     &     - 1.D0/4*upx**(-1)*z2
     &     + 212.D0/9*upx**(-1)*y
     &     + 3.D0/2*upx**(-1)*y*z2
     &     + 2*upx**(-1)*y**2
     &     - 1.D0/4*upx**(-1)*y**2*z2
     &     + 28.D0/3*Hr1(-1)
     &     + 24*Hr1(-1)*umx**(-4)
     &     - 48*Hr1(-1)*umx**(-4)*y
     &     + 24*Hr1(-1)*umx**(-4)*y**2
     &     - 48*Hr1(-1)*umx**(-3)
     &     + 96*Hr1(-1)*umx**(-3)*y
     &     - 48*Hr1(-1)*umx**(-3)*y**2
     &     + 32*Hr1(-1)*umx**(-2)
     &     - 56*Hr1(-1)*umx**(-2)*y
     &     + 16*Hr1(-1)*umx**(-2)*y**2
     &     - 8*Hr1(-1)*umx**(-1)
     &     + 8*Hr1(-1)*umx**(-1)*y
     &     + 8*Hr1(-1)*umx**(-1)*y**2
     &     + 56.D0/3*Hr1(-1)*upx**(-4)
     &     + 112.D0/3*Hr1(-1)*upx**(-4)*y
     &     + 56.D0/3*Hr1(-1)*upx**(-4)*y**2
     &     - 112.D0/3*Hr1(-1)*upx**(-3)
     &     - 224.D0/3*Hr1(-1)*upx**(-3)*y
     &     - 112.D0/3*Hr1(-1)*upx**(-3)*y**2
     &     + 32.D0/3*Hr1(-1)*upx**(-2)
     &     + 40*Hr1(-1)*upx**(-2)*y
     &     + 32.D0/3*Hr1(-1)*upx**(-2)*y**2
     &     + 8*Hr1(-1)*upx**(-1)
     &     - 8.D0/3*Hr1(-1)*upx**(-1)*y
     &     + 8*Hr1(-1)*upx**(-1)*y**2
     &     + 4*Hr1(-1)*HAr1(-1)
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-4)
     &     + 32*Hr1(-1)*HAr1(-1)*upx**(-4)*y
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-4)*y**2
     &     - 32*Hr1(-1)*HAr1(-1)*upx**(-3)
     &     - 64*Hr1(-1)*HAr1(-1)*upx**(-3)*y
     &     - 32*Hr1(-1)*HAr1(-1)*upx**(-3)*y**2
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-2)
     &     + 40*Hr1(-1)*HAr1(-1)*upx**(-2)*y
     &     + 16*Hr1(-1)*HAr1(-1)*upx**(-2)*y**2
     &     - 8*Hr1(-1)*HAr1(-1)*upx**(-1)*y
     &     + 12*Hr1(-1)*HBr1(-1)
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-4)
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-4)*y
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-4)*y**2
     &     - 32*Hr1(-1)*HBr1(-1)*upx**(-3)
     &     - 64*Hr1(-1)*HBr1(-1)*upx**(-3)*y
     &     - 32*Hr1(-1)*HBr1(-1)*upx**(-3)*y**2
     &     + 32*Hr1(-1)*HBr1(-1)*upx**(-2)
     &     + 56*Hr1(-1)*HBr1(-1)*upx**(-2)*y
     &     + 16*Hr1(-1)*HBr1(-1)*upx**(-2)*y**2
     &     - 16*Hr1(-1)*HBr1(-1)*upx**(-1)
     &     - 24*Hr1(-1)*HBr1(-1)*upx**(-1)*y
     &     + 2*HAr1(-1)
     &     - 2*HAr1(-1)*upx**(-2)*y**(-1)
     &     - 4*HAr1(-1)*upx**(-2)
     &     - 2*HAr1(-1)*upx**(-2)*y
     &     + 2*HAr1(-1)*upx**(-1)*y**(-1)
     &     + 4*HAr1(-1)*upx**(-1)
     &     + 2*HAr1(-1)*upx**(-1)*y
     &     + 2*HAr1(-1)*y**(-1)
     &     - 2*HAr1(-1)*Hr1(0)
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-4)
     &     - 16*HAr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16*HAr1(-1)*Hr1(0)*upx**(-3)
     &     + 32*HAr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16*HAr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-2)
     &     - 20*HAr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8*HAr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 4*HAr1(-1)*Hr1(0)*upx**(-1)*y
     &     + 4*HBr1(-1)
     &     - 2*HBr1(-1)*upx**(-2)*z**(-1)
     &     + 2*HBr1(-1)*upx**(-2)*y
     &     + 2*HBr1(-1)*upx**(-1)*z**(-1)
     &     - 2*HBr1(-1)*upx**(-1)*y
     &     + 2*HBr1(-1)*z**(-1)
     &     - 6*HBr1(-1)*Hr1(0)
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-4)
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-4)*y
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-4)*y**2
     &     + 16*HBr1(-1)*Hr1(0)*upx**(-3)
     &     + 32*HBr1(-1)*Hr1(0)*upx**(-3)*y
     &     + 16*HBr1(-1)*Hr1(0)*upx**(-3)*y**2
     &     - 16*HBr1(-1)*Hr1(0)*upx**(-2)
     &     - 28*HBr1(-1)*Hr1(0)*upx**(-2)*y
     &     - 8*HBr1(-1)*Hr1(0)*upx**(-2)*y**2
     &     + 8*HBr1(-1)*Hr1(0)*upx**(-1)
     &     + 12*HBr1(-1)*Hr1(0)*upx**(-1)*y
     &     - 16*Hr2(-1,-1)
     &     - 32*Hr2(-1,-1)*upx**(-4)
     &     - 64*Hr2(-1,-1)*upx**(-4)*y
     &     - 32*Hr2(-1,-1)*upx**(-4)*y**2
     &     + 64*Hr2(-1,-1)*upx**(-3)
     &     + 128*Hr2(-1,-1)*upx**(-3)*y
     &     + 64*Hr2(-1,-1)*upx**(-3)*y**2
     &     - 48*Hr2(-1,-1)*upx**(-2)
     &     - 96*Hr2(-1,-1)*upx**(-2)*y
     &     - 32*Hr2(-1,-1)*upx**(-2)*y**2
     &     + 16*Hr2(-1,-1)*upx**(-1)
     &     + 32*Hr2(-1,-1)*upx**(-1)*y
     &     + 4*HAr2(-1,-1)
     &     + 8*HAr2(-1,-1)*upx**(-2)*y
     &     - 8*HAr2(-1,-1)*upx**(-1)*y
     &     - 4*HBr2(-1,-1)
     &     - 16*HBr2(-1,-1)*upx**(-2)
     &     - 8*HBr2(-1,-1)*upx**(-2)*y
     &     + 16*HBr2(-1,-1)*upx**(-1)
     &     + 8*HBr2(-1,-1)*upx**(-1)*y
     &     + 8*Hr2(-1,0)
     &     + 16*Hr2(-1,0)*upx**(-4)
     &     + 32*Hr2(-1,0)*upx**(-4)*y
     &     + 16*Hr2(-1,0)*upx**(-4)*y**2
     &     - 32*Hr2(-1,0)*upx**(-3)
     &     - 64*Hr2(-1,0)*upx**(-3)*y
     &     - 32*Hr2(-1,0)*upx**(-3)*y**2
     &     + 24*Hr2(-1,0)*upx**(-2)
     &     + 48*Hr2(-1,0)*upx**(-2)*y
     &     + 16*Hr2(-1,0)*upx**(-2)*y**2
     &     - 8*Hr2(-1,0)*upx**(-1)
     &     - 16*Hr2(-1,0)*upx**(-1)*y
     &     - 14.D0/3*Hr1(0)
     &     - 12*Hr1(0)*umx**(-4)
     &     + 24*Hr1(0)*umx**(-4)*y
     &     - 12*Hr1(0)*umx**(-4)*y**2
     &     + 24*Hr1(0)*umx**(-3)
     &     - 48*Hr1(0)*umx**(-3)*y
     &     + 24*Hr1(0)*umx**(-3)*y**2
     &     - 16*Hr1(0)*umx**(-2)
     &     + 28*Hr1(0)*umx**(-2)*y
     &     - 8*Hr1(0)*umx**(-2)*y**2
     &     + 4*Hr1(0)*umx**(-1)
     &     - 4*Hr1(0)*umx**(-1)*y
     &     - 4*Hr1(0)*umx**(-1)*y**2
     &     - 28.D0/3*Hr1(0)*upx**(-4)
     &     - 56.D0/3*Hr1(0)*upx**(-4)*y
     &     - 28.D0/3*Hr1(0)*upx**(-4)*y**2
     &     + 56.D0/3*Hr1(0)*upx**(-3)
     &     + 112.D0/3*Hr1(0)*upx**(-3)*y
     &     + 56.D0/3*Hr1(0)*upx**(-3)*y**2
     &     - 16.D0/3*Hr1(0)*upx**(-2)
     &     - 20*Hr1(0)*upx**(-2)*y
     &     - 16.D0/3*Hr1(0)*upx**(-2)*y**2
     &     - 4*Hr1(0)*upx**(-1)
     &     + 4.D0/3*Hr1(0)*upx**(-1)*y
     &     - 4*Hr1(0)*upx**(-1)*y**2
     &     + 24*Hr2(0,-1)*umx**(-5)
     &     - 48*Hr2(0,-1)*umx**(-5)*y
     &     + 24*Hr2(0,-1)*umx**(-5)*y**2
     &     - 60*Hr2(0,-1)*umx**(-4)
     &     + 120*Hr2(0,-1)*umx**(-4)*y
     &     - 60*Hr2(0,-1)*umx**(-4)*y**2
     &     + 54*Hr2(0,-1)*umx**(-3)
     &     - 100*Hr2(0,-1)*umx**(-3)*y
     &     + 38*Hr2(0,-1)*umx**(-3)*y**2
     &     - 21*Hr2(0,-1)*umx**(-2)
     &     + 30*Hr2(0,-1)*umx**(-2)*y
     &     + 3*Hr2(0,-1)*umx**(-2)*y**2
     &     + 33.D0/2*Hr2(0,-1)*umx**(-1)
     &     - 3*Hr2(0,-1)*umx**(-1)*y
     &     + 1.D0/2*Hr2(0,-1)*umx**(-1)*y**2
     &     + 16*Hr2(0,-1)*upx**(-4)
     &     + 32*Hr2(0,-1)*upx**(-4)*y
     &     + 16*Hr2(0,-1)*upx**(-4)*y**2
     &     - 20*Hr2(0,-1)*upx**(-3)
     &     - 40*Hr2(0,-1)*upx**(-3)*y
     &     - 20*Hr2(0,-1)*upx**(-3)*y**2
     &     + 6*Hr2(0,-1)*upx**(-2)
     &     + 12*Hr2(0,-1)*upx**(-2)*y
     &     - 2*Hr2(0,-1)*upx**(-2)*y**2
     &     + 1.D0/2*Hr2(0,-1)*upx**(-1)
     &     - 3*Hr2(0,-1)*upx**(-1)*y
     &     + 1.D0/2*Hr2(0,-1)*upx**(-1)*y**2
     &     - 2*HAr2(0,-1)
     &     - 4*HAr2(0,-1)*upx**(-2)*y
     &     + 4*HAr2(0,-1)*upx**(-1)*y
     &     + 2*HBr2(0,-1)
     &     + 8*HBr2(0,-1)*upx**(-2)
     &     + 4*HBr2(0,-1)*upx**(-2)*y
     &     - 8*HBr2(0,-1)*upx**(-1)
     &     - 4*HBr2(0,-1)*upx**(-1)*y
     &     - 12*Hr2(0,0)*umx**(-5)
     &     + 24*Hr2(0,0)*umx**(-5)*y
     &     - 12*Hr2(0,0)*umx**(-5)*y**2
     &     + 30*Hr2(0,0)*umx**(-4)
     &     - 60*Hr2(0,0)*umx**(-4)*y
     &     + 30*Hr2(0,0)*umx**(-4)*y**2
     &     - 27*Hr2(0,0)*umx**(-3)
     &     + 50*Hr2(0,0)*umx**(-3)*y
     &     - 19*Hr2(0,0)*umx**(-3)*y**2
     &     + 21.D0/2*Hr2(0,0)*umx**(-2)
     &     - 15*Hr2(0,0)*umx**(-2)*y
     &     - 3.D0/2*Hr2(0,0)*umx**(-2)*y**2
     &     - 33.D0/4*Hr2(0,0)*umx**(-1)
     &     + 3.D0/2*Hr2(0,0)*umx**(-1)*y
     &     - 1.D0/4*Hr2(0,0)*umx**(-1)*y**2
     &     - 8*Hr2(0,0)*upx**(-4)
     &     - 16*Hr2(0,0)*upx**(-4)*y
     &     - 8*Hr2(0,0)*upx**(-4)*y**2
     &     + 10*Hr2(0,0)*upx**(-3)
     &     + 20*Hr2(0,0)*upx**(-3)*y
     &     + 10*Hr2(0,0)*upx**(-3)*y**2
     &     - 3*Hr2(0,0)*upx**(-2)
     &     - 6*Hr2(0,0)*upx**(-2)*y
     &     + Hr2(0,0)*upx**(-2)*y**2
     &     - 1.D0/4*Hr2(0,0)*upx**(-1)
     &     + 3.D0/2*Hr2(0,0)*upx**(-1)*y
     &     - 1.D0/4*Hr2(0,0)*upx**(-1)*y**2
     &     - 44.D0/3*DLog(tm**(-2)*w**2)
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)
     &     - 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-4)*y**2
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)
     &     + 352.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y
     &     + 176.D0/3*DLog(tm**(-2)*w**2)*upx**(-3)*y**2
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)
     &     - 88*DLog(tm**(-2)*w**2)*upx**(-2)*y
     &     - 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-2)*y**2
     &     + 88.D0/3*DLog(tm**(-2)*w**2)*upx**(-1)*y
     &     )
      end

*
**************************** 72 characters *****************************
*
*
      double precision function PoleEPm2T13(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm2T13 =  + Nc * (
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm2T13 = PoleEPm2T13 + Nc**3 * (
     &     + 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm2T13 = PoleEPm2T13 + Nc**5 * (
     &     - 16.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     )
      return
      end
*
*
*
      double precision function PoleEPm1T13(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm1T13 =  + Nc**(-1) * (
     &     + 8.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 10.D0*upy**(-1)
     &     + 8.D0*upz**(-2)
     &     - 10.D0*upz**(-1)
     &     - 2.D0*z*upy**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nl * (
     &     + 16.D0/3.D0*upy**(-2)
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nc * (
     &     - 8.D0
     &     - 160.D0/3.D0*upy**(-2)
     &     - 272.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 176.D0/3.D0*upy**(-1)
     &     - 136.D0/3.D0*upz**(-2)
     &     + 146.D0/3.D0*upz**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     - 8.D0*ypzp2**(-1)
     &     + 16.D0/3.D0*z*upy**(-1)
     &     - 8.D0*z*ypzp2**(-2)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     + 10.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nc**2*Nl * (
     &     - 32.D0/3.D0*upy**(-2)
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 40.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*ypzp2**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nc**3 * (
     &     + 8.D0
     &     + 248.D0/3.D0*upy**(-2)
     &     + 224.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 262.D0/3.D0*upy**(-1)
     &     + 112.D0/3.D0*upz**(-2)
     &     - 116.D0/3.D0*upz**(-1)
     &     + 44.D0/3.D0*ypzp2**(-2)
     &     + 136.D0/3.D0*ypzp2**(-1)
     &     - 14.D0/3.D0*z*upy**(-1)
     &     + 88.D0/3.D0*z*ypzp2**(-2)
     &     + 44.D0/3.D0*z**2*ypzp2**(-2)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nc**4*Nl * (
     &     + 16.D0/3.D0*upy**(-2)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*ypzp2**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Nc**5 * (
     &     - 112.D0/3.D0*upy**(-2)
     &     + 116.D0/3.D0*upy**(-1)
     &     - 32.D0/3.D0*ypzp2**(-2)
     &     - 112.D0/3.D0*ypzp2**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     - 64.D0/3.D0*z*ypzp2**(-2)
     &     - 32.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(-1)*Nc * (
     &     + 8.D0
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(-1)*Nc**5 * (
     &     + 16.D0*upy**(-2)
     &     - 20.D0*upy**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc**(-1) * (
     &     - 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 34.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 34.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc * (
     &     - 4.D0
     &     + 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 34.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)
     &     - 52.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     + 30.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc * (
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc**3 * (
     &     + 4.D0
     &     + 8.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 10.D0*upy**(-1)
     &     - 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)
     &     + 18.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*ypzp2**(-2)
     &     + 4.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc**3 * (
     &     + 8.D0*z*ypzp2**(-2)
     &     + 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzp2**(-2)
     &     + 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + Hr1(0)*Nc**5 * (
     &     - 8.D0*upy**(-2)
     &     + 10.D0*upy**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     - 8.D0*ypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     - 8.D0*z*ypzp2**(-2)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + HAr1(-1)*Nc * (
     &     - 16.D0
     &     + 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + HAr1(-1)*Nc**3 * (
     &     + 16.D0
     &     - 32.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + HAr1(-1)*Nc**5 * (
     &     + 16.D0*upy**(-2)
     &     - 20.D0*upy**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + HBr1(-1)*Nc * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T13 = PoleEPm1T13 + HBr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 16.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      return
      end
*
*
*
      double precision function PoleEP0T13(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,zp2,zm3,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	zp2 = 2.d0+z
	zm3 = z-3.d0
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEP0T13 =  + Nc**(-1) * (
     &     + 4.D0
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + z**(-1)
     &     + 12.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 21.D0*upy**(-1)
     &     + 12.D0*upz**(-2)
     &     - 21.D0*upz**(-1)
     &     + 4.D0*z2*upy**(-3)
     &     + 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upy**(-2)
     &     + 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**(-1) * (
     &     - 93.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z2*upy**(-1)
     &     + 4.D0*z2*upz**(-3)
     &     + 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upz**(-2)
     &     - 93.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z2*upz**(-1)
     &     + 62.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z2
     &     - 4.D0*z*upy**(-2)
     &     + z*upy**(-1)
     &     + 4.D0*z*z2*upy**(-3)
     &     - 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**(-1) * (
     &     - 3.D0*z*z2*upy**(-1)
     &     + z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-2)
     &     + y*upz**(-1)
     &     + 4.D0*y*z2*upz**(-3)
     &     - 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*y*z2*upz**(-1)
     &     + y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nl * (
     &     + 8.D0/3.D0
     &     + 8.D0/3.D0*upy**(-1)
     &     + 8.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc * (
     &     - 44.D0/3.D0
     &     - 4.D0*y**(-1)*upz**(-1)
     &     - y**(-1)*zp2**(-1)
     &     - 3.D0*y**(-1)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)
     &     - 40.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 133.D0/3.D0*upy**(-1)
     &     - 28.D0*upz**(-2)
     &     + 82.D0/3.D0*upz**(-1)
     &     + zp2**(-1)*ypzp2**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc * (
     &     - 8.D0*z2*upy**(-3)
     &     - 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 52.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z2*upy**(-1)*upz**(-1)
     &     + 107.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 37.D0*z2*upy**(-1)
     &     - 4.D0*z2*upz**(-3)
     &     - 48.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upz**(-2)
     &     - 64.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 140.D0*z2*upz**(-1)*ypzp2**(-1)
     &     + 132.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc * (
     &     - 44.D0*z2*upz**(-1)
     &     + 48.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 31.D0*z2*ypzp2**(-1)
     &     - 76.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z2
     &     + 8.D0*z*upy**(-2)
     &     - 41.D0/3.D0*z*upy**(-1)
     &     - 8.D0*z*ypzp2**(-2)
     &     + z*ypzp2**(-1)
     &     - 8.D0*z*z2*upy**(-3)
     &     + 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc * (
     &     - 11.D0*z*z2*upy**(-1)
     &     - 16.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 27.D0*z*z2*ypzp2**(-1)
     &     - 14.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     - 7.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-2)
     &     - 38.D0/3.D0*y*upz**(-1)
     &     - 4.D0*y*z2*upz**(-3)
     &     + 48.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc * (
     &     - 52.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*z2*upz**(-1)
     &     + 6.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**2*Nl * (
     &     - 8.D0/3.D0
     &     - 14.D0/3.D0*upy**(-1)
     &     - 8.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     - 4.D0/3.D0*ypzp2**(-1)
     &     - 16.D0/3.D0*z*upy**(-1)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**2*Nh * (
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     + 2.D0/3.D0*upy**(-1)
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 4.D0/3.D0*ypzp2**(-1)
     &     - 12.D0*z2*upy**(-1)*upz**(-1)
     &     + 12.D0*z2*upz**(-1)*ypzp2**(-1)
     &     + 24.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     + 32.D0/3.D0
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)*zp2**(-1)
     &     + 3.D0*y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + z**(-1)
     &     + 44.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)*zm3**(-1)
     &     - 85.D0/3.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 19.D0/3.D0*upz**(-1)
     &     - 2.D0*zp2**(-1)*ypzp2**(-1)
     &     - 8.D0*zm3**(-1)*ypzm2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     - 3.D0*ypzm2**(-1)
     &     - 20.D0/3.D0*ypzp2**(-2)
     &     + 103.D0/3.D0*ypzp2**(-1)
     &     + 4.D0*z2*upy**(-3)
     &     - 100.D0*z2*upy**(-2)
     &     - 12.D0*z2*upy**(-1)*upz**(-1)
     &     + 96.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 85.D0*z2*upy**(-1)
     &     + 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*upz**(-2)
     &     + 64.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     - 140.D0*z2*upz**(-1)*ypzp2**(-1)
     &     - 39.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 47.D0*z2*upz**(-1)
     &     - 96.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 96.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 48.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 36.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 52.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z2*ypzp2**(-2)
     &     - 2.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     - 63.D0*z2*ypzp2**(-1)
     &     + 17.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z2
     &     - 4.D0*z*upy**(-2)
     &     + 73.D0/3.D0*z*upy**(-1)
     &     - 40.D0/3.D0*z*ypzp2**(-2)
     &     + 4.D0*z*z2*upy**(-3)
     &     + 8.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 23.D0*z*z2*upy**(-1)
     &     + 18.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z*z2*ypzp2**(-2)
     &     - 96.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 23.D0*z*z2*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     + 16.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*ypzm2**(-1)
     &     - 20.D0/3.D0*z**2*ypzp2**(-2)
     &     + z**2*ypzp2**(-1)
     &     + 3.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z**2*z2*ypzp2**(-2)
     &     - 22.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 35.D0/3.D0*y*upz**(-1)
     &     - 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y*z2*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**3 * (
     &     - 9.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**4*Nl * (
     &     + 2.D0*upy**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     + 4.D0/3.D0*ypzp2**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**4*Nh * (
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 2.D0/3.D0*upy**(-1)
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 4.D0/3.D0*ypzp2**(-1)
     &     + 12.D0*z2*upy**(-1)*upz**(-1)
     &     - 12.D0*z2*upz**(-1)*ypzp2**(-1)
     &     - 24.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**5 * (
     &     - y**(-1)*zp2**(-1)
     &     - y**(-1)
     &     - 16.D0*upy**(-2)
     &     - 8.D0*upy**(-1)*zm3**(-1)
     &     + 5.D0*upy**(-1)
     &     + zp2**(-1)*ypzp2**(-1)
     &     + 8.D0*zm3**(-1)*ypzm2**(-1)
     &     + 3.D0*ypzm2**(-1)
     &     + 32.D0/3.D0*ypzp2**(-2)
     &     - 55.D0/3.D0*ypzp2**(-1)
     &     + 40.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 5.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**5 * (
     &     - 45.D0*z2*upy**(-1)
     &     + 96.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 96.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 48.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 36.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z2*ypzp2**(-2)
     &     + 10.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z2*ypzp2**(-1)
     &     - 3.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z2
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**5 * (
     &     - 35.D0/3.D0*z*upy**(-1)
     &     + 64.D0/3.D0*z*ypzp2**(-2)
     &     - z*ypzp2**(-1)
     &     - 4.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*z*z2*upy**(-1)
     &     - 2.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z*z2*ypzp2**(-2)
     &     + 16.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*z2*ypzp2**(-1)
     &     - 3.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z**2*ypzm2**(-1)
     &     + 32.D0/3.D0*z**2*ypzp2**(-2)
     &     - z**2*ypzp2**(-1)
     &     - z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Nc**5 * (
     &     + 12.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z**2*z2*ypzp2**(-2)
     &     + 6.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nl * (
     &     - 16.D0/3.D0*upy**(-2)
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nc * (
     &     + 88.D0/3.D0*upy**(-2)
     &     + 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 110.D0/3.D0*upy**(-1)
     &     + 88.D0/3.D0*upz**(-2)
     &     - 110.D0/3.D0*upz**(-1)
     &     - 22.D0/3.D0*z*upy**(-1)
     &     - 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nc**2*Nl * (
     &     + 32.D0/3.D0*upy**(-2)
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 40.D0/3.D0*upy**(-1)
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*ypzp2**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nc**3 * (
     &     - 176.D0/3.D0*upy**(-2)
     &     - 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 220.D0/3.D0*upy**(-1)
     &     - 88.D0/3.D0*upz**(-2)
     &     + 110.D0/3.D0*upz**(-1)
     &     - 44.D0/3.D0*ypzp2**(-2)
     &     - 88.D0/3.D0*ypzp2**(-1)
     &     + 44.D0/3.D0*z*upy**(-1)
     &     - 88.D0/3.D0*z*ypzp2**(-2)
     &     - 44.D0/3.D0*z**2*ypzp2**(-2)
     &     + 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nc**4*Nl * (
     &     - 16.D0/3.D0*upy**(-2)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*ypzp2**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Dlog(tm**(-2)*w**2)*Nc**5 * (
     &     + 88.D0/3.D0*upy**(-2)
     &     - 110.D0/3.D0*upy**(-1)
     &     + 44.D0/3.D0*ypzp2**(-2)
     &     + 88.D0/3.D0*ypzp2**(-1)
     &     - 22.D0/3.D0*z*upy**(-1)
     &     + 88.D0/3.D0*z*ypzp2**(-2)
     &     + 44.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*Nc * (
     &     - 16.D0
     &     + 64.D0*upz**(-1)*ypzm2**(-1)
     &     + 32.D0*upz**(-1)
     &     - 48.D0*ypzm2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 16.D0*z*ypzm2**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*Nc**3 * (
     &     + 16.D0
     &     - 96.D0*upy**(-1)*zm3**(-2)
     &     - 32.D0*upy**(-1)*zm3**(-1)
     &     - 6.D0*upy**(-1)
     &     - 64.D0*upz**(-1)*ypzm2**(-1)
     &     - 32.D0*upz**(-1)
     &     + 96.D0*zm3**(-2)*ypzm2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-2)
     &     + 32.D0*zm3**(-1)*ypzm2**(-1)
     &     + 36.D0*ypzm2**(-2)
     &     + 46.D0*ypzm2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 18.D0*z*ypzm2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*Nc**3 * (
     &     - 4.D0*z*ypzp2**(-1)
     &     + 12.D0*z**2*ypzm2**(-2)
     &     - 4.D0*z**2*ypzm2**(-1)
     &     + 4.D0*z**2*ypzp2**(-1)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*Nc**5 * (
     &     + 96.D0*upy**(-1)*zm3**(-2)
     &     + 32.D0*upy**(-1)*zm3**(-1)
     &     + 6.D0*upy**(-1)
     &     - 96.D0*zm3**(-2)*ypzm2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-2)
     &     - 32.D0*zm3**(-1)*ypzm2**(-1)
     &     - 36.D0*ypzm2**(-2)
     &     + 2.D0*ypzm2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 2.D0*z*ypzm2**(-1)
     &     - 12.D0*z*ypzp2**(-1)
     &     - 12.D0*z**2*ypzm2**(-2)
     &     + 4.D0*z**2*ypzm2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*HAr1(-1)*Nc * (
     &     + 40.D0
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 64.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*HAr1(-1)*Nc**3 * (
     &     - 32.D0
     &     + 32.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 48.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*HAr1(-1)*Nc**5 * (
     &     - 8.D0
     &     - 32.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*HBr1(-1)*Nc * (
     &     - 16.D0
     &     - 64.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(-1)*HBr1(-1)*Nc**3 * (
     &     + 16.D0
     &     + 64.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**(-1) * (
     &     - 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc * (
     &     + 8.D0
     &     + upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-1)*ypzm2**(-1)
     &     + 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-1)
     &     + 24.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 22.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)
     &     - 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzm2**(-1)
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc * (
     &     - 8.D0*z*ypzp2**(-1)
     &     - 11.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     - 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**2*Nh * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-1)*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**3 * (
     &     - 8.D0
     &     + 48.D0*upy**(-1)*zm3**(-2)
     &     + 16.D0*upy**(-1)*zm3**(-1)
     &     + 7.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upy**(-1)
     &     + 32.D0*upz**(-1)*ypzm2**(-1)
     &     - 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-1)
     &     - 48.D0*zm3**(-2)*ypzm2**(-1)
     &     - 48.D0*zm3**(-1)*ypzm2**(-2)
     &     - 16.D0*zm3**(-1)*ypzm2**(-1)
     &     - 18.D0*ypzm2**(-2)
     &     - 23.D0*ypzm2**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**3 * (
     &     - 22.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     + 12.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 9.D0*z*ypzm2**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*ypzp2**(-1)
     &     - z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z**2*ypzm2**(-2)
     &     + 2.D0*z**2*ypzm2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**3 * (
     &     - 2.D0*z**2*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**4*Nh * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-1)*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*Nc**5 * (
     &     - 48.D0*upy**(-1)*zm3**(-2)
     &     - 16.D0*upy**(-1)*zm3**(-1)
     &     - 3.D0*upy**(-1)
     &     + 48.D0*zm3**(-2)*ypzm2**(-1)
     &     + 48.D0*zm3**(-1)*ypzm2**(-2)
     &     + 16.D0*zm3**(-1)*ypzm2**(-1)
     &     + 18.D0*ypzm2**(-2)
     &     - ypzm2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - z*ypzm2**(-1)
     &     + 6.D0*z*ypzp2**(-1)
     &     + 6.D0*z**2*ypzm2**(-2)
     &     - 2.D0*z**2*ypzm2**(-1)
     &     + 2.D0*z**2*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc**(-1) * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 42.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc * (
     &     - 20.D0
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 58.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 88.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 28.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc * (
     &     + 32.D0*z*ypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     + 16.D0
     &     - 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)
     &     - 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     - 24.D0*z*ypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HAr1(-1)*Nc**5 * (
     &     + 4.D0
     &     + 16.D0*upy**(-2)
     &     - 10.D0*upy**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     - 2.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HBr1(-1)*Nc**(-1) * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 42.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HBr1(-1)*Nc * (
     &     + 8.D0
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)
     &     + 52.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HBr1(-1)*Nc * (
     &     - 4.D0*y*upz**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-2)
     &     - 10.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,-1)*Nc * (
     &     - 40.D0
     &     + 56.D0*upz**(-1)
     &     + 64.D0*ypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-1)
     &     + 24.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,-1)*Nc**3 * (
     &     + 32.D0
     &     - 20.D0*upy**(-1)
     &     - 56.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 48.D0*z*ypzp2**(-1)
     &     - 24.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,-1)*Nc**5 * (
     &     + 8.D0
     &     + 20.D0*upy**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 86.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 86.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc * (
     &     + 20.D0
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 58.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 140.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 44.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 16.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc * (
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*y*upz**(-1)
     &     + 20.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc**3 * (
     &     - 16.D0
     &     - 28.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 54.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 24.D0*z*ypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc**3 * (
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y*upz**(-1)
     &     - 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(-1,0)*Nc**5 * (
     &     - 4.D0
     &     - 10.D0*upy**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc * (
     &     + 20.D0
     &     + 128.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*upz**(-1)
     &     - 96.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     - 44.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc * (
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc**3 * (
     &     - 16.D0
     &     - 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)
     &     - 128.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*upz**(-1)
     &     + 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 72.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 104.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc**3 * (
     &     + 16.D0*ypzp2**(-1)
     &     + 38.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     - 36.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*ypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y*upz**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc**3 * (
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc**5 * (
     &     - 4.D0
     &     + 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)
     &     - 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 72.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 6.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,-1)*Nc**5 * (
     &     - 2.D0*z*upy**(-1)
     &     + 4.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**(-1) * (
     &     + 4.D0
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     - 9.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upy**(-1)
     &     - 9.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upz**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*z*upy**(-1)
     &     + 5.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y*upz**(-1)
     &     + 5.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc * (
     &     - 14.D0
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     - 5.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + upy**(-1)
     &     - 64.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*upz**(-1)*ypzp2**(-1)
     &     + 12.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)
     &     + 48.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 9.D0*ypzp2**(-1)
     &     + 30.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc * (
     &     - z*upy**(-1)
     &     - 16.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 13.D0*z*ypzp2**(-1)
     &     + z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**2*Nh * (
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upz**(-1)*ypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**3 * (
     &     + 8.D0
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     + 96.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 19.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*upy**(-1)
     &     + 64.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*upz**(-1)*ypzp2**(-1)
     &     - 3.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 11.D0*upz**(-1)
     &     - 96.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**3 * (
     &     - 52.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - ypzp2**(-1)
     &     - 25.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z*upy**(-1)
     &     + 18.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*z*ypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**3 * (
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*y*upz**(-1)
     &     + 7.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**4*Nh * (
     &     - 4.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upz**(-1)*ypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**5 * (
     &     + 2.D0
     &     - 96.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 5.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*upy**(-1)
     &     + 96.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)
     &     - 3.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(0,0)*Nc**5 * (
     &     + z*upy**(-1)
     &     - 2.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*ypzp2**(-1)
     &     - 3.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(1,0)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(1,0)*Nc * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 104.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 60.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(1,0)*Nc * (
     &     + 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + Hr2(1,0)*Nc**3 * (
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc**(-1) * (
     &     + 2.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     - y**(-2)
     &     - 2.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)
     &     - y**(-1)*z
     &     - 12.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 4.D0*upz**(-1)
     &     + 4.D0*z*upy**(-2)
     &     + 4.D0*z*upy**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc * (
     &     + 24.D0
     &     + 4.D0*y**(-2)*upz**(-1)
     &     + y**(-2)*zp2**(-1)
     &     + 3.D0*y**(-2)
     &     - 4.D0*y**(-1)*upz**(-1)
     &     - y**(-1)*zp2**(-2)
     &     + y**(-1)*zp2**(-1)
     &     + y**(-1)
     &     + 3.D0*y**(-1)*z
     &     + 40.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     - 20.D0*upz**(-1)
     &     + zp2**(-2)*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc * (
     &     - zp2**(-1)*ypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 23.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-2)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 31.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc**3 * (
     &     - 26.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     - 2.D0*y**(-2)*zp2**(-1)
     &     - 3.D0*y**(-2)
     &     + 6.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)*zp2**(-2)
     &     + 2.D0*y**(-1)*zp2**(-1)
     &     + 4.D0*y**(-1)
     &     - 3.D0*y**(-1)*z
     &     - 44.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 24.D0*upy**(-1)
     &     + 16.D0*upz**(-1)
     &     - 2.D0*zp2**(-2)*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc**3 * (
     &     - 2.D0*zp2**(-1)*ypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 14.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-2)
     &     - 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 26.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc**5 * (
     &     + y**(-2)*zp2**(-1)
     &     + y**(-2)
     &     - y**(-1)*zp2**(-2)
     &     - 3.D0*y**(-1)*zp2**(-1)
     &     - 3.D0*y**(-1)
     &     + y**(-1)*z
     &     + 16.D0*upy**(-2)
     &     - 4.D0*upy**(-1)
     &     + zp2**(-2)*ypzp2**(-1)
     &     + 3.D0*zp2**(-1)*ypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 9.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*Nc**5 * (
     &     + 5.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*HBr1(-1)*Nc * (
     &     + 8.D0
     &     - 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr1(-1)*HBr1(-1)*Nc**3 * (
     &     - 8.D0
     &     + 32.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(-1,-1)*Nc * (
     &     - 16.D0
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 4.D0*upz**(-1)
     &     + 48.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 12.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(-1,-1)*Nc**3 * (
     &     + 8.D0
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 4.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     - 12.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(-1,-1)*Nc**5 * (
     &     + 8.D0
     &     + 20.D0*upy**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(0,-1)*Nc**(-1) * (
     &     + 4.D0
     &     + 4.D0*upy**(-3)
     &     + 8.D0*upy**(-2)
     &     - 12.D0*upy**(-1)*upz**(-1)
     &     - 2.D0*upy**(-1)
     &     + 8.D0*upz**(-1)
     &     + 4.D0*z*upy**(-3)
     &     + 2.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(0,-1)*Nc * (
     &     - 8.D0*upy**(-3)
     &     - 12.D0*upy**(-2)
     &     + 28.D0*upy**(-1)*upz**(-1)
     &     - 6.D0*upy**(-1)
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     - 16.D0*upz**(-1)
     &     - 10.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-3)
     &     - 6.D0*z*upy**(-1)
     &     - 2.D0*z*ypzp2**(-1)
     &     - 12.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(0,-1)*Nc**3 * (
     &     + 4.D0*upy**(-3)
     &     + 4.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 18.D0*upy**(-1)
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 8.D0*upz**(-1)
     &     - 6.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-3)
     &     + 6.D0*z*upy**(-1)
     &     - 6.D0*z*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HAr2(0,-1)*Nc**5 * (
     &     - 4.D0
     &     - 10.D0*upy**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr1(-1)*Nc**(-1) * (
     &     + 2.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     - z**(-2)
     &     - 2.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upy**(-1)
     &     - 12.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - y*z**(-1)
     &     + 4.D0*y*upz**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr1(-1)*Nc * (
     &     - 8.D0
     &     + 4.D0*z**(-2)*upy**(-1)
     &     + 2.D0*z**(-2)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upy**(-1)
     &     + 28.D0*upz**(-2)
     &     - 24.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr1(-1)*Nc * (
     &     + 2.D0*y*z**(-1)
     &     - 4.D0*y*upz**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr1(-1)*Nc**3 * (
     &     + 6.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     - z**(-2)
     &     + 6.D0*z**(-1)*upy**(-1)
     &     + 2.D0*z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upz**(-2)
     &     + 4.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr1(-1)*Nc**3 * (
     &     - y*z**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr2(-1,-1)*Nc * (
     &     - 8.D0
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr2(-1,-1)*Nc**3 * (
     &     + 8.D0
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr2(0,-1)*Nc**(-1) * (
     &     + 4.D0
     &     - 12.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)
     &     + 4.D0*upz**(-3)
     &     + 8.D0*upz**(-2)
     &     - 2.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-3)
     &     + 2.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr2(0,-1)*Nc * (
     &     - 4.D0
     &     + 12.D0*upy**(-1)*upz**(-1)
     &     - 8.D0*upy**(-1)
     &     - 4.D0*upz**(-3)
     &     - 8.D0*upz**(-2)
     &     - 8.D0*upz**(-1)
     &     - 4.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-3)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T13 = PoleEP0T13 + HBr2(0,-1)*Nc**3 * (
     &     + 10.D0*upz**(-1)
     &     + 4.D0*ypzp2**(-1)
     &     + 4.D0*z*ypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     )
      return
      end

*
**************************** 72 characters *****************************
*
*
      double precision function PoleEPm2T34(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm2T34 =  + Nc**(-1) * (
     &     + 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm2T34 = PoleEPm2T34 + Nc * (
     &     - 8.D0
     &     - 16.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm2T34 = PoleEPm2T34 + Nc**3 * (
     &     + 8.D0
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     )
      return
      end
*
*
*
      double precision function PoleEPm1T34(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm1T34 =  + Nc**(-3) * (
     &     - 8.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 10.D0*upy**(-1)
     &     - 8.D0*upz**(-2)
     &     + 10.D0*upz**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 2.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nc**(-2)*Nl * (
     &     - 16.D0/3.D0*upy**(-2)
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nc**(-1) * (
     &     + 12.D0
     &     + 136.D0/3.D0*upy**(-2)
     &     + 224.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 146.D0/3.D0*upy**(-1)
     &     + 136.D0/3.D0*upz**(-2)
     &     - 146.D0/3.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 8.D0*ypzp2**(-1)
     &     - 10.D0/3.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 10.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nl * (
     &     + 8.D0/3.D0
     &     + 16.D0/3.D0*upy**(-2)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     + 16.D0/3.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*ypzp2**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     + 32.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 16.D0/3.D0*z**2*ypzp2**(-2)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nc * (
     &     - 44.D0/3.D0
     &     - 112.D0/3.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     + 116.D0/3.D0*upy**(-1)
     &     - 112.D0/3.D0*upz**(-2)
     &     + 116.D0/3.D0*upz**(-1)
     &     - 88.D0/3.D0*ypzp2**(-2)
     &     - 184.D0/3.D0*ypzp2**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     - 176.D0/3.D0*z*ypzp2**(-2)
     &     + 88.D0/3.D0*z*ypzp2**(-1)
     &     - 88.D0/3.D0*z**2*ypzp2**(-2)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nc**2*Nl * (
     &     - 8.D0/3.D0
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 16.D0/3.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*ypzp2**(-1)
     &     - 32.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 16.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Nc**3 * (
     &     + 8.D0/3.D0
     &     - 224.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 64.D0/3.D0*ypzp2**(-2)
     &     + 160.D0/3.D0*ypzp2**(-1)
     &     + 128.D0/3.D0*z*ypzp2**(-2)
     &     - 64.D0/3.D0*z*ypzp2**(-1)
     &     + 64.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(-1)*Nc * (
     &     + 8.D0
     &     - 16.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(-1)*Nc**3 * (
     &     - 8.D0
     &     + 16.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc**(-3) * (
     &     + 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 34.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 34.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc**(-1) * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc * (
     &     - 4.D0
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-2)
     &     - 96.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 84.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-2)
     &     + 84.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 8.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc * (
     &     - 8.D0*ypzp2**(-1)
     &     - 48.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc * (
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc**3 * (
     &     + 4.D0
     &     + 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 34.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)
     &     + 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-2)
     &     - 34.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 8.D0*ypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + Hr1(0)*Nc**3 * (
     &     - 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HAr1(-1)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HAr1(-1)*Nc * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HAr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HBr1(-1)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HBr1(-1)*Nc * (
     &     + 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T34 = PoleEPm1T34 + HBr1(-1)*Nc**3 * (
     &     - 16.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     )
      return
      end
*
*
*
      double precision function PoleEP0T34(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,zp2,zm3,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	zp2 = 2.d0+z
	zm3 = z-3.d0
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEP0T34 =  + Nc**(-3) * (
     &     - 4.D0
     &     - 2.D0*y**(-1)*upz**(-1)
     &     - y**(-1)
     &     - 2.D0*z**(-1)*upy**(-1)
     &     - z**(-1)
     &     - 12.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 21.D0*upy**(-1)
     &     - 12.D0*upz**(-2)
     &     + 21.D0*upz**(-1)
     &     - 4.D0*z2*upy**(-3)
     &     - 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-3) * (
     &     + 93.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*z2*upy**(-1)
     &     - 4.D0*z2*upz**(-3)
     &     - 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z2*upz**(-2)
     &     + 93.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*z2*upz**(-1)
     &     - 62.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z2
     &     + 4.D0*z*upy**(-2)
     &     - z*upy**(-1)
     &     - 4.D0*z*z2*upy**(-3)
     &     + 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-3) * (
     &     + 3.D0*z*z2*upy**(-1)
     &     - z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 5.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-2)
     &     - y*upz**(-1)
     &     - 4.D0*y*z2*upz**(-3)
     &     + 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y*z2*upz**(-1)
     &     - y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 5.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-2)*Nl * (
     &     - 8.D0/3.D0
     &     - 8.D0/3.D0*upy**(-1)
     &     - 8.D0/3.D0*upz**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-1) * (
     &     + 41.D0/3.D0
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + y**(-1)*zp2**(-1)
     &     + 2.D0*y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + z**(-1)*ypzp2**(-1)
     &     + 2.D0*z**(-1)
     &     + 28.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 70.D0/3.D0*upy**(-1)
     &     + 28.D0*upz**(-2)
     &     - 70.D0/3.D0*upz**(-1)
     &     - zp2**(-1)*ypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-1) * (
     &     + 26.D0*ypzp2**(-1)
     &     + 4.D0*z2*upy**(-3)
     &     - 60.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 60.D0*z2*upy**(-1)*upz**(-1)
     &     + 34.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z2*upy**(-1)
     &     + 4.D0*z2*upz**(-3)
     &     - 60.D0*z2*upz**(-2)
     &     + 34.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z2*upz**(-1)
     &     + 64.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z2*ypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-1) * (
     &     - 22.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z2
     &     - 4.D0*z*upy**(-2)
     &     + 38.D0/3.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 4.D0*z*z2*upy**(-3)
     &     + 4.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*z*z2*upy**(-1)
     &     + 128.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 19.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     + 2.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**(-1) * (
     &     + 64.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-2)
     &     + 38.D0/3.D0*y*upz**(-1)
     &     + 4.D0*y*z2*upz**(-3)
     &     + 4.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*y*z2*upz**(-1)
     &     - 13.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nl * (
     &     - 8.D0/3.D0
     &     + 2.D0*upy**(-1)
     &     + 2.D0*upz**(-1)
     &     - 16.D0/3.D0*ypzp2**(-2)
     &     + 8.D0*ypzp2**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     - 32.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 16.D0/3.D0*z**2*ypzp2**(-2)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nh * (
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 2.D0/3.D0*upy**(-1)
     &     - 2.D0/3.D0*upz**(-1)
     &     + 64.D0*ypzp2**(-2)
     &     + 8.D0/3.D0*ypzp2**(-1)
     &     + 12.D0*z2*upy**(-1)*upz**(-1)
     &     - 48.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     + 44.D0/3.D0
     &     + 2.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)*zp2**(-1)
     &     - y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)*ypzp2**(-1)
     &     - z**(-1)
     &     - 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 8.D0*upy**(-1)*zm3**(-1)
     &     + upy**(-1)
     &     - 16.D0*upz**(-2)
     &     - 8.D0*upz**(-1)*ypzm2**(-1)
     &     + upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     + 2.D0*zp2**(-1)*ypzp2**(-1)
     &     + 8.D0*zm3**(-1)*ypzm2**(-1)
     &     + 10.D0*ypzm2**(-1)
     &     + 40.D0/3.D0*ypzp2**(-2)
     &     - 86.D0*ypzp2**(-1)
     &     + 48.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 92.D0*z2*upy**(-2)
     &     + 288.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 44.D0*z2*upy**(-1)*upz**(-1)
     &     - 96.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 219.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 41.D0*z2*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     + 48.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 92.D0*z2*upz**(-2)
     &     - 96.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 16.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 219.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 41.D0*z2*upz**(-1)
     &     + 96.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 96.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 16.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     + 120.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z2*ypzp2**(-2)
     &     - 36.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z2*ypzp2**(-1)
     &     + 138.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 15.D0*z2
     &     - 35.D0/3.D0*z*upy**(-1)
     &     - 4.D0*z*ypzm2**(-1)
     &     + 80.D0/3.D0*z*ypzp2**(-2)
     &     - 52.D0/3.D0*z*ypzp2**(-1)
     &     - 48.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     - 27.D0*z*z2*upy**(-1)
     &     - 48.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 80.D0*z*z2*ypzp2**(-2)
     &     + 152.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z*z2*ypzp2**(-1)
     &     - 20.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*ypzm2**(-1)
     &     + 40.D0/3.D0*z**2*ypzp2**(-2)
     &     - 2.D0*z**2*ypzp2**(-1)
     &     + 5.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc * (
     &     - 64.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z**2*z2*ypzp2**(-2)
     &     + 44.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 35.D0/3.D0*y*upz**(-1)
     &     - 48.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 27.D0*y*z2*upz**(-1)
     &     + 20.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**2*Nl * (
     &     + 16.D0/3.D0
     &     + 2.D0/3.D0*upy**(-1)
     &     + 2.D0/3.D0*upz**(-1)
     &     + 16.D0/3.D0*ypzp2**(-2)
     &     - 8.D0*ypzp2**(-1)
     &     + 32.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 16.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**2*Nh * (
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     + 2.D0/3.D0*upy**(-1)
     &     + 2.D0/3.D0*upz**(-1)
     &     - 64.D0*ypzp2**(-2)
     &     - 8.D0/3.D0*ypzp2**(-1)
     &     - 12.D0*z2*upy**(-1)*upz**(-1)
     &     + 48.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**3 * (
     &     - 73.D0/3.D0
     &     - 2.D0*y**(-1)*upz**(-1)
     &     + y**(-1)*zp2**(-1)
     &     - 2.D0*z**(-1)*upy**(-1)
     &     + z**(-1)*ypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)*zm3**(-1)
     &     + 4.D0/3.D0*upy**(-1)
     &     + 8.D0*upz**(-1)*ypzm2**(-1)
     &     + 4.D0/3.D0*upz**(-1)
     &     - zp2**(-1)*ypzp2**(-1)
     &     - 8.D0*zm3**(-1)*ypzm2**(-1)
     &     - 10.D0*ypzm2**(-1)
     &     - 64.D0/3.D0*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**3 * (
     &     + 60.D0*ypzp2**(-1)
     &     - 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 104.D0*z2*upy**(-1)*upz**(-1)
     &     + 96.D0*z2*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 92.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z2*upy**(-1)
     &     - 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z2*upz**(-2)
     &     + 96.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 16.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**3 * (
     &     + 92.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z2*upz**(-1)
     &     - 96.D0*z2*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 96.D0*z2*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 16.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 120.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z2*ypzp2**(-2)
     &     - 12.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z2*ypzp2**(-1)
     &     - 54.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**3 * (
     &     - 20.D0*z2
     &     + 4.D0*z*ypzm2**(-1)
     &     - 128.D0/3.D0*z*ypzp2**(-2)
     &     + 76.D0/3.D0*z*ypzp2**(-1)
     &     + 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*z*z2*upy**(-1)
     &     + 48.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 80.D0*z*z2*ypzp2**(-2)
     &     - 24.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z*z2*ypzp2**(-1)
     &     + 2.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzm2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Nc**3 * (
     &     - 64.D0/3.D0*z**2*ypzp2**(-2)
     &     + 2.D0*z**2*ypzp2**(-1)
     &     - 2.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z**2*z2*ypzp2**(-2)
     &     - 12.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*y*z2*upz**(-1)
     &     - 6.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nc**(-2)*Nl * (
     &     + 16.D0/3.D0*upy**(-2)
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nc**(-1) * (
     &     - 88.D0/3.D0*upy**(-2)
     &     - 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 110.D0/3.D0*upy**(-1)
     &     - 88.D0/3.D0*upz**(-2)
     &     + 110.D0/3.D0*upz**(-1)
     &     + 22.D0/3.D0*z*upy**(-1)
     &     + 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nl * (
     &     - 8.D0/3.D0
     &     - 16.D0/3.D0*upy**(-2)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     - 16.D0/3.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*ypzp2**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     - 32.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 16.D0/3.D0*z**2*ypzp2**(-2)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nc * (
     &     + 44.D0/3.D0
     &     + 88.D0/3.D0*upy**(-2)
     &     - 110.D0/3.D0*upy**(-1)
     &     + 88.D0/3.D0*upz**(-2)
     &     - 110.D0/3.D0*upz**(-1)
     &     + 88.D0/3.D0*ypzp2**(-2)
     &     + 88.D0/3.D0*ypzp2**(-1)
     &     - 22.D0/3.D0*z*upy**(-1)
     &     + 176.D0/3.D0*z*ypzp2**(-2)
     &     - 88.D0/3.D0*z*ypzp2**(-1)
     &     + 88.D0/3.D0*z**2*ypzp2**(-2)
     &     - 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nc**2*Nl * (
     &     + 8.D0/3.D0
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 16.D0/3.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*ypzp2**(-1)
     &     + 32.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 16.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Dlog(tm**(-2)*w**2)*Nc**3 * (
     &     - 44.D0/3.D0
     &     + 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 88.D0/3.D0*ypzp2**(-2)
     &     - 88.D0/3.D0*ypzp2**(-1)
     &     - 176.D0/3.D0*z*ypzp2**(-2)
     &     + 88.D0/3.D0*z*ypzp2**(-1)
     &     - 88.D0/3.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*Nc * (
     &     - 6.D0
     &     + 96.D0*upy**(-1)*zm3**(-2)
     &     - 32.D0*upy**(-1)*zm3**(-1)
     &     - 26.D0*upy**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-2)
     &     - 32.D0*upz**(-1)*ypzm2**(-1)
     &     - 26.D0*upz**(-1)
     &     - 96.D0*zm3**(-2)*ypzm2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-2)
     &     + 32.D0*zm3**(-1)*ypzm2**(-1)
     &     - 120.D0*ypzm2**(-2)
     &     + 40.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*Nc * (
     &     + 48.D0*z*ypzm2**(-2)
     &     + 8.D0*z*ypzm2**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 24.D0*z**2*ypzm2**(-2)
     &     + 8.D0*z**2*ypzm2**(-1)
     &     - 8.D0*z**2*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*Nc**3 * (
     &     + 6.D0
     &     - 96.D0*upy**(-1)*zm3**(-2)
     &     + 32.D0*upy**(-1)*zm3**(-1)
     &     + 26.D0*upy**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-2)
     &     + 32.D0*upz**(-1)*ypzm2**(-1)
     &     + 26.D0*upz**(-1)
     &     + 96.D0*zm3**(-2)*ypzm2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-2)
     &     - 32.D0*zm3**(-1)*ypzm2**(-1)
     &     + 120.D0*ypzm2**(-2)
     &     - 40.D0*ypzm2**(-1)
     &     - 8.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*Nc**3 * (
     &     - 48.D0*z*ypzm2**(-2)
     &     - 8.D0*z*ypzm2**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 24.D0*z**2*ypzm2**(-2)
     &     - 8.D0*z**2*ypzm2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*HAr1(-1)*Nc * (
     &     - 32.D0
     &     + 32.D0*upy**(-2)
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 48.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*HAr1(-1)*Nc**3 * (
     &     + 32.D0
     &     - 32.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 48.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 16.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*HBr1(-1)*Nc * (
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(-1)*HBr1(-1)*Nc**3 * (
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**(-3) * (
     &     + 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**(-1) * (
     &     + 7.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 7.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 11.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 7.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nh * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 128.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc * (
     &     + 3.D0
     &     - 48.D0*upy**(-1)*zm3**(-2)
     &     + 16.D0*upy**(-1)*zm3**(-1)
     &     - 23.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 13.D0*upy**(-1)
     &     - 48.D0*upz**(-1)*ypzm2**(-2)
     &     + 16.D0*upz**(-1)*ypzm2**(-1)
     &     - 23.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 13.D0*upz**(-1)
     &     + 48.D0*zm3**(-2)*ypzm2**(-1)
     &     + 48.D0*zm3**(-1)*ypzm2**(-2)
     &     - 16.D0*zm3**(-1)*ypzm2**(-1)
     &     + 60.D0*ypzm2**(-2)
     &     - 20.D0*ypzm2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc * (
     &     + 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*ypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 7.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 24.D0*z*ypzm2**(-2)
     &     - 4.D0*z*ypzm2**(-1)
     &     + 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     + 35.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*ypzm2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc * (
     &     - 4.D0*z**2*ypzm2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzp2**(-1)
     &     - 7.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     + 31.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**2*Nh * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 128.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**3 * (
     &     - 3.D0
     &     + 48.D0*upy**(-1)*zm3**(-2)
     &     - 16.D0*upy**(-1)*zm3**(-1)
     &     + 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 13.D0*upy**(-1)
     &     + 48.D0*upz**(-1)*ypzm2**(-2)
     &     - 16.D0*upz**(-1)*ypzm2**(-1)
     &     + 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 13.D0*upz**(-1)
     &     - 48.D0*zm3**(-2)*ypzm2**(-1)
     &     - 48.D0*zm3**(-1)*ypzm2**(-2)
     &     + 16.D0*zm3**(-1)*ypzm2**(-1)
     &     - 60.D0*ypzm2**(-2)
     &     + 20.D0*ypzm2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**3 * (
     &     + 4.D0*ypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 24.D0*z*ypzm2**(-2)
     &     + 4.D0*z*ypzm2**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzm2**(-2)
     &     + 4.D0*z**2*ypzm2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*Nc**3 * (
     &     - 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc**(-3) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 42.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc**(-1) * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc * (
     &     + 16.D0
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-2)
     &     + 96.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 100.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upy**(-1)
     &     - 132.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 46.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc * (
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 24.D0*z*ypzp2**(-1)
     &     + 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     + 24.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     + 20.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     - 16.D0
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 42.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)
     &     + 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     - 2.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 24.D0*z*ypzp2**(-1)
     &     - 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc**(-3) * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 42.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc**(-1) * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc * (
     &     + 96.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 132.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upy**(-1)
     &     + 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)
     &     - 100.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 46.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc * (
     &     + 16.D0*z*ypzp2**(-2)
     &     + 8.D0*z*ypzp2**(-1)
     &     + 20.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     + 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)
     &     + 42.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     - 8.D0*z*ypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     - 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,-1)*Nc * (
     &     + 16.D0
     &     - 36.D0*upy**(-1)
     &     - 36.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     - 20.D0*z*upy**(-1)
     &     - 20.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,-1)*Nc**3 * (
     &     - 16.D0
     &     + 36.D0*upy**(-1)
     &     + 36.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-1)
     &     + 20.D0*z*upy**(-1)
     &     + 20.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc**(-3) * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 86.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 86.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc**(-1) * (
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 60.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 60.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc * (
     &     - 8.D0
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 232.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upy**(-1)
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 232.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 92.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*z*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc * (
     &     - 34.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*y*upz**(-1)
     &     - 34.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc**3 * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 86.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*upy**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 86.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z*upy**(-1)
     &     + 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(-1,0)*Nc**3 * (
     &     + 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*y*upz**(-1)
     &     + 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc * (
     &     - 8.D0
     &     + 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upy**(-1)
     &     + 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*upz**(-1)
     &     - 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 240.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc * (
     &     - 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*z*upy**(-1)
     &     + 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*y*upz**(-1)
     &     - 10.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc * (
     &     + 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc**3 * (
     &     + 8.D0
     &     - 192.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*upy**(-1)
     &     - 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 18.D0*upz**(-1)
     &     + 192.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 192.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 240.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc**3 * (
     &     + 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z*upy**(-1)
     &     - 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*y*upz**(-1)
     &     + 10.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,-1)*Nc**3 * (
     &     - 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**(-3) * (
     &     - 4.D0
     &     + 8.D0*upy**(-1)*upz**(-1)
     &     + 9.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*upy**(-1)
     &     + 9.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*upz**(-1)
     &     + 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z*upy**(-1)
     &     - 5.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*y*upz**(-1)
     &     - 5.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**(-1) * (
     &     - 1.D0
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     + 14.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*upy**(-1)
     &     + 14.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**(-1) * (
     &     - 7.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     + y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nh * (
     &     - 4.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc * (
     &     + 13.D0
     &     - 20.D0*upy**(-1)*upz**(-1)
     &     - 96.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 31.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + upy**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 31.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + upz**(-1)
     &     + 96.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 120.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc * (
     &     - 20.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*z*upy**(-1)
     &     - 48.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc * (
     &     - 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y*upz**(-1)
     &     + 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**2*Nh * (
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**3 * (
     &     - 8.D0
     &     + 8.D0*upy**(-1)*upz**(-1)
     &     + 96.D0*upy**(-1)*zm3**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*upy**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*upz**(-1)
     &     - 96.D0*zm3**(-2)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*zm3**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 120.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**3 * (
     &     + 20.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     + 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     + 48.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(0,0)*Nc**3 * (
     &     - 10.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(1,0)*Nc**(-3) * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(1,0)*Nc**(-1) * (
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(1,0)*Nc * (
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 192.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 168.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 168.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(1,0)*Nc * (
     &     - 16.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + Hr2(1,0)*Nc**3 * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc**(-3) * (
     &     - 2.D0
     &     + 2.D0*y**(-2)*upz**(-1)
     &     + y**(-2)
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)
     &     + y**(-1)*z
     &     + 12.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 4.D0*upz**(-1)
     &     - 4.D0*z*upy**(-2)
     &     - 4.D0*z*upy**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc**(-1) * (
     &     - 10.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     - y**(-2)*zp2**(-1)
     &     - 2.D0*y**(-2)
     &     + 6.D0*y**(-1)*upz**(-1)
     &     + y**(-1)*zp2**(-2)
     &     - y**(-1)*zp2**(-1)
     &     + y**(-1)
     &     - 2.D0*y**(-1)*z
     &     - 28.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - zp2**(-2)*ypzp2**(-1)
     &     + zp2**(-1)*ypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc**(-1) * (
     &     - 9.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-2)
     &     - z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc * (
     &     + 2.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     + 2.D0*y**(-2)*zp2**(-1)
     &     + y**(-2)
     &     - 2.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)*zp2**(-2)
     &     - 2.D0*y**(-1)*zp2**(-1)
     &     - 4.D0*y**(-1)
     &     + y**(-1)*z
     &     + 16.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upz**(-1)
     &     + 2.D0*zp2**(-2)*ypzp2**(-1)
     &     + 2.D0*zp2**(-1)*ypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc * (
     &     - 8.D0*ypzp2**(-2)
     &     + 34.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 22.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc**3 * (
     &     + 10.D0
     &     + 2.D0*y**(-2)*upz**(-1)
     &     - y**(-2)*zp2**(-1)
     &     - 6.D0*y**(-1)*upz**(-1)
     &     + y**(-1)*zp2**(-2)
     &     + 3.D0*y**(-1)*zp2**(-1)
     &     + y**(-1)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upz**(-1)
     &     - zp2**(-2)*ypzp2**(-1)
     &     - 3.D0*zp2**(-1)*ypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     - 25.D0*ypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*Nc**3 * (
     &     - 21.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*HBr1(-1)*Nc**(-1) * (
     &     + 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*HBr1(-1)*Nc * (
     &     - 64.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     - 64.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr1(-1)*HBr1(-1)*Nc**3 * (
     &     + 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(-1,-1)*Nc**(-1) * (
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 20.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(-1,-1)*Nc * (
     &     + 16.D0
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 24.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(-1,-1)*Nc**3 * (
     &     - 16.D0
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 12.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(0,-1)*Nc**(-3) * (
     &     - 4.D0
     &     - 4.D0*upy**(-3)
     &     - 8.D0*upy**(-2)
     &     + 12.D0*upy**(-1)*upz**(-1)
     &     + 2.D0*upy**(-1)
     &     - 8.D0*upz**(-1)
     &     - 4.D0*z*upy**(-3)
     &     - 2.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(0,-1)*Nc**(-1) * (
     &     + 4.D0*upy**(-3)
     &     + 4.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 6.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-3)
     &     + 4.D0*z*upy**(-1)
     &     - 2.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(0,-1)*Nc * (
     &     + 4.D0*upy**(-2)
     &     + 4.D0*upy**(-1)*upz**(-1)
     &     - 10.D0*upy**(-1)
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 16.D0*upz**(-1)
     &     + 10.D0*ypzp2**(-1)
     &     - 2.D0*z*upy**(-1)
     &     + 10.D0*z*ypzp2**(-1)
     &     + 12.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HAr2(0,-1)*Nc**3 * (
     &     + 4.D0
     &     - 8.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc**(-3) * (
     &     - 2.D0
     &     + 2.D0*z**(-2)*upy**(-1)
     &     + z**(-2)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + 2.D0*z**(-1)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upy**(-1)
     &     + 12.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + y*z**(-1)
     &     - 4.D0*y*upz**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc**(-1) * (
     &     - 11.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     - z**(-2)*ypzp2**(-1)
     &     - 2.D0*z**(-2)
     &     + 6.D0*z**(-1)*upy**(-1)
     &     - z**(-1)*ypzp2**(-1)
     &     + z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 28.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 7.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + z*ypzp2**(-1)
     &     - 2.D0*y*z**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc**(-1) * (
     &     + 4.D0*y*upz**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc * (
     &     + 16.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     + 2.D0*z**(-2)*ypzp2**(-1)
     &     + z**(-2)
     &     - 2.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)*ypzp2**(-1)
     &     - 4.D0*z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 8.D0*ypzp2**(-2)
     &     + 6.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc * (
     &     - 6.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + y*z**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr1(-1)*Nc**3 * (
     &     - 3.D0
     &     + 2.D0*z**(-2)*upy**(-1)
     &     - z**(-2)*ypzp2**(-1)
     &     - 6.D0*z**(-1)*upy**(-1)
     &     + 3.D0*z**(-1)*ypzp2**(-1)
     &     + z**(-1)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 5.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(-1,-1)*Nc**(-1) * (
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 20.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(-1,-1)*Nc * (
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 24.D0*upy**(-1)
     &     + 20.D0*upz**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(-1,-1)*Nc**3 * (
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upy**(-1)
     &     + 12.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(0,-1)*Nc**(-3) * (
     &     - 4.D0
     &     + 12.D0*upy**(-1)*upz**(-1)
     &     - 8.D0*upy**(-1)
     &     - 4.D0*upz**(-3)
     &     - 8.D0*upz**(-2)
     &     + 2.D0*upz**(-1)
     &     - 4.D0*z*upy**(-1)
     &     - 4.D0*y*upz**(-3)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(0,-1)*Nc**(-1) * (
     &     - 2.D0
     &     + 4.D0*upz**(-3)
     &     + 4.D0*upz**(-2)
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 8.D0*upz**(-1)
     &     + 10.D0*ypzp2**(-1)
     &     + 2.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-3)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(0,-1)*Nc * (
     &     + 10.D0
     &     - 12.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*upy**(-1)
     &     + 4.D0*upz**(-2)
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     - 10.D0*upz**(-1)
     &     - 10.D0*ypzp2**(-1)
     &     + 12.D0*z*upy**(-1)
     &     - 10.D0*z*ypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEP0T34 = PoleEP0T34 + HBr2(0,-1)*Nc**3 * (
     &     - 4.D0
     &     - 8.D0*upy**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     )
      return
      end

*
**************************** 72 characters *****************************
*
*

*
**************************** 72 characters *****************************
*
*
      double precision function PoleEPm2T14(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm2T14 = + Nc * (
     &     - 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm2T14 = PoleEPm2T14 + Nc**3 * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEPm2T14 = PoleEPm2T14 + Nc**5 * (
     &     - 8.D0
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     ) 
      return
      end
*
*
*
      double precision function PoleEPm1T14(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEPm1T14 =  + Nc**(-1) * (
     &     + 8.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 10.D0*upy**(-1)
     &     + 8.D0*upz**(-2)
     &     - 10.D0*upz**(-1)
     &     - 2.D0*z*upy**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nl * (
     &     + 16.D0/3.D0*upy**(-2)
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nc * (
     &     - 12.D0
     &     - 136.D0/3.D0*upy**(-2)
     &     - 272.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 146.D0/3.D0*upy**(-1)
     &     - 160.D0/3.D0*upz**(-2)
     &     + 176.D0/3.D0*upz**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     + 10.D0/3.D0*z*upy**(-1)
     &     - 8.D0*z*ypzp2**(-2)
     &     + 8.D0*z*ypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     + 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nc**2*Nl * (
     &     - 8.D0/3.D0
     &     - 16.D0/3.D0*upy**(-2)
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 32.D0/3.D0*upz**(-2)
     &     + 40.D0/3.D0*upz**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nc**3 * (
     &     + 68.D0/3.D0
     &     + 112.D0/3.D0*upy**(-2)
     &     + 224.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 116.D0/3.D0*upy**(-1)
     &     + 248.D0/3.D0*upz**(-2)
     &     - 262.D0/3.D0*upz**(-1)
     &     + 44.D0/3.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     + 88.D0/3.D0*z*ypzp2**(-2)
     &     - 88.D0/3.D0*z*ypzp2**(-1)
     &     + 44.D0/3.D0*z**2*ypzp2**(-2)
     &     - 14.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nc**4*Nl * (
     &     + 8.D0/3.D0
     &     + 16.D0/3.D0*upz**(-2)
     &     - 20.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     - 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Nc**5 * (
     &     - 32.D0/3.D0
     &     - 112.D0/3.D0*upz**(-2)
     &     + 116.D0/3.D0*upz**(-1)
     &     - 32.D0/3.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 64.D0/3.D0*z*ypzp2**(-2)
     &     + 64.D0/3.D0*z*ypzp2**(-1)
     &     - 32.D0/3.D0*z**2*ypzp2**(-2)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(-1)*Nc * (
     &     - 8.D0
     &     + 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 40.D0*upy**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(-1)*Nc**3 * (
     &     - 32.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 40.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(-1)*Nc**5 * (
     &     + 8.D0
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc**(-1) * (
     &     - 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 34.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 34.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc * (
     &     + 4.D0
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     - 52.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upy**(-1)
     &     + 8.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 34.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-1)
     &     + 22.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc * (
     &     + 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc**3 * (
     &     - 8.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-2)
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     + 18.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)
     &     + 8.D0*upz**(-2)
     &     - 10.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*ypzp2**(-2)
     &     - 12.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc**3 * (
     &     - 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-2)
     &     + 8.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzp2**(-2)
     &     + 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + Hr1(0)*Nc**5 * (
     &     - 4.D0
     &     - 8.D0*upz**(-2)
     &     + 10.D0*upz**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-2)
     &     + 8.D0*z*ypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     + 2.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + HAr1(-1)*Nc * (
     &     + 8.D0
     &     + 16.D0*upy**(-2)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + HAr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 16.D0*upy**(-2)
     &     + 20.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + HBr1(-1)*Nc * (
     &     - 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + HBr1(-1)*Nc**3 * (
     &     - 8.D0
     &     + 16.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 40.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEPm1T14 = PoleEPm1T14 + HBr1(-1)*Nc**5 * (
     &     + 8.D0
     &     + 16.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      return
      end
*
*
*
      double precision function PoleEP0T14(s,t,tm,w) 
      implicit none
      double precision Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,
     $                 HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,
     $                 HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,
     $                 en,w,tm,pi,the,s,t,u,x,xp,y,z,z2,z3,aln2,
     $                 upx,umx,upxp,umxp,upy,umy,upz,umz,zp2,zm3,
     $                 Nc,Nh,Nl,ypzp2,ypzm2,sqrypzp2,sqrypzm2
*
      complex*16 Hc1,Hc2,Hc3,Hc4
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1), 
     $          Hc4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1), 
     $          Hr4(-1:1,-1:1,-1:1,-1:1) 
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1), 
     $          Hi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HAc1,HAc2,HAc3,HAc4
      dimension HAc1(-1:1),HAc2(-1:1,-1:1),HAc3(-1:1,-1:1,-1:1), 
     $          HAc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAr1(-1:1),HAr2(-1:1,-1:1),HAr3(-1:1,-1:1,-1:1), 
     $          HAr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HAi1(-1:1),HAi2(-1:1,-1:1),HAi3(-1:1,-1:1,-1:1), 
     $          HAi4(-1:1,-1:1,-1:1,-1:1) 
*
      complex*16 HBc1,HBc2,HBc3,HBc4
      dimension HBc1(-1:1),HBc2(-1:1,-1:1),HBc3(-1:1,-1:1,-1:1), 
     $          HBc4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBr1(-1:1),HBr2(-1:1,-1:1),HBr3(-1:1,-1:1,-1:1), 
     $          HBr4(-1:1,-1:1,-1:1,-1:1) 
      dimension HBi1(-1:1),HBi2(-1:1,-1:1),HBi3(-1:1,-1:1,-1:1), 
     $          HBi4(-1:1,-1:1,-1:1,-1:1) 
*
*
      parameter ( pi   = 3.1415926535897932385d0 )
      parameter ( z2   = 1.6449340668482264365d0 )
      parameter ( z3   = 1.2020569031595942854d0 ) 
      parameter ( aln2 = 0.6931471805599453094d0 )
*
	u = 2.d0*tm**2 - s - t
*
        x  = 4.d0*tm**2/(Dsqrt(s-4.d0*tm**2)+Dsqrt(s))**2 
        y  = - t/tm**2
        z  = - u/tm**2
*
        upx = 1.d0+x
	umx = 1.d0-x
	upxp = 1.d0+xp
	umxp = 1.d0-xp
	upy = 1.d0+y
	umy = 1.d0-y
	upz = 1.d0+z
	umz = 1.d0-z
	zp2 = 2.d0+z
	zm3 = z-3.d0
	ypzp2 = y+z+2.d0
	ypzm2 = y+z-2.d0
	sqrypzp2 = Dsqrt(y+z+2.d0)
	sqrypzm2 = Dsqrt(y+z-2.d0)
*
        Nc = 3.d0
	Nl = 5.d0
	Nh = 1.d0
* 
       call hplog(x,4,Hc1,Hc2,Hc3,Hc4, 
     $                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1) 
*
       call hplog(y,4,HAc1,HAc2,HAc3,HAc4, 
     $                  HAr1,HAr2,HAr3,HAr4,HAi1,HAi2,HAi3,HAi4,-1,1) 
*
       call hplog(z,4,HBc1,HBc2,HBc3,HBc4, 
     $                  HBr1,HBr2,HBr3,HBr4,HBi1,HBi2,HBi3,HBi4,-1,1) 
*
*
      PoleEP0T14 =  + Nc**(-1) * (
     &     + 4.D0
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + z**(-1)
     &     + 12.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 21.D0*upy**(-1)
     &     + 12.D0*upz**(-2)
     &     - 21.D0*upz**(-1)
     &     + 4.D0*z2*upy**(-3)
     &     + 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upy**(-2)
     &     + 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**(-1) * (
     &     - 93.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z2*upy**(-1)
     &     + 4.D0*z2*upz**(-3)
     &     + 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upz**(-2)
     &     - 93.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z2*upz**(-1)
     &     + 62.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z2
     &     - 4.D0*z*upy**(-2)
     &     + z*upy**(-1)
     &     + 4.D0*z*z2*upy**(-3)
     &     - 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**(-1) * (
     &     - 3.D0*z*z2*upy**(-1)
     &     + z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-2)
     &     + y*upz**(-1)
     &     + 4.D0*y*z2*upz**(-3)
     &     - 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*y*z2*upz**(-1)
     &     + y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nl * (
     &     + 8.D0/3.D0
     &     + 8.D0/3.D0*upy**(-1)
     &     + 8.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*z*upy**(-1)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc * (
     &     - 53.D0/3.D0
     &     - 4.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     - z**(-1)*ypzp2**(-1)
     &     - 3.D0*z**(-1)
     &     - 28.D0*upy**(-2)
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     + 82.D0/3.D0*upy**(-1)
     &     - 40.D0*upz**(-2)
     &     + 133.D0/3.D0*upz**(-1)
     &     - 4.D0*ypzp2**(-2)
     &     - 10.D0*ypzp2**(-1)
     &     - 4.D0*z2*upy**(-3)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc * (
     &     - 48.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z2*upy**(-2)
     &     - 96.D0*z2*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 152.D0*z2*upy**(-1)*upz**(-1)
     &     - 64.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 132.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 44.D0*z2*upy**(-1)
     &     - 8.D0*z2*upz**(-3)
     &     - 24.D0*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 52.D0*z2*upz**(-2)
     &     - 140.D0*z2*upz**(-1)*ypzp2**(-1)
     &     + 107.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 37.D0*z2*upz**(-1)
     &     + 64.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc * (
     &     + 16.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 23.D0*z2*ypzp2**(-1)
     &     - 76.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 17.D0*z2
     &     + 4.D0*z*upy**(-2)
     &     - 38.D0/3.D0*z*upy**(-1)
     &     - 8.D0*z*ypzp2**(-2)
     &     + 7.D0*z*ypzp2**(-1)
     &     - 4.D0*z*z2*upy**(-3)
     &     + 48.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 52.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*z2*upy**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc * (
     &     + 16.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 27.D0*z*z2*ypzp2**(-1)
     &     - 10.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-2)
     &     - 8.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-2)
     &     - 41.D0/3.D0*y*upz**(-1)
     &     - 8.D0*y*z2*upz**(-3)
     &     + 24.D0*y*z2*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc * (
     &     - 11.D0*y*z2*upz**(-1)
     &     + 2.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 7.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**2*Nl * (
     &     - 8.D0/3.D0*upy**(-1)
     &     - 14.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     - 20.D0/3.D0*ypzp2**(-1)
     &     - 8.D0/3.D0*z*upy**(-1)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     - 16.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**2*Nh * (
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 2.D0/3.D0*upz**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     - 4.D0/3.D0*ypzp2**(-1)
     &     - 12.D0*z2*upz**(-1)*ypzp2**(-1)
     &     + 24.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**3 * (
     &     + 2.D0*y**(-1)*upz**(-1)
     &     + y**(-1)
     &     + 2.D0*z**(-1)*upy**(-1)
     &     + 2.D0*z**(-1)*ypzp2**(-1)
     &     + 3.D0*z**(-1)
     &     + 16.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 19.D0/3.D0*upy**(-1)
     &     + 44.D0*upz**(-2)
     &     + 8.D0*upz**(-1)*ypzm2**(-1)
     &     - 85.D0/3.D0*upz**(-1)
     &     - 7.D0*ypzm2**(-1)
     &     - 20.D0/3.D0*ypzp2**(-2)
     &     + 155.D0/3.D0*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**3 * (
     &     + 24.D0*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z2*upy**(-2)
     &     - 152.D0*z2*upy**(-1)*upz**(-1)
     &     + 64.D0*z2*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 39.D0*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 47.D0*z2*upy**(-1)
     &     + 4.D0*z2*upz**(-3)
     &     - 100.D0*z2*upz**(-2)
     &     + 96.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 48.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     + 140.D0*z2*upz**(-1)*ypzp2**(-1)
     &     - 9.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**3 * (
     &     + 85.D0*z2*upz**(-1)
     &     - 64.D0*z2*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 84.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 56.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z2*ypzp2**(-2)
     &     + 38.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 23.D0*z2*ypzp2**(-1)
     &     + 7.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 31.D0*z2
     &     + 35.D0/3.D0*z*upy**(-1)
     &     + 4.D0*z*ypzm2**(-1)
     &     - 40.D0/3.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**3 * (
     &     + 52.D0/3.D0*z*ypzp2**(-1)
     &     - 24.D0*z*z2*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*z*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z*z2*upy**(-1)
     &     + 48.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*z*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*z*z2*ypzp2**(-2)
     &     - 56.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 63.D0*z*z2*ypzp2**(-1)
     &     + 11.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*ypzm2**(-1)
     &     - 20.D0/3.D0*z**2*ypzp2**(-2)
     &     + z**2*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**3 * (
     &     + 3.D0*z**2*z2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z**2*z2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*z**2*z2*ypzp2**(-2)
     &     - 22.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-2)
     &     + 73.D0/3.D0*y*upz**(-1)
     &     + 4.D0*y*z2*upz**(-3)
     &     + 8.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 23.D0*y*z2*upz**(-1)
     &     - 4.D0*y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**4*Nl * (
     &     - 8.D0/3.D0
     &     + 2.D0*upz**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     + 20.D0/3.D0*ypzp2**(-1)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     + 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**4*Nh * (
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     - 2.D0/3.D0*upz**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 4.D0/3.D0*ypzp2**(-1)
     &     + 12.D0*z2*upz**(-1)*ypzp2**(-1)
     &     - 24.D0*z2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**5 * (
     &     + 41.D0/3.D0
     &     - z**(-1)*ypzp2**(-1)
     &     - z**(-1)
     &     - 16.D0*upz**(-2)
     &     - 8.D0*upz**(-1)*ypzm2**(-1)
     &     + 5.D0*upz**(-1)
     &     + 7.D0*ypzm2**(-1)
     &     + 32.D0/3.D0*ypzp2**(-2)
     &     - 125.D0/3.D0*ypzp2**(-1)
     &     + 40.D0*z2*upz**(-2)
     &     - 96.D0*z2*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     - 48.D0*z2*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*
     &    sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**5 * (
     &     - 5.D0*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 45.D0*z2*upz**(-1)
     &     + 84.D0*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z2*ypzp2**(-2)
     &     + 2.D0*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 7.D0*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*z2
     &     - 4.D0*z*ypzm2**(-1)
     &     + 64.D0/3.D0*z*ypzp2**(-2)
     &     - 73.D0/3.D0*z*ypzp2**(-1)
     &     - 48.D0*z*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*z*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*z*z2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**5 * (
     &     + 8.D0*z*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*z*z2*ypzp2**(-1)
     &     - 2.D0*z*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z**2*ypzm2**(-1)
     &     + 32.D0/3.D0*z**2*ypzp2**(-2)
     &     - z**2*ypzp2**(-1)
     &     + 12.D0*z**2*z2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*z2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*z**2*z2*ypzp2**(-2)
     &     + 6.D0*z**2*z2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 35.D0/3.D0*y*upz**(-1)
     &     - 4.D0*y*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*y*z2*upz**(-1)
     &     + y*z2*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Nc**5 * (
     &     - y**2*z2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nl * (
     &     - 16.D0/3.D0*upy**(-2)
     &     - 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 20.D0/3.D0*upy**(-1)
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     + 4.D0/3.D0*z*upy**(-1)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nc * (
     &     + 88.D0/3.D0*upy**(-2)
     &     + 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 110.D0/3.D0*upy**(-1)
     &     + 88.D0/3.D0*upz**(-2)
     &     - 110.D0/3.D0*upz**(-1)
     &     - 22.D0/3.D0*z*upy**(-1)
     &     - 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nc**2*Nl * (
     &     + 8.D0/3.D0
     &     + 16.D0/3.D0*upy**(-2)
     &     + 32.D0/3.D0*upy**(-1)*upz**(-1)
     &     - 20.D0/3.D0*upy**(-1)
     &     + 32.D0/3.D0*upz**(-2)
     &     - 40.D0/3.D0*upz**(-1)
     &     + 8.D0/3.D0*ypzp2**(-2)
     &     - 4.D0/3.D0*z*upy**(-1)
     &     + 16.D0/3.D0*z*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-1)
     &     + 8.D0/3.D0*z**2*ypzp2**(-2)
     &     - 8.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nc**3 * (
     &     - 44.D0/3.D0
     &     - 88.D0/3.D0*upy**(-2)
     &     - 176.D0/3.D0*upy**(-1)*upz**(-1)
     &     + 110.D0/3.D0*upy**(-1)
     &     - 176.D0/3.D0*upz**(-2)
     &     + 220.D0/3.D0*upz**(-1)
     &     - 44.D0/3.D0*ypzp2**(-2)
     &     + 22.D0/3.D0*z*upy**(-1)
     &     - 88.D0/3.D0*z*ypzp2**(-2)
     &     + 88.D0/3.D0*z*ypzp2**(-1)
     &     - 44.D0/3.D0*z**2*ypzp2**(-2)
     &     + 44.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nc**4*Nl * (
     &     - 8.D0/3.D0
     &     - 16.D0/3.D0*upz**(-2)
     &     + 20.D0/3.D0*upz**(-1)
     &     - 8.D0/3.D0*ypzp2**(-2)
     &     - 16.D0/3.D0*z*ypzp2**(-2)
     &     + 16.D0/3.D0*z*ypzp2**(-1)
     &     - 8.D0/3.D0*z**2*ypzp2**(-2)
     &     + 4.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Dlog(tm**(-2)*w**2)*Nc**5 * (
     &     + 44.D0/3.D0
     &     + 88.D0/3.D0*upz**(-2)
     &     - 110.D0/3.D0*upz**(-1)
     &     + 44.D0/3.D0*ypzp2**(-2)
     &     + 88.D0/3.D0*z*ypzp2**(-2)
     &     - 88.D0/3.D0*z*ypzp2**(-1)
     &     + 44.D0/3.D0*z**2*ypzp2**(-2)
     &     - 22.D0/3.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*Nc * (
     &     + 16.D0
     &     + 64.D0*upy**(-1)*zm3**(-1)
     &     + 32.D0*upy**(-1)
     &     - 64.D0*zm3**(-1)*ypzm2**(-1)
     &     - 16.D0*ypzm2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 16.D0*z*ypzm2**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*Nc**3 * (
     &     - 10.D0
     &     - 64.D0*upy**(-1)*zm3**(-1)
     &     - 32.D0*upy**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-2)
     &     - 32.D0*upz**(-1)*ypzm2**(-1)
     &     - 6.D0*upz**(-1)
     &     + 64.D0*zm3**(-1)*ypzm2**(-1)
     &     + 84.D0*ypzm2**(-2)
     &     + 42.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     - 48.D0*z*ypzm2**(-2)
     &     + 10.D0*z*ypzm2**(-1)
     &     + 20.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*Nc**3 * (
     &     + 12.D0*z**2*ypzm2**(-2)
     &     - 4.D0*z**2*ypzm2**(-1)
     &     + 4.D0*z**2*ypzp2**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*Nc**5 * (
     &     - 6.D0
     &     + 96.D0*upz**(-1)*ypzm2**(-2)
     &     + 32.D0*upz**(-1)*ypzm2**(-1)
     &     + 6.D0*upz**(-1)
     &     - 84.D0*ypzm2**(-2)
     &     - 26.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     + 48.D0*z*ypzm2**(-2)
     &     + 6.D0*z*ypzm2**(-1)
     &     - 4.D0*z*ypzp2**(-1)
     &     - 12.D0*z**2*ypzm2**(-2)
     &     + 4.D0*z**2*ypzm2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*HAr1(-1)*Nc * (
     &     - 16.D0
     &     - 64.D0*upy**(-2)
     &     + 40.D0*upy**(-1)
     &     - 32.D0*ypzp2**(-2)
     &     + 8.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z*ypzp2**(-1)
     &     - 32.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*HAr1(-1)*Nc**3 * (
     &     + 16.D0
     &     + 64.D0*upy**(-2)
     &     - 40.D0*upy**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     - 8.D0*z*upy**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-1)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*HBr1(-1)*Nc * (
     &     + 8.D0
     &     - 64.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)
     &     + 32.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)
     &     - 16.D0*z*upy**(-1)
     &     + 64.D0*z*ypzp2**(-2)
     &     + 32.D0*z**2*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*HBr1(-1)*Nc**3 * (
     &     + 64.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)
     &     + 16.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(-1)*HBr1(-1)*Nc**5 * (
     &     - 8.D0
     &     - 32.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**(-1) * (
     &     - 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc * (
     &     - 8.D0
     &     - 32.D0*upy**(-1)*zm3**(-1)
     &     + 16.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upy**(-1)
     &     + upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*zm3**(-1)*ypzm2**(-1)
     &     + 8.D0*ypzm2**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-1)
     &     - 4.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 8.D0*z*ypzm2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc * (
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)
     &     - 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 13.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**2*Nh * (
     &     - 32.D0*upz**(-1)*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**3 * (
     &     + 5.D0
     &     + 32.D0*upy**(-1)*zm3**(-1)
     &     - 8.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upy**(-1)
     &     + 48.D0*upz**(-1)*ypzm2**(-2)
     &     + 16.D0*upz**(-1)*ypzm2**(-1)
     &     + 7.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upz**(-1)
     &     - 32.D0*zm3**(-1)*ypzm2**(-1)
     &     - 42.D0*ypzm2**(-2)
     &     - 21.D0*ypzm2**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**3 * (
     &     - 4.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 24.D0*z*ypzm2**(-2)
     &     - 5.D0*z*ypzm2**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*z*ypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z**2*ypzm2**(-2)
     &     + 2.D0*z**2*ypzm2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**3 * (
     &     - 2.D0*z**2*ypzp2**(-1)
     &     - y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     + y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**4*Nh * (
     &     + 32.D0*upz**(-1)*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*Nc**5 * (
     &     + 3.D0
     &     - 48.D0*upz**(-1)*ypzm2**(-2)
     &     - 16.D0*upz**(-1)*ypzm2**(-1)
     &     - 3.D0*upz**(-1)
     &     + 42.D0*ypzm2**(-2)
     &     + 13.D0*ypzm2**(-1)
     &     - 4.D0*ypzp2**(-1)
     &     - 24.D0*z*ypzm2**(-2)
     &     - 3.D0*z*ypzm2**(-1)
     &     + 2.D0*z*ypzp2**(-1)
     &     + 6.D0*z**2*ypzm2**(-2)
     &     - 2.D0*z**2*ypzm2**(-1)
     &     + 2.D0*z**2*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HAr1(-1)*Nc**(-1) * (
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 42.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HAr1(-1)*Nc * (
     &     + 8.D0
     &     - 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-2)
     &     + 52.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 20.D0*upy**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HAr1(-1)*Nc * (
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     - 8.D0
     &     + 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 20.D0*upy**(-1)
     &     + 44.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HAr1(-1)*Nc**3 * (
     &     + 16.D0*z*ypzp2**(-1)
     &     + 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc**(-1) * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 42.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 18.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc * (
     &     - 4.D0
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 88.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 58.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 28.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     - 32.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc * (
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     + 32.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 44.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*upy**(-1)
     &     - 16.D0*upz**(-2)
     &     - 16.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc**3 * (
     &     + 8.D0*z*ypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     + 2.D0*y*upz**(-1)
     &     + 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr1(0)*HBr1(-1)*Nc**5 * (
     &     + 4.D0
     &     + 16.D0*upz**(-2)
     &     - 10.D0*upz**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-2)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,-1)*Nc * (
     &     + 24.D0
     &     + 56.D0*upy**(-1)
     &     - 64.D0*ypzp2**(-1)
     &     + 24.D0*z*upy**(-1)
     &     - 64.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,-1)*Nc**3 * (
     &     - 16.D0
     &     - 56.D0*upy**(-1)
     &     - 20.D0*upz**(-1)
     &     + 64.D0*ypzp2**(-1)
     &     - 24.D0*z*upy**(-1)
     &     + 48.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,-1)*Nc**5 * (
     &     - 8.D0
     &     + 20.D0*upz**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 86.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 86.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 36.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 6.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc * (
     &     - 12.D0
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 140.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*upy**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 58.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)
     &     + 44.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc * (
     &     + 20.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc**3 * (
     &     + 8.D0
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 54.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*upy**(-1)
     &     - 28.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     - 8.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z*upy**(-1)
     &     - 24.D0*z*ypzp2**(-1)
     &     - 6.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc**3 * (
     &     - 6.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(-1,0)*Nc**5 * (
     &     + 4.D0
     &     - 10.D0*upz**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc * (
     &     - 12.D0
     &     + 128.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 28.D0*upy**(-1)
     &     - 128.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*ypzp2**(-1)
     &     + 20.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z*upy**(-1)
     &     - 32.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc * (
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc**3 * (
     &     + 8.D0
     &     - 128.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*upy**(-1)
     &     - 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 96.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)
     &     + 128.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 168.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 112.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*ypzp2**(-1)
     &     - 6.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc**3 * (
     &     - 8.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z*upy**(-1)
     &     - 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*z*ypzp2**(-1)
     &     - 12.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y*upz**(-1)
     &     + 6.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc**3 * (
     &     - 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc**5 * (
     &     + 4.D0
     &     + 192.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 10.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 10.D0*upz**(-1)
     &     - 168.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 80.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 96.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 28.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     + 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,-1)*Nc**5 * (
     &     - 24.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**(-1) * (
     &     + 4.D0
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     - 9.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upy**(-1)
     &     - 9.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*upz**(-1)
     &     - 2.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*z*upy**(-1)
     &     + 5.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 3.D0*y*upz**(-1)
     &     + 5.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc * (
     &     - 1.D0
     &     + 16.D0*upy**(-1)*upz**(-1)
     &     - 64.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*upy**(-1)
     &     - 12.D0*upz**(-1)*ypzp2**(-1)
     &     - 5.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + upz**(-1)
     &     + 64.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 17.D0*ypzp2**(-1)
     &     - 10.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc * (
     &     - 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 13.D0*z*ypzp2**(-1)
     &     - 10.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y*upz**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**2*Nh * (
     &     + 4.D0*upz**(-1)*ypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**3 * (
     &     - 1.D0
     &     - 8.D0*upy**(-1)*upz**(-1)
     &     + 64.D0*upy**(-1)*zm3**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 11.D0*upy**(-1)
     &     + 96.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 48.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 12.D0*upz**(-1)*ypzp2**(-1)
     &     + 19.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 9.D0*upz**(-1)
     &     - 64.D0*zm3**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 84.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 56.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**3 * (
     &     - 38.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 17.D0*ypzp2**(-1)
     &     + 5.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*z*upy**(-1)
     &     + 48.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 9.D0*z*ypzp2**(-1)
     &     + 7.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 12.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**3 * (
     &     - 8.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 3.D0*y*upz**(-1)
     &     - 2.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**4*Nh * (
     &     - 4.D0*upz**(-1)*ypzp2**(-1)
     &     + 8.D0*ypzp2**(-2)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**5 * (
     &     - 2.D0
     &     - 96.D0*upz**(-1)*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*upz**(-1)*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 5.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 5.D0*upz**(-1)
     &     + 84.D0*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 40.D0*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 2.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 7.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 48.D0*z*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 14.D0*z*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*z*ypzp2**(-1)
     &     - 2.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(0,0)*Nc**5 * (
     &     + 12.D0*z**2*ypzm2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 2.D0*z**2*ypzm2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 6.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + y*upz**(-1)
     &     + y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(1,0)*Nc**(-1) * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 68.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 40.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(1,0)*Nc * (
     &     + 32.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 64.D0*upy**(-1)*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 104.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 68.D0*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 44.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(1,0)*Nc * (
     &     + 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*y*upz**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*y*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 4.D0*y**2*upz**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + Hr2(1,0)*Nc**3 * (
     &     - 16.D0*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 36.D0*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 24.D0*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*upy**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 32.D0*z*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 16.D0*z*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 4.D0*z**2*upy**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     - 16.D0*z**2*ypzp2**(-2)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     + 8.D0*z**2*ypzp2**(-1)*sqrypzm2**(-1)*sqrypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*Nc**(-1) * (
     &     + 2.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     - y**(-2)
     &     - 2.D0*y**(-1)*upz**(-1)
     &     - 2.D0*y**(-1)
     &     - y**(-1)*z
     &     - 12.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 4.D0*upz**(-1)
     &     + 4.D0*z*upy**(-2)
     &     + 4.D0*z*upy**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*Nc * (
     &     - 8.D0
     &     + 4.D0*y**(-2)*upz**(-1)
     &     + 2.D0*y**(-2)
     &     - 4.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)*z
     &     + 28.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 24.D0*upy**(-1)
     &     - 4.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-2)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-2)
     &     - 32.D0*z*ypzp2**(-2)
     &     + 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*Nc * (
     &     - 16.D0*z**2*ypzp2**(-2)
     &     + 8.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*Nc**3 * (
     &     + 6.D0
     &     - 2.D0*y**(-2)*upz**(-1)
     &     - y**(-2)
     &     + 6.D0*y**(-1)*upz**(-1)
     &     + 2.D0*y**(-1)
     &     - y**(-1)*z
     &     - 16.D0*upy**(-2)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upy**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     - 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 32.D0*z*ypzp2**(-2)
     &     - 16.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*Nc**3 * (
     &     + 16.D0*z**2*ypzp2**(-2)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*HBr1(-1)*Nc * (
     &     - 8.D0
     &     + 32.D0*upy**(-2)
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     - 32.D0*upz**(-2)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr1(-1)*HBr1(-1)*Nc**3 * (
     &     + 8.D0
     &     - 32.D0*upy**(-2)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     + 32.D0*upz**(-2)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr2(-1,-1)*Nc * (
     &     + 8.D0
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     + 20.D0*upy**(-1)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     - 48.D0*ypzp2**(-1)
     &     + 4.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr2(-1,-1)*Nc**3 * (
     &     - 8.D0
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 20.D0*upz**(-1)
     &     + 48.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr2(0,-1)*Nc**(-1) * (
     &     + 4.D0
     &     + 4.D0*upy**(-3)
     &     + 8.D0*upy**(-2)
     &     - 12.D0*upy**(-1)*upz**(-1)
     &     - 2.D0*upy**(-1)
     &     + 8.D0*upz**(-1)
     &     + 4.D0*z*upy**(-3)
     &     + 2.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr2(0,-1)*Nc * (
     &     - 8.D0
     &     - 4.D0*upy**(-3)
     &     - 8.D0*upy**(-2)
     &     + 12.D0*upy**(-1)*upz**(-1)
     &     - 8.D0*upy**(-1)
     &     - 8.D0*upz**(-1)
     &     + 4.D0*ypzp2**(-1)
     &     - 4.D0*z*upy**(-3)
     &     - 4.D0*z*upy**(-1)
     &     + 4.D0*z*ypzp2**(-1)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HAr2(0,-1)*Nc**3 * (
     &     + 4.D0
     &     + 10.D0*upy**(-1)
     &     - 4.D0*ypzp2**(-1)
     &     + 2.D0*z*upy**(-1)
     &     - 4.D0*z*ypzp2**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc**(-1) * (
     &     + 2.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     - z**(-2)
     &     - 2.D0*z**(-1)*upy**(-1)
     &     - 2.D0*z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upy**(-1)
     &     - 12.D0*upz**(-2)
     &     + 20.D0*upz**(-1)
     &     - y*z**(-1)
     &     + 4.D0*y*upz**(-2)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc * (
     &     + 9.D0
     &     + 4.D0*z**(-2)*upy**(-1)
     &     + z**(-2)*ypzp2**(-1)
     &     + 3.D0*z**(-2)
     &     - 4.D0*z**(-1)*upy**(-1)
     &     + z**(-1)*ypzp2**(-1)
     &     + z**(-1)
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 20.D0*upy**(-1)
     &     + 40.D0*upz**(-2)
     &     - 40.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-2)
     &     + 7.D0*ypzp2**(-1)
     &     - 8.D0*z*upy**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc * (
     &     + 32.D0*z*ypzp2**(-2)
     &     - z*ypzp2**(-1)
     &     + 16.D0*z**2*ypzp2**(-2)
     &     + 3.D0*y*z**(-1)
     &     - 8.D0*y*upz**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc**3 * (
     &     - 8.D0
     &     - 2.D0*z**(-2)*upy**(-1)
     &     - 2.D0*z**(-2)*ypzp2**(-1)
     &     - 3.D0*z**(-2)
     &     + 6.D0*z**(-1)*upy**(-1)
     &     + 2.D0*z**(-1)*ypzp2**(-1)
     &     + 4.D0*z**(-1)
     &     - 16.D0*upy**(-1)*upz**(-1)
     &     + 16.D0*upy**(-1)
     &     - 44.D0*upz**(-2)
     &     + 24.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     - 22.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc**3 * (
     &     - 16.D0*z*ypzp2**(-2)
     &     - 10.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     - 3.D0*y*z**(-1)
     &     + 4.D0*y*upz**(-2)
     &     - 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr1(-1)*Nc**5 * (
     &     - 3.D0
     &     + z**(-2)*ypzp2**(-1)
     &     + z**(-2)
     &     - 3.D0*z**(-1)*ypzp2**(-1)
     &     - 3.D0*z**(-1)
     &     + 16.D0*upz**(-2)
     &     - 4.D0*upz**(-1)
     &     - 8.D0*ypzp2**(-2)
     &     + 15.D0*ypzp2**(-1)
     &     - 16.D0*z*ypzp2**(-2)
     &     + 11.D0*z*ypzp2**(-1)
     &     - 8.D0*z**2*ypzp2**(-2)
     &     + y*z**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(-1,-1)*Nc * (
     &     + 32.D0*upy**(-1)*upz**(-1)
     &     - 4.D0*upy**(-1)
     &     - 64.D0*upz**(-1)*ypzp2**(-1)
     &     + 20.D0*upz**(-1)
     &     + 16.D0*ypzp2**(-1)
     &     + 12.D0*z*upy**(-1)
     &     - 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(-1,-1)*Nc**3 * (
     &     + 8.D0
     &     - 32.D0*upy**(-1)*upz**(-1)
     &     + 4.D0*upy**(-1)
     &     + 64.D0*upz**(-1)*ypzp2**(-1)
     &     - 40.D0*upz**(-1)
     &     - 16.D0*ypzp2**(-1)
     &     - 12.D0*z*upy**(-1)
     &     - 8.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(-1,-1)*Nc**5 * (
     &     - 8.D0
     &     + 20.D0*upz**(-1)
     &     + 16.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(0,-1)*Nc**(-1) * (
     &     + 4.D0
     &     - 12.D0*upy**(-1)*upz**(-1)
     &     + 8.D0*upy**(-1)
     &     + 4.D0*upz**(-3)
     &     + 8.D0*upz**(-2)
     &     - 2.D0*upz**(-1)
     &     + 4.D0*z*upy**(-1)
     &     + 4.D0*y*upz**(-3)
     &     + 2.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(0,-1)*Nc * (
     &     - 2.D0
     &     + 12.D0*upy**(-1)*upz**(-1)
     &     - 16.D0*upy**(-1)
     &     - 8.D0*upz**(-3)
     &     - 12.D0*upz**(-2)
     &     + 16.D0*upz**(-1)*ypzp2**(-1)
     &     - 6.D0*upz**(-1)
     &     - 6.D0*ypzp2**(-1)
     &     - 12.D0*z*upy**(-1)
     &     + 2.D0*z*ypzp2**(-1)
     &     - 8.D0*y*upz**(-3)
     &     - 6.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(0,-1)*Nc**3 * (
     &     - 6.D0
     &     + 8.D0*upy**(-1)
     &     + 4.D0*upz**(-3)
     &     + 4.D0*upz**(-2)
     &     - 16.D0*upz**(-1)*ypzp2**(-1)
     &     + 18.D0*upz**(-1)
     &     + 6.D0*ypzp2**(-1)
     &     + 8.D0*z*upy**(-1)
     &     + 6.D0*z*ypzp2**(-1)
     &     + 4.D0*y*upz**(-3)
     &     + 6.D0*y*upz**(-1)
     &     )
      PoleEP0T14 = PoleEP0T14 + HBr2(0,-1)*Nc**5 * (
     &     + 4.D0
     &     - 10.D0*upz**(-1)
     &     - 8.D0*z*ypzp2**(-1)
     &     - 2.D0*y*upz**(-1)
     &     )
      return
      end
