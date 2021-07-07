C     Finite part of <M0|M2>+cc according to 1312.6279
C     The normalisation is defined according to eqs. (3.6)
C
C     qqbar channel                                                   

      double precision function virtqq(beta,ctheta)
      implicit none
      double precision beta,ctheta,x,ms,qqhe,qqth,y,Pi,nf

c     Hard-coding nf in this contribution to avoid confusion,
c     since anyway cannot be changed in the soft function
      nf = 5.d0

      Pi=3.14159265359d0

      x=0.5d0*(1-beta*ctheta)
      ms=1d0/4*(1-beta**2)

      if(beta.gt.0.999d0) then
C
C     Call high-energy expansion
C
         virtqq=beta*(1-beta**2)/576/Pi*qqhe(x,ms,nf)

      elseif(beta.lt.0.0125d0) then

C     Call threshold expansion
C
         virtqq=beta*(1-beta**2)/576/Pi*qqth(beta,ctheta,nf)

      else
C
C     Call grid interpolation
C
      call interpolqq(beta,ctheta,virtqq,nf)
      endif
      return
      end

      subroutine interpolqq(x1,x2,y,nf)
      implicit none
      double precision x1,x2,x1a(1:80),x2a(1:42)
      double precision ya(1:80,1:42),y2a(1:80,1:42)
      double precision y,dy,dummy,nf
      integer n,m,i,j,nflav

      nflav = NINT(real(nf))

      m=80
      n=42

      if(nflav.eq.5) then
        open(unit=11,file='../grids/gridqq.dat',status='old')
      elseif(nflav.eq.4) then
        open(unit=11,file='../grids/gridqq_nf4.dat',status='old')
      elseif(nflav.eq.3) then
        open(unit=11,file='../grids/gridqq_nf3.dat',status='old')
      endif

      do i=1,m
      do j=1,n
         read(11,*)x1a(i),x2a(j),ya(i,j),dummy
      enddo
      enddo

      close(11)

C     Computes second derivatives
      
      call splie2(x1a,x2a,ya,m,n,y2a)

C     Performs interpolation
      
      call splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)

      return
      end

 
