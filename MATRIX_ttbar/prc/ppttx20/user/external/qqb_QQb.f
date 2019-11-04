      subroutine qqb_QQb(p,msq,iswitch) 
      implicit none

************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C      This is the four dimensional result for                         *
C      Quark antiquark annihilation in order alfa_s^2                  *
C      q(P1) + qbar(P2) --> Q(-P3) + Qbar(-P4)                         *
************************************************************************
      include 'masses.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
c      include 'breit.f'
      
      integer j,k,cs,iswitch
      logical first
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision wtqqb,wtgg,t1,t2,ro
      data first/.true./
      save first
chhhhhhhh
      double precision s12

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mt
      endif 

c---statement function
c      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      do cs=0,2
      msq_cs(cs,j,k)=0d0
      enddo
      enddo
      enddo
      call dotem(4,p,s)
ch      s12=5d0*mt**2
ch      t1=0.45d0
ch      t2=0.55d0
      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
cch      write(*,*)t1+t2
c      write(*,*)t1,t2,s(1,2),'t1'
c      t1=0.497496291880528d0
c      t2=0.502503708119500d0
c      s(1,2)=154760.880713210d0
c      ro=0.992359547138604d0
      ro=4d0*mt**2/s(1,2)
ch      ro=4d0*mt**2/s12
ch      write(*,*)s(1,2),s(1,3),s(2,3)
ch      write(*,*)gsq,'a'

c      write(*,*)p(1,4),'a1'
c      write(*,*)p(1,3),'a2'
c      write(*,*)p(1,2),'a3'
c      write(*,*)p(1,1),'a4'
c      write(*,*)p(3,4),'a5'
c      write(*,*)p(3,3),'a6'
c      write(*,*)p(3,2),'a7'
c      write(*,*)p(3,1),'a8'
c      write(*,*)-s(1,3),'s13'
c      write(*,*)-2*(p(1,4)*p(3,4)-p(1,1)*p(3,1)-p(1,2)
c     .           *p(3,2)-p(1,3)*p(3,3)),'p1p3'

      wtqqb=gsq**2*4d0/9d0*(t1**2+t2**2+ro/2d0)
c      wtqqb=4d0/9d0*(t1**2+t2**2+ro/2d0)*16d0*pi**2
c      write(*,*)wtqqb/gsq**2*16d0*pi**2
c      write(*,*)gsq,'a',dsqrt(s(3,4)+2*mt**2)

c      wtgg=gsq**2*(1d0/6d0/t1/t2-3d0/8d0)
c     . *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
      msq_cs(1,0,0)=V*xn*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t2)-4d0*t1**2-0.5d0*(ro/t2)**2)
      msq_cs(2,0,0)=V*xn*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t1)-4d0*t2**2-0.5d0*(ro/t1)**2)
      msq_cs(0,0,0)=-V/xn*2d0*gsq**2
     . *(-2d0+(1d0+ro*(1d0-0.5d0*ro))/t1/t2
     .         -0.25d0*(ro/t1)**2-0.25d0*(ro/t2)**2)
ch      write(*,*)msq_cs(0,0,0)/V/gsq**2*4d0*xn,'Bg'
ch      write(*,*)(msq_cs(1,0,0)+msq_cs(2,0,0))/V/gsq**2*4d0/xn,'Ag'
      if(iswitch.eq.1)then
       wtgg=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)
ch       wtgg=wtgg/gsq**2
ch       write(*,*),msq_cs(1,0,0),msq_cs(2,0,0),
ch     .  msq_cs(0,0,0),'l'
      elseif(iswitch.eq.0)then
       wtgg=16d0*gsq**2*(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))/2/t1/t2
      endif
ch      write(*,*)(msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0))
ch     .  *avegg
c      write(*,*)wtgg/gsq**2*16d0*pi**2,wtqqb/gsq**2*16d0*pi**2,'d'

C---fill qb-q, gg and q-qb elements
c--- the msq_cs entries for qqb and qbq are arbitrary
c--- divisions that are needed for summing up the _z contribution
c--- in virtint
      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k.eq.0)) then
          msq(j,k)=wtgg*avegg
ch          msq(j,k)=wtgg
      elseif ((j .gt. 0) .and. (k.lt.0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3d0
          msq_cs(1,j,k)=wtqqb/3d0
          msq_cs(2,j,k)=wtqqb/3d0
      elseif ((j .lt. 0) .and. (k.gt.0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3d0
          msq_cs(1,j,k)=wtqqb/3d0
          msq_cs(2,j,k)=wtqqb/3d0
      endif
      enddo

      return
      end
