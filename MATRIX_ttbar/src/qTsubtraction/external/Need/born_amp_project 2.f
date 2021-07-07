      subroutine project_amp(p)
      include 'projected_amplitudes.f'
      include 'constants.f'
      include 'masses.f'
      include 'inner_prod.f'
      double precision p(mxpart,4),dot,v34,q2,t1,u1,
     .  M0sqeps
      external dot

        q2=2d0*dot(p,3,4)+2d0*mt**2
        t1=-2d0*dot(p,1,3)
        u1=-2d0*dot(p,2,3)

        M0sq=(8d0*q2**2/u1/t1/4d0/xn**2*(
     .            (u1**2 + t1**2)/q2**2 +
     .            4d0*mt**2/q2 -
     .            4d0*mt**4/u1/t1)*16d0*pi**2)*c1gs
     .      +(t1-u1)**2/q2**2*(8d0*q2**2/u1/t1/4d0*(
     .            (u1**2 + t1**2)/q2**2 +
     .            4d0*mt**2/q2 -
     .            4d0*mt**4/u1/t1)*16d0*pi**2)*c2gs
     .      +(8d0*q2**2/u1/t1/4d0*(
     .            (u1**2 + t1**2)/q2**2 +
     .            4d0*mt**2/q2 -
     .            4d0*mt**4/u1/t1)*16d0*pi**2)*c3gs

        M0sqeps=(-4d0*(-1d0 + xn**2)*
     -    (-(q2**2*(t1**2 + q2*u1)) + 
     -      (q2**2*t1*(t1 - 3d0*u1) + q2**3*u1 - 
     -         2d0*t1**3*u1 + q2*t1*u1*(3d0*t1 + u1)
     -         )*xn**2))/(q2**2*t1*u1*xn)

        Mgs0eps(1,1)=(-4*(q2**2 - q2*u1 + u1**2))/
     -  ((q2 - u1)*u1*xn**2)
        Mgs0eps(3,3)=xn**2*Mgs0eps(1,1)
        Mgs0eps(2,2)=xn**2*(t1-u1)**2/q2**2*Mgs0eps(1,1)
        Mgs0eps(1,2)=xn*(t1-u1)/q2*Mgs0eps(1,1)
        Mgs0eps(2,3)=xn**2*(t1-u1)/q2*Mgs0eps(1,1)
        Mgs0eps(1,3)=xn*Mgs0eps(1,1)


ch      Checked !!
ch      M0sqeps=ggQQv_1/gsq**2/avegg
ch      Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs+Mgs0eps(3,3)*c3gs=Mosqeps

ch      write(*,*)Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs+
ch     .                Mgs0eps(3,3)*c3gs,'a',
ch     .          M0sqeps,'aa',q2,t1,u1,'d',
ch     .          Mgs0eps(1,1),Mgs0eps(2,2),Mgs0eps(3,3),'e',
ch     .          c1gs,c2gs,c3gs,'f'

        Mgs0(1,1)=8d0*q2**2/u1/t1/4d0/xn**2*(
     .            (u1**2 + t1**2)/q2**2 +
     .            4d0*mt**2/q2 -
     .            4d0*mt**4/u1/t1)*16d0*pi**2/M0sq

        Mgs0(3,3)=xn**2*Mgs0(1,1)
        Mgs0(2,2)=xn**2*(t1-u1)**2/q2**2*Mgs0(1,1)
        Mgs0(1,2)=xn*(t1-u1)/q2*Mgs0(1,1)
        Mgs0(2,3)=xn**2*(t1-u1)/q2*Mgs0(1,1)
        Mgs0(1,3)=xn*Mgs0(1,1)

ch      als^2*avegg*Mosq=msqc(0,0)
ch        write(*,*)M0sq/16d0/pi**2,'aa'
ch     .     Mgs0(1,1)*c1gs+Mgs0(2,2)*c2gs+Mgs0(3,3)*c3gs,'a'

      return
      end

      subroutine col_operators(p)

      include 'projected_amplitudes.f'
      include 'constants.f'
      include 'masses.f'
      include 'inner_prod.f'
      include 'born_col_correl.f'

      call project_amp(p)

      Tgg(1,1)=ca
      Tgg(2,2)=ca
      Tgg(3,3)=cf
      Tgg(4,4)=cf
      Tgg(3,4)=-cf*Mgs0(1,1)*c1gs+
     .     1d0/2d0/xn*(Mgs0(2,2)*c2gs+Mgs0(3,3)*c3gs)
      Tgg(1,4)=-0.5d0*Mgs0(1,2)*c1gs
     .          -(Mgs0(1,2)+ca/4d0*Mgs0(2,2)
     .          +(xn**2-4d0)/4d0/xn*Mgs0(2,3))*c2gs
     .          -ca/4d0*(Mgs0(2,3)+Mgs0(3,3))*c3gs
      Tgg(1,3)=0.5d0*Mgs0(1,2)*c1gs
     .          +(Mgs0(1,2)-ca/4d0*Mgs0(2,2)
     .          +(xn**2-4d0)/4d0/xn*Mgs0(2,3))*c2gs
     .          +ca/4d0*(Mgs0(2,3)-Mgs0(3,3))*c3gs
      Tgg(2,3)=Tgg(1,4)
      Tgg(2,4)=Tgg(1,3)

      Tgg4(2,3,2,3)=c1gs*(0.5d0*Mgs0(1,1)
     .                    +xn/8d0*Mgs0(1,2)
     .                    +(xn**2-4d0)/8d0/xn*Mgs0(1,3))
     .             +c2gs*(xn/4d0*Mgs0(1,2)
     .                    +(xn**2+2d0)/8d0*Mgs0(2,2)
     .                    +(xn**2-4d0)/8d0*Mgs0(2,3))
     .             +c3gs*(xn/4d0*Mgs0(1,3)
     .                    +xn**2/8d0*Mgs0(2,3)
     .                    +(xn**2-2d0)/8d0*Mgs0(3,3))
      Tgg4(1,3,1,3)=c1gs*(0.5d0*Mgs0(1,1)
     .                    -xn/8d0*Mgs0(1,2)
     .                    +(xn**2-4d0)/8d0/xn*Mgs0(1,3))
     .             +c2gs*(-xn/4d0*Mgs0(1,2)
     .                    +(xn**2+2d0)/8d0*Mgs0(2,2)
     .                    -(xn**2-4d0)/8d0*Mgs0(2,3))
     .             +c3gs*(xn/4d0*Mgs0(1,3)
     .                    -xn**2/8d0*Mgs0(2,3)
     .                    +(xn**2-2d0)/8d0*Mgs0(3,3))
      Tgg4(2,3,1,3)=c1gs*(-0.5d0*Mgs0(1,1)
     .                    +xn/8d0*Mgs0(1,2)
     .                    -(xn**2-4d0)/8d0/xn*Mgs0(1,3))
     .             +c2gs*(-xn/4d0*Mgs0(1,2)
     .                    -1/4d0*Mgs0(2,2))
     .             +c3gs*(-xn/4d0*Mgs0(1,3)
     .                    +1/4d0*Mgs0(3,3))
      Tgg4(1,3,2,3)=c1gs*(-0.5d0*Mgs0(1,1)
     .                    -xn/8d0*Mgs0(1,2)
     .                    -(xn**2-4d0)/8d0/xn*Mgs0(1,3))
     .             +c2gs*(xn/4d0*Mgs0(1,2)
     .                    -1/4d0*Mgs0(2,2))
     .             +c3gs*(-xn/4d0*Mgs0(1,3)
     .                    +1/4d0*Mgs0(3,3))
      Tgg4(2,3,3,4)=c1gs*(-1d0/xn/4d0*Mgs0(1,2))
     .             +c2gs*(cf*Mgs0(1,2)
     .                    -1/8d0*Mgs0(2,2)
     .                    -(xn**2-4d0)/8d0/xn**2*Mgs0(2,3))
     .             +c3gs*(-1/8d0*Mgs0(2,3)
     .                    -1/8d0*Mgs0(3,3))
      Tgg4(1,3,3,4)=c1gs*(1d0/xn/4d0*Mgs0(1,2))
     .             +c2gs*(-cf*Mgs0(1,2)
     .                    -1/8d0*Mgs0(2,2)
     .                    +(xn**2-4d0)/8d0/xn**2*Mgs0(2,3))
     .             +c3gs*(1/8d0*Mgs0(2,3)
     .                    -1/8d0*Mgs0(3,3))
      Tgg4(3,4,3,4)=c1gs*(cf**2*Mgs0(1,1))
     .             +c2gs*(1/4d0/xn**2*Mgs0(2,2))
     .             +c3gs*(1/4d0/xn**2*Mgs0(3,3))
      Tgg4(3,4,2,3)=c1gs*(cf/2d0*Mgs0(1,2))
     .             +c2gs*(-1d0/2/xn*Mgs0(1,2)
     .                    -1/8d0*Mgs0(2,2)
     .                    -(xn**2-4d0)/8d0/xn**2*Mgs0(2,3))
     .             +c3gs*(-1/8d0*Mgs0(2,3)
     .                    -1/8d0*Mgs0(3,3))
      Tgg4(3,4,1,3)=c1gs*(-cf/2d0*Mgs0(1,2))
     .             +c2gs*(1d0/2/xn*Mgs0(1,2)
     .                    -1/8d0*Mgs0(2,2)
     .                    +(xn**2-4d0)/8d0/xn**2*Mgs0(2,3))
     .             +c3gs*(1/8d0*Mgs0(2,3)
     .                    -1/8d0*Mgs0(3,3))
      Tgg4(1,4,1,4)=Tgg4(2,3,2,3)
      Tgg4(2,4,2,4)=Tgg4(1,3,1,3)
      Tgg4(1,4,2,4)=Tgg4(2,3,1,3)
      Tgg4(2,4,1,4)=Tgg4(1,3,2,3)
      Tgg4(1,3,2,4)=Tgg4(1,3,1,3)
      Tgg4(1,3,1,4)=Tgg4(1,3,2,3)
      Tgg4(1,4,1,3)=Tgg4(2,3,1,3)
      Tgg4(1,4,2,3)=Tgg4(2,3,2,3)
      Tgg4(2,3,1,4)=Tgg4(2,3,2,3)
      Tgg4(2,3,2,4)=Tgg4(2,3,1,3)
      Tgg4(2,4,1,3)=Tgg4(1,3,1,3)
      Tgg4(2,4,2,3)=Tgg4(1,3,2,3)
      Tgg4(2,3,3,3)=Tgg(2,3)*cf
      Tgg4(2,3,4,4)=Tgg(2,3)*cf
      Tgg4(1,3,3,3)=Tgg(1,3)*cf
      Tgg4(1,3,4,4)=Tgg(1,3)*cf
      Tgg4(3,4,3,3)=Tgg(3,4)*cf
      Tgg4(3,4,4,4)=Tgg(3,4)*cf
      Tgg4(3,4,2,4)=Tgg4(3,4,1,3)
      Tgg4(3,4,1,4)=Tgg4(3,4,2,3)
      Tgg4(1,4,3,4)=Tgg4(2,3,3,4)
      Tgg4(2,4,3,4)=Tgg4(1,3,3,4)
      Tgg4(3,3,3,4)=Tgg(3,4)*cf
      Tgg4(4,4,3,4)=Tgg(3,4)*cf
      Tgg4(3,4,1,2)=-Tgg4(3,4,1,3)-ca*Tgg(3,4)-Tgg4(3,4,1,4)
      Tgg4(1,4,3,3)=Tgg(1,4)*cf
      Tgg4(1,4,4,4)=Tgg(1,4)*cf

      Tggeps(1,1)=ca*(Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs
     .                +Mgs0eps(3,3)*c3gs)
      Tggeps(2,2)=ca*(Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs
     .                +Mgs0eps(3,3)*c3gs)
      Tggeps(3,3)=cf*(Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs
     .                +Mgs0eps(3,3)*c3gs)
      Tggeps(4,4)=cf*(Mgs0eps(1,1)*c1gs+Mgs0eps(2,2)*c2gs
     .                +Mgs0eps(3,3)*c3gs)
      Tggeps(3,4)=-cf*Mgs0eps(1,1)*c1gs+
     .     1d0/2d0/xn*(Mgs0eps(2,2)*c2gs+Mgs0eps(3,3)*c3gs)
      Tggeps(1,4)=-0.5d0*Mgs0eps(1,2)*c1gs
     .          -(Mgs0eps(1,2)+ca/4d0*Mgs0eps(2,2)
     .          +(xn**2-4d0)/4d0/xn*Mgs0eps(2,3))*c2gs
     .          -ca/4d0*(Mgs0eps(2,3)+Mgs0eps(3,3))*c3gs
      Tggeps(1,3)=0.5d0*Mgs0eps(1,2)*c1gs
     .          +(Mgs0eps(1,2)-ca/4d0*Mgs0eps(2,2)
     .          +(xn**2-4d0)/4d0/xn*Mgs0eps(2,3))*c2gs
     .          +ca/4d0*(Mgs0eps(2,3)-Mgs0eps(3,3))*c3gs
      Tggeps(2,3)=Tggeps(1,4)
      Tggeps(2,4)=Tggeps(1,3)

      Tgg4eps(2,3,2,3)=c1gs*(0.5d0*Mgs0eps(1,1)
     .                    +xn/8d0*Mgs0eps(1,2)
     .                    +(xn**2-4d0)/8d0/xn*Mgs0eps(1,3))
     .             +c2gs*(xn/4d0*Mgs0eps(1,2)
     .                    +(xn**2+2d0)/8d0*Mgs0eps(2,2)
     .                    +(xn**2-4d0)/8d0*Mgs0eps(2,3))
     .             +c3gs*(xn/4d0*Mgs0eps(1,3)
     .                    +xn**2/8d0*Mgs0eps(2,3)
     .                    +(xn**2-2d0)/8d0*Mgs0eps(3,3))
      Tgg4eps(1,3,1,3)=c1gs*(0.5d0*Mgs0eps(1,1)
     .                    -xn/8d0*Mgs0eps(1,2)
     .                    +(xn**2-4d0)/8d0/xn*Mgs0eps(1,3))
     .             +c2gs*(-xn/4d0*Mgs0eps(1,2)
     .                    +(xn**2+2d0)/8d0*Mgs0eps(2,2)
     .                    -(xn**2-4d0)/8d0*Mgs0eps(2,3))
     .             +c3gs*(xn/4d0*Mgs0eps(1,3)
     .                    -xn**2/8d0*Mgs0eps(2,3)
     .                    +(xn**2-2d0)/8d0*Mgs0eps(3,3))
      Tgg4eps(2,3,1,3)=c1gs*(-0.5d0*Mgs0eps(1,1)
     .                    +xn/8d0*Mgs0eps(1,2)
     .                    -(xn**2-4d0)/8d0/xn*Mgs0eps(1,3))
     .             +c2gs*(-xn/4d0*Mgs0eps(1,2)
     .                    -1/4d0*Mgs0eps(2,2))
     .             +c3gs*(-xn/4d0*Mgs0eps(1,3)
     .                    +1/4d0*Mgs0eps(3,3))
      Tgg4eps(1,3,2,3)=c1gs*(-0.5d0*Mgs0eps(1,1)
     .                    -xn/8d0*Mgs0eps(1,2)
     .                    -(xn**2-4d0)/8d0/xn*Mgs0eps(1,3))
     .             +c2gs*(xn/4d0*Mgs0eps(1,2)
     .                    -1/4d0*Mgs0eps(2,2))
     .             +c3gs*(-xn/4d0*Mgs0eps(1,3)
     .                    +1/4d0*Mgs0eps(3,3))
      Tgg4eps(2,3,3,4)=c1gs*(-1d0/xn/4d0*Mgs0eps(1,2))
     .             +c2gs*(cf*Mgs0eps(1,2)
     .                    -1/8d0*Mgs0eps(2,2)
     .                    -(xn**2-4d0)/8d0/xn**2*Mgs0eps(2,3))
     .             +c3gs*(-1/8d0*Mgs0eps(2,3)
     .                    -1/8d0*Mgs0eps(3,3))
      Tgg4eps(1,3,3,4)=c1gs*(1d0/xn/4d0*Mgs0eps(1,2))
     .             +c2gs*(-cf*Mgs0eps(1,2)
     .                    -1/8d0*Mgs0eps(2,2)
     .                    +(xn**2-4d0)/8d0/xn**2*Mgs0eps(2,3))
     .             +c3gs*(1/8d0*Mgs0eps(2,3)
     .                    -1/8d0*Mgs0eps(3,3))
      Tgg4eps(3,4,3,4)=c1gs*(cf**2*Mgs0eps(1,1))
     .             +c2gs*(1/4d0/xn**2*Mgs0eps(2,2))
     .             +c3gs*(1/4d0/xn**2*Mgs0eps(3,3))
      Tgg4eps(3,4,2,3)=c1gs*(cf/2d0*Mgs0eps(1,2))
     .             +c2gs*(-1d0/2/xn*Mgs0eps(1,2)
     .                    -1/8d0*Mgs0eps(2,2)
     .                    -(xn**2-4d0)/8d0/xn**2*Mgs0eps(2,3))
     .             +c3gs*(-1/8d0*Mgs0eps(2,3)
     .                    -1/8d0*Mgs0eps(3,3))
      Tgg4eps(3,4,1,3)=c1gs*(-cf/2d0*Mgs0eps(1,2))
     .             +c2gs*(1d0/2/xn*Mgs0eps(1,2)
     .                    -1/8d0*Mgs0eps(2,2)
     .                    +(xn**2-4d0)/8d0/xn**2*Mgs0eps(2,3))
     .             +c3gs*(1/8d0*Mgs0eps(2,3)
     .                    -1/8d0*Mgs0eps(3,3))
      Tgg4eps(1,4,1,4)=Tgg4eps(2,3,2,3)
      Tgg4eps(2,4,2,4)=Tgg4eps(1,3,1,3)
      Tgg4eps(1,4,2,4)=Tgg4eps(2,3,1,3)
      Tgg4eps(2,4,1,4)=Tgg4eps(1,3,2,3)
      Tgg4eps(1,3,2,4)=Tgg4eps(1,3,1,3)
      Tgg4eps(1,3,1,4)=Tgg4eps(1,3,2,3)
      Tgg4eps(1,4,1,3)=Tgg4eps(2,3,1,3)
      Tgg4eps(1,4,2,3)=Tgg4eps(2,3,2,3)
      Tgg4eps(2,3,1,4)=Tgg4eps(2,3,2,3)
      Tgg4eps(2,3,2,4)=Tgg4eps(2,3,1,3)
      Tgg4eps(2,4,1,3)=Tgg4eps(1,3,1,3)
      Tgg4eps(2,4,2,3)=Tgg4eps(1,3,2,3)
      Tgg4eps(2,3,3,3)=Tggeps(2,3)*cf
      Tgg4eps(1,3,3,3)=Tggeps(1,3)*cf
      Tgg4eps(3,4,3,3)=Tggeps(3,4)*cf
      Tgg4eps(3,4,4,4)=Tggeps(3,4)*cf
      Tgg4eps(3,4,2,4)=Tgg4eps(3,4,1,3)
      Tgg4eps(3,4,1,4)=Tgg4eps(3,4,2,3)
      Tgg4eps(1,4,3,4)=Tgg4eps(2,3,3,4)
      Tgg4eps(2,4,3,4)=Tgg4eps(1,3,3,4)
      Tgg4eps(3,3,3,4)=Tggeps(3,4)*cf
      Tgg4eps(4,4,3,4)=Tggeps(3,4)*cf
      Tgg4eps(3,4,1,2)=-Tgg4eps(3,4,1,3)-ca*Tggeps(3,4)
     .                     -Tgg4eps(3,4,1,4)
      Tgg4eps(1,4,3,3)=Tggeps(1,4)*cf
      Tgg4eps(1,4,4,4)=Tggeps(1,4)*cf

      return
      end
