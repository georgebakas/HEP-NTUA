c Here are the routines to compute <Fex,1(0)>^2
c The results are written in terms of the color
c structures T3T4^2, T3T4*CF and CF^2

      function sqFex1_CFsq(beta,cost)
      implicit none
      double precision sqFex1_CFsq,beta,cost
      double precision B
      double complex aux,li2c
      external li2c

      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)

      aux = dlog(1.d0 + B)**2.d0 + 2.d0*dlog(1.d0 + B)*li2c(-B)
     &      + li2c(-B)**2.d0

      sqFex1_CFsq = dreal(aux)

      return
      end


      function sqFex1_T3T4CF(beta,cost)
      implicit none
      double precision sqFex1_T3T4CF,beta,cost
      double precision B,v,L1_34,L1
      double complex aux,li2c
      external li2c,L1_34

      v  = (2.d0*beta)/(1.d0 + beta**2.d0)
      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
      L1 = L1_34(beta,cost)

      aux = (2.d0*L1*dlog(1.d0 + B))/v + (2.d0*L1*li2c(-B))/v + 
     &      2.d0*dlog(1.d0 + B)*li2c(-B) + 2.d0*li2c(-B)**2.d0

      sqFex1_T3T4CF = dreal(aux)

      return
      end


      function sqFex1_T3T4sq(beta,cost)
      implicit none
      double precision sqFex1_T3T4sq,beta,cost
      double precision B,v,L1_34,L1
      double complex aux,li2c
      external li2c,L1_34

      v  = (2.d0*beta)/(1.d0 + beta**2.d0)
      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
      L1 = L1_34(beta,cost)

      aux = L1**2.d0/v**2.d0 + (2.d0*L1*li2c(-B))/v + li2c(-B)**2.d0

      sqFex1_T3T4sq = dreal(aux)

      return
      end


      function L1_34(beta,cost)
      implicit none
      double precision L1_34,beta,cost
      double complex li2c,aux
      double precision v,B,c,vT,cT
      external li2c

      v  = (2.d0*beta)/(1.d0 + beta**2.d0)
      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
      c  = dsqrt((1.d0 - v)/(1.d0 + v))
      vT = (2.d0*beta*dabs(cost))/(1.d0 + beta**2.d0*cost**2.d0)
      cT = dsqrt((1.d0 - vT)/(1.d0 + vT))

      aux =  dlog(cT)**2.d0/2. 
     &  + (dlog(1.d0+B)*dlog((1.d0 + v)/(1.d0 - v)))/2. - 
     &  dlog((1.d0 + v)/(1.d0 - v))**2.d0/8. + li2c(1.d0 - c/cT) + 
     &  li2c(1.d0 - c*cT) - li2c((2.d0*v)/(1.d0 + v))

      L1_34 = dreal(aux)

      return
      end
