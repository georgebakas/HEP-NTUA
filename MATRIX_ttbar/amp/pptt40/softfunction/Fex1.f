c Result for <Fex1(0)> written in terms of the
c color factors T3T4 and CF

      function Fex1_T3T4(beta,cost)
      implicit none
      double precision Fex1_T3T4,beta,cost
      double precision v,B,L1_34,L1
      double complex aux,li2c
      external li2c,L1_34

      v  = (2.d0*beta)/(1.d0 + beta**2.d0)
      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
      L1 = L1_34(beta,cost)

      aux = L1/v + li2c(-B)

      Fex1_T3T4 = dreal(aux)

      return
      end

      function Fex1_CF(beta,cost)
      implicit none
      double precision Fex1_CF,beta,cost
      double precision B
      double complex aux,li2c
      external li2c

      B  = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)

      aux = dlog(1.d0 + B) + li2c(-B)

      Fex1_CF = dreal(aux)

      return
      end

