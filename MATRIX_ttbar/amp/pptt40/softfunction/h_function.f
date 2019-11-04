
      function h1_TT(i,j,beta,cost)
      implicit none
      double precision h1_TT,beta,cost
      integer i,j
      double precision Fex1_T3T4,Fex1_CF
      external Fex1_T3T4,Fex1_CF

      h1_TT = 0.d0

      if((i.eq.3).and.(j.eq.4))then
        h1_TT = Fex1_T3T4(beta,cost)
      endif

      if((i.eq.3).and.(j.eq.3))then
        h1_TT = Fex1_CF(beta,cost)
      endif

      return
      end


      function h2_TT(i,j,beta,cost)
      implicit none
      double precision h2_TT,beta,cost
      integer i,j
      double precision Fex2_and_11_T3T4,Fex2_TITj,Fex1_1_TITj
      double precision A,B
      double precision pi
      external Fex2_and_11_T3T4,Fex2_TITj,Fex1_1_TITj

      pi = 3.1415926535897932385d0

      h2_TT = 0.d0

      if((i.eq.3).and.(j.eq.4))then
        h2_TT = Fex2_and_11_T3T4(beta,cost)
      endif

      if( ((i.eq.3).and.(j.eq.1)).or.
     &    ((i.eq.4).and.(j.eq.2)) )then
        B = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
        A = (1.d0 - beta*cost)**2.d0/(1.d0 - beta**2.d0)
        h2_TT = Fex2_TITj(B,A)-23.d0/6.d0*Fex1_1_TITj(B,A)
      endif

      if( ((i.eq.3).and.(j.eq.2)).or.
     &    ((i.eq.4).and.(j.eq.1)) )then
        B = (beta**2.d0*(1.d0 - cost**2.d0))/(1.d0 - beta**2.d0)
        A = (1.d0 + beta*cost)**2.d0/(1.d0 - beta**2.d0) 
        h2_TT = Fex2_TITj(B,A)-23.d0/6.d0*Fex1_1_TITj(B,A)
      endif

      return
      end


      function h2_TTTT(i,j,k,l,beta,cost)
      implicit none
      double precision h2_TTTT,beta,cost
      integer i,j,k,l
      double precision Fex1sq_T3T4sq,Fex1sq_CFsq,Fex1sq_T3T4CF
      double precision sqFex1_T3T4sq,sqFex1_CFsq,sqFex1_T3T4CF 
      external Fex1sq_T3T4sq,Fex1sq_CFsq,Fex1sq_T3T4CF 
      external sqFex1_T3T4sq,sqFex1_CFsq,sqFex1_T3T4CF 

      h2_TTTT = 0.d0

      if((i.eq.3).and.(j.eq.4).and.(k.eq.3).and.(l.eq.4))then
        h2_TTTT = Fex1sq_T3T4sq(beta,cost)
     &            -1.d0/2.d0*sqFex1_T3T4sq(beta,cost)

      endif

      if((i.eq.3).and.(j.eq.4).and.(k.eq.3).and.(l.eq.3))then
        h2_TTTT = Fex1sq_T3T4CF(beta,cost)
     &            -1.d0/2.d0*sqFex1_T3T4CF(beta,cost)

      endif

      if((i.eq.3).and.(j.eq.3).and.(k.eq.3).and.(l.eq.3))then
        h2_TTTT = Fex1sq_CFsq(beta,cost)
     &            -1.d0/2.d0*sqFex1_CFsq(beta,cost)

      endif

      return
      end
