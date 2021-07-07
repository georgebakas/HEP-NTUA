C Required polylogs are Li2, Li3, Li4, Li22
C All are expected to be real
      
      function li2(x)
      implicit none
      double precision li2,x
      double complex z,HPL2
      z=dcmplx(x)
      li2=dreal(HPL2(0,1,z))
      return
      end

      function li3(x)
      implicit none
      double precision li3,x
      double complex z,HPL3
      z=dcmplx(x)
      li3=dreal(HPL3(0,0,1,z))
      return
      end

      function li4(x)
      implicit none
      double precision li4,x
      double complex z,HPL4
      z=dcmplx(x)
      li4=dreal(HPL4(0,0,0,1,z))
      return
      end


      function li22(x)
      implicit none
      double precision li22,x
      double complex z,HPL4
      z=dcmplx(x)
      li22=dreal(HPL4(0,0,1,1,z))
      return
      end
