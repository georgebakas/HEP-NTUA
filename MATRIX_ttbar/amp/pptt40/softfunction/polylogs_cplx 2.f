C Required polylogs are Li2, Li3, Li4, Li22
C All the arguments are real, but I added a small
C imaginary part to ensure the correct analytic continuation
C (i.e. consistent to what I was having in Mathematica)
      
      function li2c(x)
      implicit none
      double precision x
      double complex li2c,z,HPL2
      z=dcmplx(x,-1.d-40)
      li2c=HPL2(0,1,z)
      return
      end

      function li3c(x)
      implicit none
      double precision x
      double complex li3c,z,HPL3
      z=dcmplx(x,-1.d-40)
      li3c=HPL3(0,0,1,z)
      return
      end

      function li4c(x)
      implicit none
      double precision x
      double complex li4c,z,HPL4
      z=dcmplx(x,-1.d-40)
      li4c=HPL4(0,0,0,1,z)
      return
      end


      function li22c(x)
      implicit none
      double precision x
      double complex li22c,z,HPL4
      z=dcmplx(x,-1.d-40)
      li22c=HPL4(0,0,1,1,z)
      return
      end
