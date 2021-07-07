
      function aux_expression_1(B)
      implicit none
      double precision B,aux_expression_1
      double precision pi,z3
      double complex li2c,li3c,li4c,li22c
      double complex aux,im
      external li2c,li3c,li4c,li22c
      im = dcmplx(0.d0,1.d0);
      pi = 3.1415926535897932385d0
      z3 = 1.2020569031595942854d0


      aux =  -1.6078784942184556100263653206907164041 - 
     -  85.2672608708245054825598664345288368061195*im
     -   - (51*z3)/8. - 
     -  (7*Pi**2*dlog(1.*(1 + 1/B))**2)/24. + 
     -  dlog(1.*(1 + 1/B))**3/4. - 
     -  (7*dlog(1.*(1 + 1/B))**4)/96. - 
     -  (3*dlog(1.*(1 + 1/B))**2*(im*Pi + dlog(1./B)))/
     -   4. - dlog(1.*B)/8. + (Pi**2*dlog(1.*B))/4. + 
     -  (z3*dlog(1.*B))/2. + 
     -  (dlog(1.*(1 + 1/B))**2*dlog(1.*B))/4. + 
     -  (11*dlog(1.*B)**2)/4. - 
     -  (im*Pi*dlog(1.*B)**2)/2. + 
     -  (25*Pi**2*dlog(1.*B)**2)/48. + 
     -  dlog(1.*B)**3/4. + (25*dlog(1.*B)**4)/96. - 
     -  (2*dlog((1.*B)/(1 + B)))/(3.*Sqrt(1 + B)) + 
     -  (2*Pi**2*dlog((1.*B)/(1 + B)))/3. + 
     -  (4*Pi**2*dlog(0.5/Sqrt(1 + B)))/3. - 
     -  (73*dlog(1.*(1 + B)))/36. + 
     -  (23*Pi**2*dlog(1.*(1 + B)))/12. - 
     -  (z3*dlog(1.*(1 + B)))/2. - 
     -  (Pi**2*dlog(1.*(1 + 1/B))*dlog(1.*(1 + B)))/6. + 
     -  (dlog(1.*(1 + 1/B))**3*dlog(1.*(1 + B)))/12. - 
     -  dlog(1.*B)*dlog(1.*(1 + B)) + 
     -  im*Pi*dlog(1.*B)*dlog(1.*(1 + B)) - 
     -  (5*Pi**2*dlog(1.*B)*dlog(1.*(1 + B)))/4. - 
     -  (dlog(1.*B)**2*dlog(1.*(1 + B)))/2. - 
     -  (5*dlog(1.*B)**3*dlog(1.*(1 + B)))/12. - 
     -  (3*dlog(1.*(1 + B))**2)/32. - 
     -  (im*Pi*dlog(1.*(1 + B))**2)/8. + 
     -  (11*Pi**2*dlog(1.*(1 + B))**2)/16. + 
     -  (dlog(2.d0)*dlog(1.*(1 + B))**2)/2. + 
     -  (dlog(1.*B)*dlog(1.*(1 + B))**2)/8. - 
     -  (dlog(1.*B)**2*dlog(1.*(1 + B))**2)/8. - 
     -  (dlog((1.*B)/(1 + B))*dlog(1.*(1 + B))**2)/4. - 
     -  (dlog(0.5/Sqrt(1 + B))*dlog(1.*(1 + B))**2)/2. - 
     -  (29*dlog(1.*(1 + B))**3)/72. + 
     -  (11*dlog(1.*B)*dlog(1.*(1 + B))**3)/24. - 
     -  (17*dlog(1.*(1 + B))**4)/96. - 
     -  (4*dlog(1.*(1 - 1/Sqrt(1 + B))))/3. + 
     -  (2*dlog(1.*(1 - 1/Sqrt(1 + B))))/
     -   (3.*Sqrt(1 + B)) - 
     -  (4*Pi**2*dlog(1.*(1 - 1/Sqrt(1 + B))))/3. - 
     -  (dlog(1.*(1 + B))*dlog(1.*(1 - 1/Sqrt(1 + B))))/
     -   6. - (dlog(1.*(1 + B))*
     -     dlog(1.*(1 - 1/Sqrt(1 + B))))/(4.*Sqrt(1 + B))
     -    + (dlog(1.*(1 + B))**2*
     -     dlog(1.*(1 - 1/Sqrt(1 + B))))/2. - 
     -  (4*Pi**2*dlog(0.5*(1 + 1/Sqrt(1 + B))))/3. - 
     -  2*dlog(2.d0)*dlog(0.5*(1 + 1/Sqrt(1 + B)))**2 + 
     -  (2*dlog(1.*(1 + 1/Sqrt(1 + B))))/
     -   (3.*Sqrt(1 + B)) + 
     -  (4*Pi**2*dlog(1.*(1 + 1/Sqrt(1 + B))))/3. + 
     -  (dlog(1.*(1 + B))*dlog(1.*(1 + 1/Sqrt(1 + B))))/
     -   (4.*Sqrt(1 + B)) + 
     -  2*dlog(1.*B)*dlog(1.*(1 + B))*
     -   dlog(1.*(1 + 1/Sqrt(1 + B))) + 
     -  2*dlog(0.5*(1 + 1/Sqrt(1 + B)))*
     -   dlog(1.*(1 + 1/Sqrt(1 + B)))**2 - 
     -  (2*dlog(1.*(1 + 1/Sqrt(1 + B)))**3)/3. - 
     -  (Pi**2*dlog(1.*(1 - 1/(1 + 2*B))))/3. + 
     -  (Pi**2*dlog(1.*(1 + 1/(1 + 2*B))))/3. - 
     -  (Pi**2*dlog(2./(-1 + dsqrt(1 + B))))/4. + 
     -  (35*dlog(1.*(-1 + dsqrt(1 + B))))/24. - 
     -  (35*Pi**2*dlog(1.*(-1 + dsqrt(1 + B))))/12. + 
     -  (dlog(2.d0)*dlog(1.*(-1 + dsqrt(1 + B))))/3. - 
     -  dlog(512.d0)*dlog(1.*(-1 + dsqrt(1 + B))) - 
     -  dlog(1.*B)*dlog(1.*(-1 + dsqrt(1 + B))) + 
     -  (13*dlog((1.*B)/(1 + B))*
     -     dlog(1.*(-1 + dsqrt(1 + B))))/6. + 
     -  (13*dlog(1.*(1 + B))*
     -     dlog(1.*(-1 + dsqrt(1 + B))))/24. - 
     -  dlog(4.d0)*dlog(1.*(1 + B))*
     -   dlog(1.*(-1 + dsqrt(1 + B))) + 
     -  2*dlog(1.*(1 + B))**2*
     -   dlog(1.*(-1 + dsqrt(1 + B))) + 
     -  (13*dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -     dlog(1.*(-1 + dsqrt(1 + B))))/6. - 
     -  (13*dlog(1.*(1 + 1/Sqrt(1 + B)))*
     -     dlog(1.*(-1 + dsqrt(1 + B))))/6. - 
     -  (16*dlog(1.*(-1 + dsqrt(1 + B)))**2)/3. + 
     -  (17*Pi**2*(im*Pi + 
     -       dlog(1.*(-1 + dsqrt(1 + B)))))/24. - 
     -  (3*dlog(1.*(1 + B))*
     -     (im*Pi + dlog(1.*(-1 + dsqrt(1 + B)))))/2. + 
     -  (5*im*Pi*dlog(2.d0)*
     -     (im*Pi + dlog(1.*dsqrt(1 + B))))/2. + 
     -  (dlog(32.d0)*dlog(1.*(1 + B))*
     -     (im*Pi + dlog(1.*dsqrt(1 + B))))/4. - 
     -  (5*dlog(2.d0)*(im*Pi + dlog(1.*dsqrt(1 + B)))**2)/
     -   2. - (11*Pi**2*
     -     dlog((1.*dsqrt(1 + B))/(-1 + dsqrt(1 + B))))/
     -   12. - (2*dlog((1.*dsqrt(1 + B))/
     -        (-1 + dsqrt(1 + B)))**3)/3. - 
     -  (13*Pi**2*dlog((2.*dsqrt(1 + B))/
     -       (-1 + dsqrt(1 + B))))/12. - 
     -  2*dlog(2.d0)*dlog(1.*(1 + B))*
     -   dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B))) + 
     -  2*dlog((1.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))**2*
     -   dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B))) - 
     -  2*dlog(0.5/(-1 + dsqrt(1 + B)))*
     -   dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))**2
     -   - 2*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))**2
     -   - (4*dlog((2.*dsqrt(1 + B))/
     -        (-1 + dsqrt(1 + B)))**3)/3. + 
     -  2*dlog(2.d0)*dlog(1.*(1 + B))*
     -   dlog(0.5*(1 + dsqrt(1 + B))) + 
     -  4*dlog(2.d0)*dlog(0.5*(1 + 1/Sqrt(1 + B)))*
     -   dlog(0.5*(1 + dsqrt(1 + B))) - 
     -  (19*dlog(1.*(-1 + dsqrt(1 + B)))*
     -     dlog(0.5*(1 + dsqrt(1 + B))))/3. - 
     -  2*dlog(1.*(1 + B))*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   dlog(0.5*(1 + dsqrt(1 + B))) - 
     -  (7*dlog(0.5*(1 + dsqrt(1 + B)))**2)/6. - 
     -  2*dlog(2.d0)*dlog(0.5*(1 + dsqrt(1 + B)))**2 + 
     -  (4*dlog(0.5*(1 + dsqrt(1 + B)))**3)/3. + 
     -  dlog(1.*(1 + dsqrt(1 + B)))/8. - 
     -  (17*Pi**2*dlog(1.*(1 + dsqrt(1 + B))))/24. - 
     -  dlog(1.*B)*dlog(1.*(1 + dsqrt(1 + B))) - 
     -  (13*dlog((1.*B)/(1 + B))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/6. + 
     -  (19*dlog(1.*(1 + B))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/6. + 
     -  (11*im*Pi*dlog(1.*(1 + B))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/4. - 
     -  (3*dlog(1.*B)*dlog(1.*(1 + B))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/4. + 
     -  (5*dlog(1.*(1 + B))**2*
     -     dlog(1.*(1 + dsqrt(1 + B))))/2. + 
     -  (13*dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/2. - 
     -  dlog(1.*(1 + B))*dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -   dlog(1.*(1 + dsqrt(1 + B))) + 
     -  (13*dlog(1.*(1 + 1/Sqrt(1 + B)))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/6. + 
     -  dlog(1.*(1 + B))*dlog(1.*(1 + 1/Sqrt(1 + B)))*
     -   dlog(1.*(1 + dsqrt(1 + B))) + 
     -  (7*dlog(0.5*(-1 + dsqrt(1 + B)))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/3. - 
     -  (7*dlog(1.*(1 + B))*dlog(1.*(-1 + dsqrt(1 + B)))*
     -     dlog(1.*(1 + dsqrt(1 + B))))/2. + 
     -  5*im*Pi*(im*Pi + dlog(1.*(-1 + dsqrt(1 + B))))*
     -   dlog(1.*(1 + dsqrt(1 + B))) - 
     -  dlog(1.*(1 + dsqrt(1 + B)))**2 - 
     -  (13*dlog(1.*(1 + B))*
     -     dlog(1.*(1 + dsqrt(1 + B)))**2)/4. + 
     -  2*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   dlog(1.*(1 + dsqrt(1 + B)))**2 - 
     -  (5*(im*Pi + dlog(1.*dsqrt(1 + B)))*
     -     dlog(1.*(1 + dsqrt(1 + B)))**2)/2. + 
     -  (2*dlog(1.*(1 + dsqrt(1 + B)))**3)/3. - 
     -  (5*dlog(1.*(1 + B))*dlog(1.*(1 + dsqrt(1 + B)))*
     -     (im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -         (1 - dsqrt(1 + B)))))/2. + 
     -  5*(im*Pi + dlog(1.*dsqrt(1 + B)))*
     -   dlog(1.*(1 + dsqrt(1 + B)))*
     -   (im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -       (1 - dsqrt(1 + B)))) + 
     -  (5*im*Pi*(im*Pi + 
     -        dlog((-1.*(1 + dsqrt(1 + B)))/
     -          (1 - dsqrt(1 + B))))**2)/2. - 
     -  (8*Pi**2*dlog((1.*(1 + dsqrt(1 + B)))/
     -       (-1 + dsqrt(1 + B))))/3. - 
     -  (13*dlog(1.*B)*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -       (-1 + dsqrt(1 + B))))/2. + 
     -  dlog(1.*B)*dlog(1.*(1 + B))*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) - 
     -  (dlog(1.*(1 + B))**2*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -       (-1 + dsqrt(1 + B))))/2. + 
     -  (13*dlog(1.*(-1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -       (-1 + dsqrt(1 + B))))/2. - 
     -  dlog(1.*(1 + B))*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) - 
     -  4*dlog(1.*(-1 + dsqrt(1 + B)))**2*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) + 
     -  (13*dlog(1.*(1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -       (-1 + dsqrt(1 + B))))/2. - 
     -  dlog(1.*(1 + B))*dlog(1.*(1 + dsqrt(1 + B)))*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) + 
     -  8*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   dlog(1.*(1 + dsqrt(1 + B)))*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) - 
     -  6*dlog(1.*(1 + dsqrt(1 + B)))**2*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B))) - 
     -  dlog(1.*(1 + B))*
     -   dlog((1.*(1 + dsqrt(1 + B)))/
     -      (-1 + dsqrt(1 + B)))**2 + 
     -  (dlog(2./(-1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -        (-1 + dsqrt(1 + B)))**2)/4. + 
     -  (9*dlog(1.*(-1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -        (-1 + dsqrt(1 + B)))**2)/4. + 
     -  (9*dlog((1.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -        (-1 + dsqrt(1 + B)))**2)/4. - 
     -  (dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))*
     -     dlog((1.*(1 + dsqrt(1 + B)))/
     -        (-1 + dsqrt(1 + B)))**2)/4. + 
     -  (8*dlog((1.*(1 + dsqrt(1 + B)))/
     -        (-1 + dsqrt(1 + B)))**3)/3. - 
     -  (3*dlog(1.*(1 + 1/B))*li2c(1 + 1/B))/2. + 
     -  (3*li2c(-(1/B)))/2. + (5*Pi**2*li2c(-(1/B)))/24. + 
     -  (dlog(1.*(1 + 1/B))**2*li2c(-(1/B)))/8. + 
     -  (dlog(1.*B)*li2c(-(1/B)))/2. + 
     -  (dlog(1.*B)**2*li2c(-(1/B)))/8. + 
     -  (dlog(1.*(1 + B))*li2c(-(1/B)))/2. - 
     -  (dlog(1.*(1 + B))**2*li2c(-(1/B)))/2. + 
     -  (5*li2c(-B))/9. - dlog(1.*B)*li2c(-B) + 
     -  (7*dlog(1.*(1 + B))*li2c(-B))/12. + 
     -  ((im*Pi + dlog(1.*(-1 + dsqrt(1 + B))))*li2c(-B))/
     -   4. - (9*dlog(1.*(1 + dsqrt(1 + B)))*li2c(-B))/
     -   4. - li2c(1/(1 + B)) - 
     -  li2c(1/(1 + B))/(4.*Sqrt(1 + B)) + 
     -  (3*Pi**2*li2c(1/(1 + B)))/8. + 
     -  (dlog(1.*(1 + 1/B))*li2c(1/(1 + B)))/2. + 
     -  (3*dlog(1.*B)*li2c(1/(1 + B)))/2. - 
     -  dlog(1.*(1 + B))*li2c(1/(1 + B)) - 
     -  (3*dlog(1.*(1 + B))**2*li2c(1/(1 + B)))/8. + 
     -  dlog(1.*(-1 + dsqrt(1 + B)))*li2c(1/(1 + B)) - 
     -  dlog(1.*(1 + dsqrt(1 + B)))*li2c(1/(1 + B)) - 
     -  2*dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B)))*li2c(1/(1 + B)) + 
     -  (Pi**2*li2c(B/(1 + B)))/12. - 
     -  (dlog(1.*(1 + 1/B))**2*li2c(B/(1 + B)))/8. + 
     -  (li2c(-(1/B))*li2c(B/(1 + B)))/4. - 
     -  li2c(B/(1 + B))**2/8. + 
     -  (91*li2c(1/Sqrt(1 + B)))/12. + 
     -  li2c(1/Sqrt(1 + B))/Sqrt(1 + B) - 
     -  2*dlog((1.*B)/(1 + B))*li2c(1/Sqrt(1 + B)) - 
     -  4*dlog(0.5/Sqrt(1 + B))*li2c(1/Sqrt(1 + B)) + 
     -  2*dlog(1.*(1 + B))*li2c(1/Sqrt(1 + B)) + 
     -  4*dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -   li2c(1/Sqrt(1 + B)) - 
     -  6*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   li2c(1/Sqrt(1 + B)) + 
     -  6*dlog(1.*(1 + dsqrt(1 + B)))*
     -   li2c(1/Sqrt(1 + B)) + 
     -  4*dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B)))*li2c(1/Sqrt(1 + B)) + 
     -  4*dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))*
     -   li2c(1 - 1/Sqrt(1 + B)) - 
     -  4*dlog((2.*dsqrt(1 + B))/(-1 + dsqrt(1 + B)))*
     -   li2c(0.5 - 1/(2.*Sqrt(1 + B))) + 
     -  4*dlog(0.5*(1 + 1/Sqrt(1 + B)))*
     -   li2c((1 + 1/Sqrt(1 + B))/2.) + 
     -  (7*li2c((1 - dsqrt(1 + B))/2.))/3. - 
     -  4*dlog(0.5*(-1 + dsqrt(1 + B)))*
     -   li2c((1 - dsqrt(1 + B))/2.) - 
     -  (103*li2c(1 - dsqrt(1 + B)))/12. - 
     -  (im*Pi*li2c(1 - dsqrt(1 + B)))/2. - 
     -  (5*dlog(1.*B)*li2c(1 - dsqrt(1 + B)))/2. + 
     -  (9*dlog(1.*(1 + B))*li2c(1 - dsqrt(1 + B)))/8. + 
     -  (dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -     li2c(1 - dsqrt(1 + B)))/4. + 
     -  2*dlog(1.*(1 + 1/Sqrt(1 + B)))*
     -   li2c(1 - dsqrt(1 + B)) - 
     -  4*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   li2c(1 - dsqrt(1 + B)) - 
     -  (9*dlog((1.*(-1 + dsqrt(1 + B)))/
     -       (1 + dsqrt(1 + B)))*li2c(1 - dsqrt(1 + B)))/
     -   4. + (11*dlog(1.*(1 + dsqrt(1 + B)))*
     -     li2c(1 - dsqrt(1 + B)))/4. - 
     -  5*(im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -       (1 - dsqrt(1 + B))))*li2c(1 - dsqrt(1 + B))
     -   - 8*dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B)))*li2c(1 - dsqrt(1 + B)) - 
     -  2*li2c(-dsqrt(1 + B)) + 
     -  (11*im*Pi*li2c(-dsqrt(1 + B)))/2. + 
     -  (5*dlog(1.*B)*li2c(-dsqrt(1 + B)))/2. - 
     -  2*dlog((1.*B)/(1 + B))*li2c(-dsqrt(1 + B)) - 
     -  (dlog(1.*(1 + B))*li2c(-dsqrt(1 + B)))/2. + 
     -  2*dlog(1.*(1 - 1/Sqrt(1 + B)))*
     -   li2c(-dsqrt(1 + B)) + 
     -  2*dlog(1.*(1 + 1/Sqrt(1 + B)))*
     -   li2c(-dsqrt(1 + B)) + 
     -  3*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   li2c(-dsqrt(1 + B)) - 
     -  4*dlog(2.*dsqrt(1 + B))*li2c(-dsqrt(1 + B)) - 
     -  7*dlog(1.*(1 + dsqrt(1 + B)))*
     -   li2c(-dsqrt(1 + B)) - 
     -  4*dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B)))*li2c(-dsqrt(1 + B)) - 
     -  3*li2c(dsqrt(1 + B)) + 
     -  (3*dlog(1.*(1 + B))*li2c(dsqrt(1 + B)))/2. + 
     -  4*dlog(1.*(-1 + dsqrt(1 + B)))*
     -   li2c(1/(1 + dsqrt(1 + B))) + 
     -  4*dlog(1.*(1 + dsqrt(1 + B)))*
     -   li2c(1/(1 + dsqrt(1 + B))) - 
     -  4*dlog((1.*(1 + dsqrt(1 + B)))/
     -     (-1 + dsqrt(1 + B)))*li2c(1/(1 + dsqrt(1 + B)))
     -    - (7*li2c(2/(1 + dsqrt(1 + B))))/3. + 
     -  4*dlog(0.5*(1 + dsqrt(1 + B)))*
     -   li2c(2/(1 + dsqrt(1 + B))) + 
     -  4*dlog(0.5*(1 + 1/Sqrt(1 + B)))*
     -   li2c(dsqrt(1 + B)/(1 + dsqrt(1 + B))) + 
     -  5*(im*Pi + dlog(1.*(-1 + dsqrt(1 + B))))*
     -   li2c(1 + dsqrt(1 + B)) + 
     -  5*(im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -       (1 - dsqrt(1 + B))))*li2c(1 + dsqrt(1 + B))
     -   - 5*(im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -       (1 - dsqrt(1 + B))))*
     -   li2c((1 + dsqrt(1 + B))/(1 - dsqrt(1 + B))) + 
     -  5*(im*Pi + dlog((-1.*(1 + dsqrt(1 + B)))/
     -       (1 - dsqrt(1 + B))))*
     -   li2c((1 + dsqrt(1 + B))/(-1 + dsqrt(1 + B))) + 
     -  li22c(-(1/B))/2. + li3c(1 + 1/B)/2. - 
     -  li3c(-(1/B))/2. - (dlog(1.*B)*li3c(-(1/B)))/2. + 
     -  (dlog(1.*(1 + B))*li3c(-(1/B)))/2. + li3c(-B) - 
     -  li3c(1/(1 + B)) - 
     -  (dlog(1.*(1 + B))*li3c(1/(1 + B)))/4. + 
     -  (11*li3c(B/(1 + B)))/12. - 
     -  (dlog(1.*(1 + B))*li3c(B/(1 + B)))/2. + 
     -  8*li3c(1/Sqrt(1 + B)) + 
     -  4*li3c(1 - 1/Sqrt(1 + B)) - 
     -  4*li3c(0.5 - 1/(2.*Sqrt(1 + B))) - 
     -  4*li3c((1 + 1/Sqrt(1 + B))/2.) + 
     -  4*li3c((1 - dsqrt(1 + B))/2.) + 
     -  3*li3c(1 - dsqrt(1 + B)) + 5*li3c(-dsqrt(1 + B)) - 
     -  3*li3c(dsqrt(1 + B)) + 
     -  8*li3c(1/(1 + dsqrt(1 + B))) + 
     -  4*li3c(2/(1 + dsqrt(1 + B))) + 
     -  8*li3c((-1 + dsqrt(1 + B))/(1 + dsqrt(1 + B))) + 
     -  4*li3c(dsqrt(1 + B)/(1 + dsqrt(1 + B))) - 
     -  5*li3c(1 + dsqrt(1 + B)) - 
     -  3*li3c((1 + dsqrt(1 + B))/(1 - dsqrt(1 + B))) - 
     -  5*li3c((1 + dsqrt(1 + B))/(-1 + dsqrt(1 + B))) - 
     -  (3*li4c(1/(1 + B)))/4. - li4c(B/(1 + B))


      aux_expression_1 = dreal(aux)

      return
      end

      function Fex2_TITj(B,A)
      implicit none
      double precision B,A,Fex2_TITj
      double precision pi,z3
      double precision aux_expression_1
      double complex li2c,li3c
      external li2c,li3c

      pi = 3.1415926535897932385d0
      z3 = 1.2020569031595942854d0

      Fex2_TITj = 3.d0*aux_expression_1(B)
     &      + ((-10.d0*(56.d0 + 3.d0*Pi**2.d0)
     &      + 3.d0*(808.d0 + 33.d0*Pi**2.d0 - 756.d0*z3))*dlog(A))/432.
     &      + 5.d0*(dreal((56.d0 + 3.d0*Pi**2.d0 
     &      - 36.d0*li3c(B/(1.d0 + B)) + 60.d0*log(1.d0 + B) + 
     &      18.d0*log(1.d0 + B)**2.d0 + 6.d0*log(1.d0 + B)**3.d0 + 
     &      12.d0*li2c(-B)*(2.d0 + 3.d0*log(1.d0 + B)))/216.))

      return
      end

      function Fex1_1_TITj(B,A)
      implicit none
      double precision B,A,Fex1_1_TITj
      double precision pi
      double complex li2c,li3c,aux
      external li2c,li3c

      pi = 3.1415926535897932385d0

      aux = (-Pi**2.d0 + Pi**2.d0*dlog(A)
     &      + 12.d0*li2c(-B) - 12.d0*li3c(-B))/24.d0
      Fex1_1_TITj = dreal(aux)

      return
      end
