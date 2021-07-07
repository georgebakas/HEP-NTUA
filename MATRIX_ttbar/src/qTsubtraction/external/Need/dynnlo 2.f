      program dynnlo

      implicit none
      include 'constants.f'
      include 'gridinfo.f'
      integer itmx1,ncall1,itmx2,ncall2
      double precision integ,integ_err
      double precision p(mxpart,4),wt
      common/iterat/itmx1,ncall1,itmx2,ncall2
    

CC    Initialization

      call dyinit


CC    Warm up

      if(readin.eqv. .false.)
     & call integrate(0,itmx1,ncall1,.false.,integ,integ_err)

      
CC    Main run    

      call integrate(1,itmx2,ncall2,.true.,integ,integ_err)
    


CC    Final processing and print-out

      call hexit(integ,integ_err)
      
      stop
      end
       
