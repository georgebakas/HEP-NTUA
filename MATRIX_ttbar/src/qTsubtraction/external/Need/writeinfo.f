      subroutine writeinfo(unitno,xsec,xsec_err)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'maxwt.f'
      include 'masses.f'
      include 'facscale.f'
      include 'scale.f'
      include 'zerowidth.f'
      include 'flags.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'pdfiset.f'
      include 'dynamicscale.f'
      integer unitno
      include 'lhapdf.f'
      double precision xsec,xsec_err
      
      character*4 part
      character*30 runstring
      character *50 prefix
      logical creatent,dswhisto,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,rseed
      integer order
      double precision sqrts,Mwmin,Mwmax
      double precision Rcut
      double precision leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut
      integer lbjscheme
      logical jetsopphem
 
      common/outputflags/creatent,dswhisto      

      common/nnlo/order
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      integer nset
      common/prefix/nset,prefix
      common/mwminmax/Mwmin,Mwmax
         
    
      
      common/Rcut/Rcut
      common/makecuts/makecuts
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut,
     . lbjscheme,jetsopphem

      common/rseed/rseed

      write(unitno,*) '( Cross-section is: ',xsec,'+/-',xsec_err,')'
      write(unitno,*)
      write(unitno,*) '( Run corresponds to this input file)'
      write(unitno,*)
      write(unitno,*) '(sqrts= ',sqrts
      write(unitno,*) '(ih1= ',ih1,'  ih2= ',ih2
      write(unitno,*) '(nproc: ',nproc
      write(unitno,*) '(dynamicscale=',dynamicscale
      if(dynamicscale.eqv..false.) then
       write(unitno,*) '(muf= ',facscale
       write(unitno,*) '(mur= ',scale
      endif
      write(unitno,*) '(order= ',order
      write(unitno,*) '(part= ',part
      write(unitno,*) '(zerowidth= ',zerowidth
      write(unitno,*) '(Mwmin= ',Mwmin,'  Mwmax= ',Mwmax
      write(unitno,*) '(itmx1= ',itmx1
      write(unitno,*) '(ncall1= ',ncall1
      write(unitno,*) '(itmx2= ',itmx2
      write(unitno,*) '(ncall2= ',ncall2
      write(unitno,*) '(rnd seed= ',rseed
      write(unitno,*) '(iset=',iset,' nset=',nset
      write(unitno,*) '(PDFname=',PDFname,' PDFmember=',PDFmember
      write(unitno,*) '(runstring=',runstring
      write(unitno,*)
      write(unitno,*) '( td -b filename.top '
      write(unitno,*) 'SET DEVICE POSTSCRIPT SIDEWAYS'
      write(unitno,*) 'SET SIZE SIDEWAYS'

      return
      
      end
      
      
