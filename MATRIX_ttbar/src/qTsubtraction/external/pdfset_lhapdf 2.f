*****************
* LHAPDF version*
*****************
      subroutine pdfset
ch      implicit none
ch      include 'masses.f'
ch      include 'lhapdf.f'
ch      include 'PDFerrors.f'
ch      include 'pdlabel.f'
ch      double precision amz,alphasPDF
ch      logical validPDF
ch      character*30 oldPDFname
ch      integer i
ch      logical lhapdfs
ch      common/lhapdfs/lhapdfs

      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      include 'PDFerrors.f'
      include 'pdlabel.f'
      double precision amz,alphasPDF
      logical validPDF
      character*30 oldPDFname
      integer i, ii
      logical lhapdfs
      common/lhapdfs/lhapdfs

      common/couple/amz

      lhapdfs=.true.
      
c      if (newinput .eqv. .false.) then
c        open(unit=21,file='lhapdf.DAT',status='old',err=999)
c        call checkversion(21,'lhapdf.DAT')
c        read(21,*) PDFname
c        read(21,*) PDFmember            
c        close(21)
c      endif
      
ch     oldPDFname=PDFname
ch      validPDF=.false.
ch      i=0
ch   20 continue
ch      i=i+1    
ch      if ((oldPDFname(i:i) .eq. '.') .or.
ch     .    (oldPDFname(i:i) .eq. ' ') .or.
ch     .    (oldPDFname(i:i) .eq. '[')) then
ch        validPDF=.true.
ch        if (oldPDFname(i:i+6) .eq. '.LHgrid') then        
ch          PDFname=oldPDFname(1:i-1)//'.LHgrid'
ch        else
ch          PDFname=oldPDFname(1:i-1)//'.LHpdf'
ch        endif
ch      endif  
ch      if ((i .lt. 20) .and. (validPDF .eqv. .false.)) goto 20
      
ch      if (validPDF .eqv. .false.) then
ch        write(6,*) 'Problem with PDFname'
ch        write(6,*)
ch        stop
ch      endif
      

      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*) 'C                                                  C'
      write(6,*) 'C            DYNNLO now calling LHAPDF             C'
      write(6,*) 'C                                                  C'
      write(6,98) 'PDFname',PDFname(1:20)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) 'C                                                  C'
      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*)

      call InitPDFset('PDFsets/'//PDFname)
      amz=alphasPDF(zmass)
      if (PDFmember .lt. 0) then
        PDFerrors=.true.
ch        call numberPDF(maxPDFsets)
ch        if (maxPDFsets .gt. 50) then
ch          write(6,*) 'ERROR: Max. number of error sets is 50!'
ch          stop
ch        endif
        write(6,*)
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'        
        write(6,*) 'C        Calculating errors using      C'
        write(6,*) 'C        ',maxPDFsets,' sets of error PDFs        C'
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        call InitPDF(0)
        amz=alphasPDF(zmass)
        currentPDF=0
      else  
        call InitPDF(PDFmember)
        amz=alphasPDF(zmass)
ch        write(*,*)amz,'ee',alphasPDF(173.3d0),'rrr'
ch        write(*,*)getOrderPDF(),'aaaaaaadadad'
ch        call GetOrderAs(ii)
ch        call alphasQ(173.3d0,amz)
ch        write(*,*)ii,'ii'
      endif

c--- rename pdlabel to get sensible output name
      pdlabel=PDFname(1:7)

      return
 
   98 format(' C            ',a7,' ',a20,'          C')
   99 format(' C                ',a10,i3,'                     C')

  999 write(6,*) 'Error reading lhapdf.DAT'
      call flush(6)
      stop

      end
 

