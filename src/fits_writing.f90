subroutine writeimage(array,nvel,npix,filename,restfreq,&
& cdelt1,crpix1,cdelt2,crpix2,cdelt3,crpix3,&
& crval1,crval2,crval3)

!  Create a FITS primary array containing a 3-D image

!       input:
!       array   --- npix,npix,nvel      array of image 
!       nvel    --- number of velocity channels 
!       npix    --- number of pixels (npix1=npix2=npix)
!       filename--- output file name 
!       restfreq--- rest frequency of transition 
!       cdelt1/2--- pixel distance in degree 
!       crpix1/2--- reference pixel --- defining the origin 
!       cdelt3  --- velocity channel distance (in m/s)
!       crpix3  --- reference velocity channel 
!       crval   --- reference values (1,2 decl 3 velocity) 

      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,k,group,fpixel,nelements
      character filename*80, object*8
      logical simple,extend
      real array(npix,npix,nvel)
      real bscale,bzero,epoch,lonpole,equinox,restfreq
      real cdelt1,crpix1,crval1,cdelt2,crpix2,crval2
      real cdelt3,crpix3,crval3
      integer npix, nvel, velref


      status=0
!      filename='ATESTFILEZ.FITS'

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      call deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      call ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      call ftinit(unit,filename,blocksize,status)

!  Initialize parameters about the FITS image.
!  BITPIX = -32 means that the image pixels will consist of 32-bits
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.

!     placeholder for image size and number of velocity channels
!      npix=200
!      nvel=50
!
      simple=.true.
      bitpix=-32
      naxis=3
      naxes(1)=npix
      naxes(2)=npix
      naxes(3)=nvel
      extend=.true.

!  Write the required header keywords to the file
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Initialize the values in the image with a linear ramp function
!  this is just to fill an array

!      do k=1,naxes(3) 
!      do j=1,naxes(2)
!          do i=1,naxes(1)
!              array(i,j,k)=REAL((i - 1 +j - 1)*ABS(k-25))*0.001
!          end do
!      end do
!      end do

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'E' indicates that the array has an
!  real*4 datatype. ('I' = I*2, 'J', I*4, 'E' = Real*4, 'D' = Real*8).
! ['D' does not work]
!  The 3D array is treated as a single 1-D array with NAXIS1 * NAXIS2 * NAXIS3
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      
      call ftppre(unit,group,fpixel,nelements,array,status)
   
!  Write other keywords for CASA in the header
!  The keyword record will look like this in the FITS file:
!
! these need to be changed: cdelt1 is the pixel size in RA (always negative, now 1 arcsec)
!   crpix is the reference pixel (typically the center pixel)
!  cdelt2 is the pixel size in DEC (also 1 arcsec now)
!  cdelt3 is the velocity channel width (the cube is in velocity) in M/S
!     cpix3 the center channel

!      restfreq=345e9
!      cdelt1  =-2.7778e-4  
!      crpix1  =100e0
!      cdelt2  =2.7778e-4
!      crpix2  =100e0
!      cdelt3  =1.0e3
!      crpix3  =2.5e1

! these can remain like this (the crvals are the value of the center pixel/channel, so now the
! map is towards RA,DEC 0,0 and center velocity 0

      object = 'PORTAL'
      epoch   =2.0d3
      lonpole =1.8e2
      equinox =2.0e3
      velref  =257
!      crval1  =0.0e0
!      crval2  =0.0e0
!      crval3  =0.0e0
      bscale  =1.0e0
      bzero   =0.0e0


      call ftpkys(unit,'OBJECT',object,'',status)
      call ftpkye(unit,'EPOCH',epoch,8,'',status)
      call ftpkye(unit,'LONPOLE',lonpole,8,'',status)
      call ftpkye(unit,'EQUINOX',equinox,8,'',status)
      call ftpkys(unit,'SPECSYS','LSRK','',status)
      call ftpkye(unit,'RESTFREQ',restfreq,8,'',status)
      call ftpkyj(unit,'VELREF',velref,'',status)
      call ftpkys(unit,'CTYPE1','RA---SIN','',status)
      call ftpkye(unit,'CDELT1',cdelt1,8,'',status)
      call ftpkye(unit,'CRPIX1',crpix1,8,'',status)
      call ftpkye(unit,'CRVAL1',crval1,8,'',status)
      call ftpkys(unit,'CUNIT1','DEG','',status)
      call ftpkys(unit,'CTYPE2','DEC--SIN','',status)
      call ftpkye(unit,'CDELT2',cdelt2,8,'',status)
      call ftpkye(unit,'CRPIX2',crpix2,8,'',status)
      call ftpkye(unit,'CRVAL2',crval2,8,'',status)
      call ftpkys(unit,'CUNIT2','DEG','',status)
      call ftpkys(unit,'CTYPE3','VELO-LSR','',status)
      call ftpkye(unit,'CDELT3',cdelt3,8,'',status)
      call ftpkye(unit,'CRPIX3',crpix3,8,'',status)
      call ftpkye(unit,'CRVAL3',crval3,8,'',status)
      call ftpkys(unit,'CUNIT3','M/S','',status)
      call ftpkye(unit,'BSCALE',bscale,8,'',status)
      call ftpkye(unit,'BZERO',bzero,8,'',status)
      call ftpkys(unit,'BUNIT','K','',status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end
! *************************************************************************

      subroutine readheader(filename)

!  Print out all the header keywords in all extensions of a FITS file

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character filename*80,record*80

!  The STATUS parameter must always be initialized.
      status=0

!  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

!     name of FITS file 
!      filename='ATESTFILEZ.FITS'

!     open the FITS file, with read-only access.  The returned BLOCKSIZE
!     parameter is obsolete and should be ignored. 
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

      j = 0
100   continue
      j = j + 1

      print *,'Header listing for HDU', j

!  The FTGHSP subroutine returns the number of existing keywords in the
!  current header data unit (CHDU), not counting the required END keyword,
      call ftghsp(unit,nkeys,nspace,status)

!  Read each 80-character keyword record, and print it out.
      do i = 1, nkeys
          call ftgrec(unit,i,record,status)
          print *,record
      end do

!  Print out an END record, and a blank line to mark the end of the header.
      if (status .eq. 0)then
          print *,'END'
          print *,' '
      end if

!  Try moving to the next extension in the FITS file, if it exists.
!  The FTMRHD subroutine attempts to move to the next HDU, as specified by
!  the second parameter.   This subroutine moves by a relative number of
!  HDUs from the current HDU.  The related FTMAHD routine may be used to
!  move to an absolute HDU number in the FITS file.  If the end-of-file is
!  encountered when trying to move to the specified extension, then a
!  status = 107 is returned.
      call ftmrhd(unit,1,hdutype,status)

      if (status .eq. 0)then
!         success, so jump back and print out keywords in this extension
          go to 100

      else if (status .eq. 107)then
!         hit end of file, so quit
          status=0
      end if

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end

! *************************************************************************
      subroutine printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.

      integer status
      character errtext*30,errmessage*80

!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return

!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
! *************************************************************************
      subroutine deletefile(filename,status)

!  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!  Simply return if status is greater than zero
      if (status .gt. 0)return

!  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
!         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!  Free the unit number for later reuse
      call ftfiou(unit, status)
      end
