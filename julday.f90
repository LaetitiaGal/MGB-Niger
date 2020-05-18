      FUNCTION julday(mm,id,iyyy)

!---------------------------------------------------------------------------------
! julday.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This function determines the day of the Julian Calendar corresponding to the given day, mounth and year
! Does the opposite of caldat.f90
!
! Usage:
!
! CALL julday(mm,id,iyyy)
!
! uses modules, functions, and subroutines
!
! * no modules, functions, or subroutines are used in this funcation
!
! opens
!
! * no files are opened in this routine
!
! reads
!
! * no files are read in this routine
!
! creates
!
! * no files are created in this routine
!
!---------------------------------------------------------------------------------
! Licensing:
!
! This code is distributed under the...
!
! Version/Modified:
!
! 2014.11.12 - 12 November 2014
! By: Fernando Mainardi Fan
!
! Authors:
!
! Original fortran version by Walter Collischonn
! Present fortran version by:
! * Walter Collischonn
! * Rodrigo Cauduro Dias de Paiva
! * Diogo da Costa Buarque
! * Paulo Pontes Rógenes
! * Mino Viana Sorribas
! * Fernando Mainardi Fan
! * Juan Martin Bravo
!
! Main References:
! COLLISCHONN, W. ; ALLASIA, D. G. ; SILVA, B. C. ; TUCCI, C. E. M. The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
! COLLISCHONN, W., TUCCI, C. E. M. Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n. 2, 2001.
! COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------
	
!Function:

      INTEGER julday,id,iyyy,mm,IGREG
      PARAMETER (IGREG=15+31*(10+12*1582))
      INTEGER ja,jm,jy
      jy=iyyy
      if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
        jm=mm+1
      else
        jy=jy-1
        jm=mm+13
      endif
      julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
      if (id+31*(mm+12*iyyy).ge.IGREG) then
        ja=int(0.01*jy)
        julday=julday+2-ja+int(0.25*ja)
      endif

      return
      END
      !End of function
