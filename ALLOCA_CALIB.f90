	SUBROUTINE ALLOCA_CALIB(IOP)
!---------------------------------------------------------------------------------
! ALLOCA_CALIB.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This routine allocate or deallocate the MGB-IPH automatic calibration procedures variables.
!
! Usage:
!
! CALL ALLOCA_CALIB(IOP)
!
! uses modules, functions, and subroutines
!
! *Module VARS_CALIB !module with the calibration variables
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
! 2014.11.11 - 11 November 2014
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

!Variables and modules

	USE VARS_CALIB
	IMPLICIT NONE
	INTEGER IOP
	SAVE
	
!First we verify if we are doing allocation or deallocation. If IOP =0 then allocate. If IOP =0 then deallocate (clean).

	ALLOC_CASE: SELECT CASE (IOP) !Verify if we are doing allocation or deallocation
	CASE (0) !allocation
		ALLOCATE (IRUIM(NS))
		ALLOCATE (FO(NS,NF)) !Objective functions
		ALLOCATE (IPARET(NS))
		ALLOCATE (XPAR(NPAR))
		ALLOCATE (SOMAPAR(NPAR),REFLEX(NPAR),CONTRA(NPAR)) !Sums of COORD., reflection point, point of contraction
		ALLOCATE (MEDIA(NS)) !Average of objective functions of the population
		ALLOCATE (PPAR(NPAR,NS)) !relative value of the parameter at the validity range
		ALLOCATE (PAR(NPAR,NS)) !parameter value 
!		ALLOCATE (PMIN(NPAR),PMAX(NPAR)) !validity range of the parameters
		ALLOCATE (FOLD(NS)) !former value of the function at each point of the sample
		ALLOCATE (PROB(NS)) !probability of choice of point 1..NS
		ALLOCATE (VMIN(NF)) !minimum value of objective functions
		ALLOCATE (PARX(NPAR)) !parameter value
		
	CASE (1) !deallocation
		DEALLOCATE (IRUIM,FO,IPARET,XPAR)
		DEALLOCATE (SOMAPAR,REFLEX,CONTRA) !Sums of COORD., reflection point, point of contraction
		DEALLOCATE (MEDIA) !Average of objective functions of the population
		DEALLOCATE (PPAR) !relative value of the parameter at the validity range
		DEALLOCATE (PAR) !parameter value
		DEALLOCATE (PMIN,PMAX) !validity range of the parameters
		DEALLOCATE (FOLD) !former value of the function at each point of the sample
		DEALLOCATE (PROB) !probability of choice of point 1..NS
		DEALLOCATE (VMIN) !minimum value of objective functions
	CASE DEFAULT
		STOP ' ERROR: UNKNOWN VALUE OF IOP AT SUBROUTINE ALLOCA_CALIB!!!'
	END SELECT ALLOC_CASE

	!End the routine
	RETURN
	END
