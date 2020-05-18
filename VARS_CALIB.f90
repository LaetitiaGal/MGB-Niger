	MODULE VARS_CALIB
	
!---------------------------------------------------------------------------------
! VARS_CALIB.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This module contains the variables of MGB-IPH model used in the automatic calibration procedures
!
! Usage:
!
! USE VARS_CALIB
!
! uses modules, functions, and subroutines
!
! * no modules, functions, or subroutines are used in this module
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
	
	!Variables declaration
	
	IMPLICIT NONE
	SAVE
	INTEGER ISEED
	INTEGER KOUNTF,ISHUFFLE,ISOMA,IDOMIN,IDD,NPC
	REAL FMIN,SPARET,RMAX,XRAND,PINI,DIFPAR
	INTEGER IREJECT,LRUIM
	REAL,ALLOCATABLE:: FO(:,:)
	INTEGER, ALLOCATABLE:: IPARET(:),IRUIM(:)
	REAL, ALLOCATABLE:: PPAR(:,:) 
	REAL, ALLOCATABLE:: PAR(:,:),PARX(:)
	REAL, ALLOCATABLE:: PMIN(:),PMAX(:)
	REAL, ALLOCATABLE:: FOLD(:)
	REAL, ALLOCATABLE:: MEDIA(:)
	REAL,ALLOCATABLE::XPAR(:)
	REAL,ALLOCATABLE:: SPAR(:,:,:),FPLEX(:,:,:)
	REAL,ALLOCATABLE:: SOMAPAR(:),REFLEX(:),CONTRA(:)
	REAL,ALLOCATABLE:: PROB(:)
	REAL,ALLOCATABLE:: VMIN(:) !Minimum value of the objective functions
	INTEGER:: NPAR !Number of parameters of the function to be optimized
	INTEGER:: NS !Number of points in the initial sample
	INTEGER:: NF !Number of objective functions to be optimized
	
	INTEGER,ALLOCATABLE:: p_Calib(:,:),iFO(:),IBCONGEL(:)
	INTEGER:: iMaxGen,NCONGEL

	!End of module
	END MODULE VARS_CALIB
