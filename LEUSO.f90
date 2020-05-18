    !*********************************************************************************
    !
    !  SUBROUTINE LEUSO reads  file with calibration parameters (File: PARUSO.CAL)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads file with calibration parameters (File: PARUSO.CAL)

	!
	!	 LEUSO is called inside 1main.
	!
	!	 Saves discharge time series: 
	!     	WM(NB,NU),B(NB,NU),KINT(NB,NU),KBAS(NB,NU),PLAM(NB,NU),CAP(NB,NU),WC(NB,NU) = Parameters related to HRUs: 
	!       CS(NB),CI(NB),CB(NB),CAUX(NB) = Parameters for catchment routing - linear reservoirs 
	!       QESP(NB) = Specific base flow (M3/S/KM2) 
	!
	!
	!
    !
    !  	Usage:
    !
    !    * no subroutine is called in this subroutine
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90
	!
    !	 opens
    !
    !    * file with calibration parameters (File: PARUSO.CAL). 
    !
    !    reads
    !
    !    * file with calibration parameters (File: PARUSO.CAL)
    !
    !    creates
    !
    !    * Does not create files
    !    
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the ...
    !
    !  Version/Modified: 
    !
    !    2015.27.05 - 27 MAY 2015 (By: Rodrigo Paiva)    
    !
    !  Authors:
    !
    !    Original fortran version by Walter Collischonn
    !    Present fortran version by:
    !    * Walter Collischonn
    !    * Rodrigo Cauduro Dias de Paiva
    !    * Diogo da Costa Buarque
    !    * Paulo Pontes Rógenes
    !    * Mino  Viana Sorribas
    !    * Fernando Mainardi Fan
    !    * Juan Martin Bravo 
    !
    !  Main Reference:
    !
    !    Walter Collischonn,
    !    Modelo de Grandes Bacias - Thesis
    !    Porto Alegre, 2001
    !    ISBN: XXXXXXXXXXX,
    !
    !---------------------------------------------------------------------------------
    !  Variables and Parameters:
    !
    !   *Variables declarations and routines calls are all commented below.
	!	* All variables are global!?
    !
    !---------------------------------------------------------------------------------		

	SUBROUTINE LEUSO
	! Reads file with calibration parameters (File: PARUSO.CAL)
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER KB,K !,IU
	!Parameters related to HRUs: WM(NB,NU),B(NB,NU),KINT(NB,NU),KBAS(NB,NU),PLAM(NB,NU),CAP(NB,NU),WC(NB,NU)
	REAL AUX(NB,NU),CAUX(NB)
	! Parameters for catchment routing - linear reservoirs : CS(NB),CI(NB),CB(NB),CAUX(NB)
	! Specific base flow (M3/S/KM2): QESP(NB)
	WM=0.0
	B=0.0
	KINS=0.0
	KBAS=0.0
	PLAM=0.0
	AUX=0.0
	CS=0.0
	CAUX=0.0
	CB=0.0
	CI=0.0
	OPEN(FILUSO,FILE=INPUT_DIRECTORY // 'PARUSO.cal',STATUS='OLD')
	DO KB=1,NB
		READ(FILUSO,72)CABE(1)
		READ(FILUSO,72)(CABE(K),K=1,5)
		!write(*,*) kb
		DO IU=1,NU
			READ(FILUSO,71)AUSO(IU),WM(KB,IU),B(KB,IU),KBAS(KB,IU),KINS(KB,IU),PLAM(KB,IU),CAP(KB,IU),WC(KB,IU)

		ENDDO
		READ(FILUSO,73)CABE(1),CS(KB)
		READ(FILUSO,73)CABE(2),CI(KB)
		READ(FILUSO,73)CABE(3),CB(KB)
		READ(FILUSO,74)CABE(4),QESP(KB)
		!write(*,*) CS(KB),CI(KB),CB(KB),QESP(KB)
	ENDDO
	CLOSE (FILUSO)
71	FORMAT(A10,7F8.3)
72	FORMAT(7A10)
73	FORMAT(A10,F8.1)
74	FORMAT(A10,F8.6)
	RETURN
	END
