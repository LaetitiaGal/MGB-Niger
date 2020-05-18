    !*********************************************************************************
    !
    !  SUBROUTINE LEQOBS reads the file observed streamflow data (File with extension .QOB)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads the file observed streamflow data (File with extension .QOB)

	!
	!	 LEQOBS is called inside 1main.
	!
	!	 Saves discharge time series: 
	!     	QOBS(.,.)
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
    !    * File with extension .QOB containing time series of observed stream flow. 
    !
    !    reads
    !
    !    * File with extension .QOB containing time series of observed stream flow
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
    !    2014.26.11 - 25 November 2014 (By: Rodrigo Paiva)    
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

	SUBROUTINE LEQOBS
		USE VARS_MAIN
	IMPLICIT NONE
	INTEGER I,J,K,L
	OPEN(FILOBS,FILE=INPUT_DIRECTORY // ''//ARQOBS,STATUS='OLD')
	READ(FILOBS,701)(CABE(K),K=1,NOBS)
	! write(*,*)'NOBS = ', NOBS
	! write(*,*)'NT = ', NT
	! pause
	DO IT=1,NT
		! READ(FILOBS,702)CABE(1),(QOBS(K,IT),K=1,NOBS)
		READ(FILOBS,*)I,J,L,(QOBS(K,IT),K=1,NOBS)
		! write(*,*)'QOBS=', (QOBS(K,IT),K=1,NOBS)
		! pause
		! IF (IT == NT) then
		! WRITE(*,*)'QOBS = ', QOBS(K,IT)
		! PAUSE
		! ENDIF
	ENDDO
	CLOSE (FILOBS)
	701	FORMAT(<NOBS>A10)
	!702 FORMAT(A15,F9.2)
	702	FORMAT(2I7,I6,258F16.6)
	! END READ OBS DISCHARGE DATA --------------------!
	RETURN
	END
