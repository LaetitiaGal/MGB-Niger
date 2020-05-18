    !*********************************************************************************
    !
    !  SUBROUTINE LECLIMED reads meteorological data
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads meteorological data

	!
	!	 LECLIMED is called inside subroutine MGB_Inercial.
	!
	!	 Saves global variable: 
	!     	XYC = coordinates of stations
	!		TAMM = monthly climatology of near surface temperature
	!       URMM = monthly climatology of near surface relative humidity (%)  
	!       SOLMM = monthly climatology of sunshine hours (hours/day)  
	!       VVMM = monthly climatology of near surface wind speed)  
	!       URMM = monthly climatology of near surface relative humidity (%)  
	!       PAMM = monthly climatology of near surface atmospheric pressure (kPa) 
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
    !    * Opens ACLIMED file  containing monthly climatological data.
    !
    !    reads
    !
    !    * Reads ACLIMED file containing monthly climatological data.
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
    !    * Paulo Pontes RÃ³genes
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



	SUBROUTINE LECLIMED
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER IPCL,K

	! Text:
	CHARACTER (20) ANOME
	CHARACTER (10) AMES(12)

	OPEN(FILMED,FILE=INPUT_DIRECTORY // ''//ACLIMED,STATUS='OLD')

	! Geographic coordinates of stations
	DO IPCL=1,NCLI
!		READ(FILMED,710)XYC(IPCL,1),XYC(IPCL,2)
		READ(FILMED,*)XYC(IPCL,1),XYC(IPCL,2) ! RP formato livre
	ENDDO

	! near surface air temperature (oC)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
	DO IPCL=1,NCLI        
		READ(FILMED,713)ANOME,(TAMM(IPCL,K),K=1,12)        
	ENDDO
	! near surface relative humidity (%)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(URMM(IPCL,K),K=1,12)
	ENDDO
	! Sunshine hours (horas/dia)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
 	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(SOLMM(IPCL,K),K=1,12)
        !write(*,*)(SOLMM(IPCL,K),K=1,12)
	ENDDO
	! Wind speed (m/s)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
  	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(VVMM(IPCL,K),K=1,12)
        
	ENDDO
	! Atmospheric pressure (kPa)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
  	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(PAMM(IPCL,K),K=1,12)        
	ENDDO

	CLOSE(FILMED)

710	FORMAT(2F10.1)
711	FORMAT(A20)
712	FORMAT(A20,12A10)
713	FORMAT(A20,12F10.2)

	RETURN
	END
