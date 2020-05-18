    !*********************************************************************************
    !
    !  SUBROUTINE LECHUVA reads precipitation data
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads precipitation data

	!
	!	 LECHUVA is called inside subroutine Modelo.
	!
	!	 Saves global variable: 
	!     	P = precipation at time t for all catchments (mm/dt)
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
    !    * Does not open files
    !
    !    reads
    !
    !    * CHUVAbin.pbi, the precipitation binary file
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

	SUBROUTINE LECHUVA
	USE VARS_MAIN
    USE VARS_INERC
	IMPLICIT NONE
	INTEGER KC
	! Reads precipitation data at current time interval (ITCHUVA)
	READ(FILPLU,REC=ITCHUVA)(P(KC),KC=1,NC)
	!WRITE(*,*)p(1)
	! PAUSE
    
    !DISABLE PRECIPITATION IN DELTA
    !DO KC=1,NC
    !IF(INLAND_DELTA(KC)==1) P(KC)=0
    !ENDDO
    
	RETURN
	END
