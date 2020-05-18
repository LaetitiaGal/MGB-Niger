    !*********************************************************************************
    !
    !  SUBROUTINE PARCEL is the routine that calculates additional catchment and river
	!			parameters necessary for routing and other operations
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine sets catchment additional parameters, mainly width and depth given by
	!		geomorphological relations. Also reference discharge for Muskingum-Cunge method
	!		and concentration time for catchment routing are calculated.
	!
	!	 PARCEL is called in MGB-IPH Main routine.
	!	
    !	 It calls subroutine INTECLIM to assing nearest climate station for each catchment
	!	 	also it calls subroutine REGION to get hydraulic parameters.
	!
    !
    !  	Usage:
    !
    !    CALL PARCEL
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90
    !    * module     VARS_INERC  in      VARS_INERC.f90  
	!    * subroutine INTECLIM    in      INTECLIM.f90 
	!    * subroutine REGION      in      REGION.f90 
    !
    !	 opens
    !
    !    * Does not open files
    !
    !    reads
    !
    !    * Does not read files
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
    !    2014.25.11 - 25 November 2014 (By: Mino V. Sorribas)    
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
    !
    !---------------------------------------------------------------------------------	
	SUBROUTINE PARCEL

	USE VARS_MAIN
	USE VARS_INERC
	IMPLICIT NONE
	INTEGER IW
	INTEGER AUXRG

	! Calculates Time of Concentration (i.e. time of travel) for Catchment
	DO IW=1,NC
		!XLADO=ACEL(IW)**0.5 	! catchment length estimative (old)
		XLADO=LCEL(IW) 			! main river length in catchment (in km) !RP
		!DH=HCEL(IW)-LCEL(IW)
		DH=HCEL(IW)*LCEL(IW) 	! altitude variation in main river (in meters) !RP

		TIND(IW)=3600.*((0.868*XLADO**3.0)/DH)**0.385	!Concentration Time by Kirpich eq. (ni seconds)
	ENDDO

	! Assign nearest climate station
	CALL INTECLIM 				

	!--------------------------------------------------------------------------------------- !PR 15SET14
	! Code for river width and depth based on Geomorphological relations assigned in REGION subroutine	
	if(1==0)then
	! Calculates width, depth and reference discharge for MC-Method
	DO IC=1,NC
		ARX=ACUR(IC)  	! Cumulative Drainage Area
		CALL REGION		! Subroutine calculates river width and full bank depth
		BRIO(IC)=BRX	 ! store river width
		HRIO(IC)=HRX !	 ! store full bank depth
		QREF(IC)=QRX 	 ! store reference discharge for Muskingum-Cunge calcs.
	ENDDO	
	endif
	
	!--------------------------------------------------------------------------------------- !PR 15SET14
	! Code for river width and depth based on Geomorphological relations read in external file
	if(1==0)then	
	OPEN(FILRG,FILE=INPUT_DIRECTORY // 'GEOR.RG',STATUS='OLD')
	DO IW=1,NC
        READ(FILRG,'(I12,3F16.3)')AUXRG,BRIO(IW),HRIO(IW),QREF(IW)
        !WRITE(*,*)IW,AUXRG
    ENDDO
	CLOSE(FILRG)
	!--------------------------------------------------------------------------------------- !PR 15SET14
	endif

	RETURN
	END
