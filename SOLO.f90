    !*********************************************************************************
    !
    !  SUBROUTINE SOLO is the routine that calculates vertical water balance at URH
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine calculates URH flow generation for surface, subsurface, groundwater
	!		pathways and also capilarity.
	!		Evaporation can be reduced to account for drying.
	!
	!	 SOLO is called inside subroutine CELULA.
	!
	!    In general the main variables saved are DSUP, DINT and DBAS
	!
	!	 *Correction from Fernando & Paulo to Parana River Basin!
    !
    !  	Usage:
    !
    !    CALL SOLO(PX,EX,WX,WMX,BX,KIX,KBX,XL,DSUP,DINT,DBAS,CAPX,WCX,DTP)
    !
    !    where
    !
    !    * PX   is the precipitation (input)
	!	 * EX   is the evapotranspiration (input/output)
	!	 * WMX  is the in-soil water storage (input/parameter)
	!	 * BX   is the b parameter for saturation flow in ARNO model (input/parameter)
	!	 * KIX  is the subsurface hydraulic conductivity (input/parameter)
	!	 * KBX  is the groundwater hydraulic conductivity (input/parameter)
	!	 * XL   is the lambda parameter, a porosity index for subsurface flow (input/parameter)
	!	 * DSUP is the water-depth of surface water generated (output)
	!	 * DINT is the water-depth of subsurface water generated (output)
	!	 * DBAS is the water-depth of groundwater generated (output)
	!	 * CAPX is the water-depth of water ascension by capilarity (output)
	!	 * WCX	is the moisture threshold to start capilarity ascension (input/parameter)
	!	 * DTP	is the time step in seconds (input)
    !
    !    uses modules and functions
    !
    !    * Does not use modules,functions and subroutines
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
    !   *Variables declarations and routines calls are all commented below or above.	
    !
    !---------------------------------------------------------------------------------	
	SUBROUTINE SOLO(PX,EX,WX,WMX,BX,KIX,KBX,XL,DSUP,DINT,DBAS,CAPX,WCX,DTP,IC,ACELX)
    USE VARS_INERC
	IMPLICIT NONE
    
	REAL WX				! In-soil water storage	
	REAL BX,WMX 		! Parameters for storage
	REAL KIX,KBX 		! Hydraulic Conductivity
	REAL WZ,WB 			! Limits for subsurface and gw drainage
	REAL XL 			! shape parameter for subsurface flow
	REAL CAPX,WCX,WC 	! ascending fluxes			
	REAL PX,EX			! precipitation and evaporation	
	REAL DSUP,DINT,DBAS	! surface, subsurface and groundwater drainage	
	REAL VTEMP			! auxiliar
	REAL WXOLD			! auxiliar
	REAL DTP 			! time step in seconds
    INTEGER IC
    REAL ACELX

	! Some notes
	!XL controls subsurface drainage and soil moisture according to Brooks and Corey (see Maidment)	
	!XL is the pore size index w values between 0.165 (clay) and 0.694 (sand) , according to Rawls.
	
	
	!limiters	
	WB=WCX*WMX  ! WB is the limit to groundwater generation
				! WCX is a temporary variable with values equal to WC(IB,IU) defined in parameter file.
	
	WZ=WB		! WZ is the limit to subsurface water generation
				! Assumption: equals to WB
	
	WC=WCX*WMX	! WC is the limit to start ascending capilar flux

	! Calculates direct flow (surface waters)
	IF(PX.GT.0.0) THEN
	
		VTEMP=(1.-WX/WMX)**(1./(BX+1.))-(PX/((BX+1.)*WMX))
		
		IF(VTEMP.LE.0.0) THEN ! Rainfall saturates soil
		
				!DSUP=PX-(WMX-WX)
				DSUP=MAX(PX-(WMX-WX),0.0) !Fernando & Paulo Correction 25/01/2013)
				
		ELSE                  ! Rainfall doesnt saturates soil
		
				!DSUP=PX-(WMX-WX)+WMX*(VTEMP)**(BX+1.)
				DSUP=MAX(PX-(WMX-WX)+WMX*(VTEMP)**(BX+1.),0.0)			!Fernando & Paulo Correction 25/01/2013)

		ENDIF
		
	ELSE
	
		DSUP=0.0 ! if it doesnt rain, there is no surface flow
		
	ENDIF
	

	! Calculates groundwater (slow flows)
	IF(WX.GE.WB) THEN
		DBAS=KBX*(WX-WB)/(WMX-WB)	!simple linear approach
	ELSE
		DBAS=0.0
    ENDIF
	
    !Calculates infiltration from floodplain (Inertial module only)
    IF (HDFLAG(iC)==0) THEN   
        DINFILTX=0 !Inertial module is disabled
    ELSE
        DINFILTX=(Area2(ic)/ACELX)*20*(1-WX/WMX)**1!Inertial module is enabled !in mm/day (Kinfilt=5 mm/d)
       ! DINFILTX=(Area2(ic)/ACELX)*5*(1-WX/WMX)**1!Inertial module is enabled !in mm/day (Kinfilt=5 mm/d) ! test KINF Cécile
        !DINFILTX=(Area2(ic)/ACELX)*40*(1-WX/WMX)**1!Inertial module is enabled !in mm/day (Kinfilt=5 mm/d) ! test KINF Cécile

    ENDIF
    DINFILTX=MAX(DINFILTX,0.0)
    !if (ic>4219.and.ic<4274)  DINFILTX=0  !disable infiltration between Diré and Ansongo
    !DINFILTX=0 !disable infiltration

    ! Calculates capilarity
	CAPX=MAX(CAPX*(WC-WX)/WC,0.0) 	! Bremicker.	

	! Calculates subsurface flow (intermediate)
	! uses Brooks & Corey relation (5.6 from Handbook of Hydrology)
	IF(WX.LE.WZ) THEN
		DINT=0.0 		! No drainage in low moisture condition
	ELSE
		DINT=KIX*((WX-WZ)/(WMX-WZ))**(3.+2./XL)	
	ENDIF

	! Limits drainage to avoid negative storage
	DINT=MIN(DINT,WX/10.0) ! this assumptions are ok cause
	DINT=MAX(DINT,0.0)     ! hydr. conductivy can drop fast
	
	! Updates in-soil stored water for calibration
	WXOLD=WX
	
	! Convert to time-step when necessary (smaller than daily)
	DSUP=DSUP 			
	CAPX=CAPX*DTP/86400.	
	DINT=DINT*DTP/86400.
	DBAS=DBAS*DTP/86400.
	DINFILTX=DINFILTX*DTP/86400.
    
	! Starts vertical balance loop
	!  ends only if calculated storage is not negative
	DO 		
			
		!WC DTP WX=WX+PX-EX-DSUP-DINT-DBAS+CAPX  !water balance
		WX=WX+PX+(DINFILTX+CAPX-EX-DSUP-DINT-DBAS)  		 !updates water balance

!		! Older adjustment (Walter)
!		WX=MIN(WX,WMX) !check rounding error
!		IF(WX.LT.0.0)THEN
!			WRITE(*,*) 'Soil went dry'
!			! reduce evaporation cause soil is drying
!			EX=WXOLD+PX-DSUP-DINT-DBAS 
!			WX=WXOLD
!			CYCLE
!		ELSE
!			EXIT
!		ENDIF

		! New adjustment (Paiva)
		! Considers adjusting evapotranspiration and drainage
		if (WX<WMX+0.01.and.WX>WMX) WX=WMX
		if (WX>-0.01.and.WX<0.0) WX=0.0
		if(WX>WMX)then
			! if in-soil water storage is higher than maximum, increase surface drainage			
			DSUP=DSUP+WX-WMX
			WX=WXOLD
!write(*,*) '4'

		elseif (WX<0.0) then

			!write(*,*) '4.5'
			! if in-soil water in less than zero, decrease evapotranspiration
			EX=EX-(-WX)
			if (EX>=0.0) then
				WX=WXOLD
				cycle				
			endif			
			! if evaporation is zero and water balance is not closed
			WX=EX 		!accept storage as the residual
			EX=0.0		!then zero evaporation

			!write(*,*) '5'
			! reduces percolation
			DBAS=DBAS-(-WX)
			if (DBAS>=0.0) then
				WX=WXOLD
				cycle
			endif
			! percolation is zero but water balance is not closed
			WX=DBAS    !accept storage as the residual
			DBAS=0.0	!zero groundwater flux

			
			!write(*,*) '6'
			! reduces subsurface
			DINT=DINT-(-WX)
			if (DINT>=0.0) then
				WX=WXOLD
				cycle
			endif			
			! subsurface is zero but water balance is not closed
			WX=DINT 	!accept storage as the residual
			DINT=0.0			
			
			!write(*,*) '7'
			! reduces surface water
			DSUP=DSUP-(-WX)
			if (DSUP>=0.0) then
				WX=WXOLD
				cycle
			endif
			DSUP=0.0
			
			!write(*,*) '8'
			! surface water is zero but water balance didnt close
			! nothing to do anymore.
			write(*,*) 'Soil went dry'
			WX=0.001
			exit
		else
			exit
		endif
	
!		write(*,*) '9'

	ENDDO


	RETURN
	END
