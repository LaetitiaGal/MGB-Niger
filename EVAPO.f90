    !*********************************************************************************
    !
    !  SUBROUTINE EVAPO calculates evaporation and evapotranspiration
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
	!    This routine calculates evapotranspiration for a given URH
	!	 It is called inside sobroutine CELULA and updates evaporation variables.
	! 		EX: evaporation through soil+plant (EX)
	! 		ETX: total evaporation (interceptation+ soil/plant)
	! 		EIX : potential interceptation
	!               EIXP : potential evaporation
	!
	!	* Maximum Interception  is adapted for Parana River Basin! multipler is 0.4*LAI
	!	
	!
    !  Usage:
    !
    !    CALL EVAPO
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

	SUBROUTINE EVAPO
	USE VARS_MAIN
	IMPLICIT NONE


	!REAL WX,WMX,PAX,TAX,VVX,RIAFX,DTP
	!REAL RLX 							!Net Radiation (MJ/m2/dia)
	!REAL CLATE 						!Latent Heat (vap)
	!REAL ED,ES 						!Vapor Pressure (real,saturation)
	!REAL D 							!Vapor pressure deficit
	!REAL DEL 							!Derivative of ES x T function
	!REAL GAMA 							!Psicrometric constant
	!REAL MSPA 							!Air specific mass
	!REAL ZX,RSX 						!Height and Surface Resistance from Vegetation
	!REAL RUG 							!Surface roughness
	!REAL SIX 							!Water depth in Interceptation
	!REAL SF 							!Water depth in Interceptation
	!REAL SIL 							!Maximum water depth in interceptation
	!REAL EIX 							!Water depth evaporated from interceptation (potential)
	!REAL REIX 							!Water depth evaporated from interceptation (real)
	!REAL FDE 							!Fraction of evaporative demand (in soil/plant)
	!REAL PX 							!Precipitation
	!REAL,PARAMETER:: MESP=1000.0 	    !Water Specific Weigth (KG/M3)
	!REAL,PARAMETER:: CSPA=0.001013 	!Specific Heat of Wet Air  (MJ/kg/C)
	!REAL RA 							!Aerodynamic resistance
	!REAL WL 							!Moisture limit above which soil moisture wont limit evapotranspiration (Shuttleworth)
	!REAL WPM 							!Wilting point
	!REAL F4 							!Correction factor for surface resistance in function of soil moisutre
	!REAL EX 							!Evaporation in soil-plant
	!REAL ETX 							!Total (real) Evapotranspiration (water intercepted + soil + plant)

	D=ES-ED 				 			! Vapor pressure deficit in KPa
	GAMA=0.0016286*PAX/CLATE 			! KPa/C (EQ. 4.2.28)
	DEL=(4098.0*ES)/((237.3+TAX)**2.) 	! (EQ. 4.2.3)
	MSPA=3.486*(PAX/(275.0+TAX)) 		! KG/M3 (EQ. 4.2.4)

	RLX=RLX/(3600.0*24.0) 				! Net Radiation in MJ/m2/s

	! Aerodynamic roughness of vegetation
	RUG=0.1*ZX


	! Aerodynamic resistance is based on Bremicker work		
	!Fix wind speed at 2m to 10m above the ground
	
    !IF (flagaclimed == 1) then
    !IF (CRU_INDEX == 1) then
		! If use CRU data, then wind is already at 10m
	!	VVX=VVX
	!else
		!VVX=VVX*((LOG(10.0/RUG))/LOG(2.0/RUG)) !EQ. 3.46		!why it is commented ??? #mv
!		write(*,*)'VVX,RUG = ',VVX,RUG
!		pause
	!endif	

	VVX=VVX+0.1 !	Avoid zero value issues    
    VVX2 = VVX/((LOG(10.0/RUG))/LOG(2.0/RUG)) !wind speed at 2m for calculation of RA and potential evapotranspiration


	!Aerodynamic resistance of LARSIM model
	!IF(ZX.LT.10.0) THEN
	!	RA=(6.25/VVX)*(LOG(10.0/RUG))**2.0 !EQ. 3.25
	!ELSE
	!	RA=94.0/VVX !EQ. 3.26
	!
	!ENDIF


	! Aerodynamic resistance according to Shuttleworth
	ZX=MIN(ZX,12.5) !limits vegetation height to avoid problems
	RA=((LOG((10.-ZX*0.67)/(ZX*0.123)))*(LOG((10.-ZX*0.67)/(ZX*0.0123))))/(0.41*0.41*VVX)
	RAEP=(4.72*(LOG(2./0.00137))**2.)/(1.+0.536*VVX2)   !RA for open water evaporation ! (EQ.4.2.29) from Shuttleworth
!    write(*,*)'RA-evapo,zx,vvx = ',RA,zx,vvx
!    pause
	
	! Calculate Surface Resistance (RSX) for soil moisture conditions
	WL= 0.8*WMX 		! Limit from which soil moisture dont limit evapotranspiration  ! MODIF POST ETP Cécile
	!WL= 0.5*WMX 		! Limit from which soil moisture dont limit evapotranspiration
	WPM=0.1*WMX  		! Wilting Point
	!RSX=RSX 			! RS is read in fixed parameters file

	! Surface Resistance Variation with Soil moisture (Wigmosta et al. 1994 - eq. 16)
	IF(WMX<0.1)THEN 		! Open water has no surface resistance
		RSX=0.0
	ELSE
		IF(WX.LE.WPM) THEN
			RSX=999999999.0 ! Dry soil, assume superhigh surface resistance
		ELSEIF(WX.LT.WL)THEN
			F4=(WX-WPM)/(WL-WPM)
			F4=1./F4		! Moisture deficit
		ELSE
			F4=1.0
		ENDIF
		RSX=F4*RSX			!if F4>1 increases surface resistance
	ENDIF

	! Calculates evaporation from intercepted water
	! Maximum Interception
	!SIL=0.2*RIAFX !*(DTP/86400.)	!SIL=0.2*RIAFX !According to Wigmosta et al. 1994, also used by Bremicker
    ! Maximum Interception ! Modified by Fernando & Paulo (25/01/2013) to land-use in Parana Basin
	SIL=0.4*RIAFX !*(DTP/86400.)	!SIL=0.2*RIAFX !According to Wigmosta et al. 1994, also used by Bremicker

	! Adds Rainfall to water that can be evaporated.
	SF=MIN(SIX+PX,SIL)
	! Updates Rainfall
	PX=PX-(SF-SIX)

	! Potential Evaporation from interceptation (note: rsx = 0.0)	
	EIX=((DEL*RLX+MSPA*CSPA*D/RA)/(DEL+GAMA*(1.0+0.0/RA)))*((1000.0*DTP)/(CLATE*MESP)) 		!mm/dtp
	EIXP=((DEL*RLX+MSPA*CSPA*D/RAEP)/(DEL+GAMA*(1.0+0.0/RAEP)))*((1000.0*DTP)/(CLATE*MESP)) !mm/dtp !open water (Shuttleworth)
    IF (IC==1760) ETPOTENCIAL(IT)=EIXP !PARA AVALIAÇÃO DE ET POTENCIAL
    
	! Real Evaporation (cant be larger than interceptation)
	REIX=MIN(EIX,SF) 			! evaporation of intercepted water
	! Takes from interceptation reservoir
	SIX=SF-REIX	 
	
	! Calculates evaporative demand fraction after evaporation
	FDE=(EIX-REIX)/EIX

	! Calculates Penman in mm/dtpCALCULA PENMANN DIRETO EM MM/DIA
	!EX=((DEL*RLX+MSPA*CSPA*D/RA)/(DEL+GAMA*(1.0+RSX/RA)))*((1000.0*3600.0*24.)/(CLATE*MESP)) !mm/day
	EX=((DEL*RLX+MSPA*CSPA*D/RA)/(DEL+GAMA*(1.0+RSX/RA)))*((1000.0*DTP)/(CLATE*MESP)) 		  !mm/dtp
	EX=FDE*EX   																	! evaporation through soil+plant 
	ETX=EX+REIX 																	! total evaporation (interceptation+ soil/plant)

	EX=EX*1.0
	
	RETURN
	END
