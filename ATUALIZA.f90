	SUBROUTINE ATUALIZA
!---------------------------------------------------------------------------------
! ATUALIZA.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This subroutine update the models state variables (data assimilation) when the model runs in forecast mode.
! The method used is based in the comparison between simulated and forecasted flows.
! The method description can be found in:
! FAN, F. M. ; MELLER, A. ; COLLISCHONN, W. . Incorporação de filtro numérico de separação de escoamento na assimilação de dados para previsão de vazões utilizando modelagem hidrológica. Revista Brasileira de Recursos Hídricos, 2015. 
! PAZ AR, COLLISCHONN W, TUCCI C, CLARKE R, ALLASIA D. Data Assimilation in a Large-scale Distributed Hydrological Model for Medium Range Flow Forecasts. IAHS Press, Wallingford, UK, IAHS Publication, n. 313, p. 471–478, 2007.
!
! Usage:
!
! CALL ATUALIZA
!
! uses modules, functions, and subroutines
!
! * VARS_MAIN.f90
!
! opens
!
! * No files are open
!
! reads
!
! * No files are read
!
! creates
!
! * No files are created
!
!---------------------------------------------------------------------------------
! Licensing:
!
! This code is distributed under the...
!
! Version/Modified:
!
! 2014.11.12 - 22 November 2014
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
! FAN, F. M. ; MELLER, A. ; COLLISCHONN, W. . Incorporação de filtro numérico de separação de escoamento na assimilação de dados para previsão de vazões utilizando modelagem hidrológica. Revista Brasileira de Recursos Hídricos, 2015. 
! PAZ AR, COLLISCHONN W, TUCCI C, CLARKE R, ALLASIA D. Data Assimilation in a Large-scale Distributed Hydrological Model for Medium Range Flow Forecasts. IAHS Press, Wallingford, UK, IAHS Publication, n. 313, p. 471–478, 2007.
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------
	
	!Call used modules
	USE VARS_MAIN
	
	!Variables declaration
	IMPLICIT NONE
	REAL PESOBAC
	INTEGER ITR,NTRECH,NTMUSK,LT,I,K,KHID
	INTEGER POSTORIO(NB) !Indicated the gauge that is going to be used to update the streamflow
	INTEGER POSTOBAS(NB) !Indicated the gauge that is going to be used to update the groundwater reservoir
	INTEGER POSTOINT(NB) !Indicated the gauge that is going to be used to update the sub-superficial reservoir
	INTEGER POSTOSUP(NB) !Indicated the gauge that is going to be used to update the superface water reservoir
	REAL EXPOENTEBAC(NB) !Indicated the expoent used to weight the bsin area
	REAL EXPOBAS(NB) !Indicated the expoent Iused to updated the groundwater reservoir: 1 = fast 0.1 = slow
	REAL EXPOINT(NB) !Indicated the expoent Iused to updated the sub-superficial reservoir: 1 = fast 0.1 = slow
	REAL EXPOSUP(NB) !Indicated the expoent Iused to updated the surface water reservoir: 1 = fast 0.1 = slow
	REAL PLIMBAS(NB) !Limit proportion to update the groundwater reservoir
	REAL PLIMINT(NB) !Limit proportion to update the sub-superficial reservoir
	REAL PLIMSUP(NB) !Limit proportion to update the superface water reservoir
	INTEGER NDIASFSUP(NB),NDIASFINT(NB),NDIASFBAS(NB) !Number of days used in the erros estimative
	REAL SOMACAL,SOMAOBS
	REAL FCORSUP(NB),FCORINT(NB),FCORBAS(NB) !Correction factors for each basin and for each portion of the streamflow or area
	INTEGER IWW
    INTEGER J
    !REAL AFILTRO,BFIMAX !Parameters activated when using the basefow filter of ECKHARDT !Activate this one to use baseflow filter
    !REAL PFILTROBASE(0:NOBS) !Parameters activated when using the basefow filter of ECKHARDT !Activate this one to use baseflow filter
        
    POSTORIO=0
    POSTOBAS=0
    POSTOINT=0
    POSTOSUP=0
    
	!- Start the sequence of definitions of what subbasin will be update by what gauge
	
	!Indicated what is the gauge that is going to be used to correct the streamflow at each sub-basin of the model
	!POSTORIO(SUB-BASIN)=GAUGE ORDER IN THE PARHIG FILE
	POSTORIO(01)=1 !Ex.: PCH BRECHA updtates with PONTE NOVA - JUSANTE

	!Indicated what is the gauge that is going to be used to correct the groundwater at each sub-basin of the model
	!POSTOBAS(SUB-BASIN)=GAUGE ORDER IN THE PARHIG FILE
    POSTOBAS(01)=1 !Ex.: PCH BRECHA updtates with PONTE NOVA - JUSANTE
	
	!Indicated what is the gauge that is going to be used to correct the sub-superficial water at each sub-basin of the model
	!POSTOINT(SUB-BASIN)=GAUGE ORDER IN THE PARHIG FILE
	POSTOINT(01)=1 !Ex: PCH BRECHA updates with PONTE NOVA - JUSANTE

	!Indicated what is the gauge that is going to be used to correct the surface water at each sub-basin of the model
	!POSTOSUP(SUB-BASIN)=GAUGE ORDER IN THE PARHIG FILE	
	POSTOSUP(01)=1 !Ex.: PCH BRECHA ATUALIZA COM PONTE NOVA - JUSANTE
	
	!- End the sequence of definitions of what subbasin will be update by what gauge
	
	!Indicates the exponent used to weight the basin area
	EXPOENTEBAC=1.0
	!Indicates the exponent for groundwater correction; EXPOBAS = 1 updates quickly; EXPOBAS= 0.1 updates slowly;
    EXPOBAS=0.5
	!Indicates the exponent for water correction subsurface reservoir; EXPOBAS = 1 updates quickly; EXPOBAS= 0.1 updates slowly;
    EXPOINT=0.5
	!Indicates the exponent for surface water correction; EXPOBAS = 1 updates quickly; EXPOBAS= 0.1 updates slowly;
    EXPOSUP=0.8 !1.0
	!Indicates the proportion of groundwater in the total flow to get the correct underground storage
    PLIMBAS=0.5
	!Indicates the proportion of sub-surface water in the total flow to get the correct storage
    PLIMINT=0.2
	!Indicates the proportion of surface water in the total flow to get the correct storage
    PLIMSUP=0.2
    !Indicates the number of days to do the mean and calculate the error and the correction factor
    NDIASFBAS=5
    NDIASFINT=2
    NDIASFSUP=1

	NTMUSK=DTP/DT(NC) !Considers that the number of Cunge musk time intervals is the same in all cells
    
    !Correction factor from  QOBS
	FCORBAS=1.0
	FCORSUP=1.0
	FCORINT=1.0
	
	DO I=1,NOBS             !I is number of counter positions observed, but not sub-basins as expected in  NDIASFSUP(I)

        K= IBAC(IQOBS(I)) !K is the number of sub-basin (IBAC)
        
		!Calculates factor for the surface flow correction
		SOMACAL=0.0
		SOMAOBS=0.0
		
		DO IWW=1,NDIASFSUP(I) !Uses the indicated number of days
			IF(QOBS(I,IT+1-IWW)>0.0)THEN
				SOMAOBS=SOMAOBS+QOBS(I,IT+1-IWW)
				SOMACAL=SOMACAL+QR(I,IT+1-IWW)
			ENDIF
		ENDDO
		IF(SOMACAL>0.0)THEN
			FCORSUP(I)=SOMAOBS/SOMACAL
		ELSE
			FCORSUP(I)=1.0
		ENDIF

		!Calculates factor for subsurface flow correction
		SOMACAL=0.0
		SOMAOBS=0.0
		DO IWW=1,NDIASFINT(I) !Uses the indicated number of days
			IF(QOBS(I,IT+1-IWW)>0.0)THEN
				SOMAOBS=SOMAOBS+QOBS(I,IT+1-IWW)
				SOMACAL=SOMACAL+QR(I,IT+1-IWW)
			ENDIF
		ENDDO
		IF(SOMACAL>0.0)THEN
			FCORINT(I)=SOMAOBS/SOMACAL
		ELSE
			FCORINT(I)=1.0
		ENDIF

		!Calculates factor for groundwater flow correction
		SOMACAL=0.0
		SOMAOBS=0.0
		DO IWW=1,NDIASFBAS(I) !Uses the indicated number of days
			IF(QOBS(I,IT+1-IWW)>0.0)THEN
				SOMAOBS=SOMAOBS+QOBS(I,IT+1-IWW)
				SOMACAL=SOMACAL+QR(I,IT+1-IWW)
			ENDIF
		ENDDO
		IF(SOMACAL>0.0)THEN
			FCORBAS(I)=SOMAOBS/SOMACAL
		ELSE
			FCORBAS(I)=1.0
		ENDIF
		
		!ACTIVATE THE FOLLOWING LOOPS IF YOU WANT TO USE ECKHARDT BASEFLOW FILTER FOR DATA ASSIMILATION
		!FOR MORE INFORMATION, SEE: FAN, F. M. ; MELLER, A. ; COLLISCHONN, W. . Incorporação de filtro numérico de separação de escoamento na assimilação de dados para previsão de vazões utilizando modelagem hidrológica. Revista Brasileira de Recursos Hídricos, 2015.
		
!		!Calculates the post base flow using filter eckhardt
!		AFILTRO=0.999661141
!		BFIMAX=0.578
!		QBASE(I,1)=QBASE(I,2)
!		IF(IT==7939)THEN
!		    WRITE(*,*)
!		ENDIF
!		IF(QOBS(I,IT)>0.0)THEN
!		    QBASE(I,2)=((1.-BFIMAX)*AFILTRO*QBASE(I,1)+(1.-AFILTRO)*BFIMAX*QOBS(I,IT))/(1.-AFILTRO*BFIMAX)
!		    QBASE(I,2)=MIN(QBASE(I,2),QOBS(I,IT))
!		    PFILTROBASE(I)=QBASE(I,2)/QOBS(I,IT)
!		ELSE
!		    QBASE(I,2)=QBASE(I,1) !If has failed in the observed data, flow keeps constant 
!		    PFILTROBASE(I)=0.0
!		ENDIF
		
		FCORBAS(K)=MIN(FCORBAS(K),1.2)
		FCORBAS(K)=MAX(FCORBAS(K),0.8)
		FCORSUP(K)=MIN(FCORSUP(K),5.0)
		FCORSUP(K)=MAX(FCORSUP(K),0.2)
		FCORINT(K)=MIN(FCORINT(K),1.2)
		FCORINT(K)=MAX(FCORINT(K),0.8)
		
	ENDDO
	
	DO I=1,NOBS
        !Limits the maximum and the minimum of correction for not overdo the correction
	    FCORBAS(I)=MIN(FCORBAS(I),1.2)
		FCORBAS(I)=MAX(FCORBAS(I),0.8)
		FCORSUP(I)=MIN(FCORSUP(I),5.0)
		FCORSUP(I)=MAX(FCORSUP(I),0.2)
		FCORINT(I)=MIN(FCORINT(I),1.2)
		FCORINT(I)=MAX(FCORINT(I),0.8)
    END DO
		
	!Start the corrections
	DO IC=1,NC
        !WRITE(*,*)'Apdating Minibasin: ',IC
		NTRECH=NSUBT(IC)
		I=IBAC(IC)		
		
		!Corrects groundwater storage
		IF(POSTOBAS(I)>0.AND.PJB2(IC)>PLIMBAS(I))THEN !Deactivate this one to use baseflow filter
		!IF(POSTOBAS(I)>0.AND.PFILTROBASE(POSTOBAS(I))>PLIMBAS(I))THEN !Activate this one to use baseflow filter
		    !WRITE(*,*)'BACIA',I,' POSTO BAS',POSTOBAS(I)		
            !WRITE(*,*)FCORBAS(POSTOBAS(I))
			!WRITE(*,*)IT,I,IC,VBAS(IC)
			VBAS(IC)=VBAS(IC)*(FCORBAS(POSTOBAS(I))**EXPOBAS(I))*PJB2(IC)+VBAS(IC)*(1.-PJB2(IC)) !Corrects the reservoir volume 
			VBASPREV(IC)=VBASPREV(IC)*(FCORBAS(POSTOBAS(I))**EXPOBAS(I))*PJB2(IC)+VBASPREV(IC)*(1.-PJB2(IC)) !Corrects the reservoir volume 
		ENDIF
		!Corrects subsurface water
		IF(POSTOINT(I)>0.AND.PJI2(IC)>PLIMINT(I))THEN
		    !WRITE(*,*)'BACIA',I,' POSTO INT',POSTOINT(I)
		    !WRITE(*,*)FCORINT(POSTOINT(I))
			VINT(IC)=VINT(IC)*(FCORINT(POSTOINT(I))**EXPOINT(I))*PJI2(IC)+VINT(IC)*(1.-PJI2(IC)) !Corrects the reservoir volume 
			VINTPREV(IC)=VINTPREV(IC)*(FCORINT(POSTOINT(I))**EXPOINT(I))*PJI2(IC)+VINTPREV(IC)*(1.-PJI2(IC)) !Corrects the reservoir volume 
		ENDIF
		!Corrects surface water
		IF(POSTOSUP(I)>0.AND.PJS2(IC)>PLIMSUP(I))THEN
		!IF(POSTOSUP(I)>0.AND.PFILTROBASE(POSTOBAS(I))<0.6)THEN
		    !WRITE(*,*)'BACIA',I,' POSTO SUP',POSTOSUP(I)
            !WRITE(*,*)FCORSUP(POSTOSUP(I))        		
			VSUP(IC)=VSUP(IC)*(FCORSUP(POSTOSUP(I))**EXPOSUP(I))*PJS2(IC)+VSUP(IC)*(1.-PJS2(IC)) !Corrects the reservoir volume 
			VSUPPREV(IC)=VSUPPREV(IC)*(FCORSUP(POSTOSUP(I))**EXPOSUP(I))*PJS2(IC)+VSUPPREV(IC)*(1.-PJS2(IC)) !Corrects the reservoir volume 			
		ENDIF

		VSUP(IC)=MAX(VSUP(IC),0.0)
		VINT(IC)=MAX(VINT(IC),0.0)
		VBAS(IC)=MAX(VBAS(IC),0.0)
		
		VSUPPREV(IC)=MAX(VSUPPREV(IC),0.0)
		VINTPREV(IC)=MAX(VINTPREV(IC),0.0)
		VBASPREV(IC)=MAX(VBASPREV(IC),0.0)

		IF(POSTORIO(I)>0.AND.OD(IC)>1)THEN  
		    !WRITE(*,*)'BACIA',I,' POSTO SUP',POSTORIO(I)
            !WRITE(*,*)FCORSUP(POSTORIO(I))        		
		    K=IQOBS(POSTORIO(I))
    		PESOBAC=(ACUR(IC)/ACUR(K))**EXPOENTEBAC(I) !Weight for correction in other cells
			QJ2PREV(IC)=(QJ2(IC)*FCORSUP(POSTORIO(I)))*PESOBAC+QJ2(IC)*(1.-PESOBAC)
			QJ2(IC)=QJ2PREV(IC)
			DO ITR=1,NTRECH+1 !Updated initial conditions of Muskingum-Cunge
				QRIOINI(IC,ITR)=QRIOINI(IC,ITR)*FCORSUP(POSTORIO(I))*PESOBAC+QRIOINI(IC,ITR)*(1.-PESOBAC)
				QRIOINIPREV(IC,ITR)=QRIOINI(IC,ITR)
			ENDDO
			DO LT=1,NTMUSK+1
				QCONTORM(IC,LT)=(QCONTORM(IC,LT)*FCORSUP(POSTORIO(I)))*PESOBAC+QCONTORM(IC,LT)*(1.-PESOBAC)
			ENDDO
		ENDIF
		
	ENDDO

    !WRITE(*,*)'End of Subroutine ATUALIZA'
	RETURN
	END
