    !*********************************************************************************
    !
    !  SUBROUTINE MODELO controls the main time loop in MGB-IPH and call routines
	!					that realize catchment and river routing
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine has the main loop of the MGB-IPH, the time loop, from iT=1,nT
    !     where iT is the time interval and nT is the number of time steps.
	!
	!	 For each time interval date info is set, then rainfall and climate data are
	!	   loaded through subroutines LECHUVA and LECLIMA. Afterwards catchment flow
	!	   generation is done using subroutine CELULA. Finally, river routing is 
	!	   achieved by calling (i) REDE, for Muskingum-Cunge Method or
	!	   (ii) flood_inertial, for a 1D Hydrodynamic Model (w/ inertia and pressure)
	!	
	!	At the end discharge time series are stored in :
	!		QRB: calculated discharge time series in all subbasins
	!		QRG: calculated discharge time series in catchments pre-defined by user
	!		QR:  calculated discharge time series in subbasin outlet w observed data
	!	Those are recorded in files at the end when returns to SIMULA
	!
	!
	!	 * iT value should not be changed inside subroutines!
	!
	!    * AUX_MOD module from full hydrodynamic is deactivated (commented)!
	!    * Hidrodinamico2 subroutine from full hydrodynamic is deactivated (commented)!
    !
    !  Usage:
    !
    !    CALL MODELO
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
	!	 * function CALDAT	      		in	  caldat.f90
    !    * module     VARS_MAIN   		in    VARS_MAIN.f90
    !    * module     VARS_INERC  		in    VARSINERC.f90  !Quais? 
 	!	 * subroutine LECHUVA	  		in	  LECHUVA.f90
	!	 * subroutine LECLIMA	  		in	  LECLIMA.f90
	!	 * subroutine CELULA	  		in	  CELULA.f90
	!	 * subroutine REDE	      		in	  REDE.f90	
	!	 * subroutine flood_inercial	in	  flood_inercial.f90
	!
	!    Deactivated:	
	!	 * module	  AUX_MOD        	in	  HD_AUX_MOD.f90
	!	 * subroutine Hidrodinamico2	in	  Hidrodinamico2.f90
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
    !   * Variables declarations and routines calls are all commented below.
    !
    !---------------------------------------------------------------------------------
	SUBROUTINE MODELO

	USE VARS_MAIN
	USE VARS_INERC
	!USE AUX_MOD
	IMPLICIT NONE
	INTEGER K,KHID 					!indexes and counters
	INTEGER KC,JB,JC,KB,JCB,MWP 	!... same
    INTEGER count_daniel,dt_daniel
	INTEGER iTwrite


	! Initialize
	IT=0
    Dt_daniel=INT((NT-IT)/10)

	!WRITE(*,*)' Do you want to avoid simulation in any sub-basin? (0) No, (>=1) How many?
	!READ(*,*)NCONGEL
	!WRITE(*,*)' Which one? List sub-basin numbers with space separator.'
	!READ(*,*)(IBCONGEL(KB),KB=1,NCONGEL)

	!Write(*,*)'Percentage of simulation complete'
    !Write(*,*)'10------50--------100'
 
    !CALCULA PREC E ET MENSAL PARA INLAND DELTA
    PRECTUDO_NIGER_MES=0
    EVAPTUDO_NIGER_MES=0
    CONT_MES_NIGER=1
    PRECTUDO_NIGER_ANO=0
    EVAPTUDO_NIGER_ANO=0
    CONT_ANO_NIGER=1
    
    STUDO_DELTA=0.0 !variável do armazenamento de água no solo para todos IC e IT
    
    
     ! Time Loop     
    DO WHILE (IT<NT)
		
		IT=IT+1
		!write(*,*)iT
		
       if(mod(it,100)==0)  write(*,*)100*IT/NT,' %'
          
!		! Save individual records of Discharge, water level and depth for all catchments
		! only at first time-step ?!

!!ATTENTION!! TURNED OFF!!! FMF

!!		if(iT==1)then
!!		    call flood_write
!!		endif
!!		
		if (icalib==0.and.itWrite<iT) then
			!WRITE(*,*) IT,maxval(QJ2),maxloc(QJ2),minval(QJ2),minloc(QJ2),QJ2(nC) !QJ2(6848)! RP
!			write(7000,*) 'iT',iT
!			pause
			itWrite=itWrite+1
        endif
        
		if(it==count_daniel)then
	!	   WRITE(*,701)'**'
			701		FORMAT(A2,$)
			count_daniel=Count_daniel+Dt_daniel
        endif


		JDIA=IDINI+INT((IT+HORAINI-1)/(86400./DTP)) 	! Gets calendar day (big number)
		CALL CALDAT(JDIA,IMES,IDIA,IANO)

		DIAH(IT)=IDIA !stores day corresponding to iT time interval
		MESH(IT)=IMES !stores month corresponding to iT time interval
		ANOH(IT)=IANO !stores year corresponding to iT time interval
		HORAH(IT)=MOD(IT,24)+HORAINI 	!hour of day corresponding to iT time interval

        !CALCULA ET E PREC, DIARIOS E MENSAIS, PARA INLAND DELTA
        !WHOLE DELTA
        EVAPTUDO_NIGER(IT)=0 !inicializa vetor de evaporação no inland delta como zero
        EVAPTUDO_NIGER_INT(IT)=0 !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        EVAPTUDO_NIGER_SOLO(IT)=0  !ET FOR INLAND DELTA NIGER ON SOIL
        EVAPTUDO_NIGER_FLOOD(IT)=0  !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        PRECTUDO_NIGER(IT)=0 !PREC FOR INLAND DELTA NIGER 
        
        !SOUTHERN DELTA
        S_EVAPTUDO_NIGER(IT)=0 !inicializa vetor de evaporação no inland delta como zero
        S_EVAPTUDO_NIGER_INT(IT)=0 !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        S_EVAPTUDO_NIGER_SOLO(IT)=0  !ET FOR INLAND DELTA NIGER ON SOIL
        S_EVAPTUDO_NIGER_FLOOD(IT)=0  !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        S_PRECTUDO_NIGER(IT)=0 !PREC FOR INLAND DELTA NIGER
        
        !NORTHERN DELTA
        N_EVAPTUDO_NIGER(IT)=0 !inicializa vetor de evaporação no inland delta como zero
        N_EVAPTUDO_NIGER_INT(IT)=0 !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        N_EVAPTUDO_NIGER_SOLO(IT)=0  !ET FOR INLAND DELTA NIGER ON SOIL
        N_EVAPTUDO_NIGER_FLOOD(IT)=0  !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        N_PRECTUDO_NIGER(IT)=0 !PREC FOR INLAND DELTA NIGER
        
        IF (IT>1) THEN
        IF(MESH(IT).NE.MESH(IT-1)) CONT_MES_NIGER=CONT_MES_NIGER+1 !CONTADOR DE NUMERO DE MESES PARA PREC_NIGER_MES
        IF(ANOH(IT).NE.ANOH(IT-1)) CONT_ANO_NIGER=CONT_ANO_NIGER+1 !CONTADOR DE NUMERO DE MESES PARA PREC_NIGER_MES
        ENDIF
        
		! Read rainfall data for all catchments in iT time-step
		ITCHUVA=IT
		CALL LECHUVA
        
	    
		DO KC=1,NC
			PM(KC)=PM(KC)+P(KC) !Cumulative rainfall (used later for average)
		ENDDO		
	
		! Reads climate data
		TONTEM=TA
		CALL LECLIMA
	
		! Calls catchment routing/discharge for lateral inflow
        DINFILT=0.0
		CALL CELULA		
	
		! Saves detailed info for soil state variables in file NOSOLO.HIG.
		! Uses JC index for catchment and JB index for URH of interest.
		IF(ICALIB.EQ.0)THEN !only if not calibrating.
			JB=1
			JC=1
			!JC=57
			WRITE(FILSOL,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)
			JB=2
			JC=2
			WRITE(FILSOL2,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)			

!			write(971,66) (E0agua(iC),iC=1,nC)
!			write(972,66) (E0topo(iC),iC=1,nC)
!			write(973,66) (E0sup(iC),iC=1,nC)

		ENDIF

	
		! Call main river routing routine using Muskingum-Cunge Method
		if(hdFLAG0==0)then
		    CALL REDE
		endif
		
		! Calculates lateral inflow 
        do ic=1,nc
		    QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC) !sums surface, subsurface and 
		enddo

        ! Calls river routing routine using Inertial Method
        if(hdFLAG0>0)then
            CALL flood_inercial
        endif        

		! Stores calculated discharges in sub-basins with observed data - for calibration and assessment.
		!DO K=1,NOBS
		!	KHID=IQOBS(K) 		! outlet catchment id
		!	QR(K,IT)=QJ2(KHID)  ! saves on QR(K,IT), that will be compared to QOBS in FOBJ routine.
  !      ENDDO
        !do ic=1,nc
        !    QTUDO(IC,IT)=QJ2(IC) !stores discharge at all time steps 
            QTUDO(:,IT)=QJ2 !stores discharge at all time steps 
        !    QITUDO(IC,IT)=QCEL2(IC) !stores discharge at all time steps 
        !    YTUDO(IC,IT)=Yfl(IC) !stores water level at all time steps 
            YTUDO(:,IT)=Yfl !stores water level at all time steps 
            !HTUDO(IC,IT)=Hfl(IC) !stores water depth at all time steps 
            HTUDO(:,IT)=Hfl !stores water depth at all time steps 
        !    QVIZTUDO(IC,IT)=Q2viz(IC) !stores lateral connection discharge at all time steps 
        !enddo

		! Stores discharges for file ouput when in simulation model
		IF(ICALIB.EQ.0)THEN 	! only if it is not calibration
			
			DO KB=1,NB				! store discharge in sub-basin
				JCB=IEXUT(KB) 		! outlet catchment of KB-th subbasin
				QRB(KB,IT)=QJ2(JCB)
			ENDDO
	
			DO K=1,NUMHIDG 				! store discharge in catchments defined in PARHIG.HIG
				KHID=IHIDG(K) 			! indexed catchment
				QRG(K,IT)=QJ2(KHID) 	! stored in QRG
				
			if(hdFLAG0>0)then !Stores the water level for inertial model 	           
                HRG(K,IT)=Hfl(KHID) ! stored in HRG                  
                !HRG(K,IT)=Yfl(KHID) ! stored in HRG                  
            endif     
				
			ENDDO
	
			! Store discharge by water origin (i.e. surface, subsurface, groundwater) in a specified catchment
			MWP=1 						!catchment (this is manual)
			QB(IT)=QJ2(MWP)*PJB2(MWP)
			QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
			QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)
		ENDIF
        
	
	ENDDO !End time loop
	

	

75	FORMAT(I6,8F10.4)

66	FORMAT(<nC>F10.4)

     
	RETURN
	END
