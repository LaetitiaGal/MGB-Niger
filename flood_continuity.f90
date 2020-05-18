	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the level and depth for each catchment from Continuity equation   (Inertial version).
    !
    !
    ! Usage:
    !
    ! *
    !
    ! uses modules, functions, and subroutines
    !
    ! * USE VARS_MAIN
    ! * USE VARS_INERC (only to Inertial version)
    !
    ! opens
    !
    ! * no files are created in this routine
    !
    ! reads
    !
    ! * no files are created in this routine
    !
    ! creates
    !
    ! * no files are created in this routine
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. - VER ISSO.
    !
    !  Version/Modified: 
    !
    ! 2015.07.06 - 07 July 2015 (by Paulo Pontes)
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
    ! Variables and Parameters:
    ! *Variables declarations and routines calls are all commented below.
    !---------------------------------------------------------------------------------
    ! End of header
    !---------------------------------------------------------------------------------

	subroutine flood_continuity

	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use VARS_INERC
	use VARS_MAIN
	implicit none
	
	real*8 areajtab, yjtab
    
    !$OMP PARALLEL NUM_THREADS(8)
    !$OMP DO PRIVATE(Nentradas,SumQup,Areajtab,yjtab)
    !-------------------------------------------------------------------------------------
    do iC=1,nC
    
        IB=IBAC(IC) 
		!IF((IB<SUBini).OR.(IB>SUBfim))CYCLE ! Manual controls for especialized calibration
    
        !SIMULATE ONLY RED FLOODS
        !IF(IBAC(IC)<5.AND.IBAC(IC)>7) cycle
        !IF(IBAC(IC).NE.4) cycle
        !IF(IBAC(IC)>4) cycle
        !IF((IBAC(IC)>4).AND.(IBAC(IC).NE.9)) cycle
        !IF(IBAC(IC)>9) cycle
        !IF(IBAC(IC).NE.9) cycle
        !IF(IBAC(IC).NE.5) cycle
        !IF(IBAC(IC)<10.OR.IBAC(IC)>14) cycle   !modif cécile
!        IF((IBAC(IC)<10.OR.IBAC(IC)>14).AND.(IBAC(IC).NE.20)) cycle   !modif cécile

        
            !Number of upstream catchments of IC
            Nentradas = MINIMONT(iC,1)
            
            !Sum of downstream IC flows
            if(Nentradas==0)then
                SumQup=0.0
            else
                SumQup = SUM(Q2fl(MINIMONT(iC,2:1+Nentradas)))
            endif
            
            !>>>>>>>>>>>>>>>>>>>>>FOR FLOW SUBSTITUTION (E.G. ANSONGO)
            !if(iC==4278) then
            if(iC==ISUBST(1)) then
             Q2fl(iC)=QLIDO(1,IT)
                Vol2(iC) = Vol1(iC)+dtflood*(QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
            else
                Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            endif
            
            
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>WORKING LINE:            
                !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END WORKING LINE
            
            
            !>>>>>>>>>>>>OFFICE DU NIGER
            !if (inland_delta(ic)==1.and.mesh(it)<6) then
                !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((1*E0agua(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !if (ic==4127) then !Office du Niger     
            !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC)-100)-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !else
            !if (Hfl(IC)<=0.1) E0agua(IC)=0     
            !E0agua(IC)=0     
                !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))
            !else
                !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((E0agua(IC)+DINFILT(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !    endif
           ! elseif (ic==4127) then !Office du Niger 
            !    Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC)-100)-((E0agua(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !else
            !    Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-((1*E0agua(IC))*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))    
            !endif
            !>>>>>>>>>>>>>>>>>>>
            
            Vol2(iC) = max(Vol2(iC),0.0)
            
            !Interpolates the Area and Level from Volume
            CALL hunt(VTAB(:,iC),Vol2(iC),jtab(iC),ATAB(:,iC),Areajtab,ZTAB(:,IC),yjtab,NPFL(IC)+2)
            Area2(iC) = min(Areajtab,ACEL(iC))
            
            !Updates variables:
            y2_fl=max(yjtab,ZTAB(1,IC))
            y2_fl=y2_fl+0.001
            
            !Calculates depth
            Hfl(iC)=y2_fl-ZTAB(1,IC)
            Yfl(iC)=y2_fl       
            
            !Updates the volume
            Vol1(iC)=Vol2(iC)      
            
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR WHOLE INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IF (INLAND_DELTA(IC)==1) THEN ! calcula ET do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)  
        EVAPTUDO_NIGER(IT)=EVAPTUDO_NIGER(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        EVAPTUDO_NIGER_FLOOD(IT)=EVAPTUDO_NIGER_FLOOD(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) 
        PRECTUDO_NIGER(IT)=PRECTUDO_NIGER(IT)+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
        EVAPTUDO_NIGER_MES(CONT_MES_NIGER)=EVAPTUDO_NIGER_MES(CONT_MES_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)=EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
    ENDIF

    IF (INLAND_DELTA(IC)==1.AND.tflood/dtflood0==1) then ! calcula área inundada total do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)   
       AFLTUDO(1,IT)=AFLTUDO(1,IT)+Area2(IC) !Em km²
    ENDIF

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR SOUTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF (SOUTHERN_DELTA(IC)==1) THEN ! calcula ET do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)  
        S_EVAPTUDO_NIGER(IT)=S_EVAPTUDO_NIGER(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        S_EVAPTUDO_NIGER_FLOOD(IT)=S_EVAPTUDO_NIGER_FLOOD(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))     
        S_PRECTUDO_NIGER(IT)=S_PRECTUDO_NIGER(IT)+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
        S_EVAPTUDO_NIGER_MES(CONT_MES_NIGER)=S_EVAPTUDO_NIGER_MES(CONT_MES_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        S_EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)=S_EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
    ENDIF

    IF (SOUTHERN_DELTA(IC)==1.AND.tflood/dtflood0==1) then ! calcula área inundada total do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)   
       S_AFLTUDO(1,IT)=S_AFLTUDO(1,IT)+Area2(IC) !Em km²
    ENDIF
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR NORTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (NORTHERN_DELTA(IC)==1) THEN ! calcula ET do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)  
        N_EVAPTUDO_NIGER(IT)=N_EVAPTUDO_NIGER(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        N_EVAPTUDO_NIGER_FLOOD(IT)=N_EVAPTUDO_NIGER_FLOOD(IT)+(1*E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))     
        N_PRECTUDO_NIGER(IT)=N_PRECTUDO_NIGER(IT)+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
        N_EVAPTUDO_NIGER_MES(CONT_MES_NIGER)=N_EVAPTUDO_NIGER_MES(CONT_MES_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
        N_EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)=N_EVAPTUDO_NIGER_ANO(CONT_ANO_NIGER)+(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.)) !em m³
    ENDIF

    IF (NORTHERN_DELTA(IC)==1.AND.tflood/dtflood0==1) then ! calcula área inundada total do INLAND DELTA (FLAG INLAND_DELTA(IC)=1)   
       N_AFLTUDO(1,IT)=N_AFLTUDO(1,IT)+Area2(IC) !Em km²
    ENDIF

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR NORTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    
    enddo
    
    !$OMP END DO
    !$OMP END PARALLEL
	
	endsubroutine