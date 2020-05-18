    !*********************************************************************************
    !
    !  SUBROUTINE SIMULA is the main routine for hydrological simulation
	!                    (wo calibration)
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
    !    CALL SIMULA
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
	!	 * module     IFPORT  			from (visual studio 2008)
    !    * module     VARS_MAIN   		in    VARS_MAIN.f90
	!    * module     VARS_CALIB   		in    VARS_CALIB.f90
    !    * module     VARS_INERC  		in    VARSINERC.f90  !Quais? 
 	!	 * subroutine CONDINIC	  		in	  CONDINIC.f90
	!	 * subroutine MODELO	  		in	  MODELO.f90
	!	 * subroutine FOBJ	  			in	  FOBJ.f90	
	!
    !
    !	 opens
    !
    !    * Opens QPROP.PRP  	   output  ascii  file with calculated discharge by origin (i.e. surface, subsurface, gw)
	!	 * Opens VAZAO.QCL 		   output  ascii  file with calculated discharge in subbasins defined in setup file
	!	 * Opens QTUDO.QBI 		   output  binary file with calculated discharge in all catchments
	!	 * Opens Qmes90.TXT 	   output  ascii  file with calculated Q90 in all catchments
	!	 * Opens Qmes95.TXT 	   output  ascii  file with calculated Q95 in all catchments
	!	 * Opens RESUMO_SIAQUA.TXT output  ascii  file for SIAQUA	
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
    ! 2015.06.21 - 21 June 2015 (By: Fernando Mainardi Fan)   
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

	SUBROUTINE SIMULA	
	USE VARS_MAIN
	USE VARS_CALIB
	USE VARS_INERC
	USE IFPORT

	IMPLICIT NONE
	INTEGER KB,K					!counters
	CHARACTER*10 ALFAQ(NC)			!auxiliar string
	REAL,ALLOCATABLE:: QSOR(:) 		!observed discharges (without "missing data") for sorting and exceedance probaility curve
	CHARACTER (10):: strIHIDG
	CHARACTER (50):: strARQUIVO	
	INTEGER MES,KSOR,IPOS90,IPOS95 				  !indexes and counter
	REAL QMES90(NC,12),QMES95(NC,12),QMLT(NC,12)  !Q90, Q95 and Long-term Mean Discharge for each month in each catchment
	INTEGER CONTAMES(12)						  !counts number of existent data for mean calculation
	INTEGER IAUX,IPOS,IPERMA
	REAL QPERMA(NC,100) 			! Discharge values in Exceedance Probability Curve for each catchment
	REAL HPERMA(NC,100) 			! Water Depth values in Exceedance Probability Curve for each catchment
	REAL VPERMA(NC,100) 			! Velocity values in Exceedance Probability Curve for each catchment
	!REAL RUGMAN						! Manning bed roughness coeficient !RUGMAN IS DECLARED IN VARS_MAIN AND ALLOCA_VARS SUBROUTINES	!FMF 09/09/2015 
    INTEGER ITini
    REAL,ALLOCATABLE:: Qmin(:),Qmax(:),QmaxTR(:,:),Qmax_aux(:) 
    INTEGER,ALLOCATABLE:: ANOaux(:)  ! RP abril/2013
    INTEGER:: i,j
    REAL:: gumbel_v,gumbel_b,gumbel_med,gumbel_var,gumbel_T
    REAL AREA_SUBBACIA !AREA TOTAL DA SUB-BACIA DO INLAND DELTA
	REAL std,Qmean
	!Manning bed roughness coeficient
    !RUGMAN=0.03 !RUGMAN can be set to 0,03 in all rivers if desired !FMF 09/09/2015   
	
	!Activate all basins to simulate
    !ATIVABACIA=1 

	! Main Routines
	WRITE(FILLOG,*)'Calling initial conditions...'
	CALL CONDINIC  ! Initial Conditions
	WRITE(FILLOG,*)'Running main model routine...'
	CALL MODELO    ! Time Loop
	WRITE(FILLOG,*)'Calculating objective functions...'
	!CALL FOBJ 	   ! Error/objective Functions Evaluation

	! Outputs
	WRITE(FILLOG,*)'Writting outputs...'
	OPEN(FILPRP,FILE=OUTPUT_DIRECTORY // 'QPROP.PRP',STATUS='UNKNOWN')
	OPEN(FILHID,FILE=OUTPUT_DIRECTORY // 'VAZAO.QCL',STATUS='UNKNOWN')
	OPEN(FILTUD,FILE=OUTPUT_DIRECTORY // 'QTUDO.QBI',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')
	OPEN(1960,FILE=OUTPUT_DIRECTORY // 'Qmes90.TXT',STATUS='UNKNOWN',ACTION='WRITE')
	OPEN(1961,FILE=OUTPUT_DIRECTORY // 'Qmes95.TXT',STATUS='UNKNOWN',ACTION='WRITE')
	OPEN(1971,FILE=OUTPUT_DIRECTORY // 'RESUMO_SIAQUA.TXT',STATUS='UNKNOWN',ACTION='WRITE')
	OPEN(1972,FILE=OUTPUT_DIRECTORY // 'QITUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') 
    OPEN(1979,FILE=OUTPUT_DIRECTORY // 'YTUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')   !modif Cécile - bug 1972 ?
!    OPEN(1972,FILE=OUTPUT_DIRECTORY // 'YTUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')   !modif Cécile
	OPEN(1973,FILE=OUTPUT_DIRECTORY // 'QBASTUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') 
	OPEN(1974,FILE=OUTPUT_DIRECTORY // 'EVAPTUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')
	OPEN(1975,FILE=OUTPUT_DIRECTORY // 'WTUDO_bin.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')

    ! Record header in RESUMO_SIAQUA.TXT
    write(1971,'(A240)')'      MINI       Q5       Q10       Q30       Q50       Q70       Q90       Q95       V10       V30       V50       V70       V90   LARGURA     Qmean     Qmax     Qmin    Qmax05ano    Qmax10ano     Qmax25ano' ! RP abril/2013
!
!	! Store string for catchment record (?!) - na real nao sei bem onde usa isso
!	DO IC=1,NC
!		IF(IC<NC)THEN
!			WRITE(ALFAQ(IC),'(I9,A1)')IC,';'
!		ELSE 	!if it is the last column, doesnt rec ";"
!			WRITE(ALFAQ(IC),'(I9)')IC
!		ENDIF
!	ENDDO
!
!	!WRITE(FILTUD,'(A11,<NC>A10)')'DD;MM; ANO;',(ALFAQ(IC),IC=1,NC)
!	
!	!-------------------------------------------------------------
!	! Record outputs in files
!	
	DO IT=1,NT
		! record discharge time series in nb predefined catchments outlets 
		write(filhid,71)it,(qrg(kb,it),kb=1,numhidg)
		
		! record special characteristics for a defined (below) sub-basin
		ib=1 		!subbasin
		write(filbac,'(2f10.2,8f10.3)')wbm(ib,it),pm2(ib,it),sim(ib,it),em(ib,it),dsum(ib,it),dinm(ib,it),dbam(ib,it)	! water balance components time series 
		write(filprp,71)it,qb(it),qbi(it),qbis(it)																		! discharge by origin (i.e. surface, subs, gw) time series
		
		! record binary file with discharge time series for all catchments
		write(filtud,rec=it)(qtudo(ic,it),ic=1,nc)
		write(1972,rec=it)(qitudo(ic,it),ic=1,nc)
        
       ! write(1972,rec=it)(htudo(ic,it),ic=1,nc) ! modif Cécile - lignes commentées 
       ! write(1972,rec=it)(ytudo(ic,it)-zfundofl(ic),ic=1,nc) : 
        
       write(1979,rec=it)(ytudo(ic,it)-zfundofl(ic),ic=1,nc)  ! cécile test 1 (écriture de cette ligne avec changement de num de fichier)
        !write(1979,rec=it)(htudo(ic,it),ic=1,nc) ! Cécile test2 (écriture de cette ligne avec changement de num de fichier)

        write(1973,rec=it)(qbastudo(ic,it),ic=1,nc)
		write(1974,rec=it)(evaptudo(ic,it),ic=1,nc)
		write(1975,rec=it)(wtudo(ic,it),ic=1,nc)

!!!ATENÇÃO, VERIFICAR AS LINHAS DE CÓDIGO ABAIXO COM A GALERA
!!!		DO IC=1,NC
!!!			exit	! ? what?
!!!
!!!			IF(IC<NC)THEN
!!!				WRITE(ALFAQ(IC),'(F9.2,A1)')QTUDO(IC,IT),';'
!!!			ELSE 
!!!				WRITE(ALFAQ(IC),'(F9.2)')QTUDO(IC,IT)
!!!			ENDIF
!!!			
!!!		ENDDO
!
!		!WRITE(*,'(10A10)')(ALFAQ(IC),IC=1,10)
!		JDIA=IDINI+IT-1 !VERIFICA QUAL É O DIA DO CALENDÁRIO 
!		CALL CALDAT(JDIA,IMES,IDIA,IANO)
!		!WRITE(FILTUD,'(I2,A1,I2,A1,I4,A1,<NC>A10)')IDIA,';',IMES,';',IANO,';',(ALFAQ(IC),IC=1,NC)
	ENDDO
	!
	!!Close file units
	!CLOSE (FILBAC)
	CLOSE (FILHID)
	CLOSE (FILPRP)
	CLOSE (FILTUD)
	CLOSE (1972)
	CLOSE (1973)
	CLOSE (1974)
	CLOSE (1975)
    
	!
	!-------------------------------------------------------------
	! Calculates mean annual precipitation by basin (long-term)
	PMB=0.0
	KPM=0
	PM=(PM/NT)*(365.*24.*3600./DTP)

	DO IC=1,NC
		IB=IBAC(IC)
		PMB(IB)=PMB(IB)+PM(IC)
		KPM(IB)=KPM(IB)+1
	ENDDO
	DO IB=1,NB
		PMB(IB)=PMB(IB)/KPM(IB)
		WRITE(*,*)'Sub-basin ',IB,'   Annual Rainfall',PMB(IB)
		WRITE(FILLOG,*)'Sub-basin ',IB,'   Annual Rainfall',PMB(IB)
	ENDDO
	
	!-------------------------------------------------------------
	! Record Objective-Function for each station
	!WRITE(*,'(A9,3A8)')'Sub-basin','Nash','Nashlog','Evol' !SCREEN
	!DO K=1,NOBS
	!	WRITE(FILAJU,'(3F10.3)')R2(K),R2L(K),ERRV(K) !FILE
 !       WRITE(*,'(I5,3F10.3)')K,R2(K),R2L(K),ERRV(K) !SCREEN
	!ENDDO
	!CLOSE (FILAJ)
 !   !-------------------------------------------------------------
 !   ! Calculates mean, min and max discharge
 !   allocate(Qmean(nC),Qmin(nC),Qmax(nC))
 !   ITini= 1 !367 !Should be 365 days due to intial conditions	
	!KSOR=NT-ITini+1 ! RP abril/2013
 !   Qmean=sum(QTUDO(:,ITini:NT),DIM=2)/KSOR
 !   Qmin=minval(QTUDO(:,ITini:NT),DIM=2)
 !   Qmax=maxval(QTUDO(:,ITini:NT),DIM=2)
 !   
				
  !  ! Calculates long-term mean discharge for each month (climatology)
  !  QMLT=0.0
  !  CONTAMES=0
  !  DO IT=1,NT
		!JDIA=IDINI+IT-1
		!CALL CALDAT(JDIA,IMES,IDIA,IANO)
  !      CONTAMES(IMES)=CONTAMES(IMES)+1        
  !      DO IC=1,NC
  !          QMLT(IC,IMES)=QMLT(IC,IMES)+QTUDO(IC,IT)
  !      ENDDO
  !  ENDDO
  !  DO IMES=1,12
  !      DO IC=1,NC
  !          QMLT(IC,IMES)=QMLT(IC,IMES)/CONTAMES(IMES)
  !      ENDDO
  !  ENDDO
  ! 
    
    ! Open files and record individual catchment discharge time series
    DO K=1,NUMHIDG
        WRITE(strIHIDG,'(I10)')IHIDG(K)
        !strARQUIVO = 'SIM_'//TRIM(adjustl((strIHIDG)))//'.TXT'
        if(hdFLAG0>1)then 
            !strARQUIVO = 'SIM_INERC_'//TRIM(adjustl((strIHIDG)))//'_teste2.TXT'
            strARQUIVO = 'SIM_INERC_'//TRIM(adjustl((strIHIDG)))//'.TXT'
        elseif(hdFLAG0==0)then
            strARQUIVO = 'SIM_MC_'//TRIM(adjustl((strIHIDG)))//'.TXT'
        endif
        !OPEN(801,FILE=OUTPUT_DIRECTORY // ''//'SIM_'//ARQHID(K),STATUS='UNKNOWN',ACTION='WRITE')
        OPEN(801,FILE=OUTPUT_DIRECTORY // ''//strARQUIVO,STATUS='UNKNOWN',ACTION='WRITE')
            !WRITE(801,*)'   DIA   MES   ANO  HORA        Q_(m3/s)'
            DO IT=1,NT
                WRITE(801,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),QRG(K,IT)
            ENDDO
        CLOSE(801)
    ENDDO
    
    if(hdFLAG0>1)then 
    ! Open files and record individual catchment water level time series
    ITini= 500 !367 !Should be 365 days due to intial conditions	
	KSOR=NT-ITini+1 ! RP abril/2013
    DO K=1,NUMHIDG
        Qmean=sum(HRG(K,ITini:NT))/KSOR
        !std = sqrt((sum(HRG(K,ITini:NT)**2)-sum(HRG(K,ITini:NT))**2/size(HRG(K,ITini:NT)))/(size(HRG(K,ITini:NT))-1))     
        WRITE(strIHIDG,'(I10)')IHIDG(K)
        strARQUIVO = 'SIM_INERC_HflAnomaly_'//TRIM(adjustl((strIHIDG)))//'.TXT'
        OPEN(801,FILE=OUTPUT_DIRECTORY // ''//strARQUIVO,STATUS='UNKNOWN',ACTION='WRITE')
            DO IT=1,NT
                !WRITE(801,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),(HRG(K,IT)-Qmean)/std
                WRITE(801,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),(HRG(K,IT)-Qmean)
            ENDDO
        CLOSE(801)
    ENDDO
    endif
    
     if(hdFLAG0>1)then 
    ! Open files and record individual catchment water level time series
    DO K=1,NUMHIDG
        WRITE(strIHIDG,'(I10)')IHIDG(K)
        strARQUIVO = 'SIM_INERC_Hfl_'//TRIM(adjustl((strIHIDG)))//'.TXT'
        OPEN(801,FILE=OUTPUT_DIRECTORY // ''//strARQUIVO,STATUS='UNKNOWN',ACTION='WRITE')
            DO IT=1,NT
                WRITE(801,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),HRG(K,IT)
                !write(*,*)Qmean,std,HRG(K,IT),KSOR
            ENDDO
        CLOSE(801)
    ENDDO
     endif
    !
    !if(hdFLAG0>1)then 
    !! Open files and record individual connection discharges
    !DO K=1,NUMHIDG
    !    WRITE(strIHIDG,'(I10)')IHIDG(K)
    !    strARQUIVO = 'SIM_INERC_Connection_'//TRIM(adjustl((strIHIDG)))//'.TXT'
    !    OPEN(801,FILE=OUTPUT_DIRECTORY // ''//strARQUIVO,STATUS='UNKNOWN',ACTION='WRITE')
    !        DO IT=1,NT
    !            WRITE(801,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),QRG_viz(K,IT)
    !        ENDDO
    !    CLOSE(801)
    !ENDDO
    !endif

     
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR WHOLE INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !inland delta flooded area at each IT
    If(hdFLAG0>1)then 
    OPEN(19751,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_FA.MGB',STATUS='UNKNOWN')
    DO IT=1,NT
		! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19751,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),AFLTUDO(1,IT)
    ENDDO
    CLOSE(19751)
    ENDIF
    
    !maximum absolute value of discharge, water level, water depth and lateral connection for each catchment
    OPEN(19760,FILE=OUTPUT_DIRECTORY // 'MAX LEVEL DISCHARGE.MGB',STATUS='UNKNOWN')
    DO IC=1,NC
        if (abs(maxval(QTUDO(IC,500:NT)))>abs(minval(QTUDO(IC,500:NT)))) then
            if (abs(maxval(QVIZTUDO(IC,500:NT)))>abs(minval(QVIZTUDO(IC,500:NT)))) then
            WRITE(19760,'(I5,4F16.6)') IC,MAXval(QTUDO(IC,500:NT)), MAXval(QVIZTUDO(IC,500:NT)), MAXval(HTUDO(IC,500:NT))-HRIO(IC), MAXval(YTUDO(IC,500:NT))
            else
            WRITE(19760,'(I5,4F16.6)') IC,MAXval(QTUDO(IC,500:NT)), MINval(QVIZTUDO(IC,500:NT)), MAXval(HTUDO(IC,500:NT))-HRIO(IC), MAXval(YTUDO(IC,500:NT))        
            endif
        else
            if (abs(maxval(QVIZTUDO(IC,500:NT)))>abs(minval(QVIZTUDO(IC,500:NT)))) then
            WRITE(19760,'(I5,4F16.6)') IC,MINval(QTUDO(IC,500:NT)), MAXval(QVIZTUDO(IC,500:NT)), MAXval(HTUDO(IC,500:NT))-HRIO(IC), MAXval(YTUDO(IC,500:NT))
            else
            WRITE(19760,'(I5,4F16.6)') IC,MINval(QTUDO(IC,500:NT)), MINval(QVIZTUDO(IC,500:NT)), MAXval(HTUDO(IC,500:NT))-HRIO(IC), MAXval(YTUDO(IC,500:NT))        
            endif
        endif
    ENDDO
    CLOSE(19760)
    
    !!evaluation of potential ET
    !OPEN(19761,FILE=OUTPUT_DIRECTORY // 'POTENTIAL_ET.MGB',STATUS='UNKNOWN')
    ! DO IT=1,NT
    !     WRITE(19761,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),ETPOTENCIAL(IT)
    !     ENDDO
    !CLOSE(19761)
    
    !inland delta ET at each IT or Month
    If(hdFLAG0>1)then 
    OPEN(19752,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET.MGB',STATUS='UNKNOWN')
    OPEN(19753,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET_INT.MGB',STATUS='UNKNOWN')
    OPEN(19754,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET_SOIL.MGB',STATUS='UNKNOWN')
    OPEN(19755,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET_FLOOD.MGB',STATUS='UNKNOWN')
    OPEN(19756,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_PREC.MGB',STATUS='UNKNOWN')
    OPEN(19757,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_PREC_month.MGB',STATUS='UNKNOWN')
    OPEN(19758,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET_month.MGB',STATUS='UNKNOWN')
    OPEN(19759,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_PREC_year.MGB',STATUS='UNKNOWN')
    OPEN(19760,FILE=OUTPUT_DIRECTORY // 'INLAND_DELTA_ET_year.MGB',STATUS='UNKNOWN')
    
    AREA_SUBBACIA=0 !CALCULO DA AREA DA SUBBACIA DO INLAND DELTA   (FLAG INLAND_DELTA(IC)=1)  
    DO IC=1,NC
    IF (INLAND_DELTA(IC)==1)    AREA_SUBBACIA=AREA_SUBBACIA+ACEL(IC)
    ENDDO
    
    DO IT=1,NT
        EVAPTUDO_NIGER(IT)=1000*EVAPTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !DIVIDE VOLUME DE ET PELA AREA, PARA UM DADO IT. OBTÉM ET EM MM/DTP
        EVAPTUDO_NIGER_INT(IT)=1000*EVAPTUDO_NIGER_INT(IT)/(AREA_SUBBACIA*1000000)    !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        EVAPTUDO_NIGER_SOLO(IT)=1000*EVAPTUDO_NIGER_SOLO(IT)/(AREA_SUBBACIA*1000000)   !ET FOR INLAND DELTA NIGER ON SOIL
        EVAPTUDO_NIGER_FLOOD(IT)=1000*EVAPTUDO_NIGER_FLOOD(IT)/(AREA_SUBBACIA*1000000) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        PRECTUDO_NIGER(IT)=1000*PRECTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !PREC FOR INLAND DELTA

        ! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19752,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),EVAPTUDO_NIGER(IT)
        WRITE(19753,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),EVAPTUDO_NIGER_INT(IT)
        WRITE(19754,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),EVAPTUDO_NIGER_SOLO(IT)
        WRITE(19755,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),EVAPTUDO_NIGER_FLOOD(IT)
        WRITE(19756,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),PRECTUDO_NIGER(IT)
    ENDDO
    
    DO IT=1,CONT_MES_NIGER
        PRECTUDO_NIGER_MES(IT) =1000*PRECTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        EVAPTUDO_NIGER_MES(IT) =1000*EVAPTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        !PRECTUDO_NIGER_MES(IT) =PRECTUDO_NIGER_MES(IT) /1000000
        !EVAPTUDO_NIGER_MES(IT) =EVAPTUDO_NIGER_MES(IT)/1000000
        
        WRITE(19757,'(I6,F16.6)')IT,PRECTUDO_NIGER_MES(IT) 
        WRITE(19758,'(I6,F16.6)')IT,EVAPTUDO_NIGER_MES(IT) 
    ENDDO
    
    DO IT=1,CONT_ANO_NIGER
        PRECTUDO_NIGER_ANO(IT) =1000*PRECTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        EVAPTUDO_NIGER_ANO(IT) =1000*EVAPTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        WRITE(19759,'(I6,F16.6)')IT,PRECTUDO_NIGER_ANO(IT) 
        WRITE(19760,'(I6,F16.6)')IT,EVAPTUDO_NIGER_ANO(IT) 
    ENDDO
         
    CLOSE(19751)
    CLOSE(19752)
    CLOSE(19753)
    CLOSE(19754)
    CLOSE(19755)
    CLOSE(19756)
    CLOSE(19757)
    CLOSE(19758)
    CLOSE(19759)
    CLOSE(19760)
    ENDIF
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR SOUTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !inland delta flooded area at each IT
    If(hdFLAG0>1)then 
    OPEN(19751,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_FA.MGB',STATUS='UNKNOWN')
    DO IT=1,NT
		! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19751,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_AFLTUDO(1,IT)
    ENDDO
    CLOSE(19751)
    ENDIF
    
    !inland delta ET at each IT or Month
    If(hdFLAG0>1)then 
    OPEN(19750,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ETkm3.MGB',STATUS='UNKNOWN')
    OPEN(19752,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ETmm.MGB',STATUS='UNKNOWN')
    OPEN(19753,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ET_INT.MGB',STATUS='UNKNOWN')
    OPEN(19754,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ET_SOIL.MGB',STATUS='UNKNOWN')
    OPEN(19755,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ET_FLOOD.MGB',STATUS='UNKNOWN')
    OPEN(19756,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_PREC.MGB',STATUS='UNKNOWN')
    OPEN(19757,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_PREC_month.MGB',STATUS='UNKNOWN')
    OPEN(19758,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ET_month.MGB',STATUS='UNKNOWN')
    OPEN(19759,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_PREC_year.MGB',STATUS='UNKNOWN')
    OPEN(19760,FILE=OUTPUT_DIRECTORY // 'S_INLAND_DELTA_ET_year.MGB',STATUS='UNKNOWN')
    
    AREA_SUBBACIA=0 !CALCCULO DA AREA DA SUBBACIA DO INLAND DELTA   (FLAG INLAND_DELTA(IC)=1)  
    DO IC=1,NC
    IF (SOUTHERN_DELTA(IC)==1)    AREA_SUBBACIA=AREA_SUBBACIA+ACEL(IC)
    ENDDO

    DO IT=1,NT
        WRITE(19750,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_EVAPTUDO_NIGER(IT)/1000000
        S_EVAPTUDO_NIGER(IT)=1000*S_EVAPTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !DIVIDE VOLUME DE ET PELA AREA, PARA UM DADO IT. OBTÉM ET EM MM/DTP
        S_EVAPTUDO_NIGER_INT(IT)=1000*S_EVAPTUDO_NIGER_INT(IT)/(AREA_SUBBACIA*1000000)    !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        S_EVAPTUDO_NIGER_SOLO(IT)=1000*S_EVAPTUDO_NIGER_SOLO(IT)/(AREA_SUBBACIA*1000000)   !ET FOR INLAND DELTA NIGER ON SOIL
        S_EVAPTUDO_NIGER_FLOOD(IT)=1000*S_EVAPTUDO_NIGER_FLOOD(IT)/(AREA_SUBBACIA*1000000) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        S_PRECTUDO_NIGER(IT)=1000*S_PRECTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !PREC FOR INLAND DELTA

        ! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19752,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_EVAPTUDO_NIGER(IT)
        WRITE(19753,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_EVAPTUDO_NIGER_INT(IT)
        WRITE(19754,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_EVAPTUDO_NIGER_SOLO(IT)
        WRITE(19755,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_EVAPTUDO_NIGER_FLOOD(IT)
        WRITE(19756,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),S_PRECTUDO_NIGER(IT)
    ENDDO
    
    DO IT=1,CONT_MES_NIGER
        S_PRECTUDO_NIGER_MES(IT) =1000*S_PRECTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        S_EVAPTUDO_NIGER_MES(IT) =1000*S_EVAPTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        WRITE(19757,'(I6,F16.6)')IT,S_PRECTUDO_NIGER_MES(IT) 
        WRITE(19758,'(I6,F16.6)')IT,S_EVAPTUDO_NIGER_MES(IT) 
    ENDDO
    
    DO IT=1,CONT_ANO_NIGER
        S_PRECTUDO_NIGER_ANO(IT) =1000*S_PRECTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        S_EVAPTUDO_NIGER_ANO(IT) =1000*S_EVAPTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        WRITE(19759,'(I6,F16.6)')IT,S_PRECTUDO_NIGER_ANO(IT) 
        WRITE(19760,'(I6,F16.6)')IT,S_EVAPTUDO_NIGER_ANO(IT) 
    ENDDO
         
    CLOSE(19750)
    CLOSE(19751)
    CLOSE(19752)
    CLOSE(19753)
    CLOSE(19754)
    CLOSE(19755)
    CLOSE(19756)
    CLOSE(19757)
    CLOSE(19758)
    CLOSE(19759)
    CLOSE(19760)
    ENDIF
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>COMPUTATION FOR NORTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !inland delta flooded area at each IT
    If(hdFLAG0>1)then 
    OPEN(19751,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_FA.MGB',STATUS='UNKNOWN')
    DO IT=1,NT
		! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19751,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_AFLTUDO(1,IT)
    ENDDO
    CLOSE(19751)
    ENDIF
    
    !inland delta ET at each IT or Month
    If(hdFLAG0>1)then 
    OPEN(19750,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ETkm3.MGB',STATUS='UNKNOWN')
    OPEN(19752,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ETmm.MGB',STATUS='UNKNOWN')
    OPEN(19753,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ET_INT.MGB',STATUS='UNKNOWN')
    OPEN(19754,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ET_SOIL.MGB',STATUS='UNKNOWN')
    OPEN(19755,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ET_FLOOD.MGB',STATUS='UNKNOWN')
    OPEN(19756,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_PREC.MGB',STATUS='UNKNOWN')
    OPEN(19757,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_PREC_month.MGB',STATUS='UNKNOWN')
    OPEN(19758,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ET_month.MGB',STATUS='UNKNOWN')
    OPEN(19759,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_PREC_year.MGB',STATUS='UNKNOWN')
    OPEN(19760,FILE=OUTPUT_DIRECTORY // 'N_INLAND_DELTA_ET_year.MGB',STATUS='UNKNOWN')
    
    AREA_SUBBACIA=0 !CALCCULO DA AREA DA SUBBACIA DO INLAND DELTA   (FLAG INLAND_DELTA(IC)=1)  
    DO IC=1,NC
    IF (NORTHERN_DELTA(IC)==1)    AREA_SUBBACIA=AREA_SUBBACIA+ACEL(IC)
    ENDDO
    
    DO IT=1,NT
        WRITE(19750,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_EVAPTUDO_NIGER(IT)/1000000 !km³
        N_EVAPTUDO_NIGER(IT)=1000*N_EVAPTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !DIVIDE VOLUME DE ET PELA AREA, PARA UM DADO IT. OBTÉM ET EM MM/DTP
        N_EVAPTUDO_NIGER_INT(IT)=1000*N_EVAPTUDO_NIGER_INT(IT)/(AREA_SUBBACIA*1000000)    !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        N_EVAPTUDO_NIGER_SOLO(IT)=1000*N_EVAPTUDO_NIGER_SOLO(IT)/(AREA_SUBBACIA*1000000)   !ET FOR INLAND DELTA NIGER ON SOIL
        N_EVAPTUDO_NIGER_FLOOD(IT)=1000*N_EVAPTUDO_NIGER_FLOOD(IT)/(AREA_SUBBACIA*1000000) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        N_PRECTUDO_NIGER(IT)=1000*N_PRECTUDO_NIGER(IT)/(AREA_SUBBACIA*1000000) !PREC FOR INLAND DELTA

        ! Record Discharge Time Series in NB predefined catchments outlets 
		WRITE(19752,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_EVAPTUDO_NIGER(IT)
        WRITE(19753,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_EVAPTUDO_NIGER_INT(IT)
        WRITE(19754,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_EVAPTUDO_NIGER_SOLO(IT)
        WRITE(19755,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_EVAPTUDO_NIGER_FLOOD(IT)
        WRITE(19756,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),N_PRECTUDO_NIGER(IT)
    ENDDO
    
    DO IT=1,CONT_MES_NIGER
        N_PRECTUDO_NIGER_MES(IT) =1000*N_PRECTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        N_EVAPTUDO_NIGER_MES(IT) =1000*N_EVAPTUDO_NIGER_MES(IT) /(AREA_SUBBACIA*1000000) 
        WRITE(19757,'(I6,F16.6)')IT,N_PRECTUDO_NIGER_MES(IT) 
        WRITE(19758,'(I6,F16.6)')IT,N_EVAPTUDO_NIGER_MES(IT) 
    ENDDO
    
    DO IT=1,CONT_ANO_NIGER
        N_PRECTUDO_NIGER_ANO(IT) =1000*N_PRECTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        N_EVAPTUDO_NIGER_ANO(IT) =1000*N_EVAPTUDO_NIGER_ANO(IT) /(AREA_SUBBACIA*1000000) 
        WRITE(19759,'(I6,F16.6)')IT,N_PRECTUDO_NIGER_ANO(IT) 
        WRITE(19760,'(I6,F16.6)')IT,N_EVAPTUDO_NIGER_ANO(IT) 
    ENDDO
         
    CLOSE(19750)
    CLOSE(19751)
    CLOSE(19752)
    CLOSE(19753)
    CLOSE(19754)
    CLOSE(19755)
    CLOSE(19756)
    CLOSE(19757)
    CLOSE(19758)
    CLOSE(19759)
    CLOSE(19760)
    ENDIF
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END COMPUTATION FOR NORTHERN INLAND DELTA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
    
    
    wwm_mean=0.0
    DO IT=1,NT         
        AREA_SUBBACIA=0.0
        DO IC=1,NC !computes w/wm values of each upstream catchment, pondered by area
            if (NORTHERN_DELTA(IC)==1) then
                wwm_mean(IT)=wwm_mean(IT)+STUDO_DELTA(IC,IT)*ACEL(IC) !computes first catchment w/wm value
                AREA_SUBBACIA=AREA_SUBBACIA+ACEL(IC)
            endif
        ENDDO
        wwm_mean(IT)=wwm_mean(IT)/AREA_SUBBACIA
    ENDDO
    
    OPEN(19750,FILE=OUTPUT_DIRECTORY // 'STUDO_N_INLAND_DELTA.MGB',STATUS='UNKNOWN')
    DO IT=1,NT
        WRITE(19750,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),wwm_mean(IT)
    ENDDO
    CLOSE(17950)
    
    wwm_mean=0.0
    DO IT=1,NT    
        AREA_SUBBACIA=0.0
        DO IC=1,NC !computes w/wm values of each upstream catchment, pondered by area
            if (SOUTHERN_DELTA(IC)==1) then
                wwm_mean(IT)=wwm_mean(IT)+STUDO_DELTA(IC,IT)*ACEL(IC) !computes first catchment w/wm value
                AREA_SUBBACIA=AREA_SUBBACIA+ACEL(IC)
            endif
        ENDDO
        wwm_mean(IT)=wwm_mean(IT)/AREA_SUBBACIA
    ENDDO
    
    OPEN(19750,FILE=OUTPUT_DIRECTORY // 'STUDO_S_INLAND_DELTA.MGB',STATUS='UNKNOWN')
    DO IT=1,NT
        WRITE(19750,'(3I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),wwm_mean(IT)
    ENDDO
    CLOSE(17950)
          
	RETURN

71	FORMAT(I10,<NUMHIDG>F15.3)
72	FORMAT(<NUMHIDG>F10.3)
	END
