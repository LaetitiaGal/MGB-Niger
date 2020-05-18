	SUBROUTINE PREVISAO

!---------------------------------------------------------------------------------
! PREVISAO.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This subroutine contains the routine that runs the model in forecast mode.
!
! Usage:
!
! CALL PREVISAO
!
! uses modules, functions, and subroutines
!
! * VARS_MAIN.f90
! * CELULA.f90
! * REDE.f90
! * CONDINIC.f90
! * ATUALIZA.f90
! * LECHUVA.f90
! * LECLIMA.f90
!
! opens
!
! * QTUDO_ATUAL.MGB
! * QTUDO_PREV.MGB   
! * QITUDO_ATUAL.MGB    
! * '//'Previsao_'//ARQHID(K)
! * Previsao.meta
! * VAZAO.HIG
! *	QPROP.HIG
! * '//ARQHID(K)	
!
! reads
!
! * Previsao.meta
!
! creates
!
! * QTUDO_ATUAL.MGB
! * QTUDO_PREV.MGB   
! * QITUDO_ATUAL.MGB    
! * '//'Previsao_'//ARQHID(K)
! * VAZAO.HIG
! *	QPROP.HIG
! * '//ARQHID(K)	
!
!---------------------------------------------------------------------------------
! Licensing:
!
! This code is distributed under the...
!
! Version/Modified:
!
! 2014.11.12 - 13 November 2014
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
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------
	
	!Variables declaration

	USE VARS_MAIN
	IMPLICIT NONE
	
	INTEGER K,KHID !Counters
	INTEGER KB,JCB,MWP 	!Counters
	INTEGER IANTE,ITP,IHORIZ,INIPREV,IFIPREV,LINIARQ
	CHARACTER (15) ARQPREV !Forecasted rainfall filename
	INTEGER ITEMP
	REAL ERROANT1,ERROANT2,ERROCORR
	REAL RIANTE,RIHORIZ
	INTEGER IOPCHU
	CHARACTER (4) ABOBO
	INTEGER IDIAPREV,IANOPREV,IMESPREV,IHORAPREV,IP
	INTEGER DIA_ANTE,MES_ANTE,ANO_ANTE,HORA_ANTE

	INTEGER MDIAS,CONTADOR
	INTEGER ITINI,IBOBO
	CHARACTER (1) BARRA
	INTEGER ITR,NTRECH,NTMUSK,LT,I
	INTEGER IANTEMES
	
    CHARACTER(40) strTrash
    INTEGER J,L,IH,ITAUX
    
    !Call routine that sets the initial conditions of a simulation with the model
	CALL CONDINIC

	!Activate all sub-basins during the forecasting loop
    ATIVABACIA=1 
	
	!Sets an initial value for used variables
	ITINI=1
    QPREV=-1.
    ISUBSTAUX=ISUBST
    QTUDO=-1.
    QITUDO=-1.
    QRG=-1.      !All streamflows are settled to start as -1 as a standard (NoData code).
    
    !Open the output files with simulated and forecasted streamflow: Simulated (until t0 of the forecast), forecasted, Simulated incremental flows
    OPEN(FILTUD,FILE='c:\MGB-hora\output\QTUDO_ATUAL.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')
    OPEN(FILTUD2,FILE='c:\MGB-hora\output\QTUDO_PREV.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')   
    OPEN(FILTUD3,FILE='c:\MGB-hora\output\QITUDO_ATUAL.MGB',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')    
 
    !Open the output files with individual forecasted hydrographs
    DO K=1,NUMHIDG
        OPEN(801+K,FILE='c:\MGB-hora\output\'//'Previsao_'//ARQHID(K),STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(801+K,*)'   DIA   MES   ANO  HORA        Q_(m3/s)'
    ENDDO
  
    !Open the output file with the forecast description (dates, horizon, start time, and etc...)
    IF (1==1) THEN
         open(41,FILE='c:\MGB-hora\input\Previsao.meta',STATUS='OLD',ACTION='READ')
         do i=1,7
           READ(41,*)strTrash
         end do
         WRITE(*,*)'Start forecast at: Year, Month, Day, Hour, Horizon'
         READ(41,*)IANOPREV,IMESPREV,IDIAPREV,IHORAPREV,IHORIZ
         WRITE(*,'(5I6)')IANOPREV,IMESPREV,IDIAPREV,IHORAPREV,IHORIZ
         IOPCHU=1
    END IF
    
    WRITE(*,*)'Simulating and performing data assimilation...'
    
	!Starts the time loop

	IT=ITINI-1
		
	DO WHILE (IT<NT-IHORIZ+1) !Runs the simulation time loop until the first day of the forecast
		
		IT=IT+1
		WRITE(*,72)IT,100*IT/NT
		72 FORMAT('+',' INTERVALO DE TEMPO:',I10,I10,'%')

		!-- This IF checks for the model states
		IF(IT==ITINI)THEN !At the first time step initial conditions ate settled      

           	!INITIAL CONDITIONS SETTING
			!In the interceptation
			SI=0.0 !Water amount is zero
			!In the soil
			DO IC=1,NC
				IB=IBAC(IC)
				DO IU=1,NU
					W(IC,IU)=WM(IB,IU)*0.40 !Have a % of water
				ENDDO
			ENDDO
			!From the reservoirs
			DO IC=1,NC
				IB=IBAC(IC)
				VBAS(IC)=QESP(IB)*ACEL(IC)*CB(IB)*3600.
				VINT(IC)=0.0
				VSUP(IC)=0.0
			ENDDO
			!In the temperature
			TA=20.0
			!Forecasted streamflows
			QPREV=-1.0 !Use NoData code
            WRITE(*,*)'OK'
			!END OF INITIAL CONDITIONS SETTINGS
		
		ELSE !At others time steps, we use values saved during the simulation comming from the last time step
			!In the soil
			DO IC=1,NC
				IB=IBAC(IC)
				DO IU=1,NU
					W(IC,IU)=WPREV(IC,IU)
				ENDDO
			ENDDO
			!In the reservoirs
			VBAS=VBASPREV
			VINT=VINTPREV
			VSUP=VSUPPREV
			TA=TAPREV
			QM2=QM2PREV
			QJ2=QJ2PREV
			SI=SIPREV
			QJ1=QJ1PREV
			QRIOINI=QRIOINIPREV
			QCONTORM=QCONTORMPREV
			QCONTORJ=QCONTORJPREV
			QCEL1=QCEL1PREV
			QCEL2=QCEL2PREV
			PMB2=PMB2PREV
			PMI2=PMI2PREV
			PMS2=PMS2PREV
			PJB2=PJB2PREV
			PJI2=PJI2PREV
			PJS2=PJS2PREV
		ENDIF
		!-- End of checking for the model states
		
		JDIA=IDINI+INT((IT+HORAINI)/(86400./DTP)) !Julian calendar day verification
		CALL CALDAT(JDIA,IMES,IDIA,IANO)
        DIAH(IT)=IDIA
        MESH(IT)=IMES
        ANOH(IT)=IANO
        HORAH(IT)=MOD(IT,24)+HORAINI !Hour of the day corresponding to the interval IT
      
		!Subroutine that reads and prepares rainfall data
		ITCHUVA=IT
		CALL LECHUVA
        
		!Subroutine that reads and prepare climatological
		TONTEM=TA
		CALL LECLIMA
	
		!Subroutine that calculate the vertical and horizontal fluxes at the minibasins (hydrology)
		CALL CELULA
	
		!Subroutine that does the flow routing (hydraulics)
		ISUBST=ISUBSTAUX
		CALL REDE

		!Save the streamflow at the minibasins where we have observed streamflow
		DO K=1,NOBS
			KHID=IQOBS(K) !Minibasins that correspoounds to the gauge
			QR(K,IT)=QJ2(KHID) !QR(K,IT) will be compared with QOBS(K,IT) in the routine FOBJ
		ENDDO
		
		!Save SubBasins data
		DO KB=1,NB	
			JCB=IEXUT(KB) !Minibasin in the outlet of the Sub-basin
			QRB(KB,IT)=QJ2(JCB)
		ENDDO
	
		!Save the hydrographs in the places defined by the file PARHIG
		DO K=1,NUMHIDG
			KHID=IHIDG(K) !IHIDG(K) is the number of the minibasin which we want to save the hydrographs
			QRG(K,IT)=QJ2(KHID) !QRG is where the hydrographs are saved		
		ENDDO
	
		!Save the streaflows according to the origin
		MWP=NC !Minibasin where they are saved
		QB(IT)=QJ2(MWP)*PJB2(MWP)
		QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
		QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)

		!Salve the model most important states before the next cycle
		WPREV=W
		VBASPREV=VBAS
		VINTPREV=VINT
		VSUPPREV=VSUP
		TAPREV=TA
		QM2PREV=QM2
		QJ2PREV=QJ2
		SIPREV=SI
		QJ1PREV=QJ1
		QRIOINIPREV=QRIOINI
		QCONTORMPREV=QCONTORM
		QCONTORJPREV=QCONTORJ
		QCEL1PREV=QCEL1
		QCEL2PREV=QCEL2
		PMB2PREV=PMB2
		PMI2PREV=PMI2
		PMS2PREV=PMS2
		PJB2PREV=PJB2
		PJI2PREV=PJI2
		PJS2PREV=PJS2

		IFIPREV=NT !Time step when the forecast should end

	    !Defines the time step when we should start to do data assimilation. To define this it is usefull to check ATUALIZA routine prior, due to conflicts.
	    IF(IT>48)THEN 
		    CALL ATUALIZA !Call the data assimilation subroutine
		ENDIF
		
		ITEMP=IT
		JDIA=IDINI+INT((ITEMP+HORAINI)/(86400./DTP)) !Verify Julian day
		CALL CALDAT(JDIA,IMES,IDIA,IANO) !Verify dates
		IDIA=IDIA !Day corresponding to the IT
		IMES=IMES !Mounth corresponding to the IT
		IANO=IANO !Year corresponding to the IT
		IHORA=MOD(IT,24)+HORAINI !Hour corresponding to the IT	

		!Here we check if we should start the forecasting loop
		!If we should start, then we start the forecasting procedure
		IF(IDIA==IDIAPREV.AND.IMES==IMESPREV.AND.IANO==IANOPREV.AND.IHORA==IHORAPREV)THEN		

		    WRITE(*,*)'Realizando a previsao...'
		    
			ITEMP=IT !Save the IT which the forecasting starts	
            
            !Starts the forecasting cycle
		    DO ITP=IT+1,IT+IHORIZ !Runs the forecast between the first time interval and the final one
    	
    			IANTE=ITP-ITEMP !IANTE is the lead-time of the forecast    			
		    
        		JDIA=IDINI+INT((ITP+HORAINI)/(86400./DTP)) !Verify Julian day
		        CALL CALDAT(JDIA,IMES,IDIA,IANO)
		        DIA_ANTE=IDIA
		        MES_ANTE=IMES
		        ANO_ANTE=IANO
		        HORA_ANTE=MOD(ITP,24)+HORAINI !Hour corresponding to IT

    			!IOPCHU=1 !Keep reading the one single file	of rainfall		
    		    !IOPCHU=2 !Use a separete file of forecasted rainfall
    			!IOPCHU=3 !Use zero rainfall
    		   			    
    			SELECT CASE (IOPCHU)
    				CASE (1) !Keep reading the single file    				    
	    				ITCHUVA=ITP
		    			CALL LECHUVA !Continue reading the single file
			    	CASE (2) !Separate file
				    	ITCHUVA=ITP-ITEMP !Reads the separate file
				    	WRITE(*,*)'LOOP DA PREVISAO',IT,ITP,ITCHUVA,ITEMP	
					    CALL LEPREV
    				CASE(3) !Null rainfall
	    			    P=0.0
		    		CASE DEFAULT
			    		STOP 'OPCAO DE CHUVA DESCONHECIDA'
				ENDSELECT

		    	IT=ITP
				
    			!Subroutine that reads and prepare climatological
	    		TONTEM=TA
		    	CALL LECLIMA
				
			    !Subroutine that calculate the vertical and horizontal fluxes at the minibasins (hydrology)
    			CALL CELULA
				
	    		!Subroutine that does the flow routing (hydraulics)
		    	!ISUBST=0 !Be careful here: it do or do not do data substitution in the future
			    CALL REDE
	
    			!Save the streamflow at the minibasins where we have observed streamflow
	    		DO K=1,NOBS
		    		KHID=IQOBS(K) !Minibasin that corresponds to the gauge
			    	QR(K,IT)=QJ2(KHID) !QR(K,IT) can be compared to QOBS(K,IT) in the subroutine FOBJ
    			ENDDO    			        

    			DO K=1,NUMHIDG !Save the hydrographs in the places defined by the file PARHIG 
	    			KHID=IHIDG(K) !IHIDG(K) is the number of the minibasins where we want to save the hydrographs
		    		QPREV(K,IANTE)=QJ2(KHID) !QPREV save the forecasted streamflows at the desired places
			    ENDDO
			    
                !Save hydrographs of forecasted stramflow
                DO K=1,NUMHIDG
                    WRITE(801+K,'(4I6,F16.6)')DIA_ANTE,MES_ANTE,ANO_ANTE,HORA_ANTE,QPREV(K,IANTE) !!GRAVA RSULTADOS EM PREVISAO_POSTO.TXT
                ENDDO
                			    
		    ENDDO
		    !End of the forecasting loop
                 
            !Write -1 at the end of the files
            DO K=1,NUMHIDG
                WRITE(801+K,'(4I6,F16.6)')DIA_ANTE,MES_ANTE,ANO_ANTE,HORA_ANTE,-1.0
            ENDDO
            
		    IT=ITEMP !Return the value of IT at the end of the forecast (used on batch forecasting)    
       	    
			!For a single forecast, exit after the forecast loop:
			EXIT
			
		ENDIF !End of the forecast
        
	ENDDO !End of the time loop
	
	!Close all files
    DO K=1,NUMHIDG
        CLOSE(801+K)
    ENDDO
    
	OPEN(FILHID,FILE='c:\MGB-hora\output\VAZAO.HIG',STATUS='UNKNOWN') !File with simulated hydrographs
	OPEN(FILPRP,FILE='c:\MGB-hora\output\QPROP.HIG',STATUS='UNKNOWN') !File with streamflows according to the origin
	
	DO IT=1,NT
		!Save the hydrographs at the NB desired outlets
		WRITE(FILHID,71)IT,(QRG(KB,IT),KB=1,NUMHIDG)
		!Save the file with streamflows according to the origin
		WRITE(FILPRP,71)IT,QB(IT),QBI(IT),QBIS(IT)
		
		IF (IT<=ITEMP) THEN !Save a binary with assimilated streamflows at all the minibasins	    

		    WRITE(FILTUD,REC=IT)(QTUDO(IC,IT),IC=1,NC)		     
		    WRITE(FILTUD3,REC=IT)(QITUDO(IC,IT),IC=1,NC)
		    
		ELSE !Save a binary file woth the forecasts at all mini-basins at all lead-times
		     
		    WRITE(FILTUD2,REC=IT-ITEMP)(QTUDO(IC,IT),IC=1,NC) !REC=IT-ITEMP  - CAUTIN: the first line must be the streamflow at the start of the forecast A

		END IF		
	ENDDO
	
	CLOSE (FILHID)
	CLOSE (FILPRP)
	
    !Open the files and save the invidual hydrographs
    DO K=1,NUMHIDG
        OPEN(801,FILE='c:\MGB-hora\output\'//ARQHID(K),STATUS='UNKNOWN',ACTION='WRITE')
        
            WRITE(801,*)'   DIA   MES   ANO  HORA        Q_(m3/s)'
            DO IT=1,NT
                WRITE(801,'(4I6,F16.6)')DIAH(IT),MESH(IT),ANOH(IT),HORAH(IT),QRG(K,IT)
            ENDDO
        CLOSE(801)
    ENDDO

71	FORMAT(I10,<NUMHIDG>F15.3)
75	FORMAT(I6,24F7.0)
     
	RETURN
	END
	!End the routine
