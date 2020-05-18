    !*********************************************************************************
    !
    !  SUBROUTINE FOBJ evaluates objective functions for model assesment and calibration
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
	!    This routine calculates error functions for observed x calculated discharges
	!       time series. There are three error functions:
	!		- Nash-Sutcliff effiency
	!		- Nash-Sutcliff effiency of log-tranformed values
	!		- Bias on Volume (i.e. sum[Vcalc(i)]/sum[Vobs(i)]-1. )
	!
	!	For calibration, each error function can be weighted by a a number of subbasin(or stations)
	!		as defined in calibration input file.
	!
	!	Errors are stored in variable VFO
	!
    !  Usage:
    !
    !    CALL FOBJ
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90    
	!    * module     VARS_CALIB  in      VARS_CALIB.f90  
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
	
	SUBROUTINE FOBJ
	USE VARS_MAIN
    USE VARS_INERC
	USE VARS_CALIB !RP
	IMPLICIT NONE
	
	INTEGER NANO
	INTEGER KB				! loop counter
	REAL SR,ERRVT,SRL		! for summation of errors
	REAL POND(NOBS,3)		! sub-basin weight for each objective-function (i.e. 3 function)
!	REAL POND(NOBS,nf)

	INTEGER INIAJU 			! time interval for start of comparisons 
	REAL SOMAOBS 			! sum of observed values
	REAL SOMACAL 			! sum of calculated values
	REAL SOMALOGS 			! sum of log-tranformed observed values
	REAL SOMALCAL 			! sum of log-tranformed calculated values
	INTEGER ITVAL 			! counts number of intervals w observed data
	REAL XMOBS,XLOGS 		! average of observed and log-transformed observed
	REAL SQDES(NOBS) 		! sum of squared deviations for observed
	REAL SOMQ 				! sum of squared deviations between observed and its average
	REAL SQLOG 				! sum of squared deviations for log-tranformed observed
	REAL SOLOG 				! sum of squared deviations between log-transformed observed and its average
	INTEGER i,j

	!QR(NOBS,NT),QOBS(NOBS,NT)  ! discharge on sub-basin outlets: calc and observed
	! NOBS,IQOBS(NOBS) 			! number of days in observed data

    iniaju=360		!starts only after 360 intervals (1st year for daily data)
	SQDES=0.0
	R2=-100.0
	R2L=-100.0
	ERRV=1000.0

	DO KB=1,NOBS
		IB=KB
		
		SOMAOBS=0.0
		SOMACAL=0.0
		SOMALOGS=0.0
		SOMALCAL=0.0
		SOMQ=0.0
		SQLOG=0.0
		SOLOG=0.0
		ITVAL=0
		IF(INIAJU.LT.NT)THEN
			DO IT=INIAJU,NT
				!WRITE(*,*)IB,IT,QOBS(IB,IT),QR(IB,IT)
				IF(QOBS(IB,IT).GE.0.0) THEN
					SOMAOBS=SOMAOBS+QOBS(IB,IT)
					SOMACAL=SOMACAL+QR(IB,IT)
					SOMALOGS=SOMALOGS+LOG(QOBS(IB,IT)+0.0001) !Sums a small value to avoid undefinition in log(0.0)
                    SOMALCAL=SOMALCAL+LOG(QR(IB,IT)+0.0001)
					ITVAL=ITVAL+1
				ENDIF
			ENDDO
            
			XMOBS=SOMAOBS/ITVAL
			XLOGS=SOMALOGS/ITVAL
			IMES=0
			NANO=0
			DO IT=INIAJU,NT
				!WRITE(*,*)IB,IT,QOBS(IB,IT),QR(IB,IT)
				IF(QOBS(IB,IT).GE.0.0) THEN
					SQDES(IB)=SQDES(IB)+(QOBS(IB,IT)-QR(IB,IT))**2.0
					SOMQ=SOMQ+(QOBS(IB,IT)-XMOBS)**2.0
					SQLOG=SQLOG+(LOG(QOBS(IB,IT)+0.0001)-LOG(QR(IB,IT)+0.0001))**2.0
					SOLOG=SOLOG+(LOG(QOBS(IB,IT)+0.0001)-XLOGS)**2.0
				ENDIF
			ENDDO
			
			if (ITVAL>0) then
				R2(IB)=1.-SQDES(IB)/SOMQ
				R2L(IB)=1.-SQLOG/SOLOG
				ERRV(IB)=(SOMACAL/SOMAOBS-1.0)*100.
			endif
		ENDIF
	ENDDO



	
	SR=0.0
	SRL=0.0
	ERRVT=0.0




	! Considers objective function as weighted-average from observed sub-basin data (or station) that
	! were selected in calibration. Only for calibration!
    if (icalib==0) then
		calibFLAG=1.0
	endif
	
	POND=0.0
	do kb=1,nobs
		do j=z1,nf
			if (calibFLAG(kb,j)>0.0) POND(kb,j)=calibFLAG(kb,j)/sum(calibFLAG(:,j))
        pause
		enddo
    enddo
    
	VFO=0.0
	do kb=1,nobs
!	write(*,*) kb,VFO(1)
!	write(*,*) R2(kb)
!	write(*,*) POND(kb,1)
		if (POND(kb,1)>0.0)  VFO(1)=VFO(1)+(1.-R2(kb))*POND(kb,1)
		if (POND(kb,2)>0.0) VFO(2)=VFO(2)+(1.-R2L(kb))*POND(kb,2)
		if (POND(kb,3)>0.0) VFO(3)=VFO(3)+ABS(ERRV(kb))*POND(kb,3)
	enddo
	

!	WRITE(*,*)(R2(IB),IB=1,NOBS)
!	WRITE(*,*)(R2L(IB),IB=1,NOBS)
!	WRITE(*,*)(ERRV(IB),IB=1,NOBS)


	write(*,*) VFO(1),VFO(2),VFO(3)


	RETURN
	END
