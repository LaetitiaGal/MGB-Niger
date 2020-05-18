    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This Sub-Routine calculates Flows by Muskingum-Cunge method.
    !
    !
    ! Usage:
    !
    !
    ! uses modules, functions, and subroutines
    !
    ! * USE VARS_MAIN
    !
    ! opens
    !
    ! * no files are opened in this routine
    !
    ! reads
    !
    ! * no files are read in this routine
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
    !    2014.09.001 - 09 September 2014
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
    
	SUBROUTINE MUSK(QMX1,QMX2,QJX1,QJX2,NTRECH,NTC,C1,C2,C3)

    ! Variables and Parameters
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER ITC,ITR
	
	!Muskingum-Cunge Coefficients 
	REAL C1,C2,C3
	!Auxiliary variables
	REAL QMX1,QMX2,QJX1,QJX2
	!Sub-Rivers number / Number of intervals in a day
	INTEGER NTRECH,NTC
	
	REAL,ALLOCATABLE:: QC(:,:)
	REAL QPERDAS
	ALLOCATE (QC(NTRECH+1,NTC+1))


	IF(NTRECH>NUMUSK) STOP 'MUITOS SUBTRECHOS NA ROTINA MUSK'

	QC=0.0 
	QCONTORJ(IC,1)=QJX1


	QPERDAS=0.0


	!Upstream boundary conditions	
	DO ITC=1,NTC+1	  
	    !Modified to remove evapotranspiration of water surfaces
		QC(1,ITC)=MAX(QCONTORM(IC,ITC)-EVQ(IC)-QPERDAS,0.0) 
	ENDDO

	

	!Initial conditions
	DO ITR=1,NTRECH+1
		QC(ITR,1)=QRIOINI(IC,ITR)
	ENDDO

	DO ITC=1,NTC
		DO ITR=1,NTRECH
			QC(ITR+1,ITC+1)=C1*QC(ITR,ITC)+C2*QC(ITR,ITC+1)+C3*QC(ITR+1,ITC)
			!Avoid negative flows
			QC(ITR+1,ITC+1)=MAX(QC(ITR+1,ITC+1),0.0)
		ENDDO
		!Save output flows (all stretches and times)
		QCONTORJ(IC,ITC+1)=QC(NTRECH+1,ITC+1)
	ENDDO
	
    !Save Initial Condition values to next Muskingum-Cunge simulation
	QRIOINI(IC,1)=QCONTORM(IC,NTC+1)
	DO ITR=2,NTRECH+1
		QRIOINI(IC,ITR)=QC(ITR,NTC+1)
	ENDDO

	!Save last downstream value. This value is returned to REDE SUBROUTINE.
	QJX2=QC(NTRECH+1,NTC+1) 

	DEALLOCATE (QC)
	RETURN
	END
