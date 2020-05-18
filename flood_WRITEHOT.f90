    SUBROUTINE flood_WRITEHOT
	!Esta subrotina escreve o arquivo de condições iniciais completo
    ! Modificado por Rodrigo Paiva 18/10/2012 ! RP
    ! Modificado por Vinícius Siqueira em 23/01/2017

    !DECLARAÇÃO DE VARIÁVEIS
	USE VARS_MAIN
    USE VARS_INERC
	
  	!DEFINE QUE NÃO PODEMOS TER VARIÁVEIS IMPLÍCITAS
	IMPLICIT NONE
    integer:: i,nStatVar,j
    real,allocatable:: WB_stateVec(:)	 
    real*8,allocatable:: h_stateVec(:)	 
		
	nStatVar=(nC*nU)+(6*nC)+(nC*nU)+(sum(NSUBT)+nC)+nC !
	allocate(WB_stateVec(nStatVar))
	   
    write(*,*)'Writing Water Balance State Variables...'
	open(FILHOTSTART_WB,FILE=OUTPUT_DIRECTORY // 'StateVars_WB.hot',STATUS='unknown',RECL=nStatVar,FORM='UNFORMATTED',ACCESS='DIRECT')

    !******************************************************************************
	! Writing state variables:
	i=0
	do iC=1,nC
		do iU=1,nU
			i=i+1
			WB_StateVec(i)=W(iC,iU) !Soil Water for each HRU
		enddo
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=VBAS(iC) !Reservoir Volume for baseflow
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=VINT(iC) !Reservoir Volume for subsurface flow
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=VSUP(iC) !Reservoir Volume for surface flow
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=TA(iC) !Temperature???
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=QM2(iC) !Streamflow at upstream of catchments
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=QJ2(iC) !Streamflow at downstream of catchments
	enddo
	do iC=1,nC
		do iU=1,NU
			i=i+1
			WB_StateVec(i)=SI(iC,iU) !Intercepted water stored on canopy reservoir
		enddo
	enddo
	do iC=1,nC
		do j=1,NSUBT(IC)+1
			i=i+1
			WB_StateVec(i)=QRIOINI(iC,j)
		enddo
	enddo
	do iC=1,nC
			i=i+1
			WB_StateVec(i)=QCEL2(iC) !Runoff generated at catchment
	enddo
		
	write(FILHOTSTART_WB,REC=1) (WB_StateVec(i),i=1,nStatVar)
    deallocate(WB_stateVec)	
    
	close(FILHOTSTART_WB)
	!************************************************************
	    
    nStatVar=(9*nC)! Variables from inertial routing
    
    write(*,*)'Writing Hydraulic State Variables...'
    !All variables below are in double precision (Real*8)
	open(FILHOTSTART_Inert,FILE=OUTPUT_DIRECTORY // 'StateVars_inertial.hot',STATUS='unknown',RECL=(nStatVar*2),FORM='UNFORMATTED',ACCESS='DIRECT')
    allocate(H_stateVec(nStatVar))
    
    i=0
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Q2fl(iC) !Streamflow at catchment (iC)
    enddo    
            
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Hfl(iC) !Water depth at catchment (iC)
	enddo    
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Yfl(iC) !Water elevation at catchment (iC)
    enddo        
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Vol1(iC) !Previous Water storage at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Vol2(iC) !Actual Water storage at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Area2(iC) !Flooded Area at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=dble(jtab(iC)) !index of Volume level at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Q2face(iC) !Streamflow coming (or leaving) the catchment through connections 
    enddo
    
    do iC=1,nC
			i=i+1
			H_StateVec(i)=Q2viz(iC) !Updated flow for interconnections (signal can be positive or negative)
    enddo
    
    write(FILHOTSTART_Inert,REC=1) (H_StateVec(i),i=1,nStatVar)
    
    deallocate(H_stateVec)	
    close(FILHOTSTART_Inert)
    
	RETURN
	END