    SUBROUTINE flood_READHOT !chili peppers
	!Esta subrotina le o arquivo de condições iniciais completo
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
	
    write(*,*)
    write(*,*)'Reading Water Balance State Variables...'
	open(FILHOTSTART_WB,FILE=INPUT_DIRECTORY // 'StateVars_WB.hot',STATUS='old',RECL=nStatVar,FORM='UNFORMATTED',ACCESS='DIRECT')
	
	!******************************************************************************
	read(FILHOTSTART_WB,REC=1) (WB_StateVec(i),i=1,nStatVar) ! Reading state variables:
	
    i=0
	do iC=1,nC
		do iU=1,nU
			i=i+1
			W(iC,iU)=WB_StateVec(i) !Reading Soil Moisture
		enddo
	enddo
	do iC=1,nC
			i=i+1
			VBAS(iC)=WB_StateVec(i) !Reading water in baseflow reservoir
	enddo
	do iC=1,nC
			i=i+1
			VINT(iC)=WB_StateVec(i) !Reading water in subsurface reservoir
	enddo	
	do iC=1,nC
			i=i+1
			VSUP(iC)=WB_StateVec(i) !Reading water surface volume
	enddo	
	do iC=1,nC
			i=i+1
			TA(iC)=WB_StateVec(i)  !Reading Temperature
	enddo	
	do iC=1,nC
			i=i+1
			QM2(iC)=WB_StateVec(i) !Reading Flow entering into the catchment
	enddo
	do iC=1,nC
			i=i+1
			QJ2(iC)=WB_StateVec(i) !Reading Flow going out of catchment
	enddo
	do iC=1,nC
		do iU=1,NU
			i=i+1
			SI(iC,iU)=WB_StateVec(i)  !Reading Canopy interception volume 
		enddo
	enddo
	do iC=1,nC
		do j=1,NSUBT(IC)+1
			i=i+1
			QRIOINI(iC,j)=WB_StateVec(i)  !Reading Canopy interception volume
		enddo
	enddo
	do iC=1,nC
			i=i+1
			QCEL2(iC)=WB_StateVec(i)   !Reading Runoff generation at each cell
	enddo
		
	close(FILHOTSTART_WB)	  
    deallocate(WB_stateVec)	

    !************************************************************
    nStatVar=(9*nC)! Variables from inertial routing
    allocate(H_StateVec(nStatVar))
    
    write(*,*)'Reading Hydraulic State Variables...'
    write(*,*)
    
    open(FILHOTSTART_Inert,FILE=INPUT_DIRECTORY // 'StateVars_inertial.hot',STATUS='old',RECL=(nStatVar*2),FORM='UNFORMATTED',ACCESS='DIRECT')    
    read(FILHOTSTART_Inert,REC=1) (H_StateVec(i),i=1,nStatVar) ! Reading state variables:
    
    i=0
    do iC=1,nC
			i=i+1
			Q2fl(iC)=H_StateVec(i) !Streamflow at catchment (iC)
    enddo    
            
    do iC=1,nC
			i=i+1
			Hfl(iC)=H_StateVec(i) !Water depth at catchment (iC)
	enddo    
    
    do iC=1,nC
			i=i+1
			Yfl(iC)=H_StateVec(i) !Water elevation at catchment (iC)
    enddo        
    
    do iC=1,nC
			i=i+1
			Vol1(iC)=H_StateVec(i) !Previous Water storage at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			Vol2(iC)=H_StateVec(i) !Actual Water storage at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			Area2(iC)=H_StateVec(i) !Flooded Area at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			jtab(iC)=int(H_StateVec(i)) !Index of Volume level at catchment (iC)
    enddo
    
    do iC=1,nC
			i=i+1
			Q2face(iC)=H_StateVec(i) !Streamflow coming (or leaving) the catchment through connections 
    enddo
    
    do iC=1,nC
			i=i+1
			Q2viz(iC)=H_StateVec(i) !Updated flow for interconnections (signal can be positive or negative)
    enddo
       
    deallocate(H_StateVec)	
    close(FILHOTSTART_Inert)
    
	RETURN
	END