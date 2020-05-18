# MGB-Niger

Niger Source code with Write and/or Read initial conditions

! Processing

	- To read initial condition, in SIMULA.f90
	"if (1===1) then call flood_READHOT" 
	"if (0===1) then call flood_WHRITEHOT"
	
	- To write initial condition 	
	if (1===1) then call flood_WHRITEHOT" 
	if (0===1) then call flood_READHOT"
	
	- To write and read initial condition 
	"if (1===1) then call flood_WHRITEHOT" 
	"if (1===1) then call flood_READHOT

! Modification of SIMULA.f90 

	WRITE(FILLOG,*)'Calling initial conditions...'
	CALL CONDINIC  ! Initial Conditions
	
    if (0==1) then
            call flood_READHOT  
    end if
    WRITE(FILLOG,*)'Running main model routine...'
    CALL MODELO    ! Time Loop
    
	if (1==1) then
            call flood_WRITEHOT  
    end if
	WRITE(FILLOG,*)'Running main model routine...'
	CALL FOBJ 	   ! Error/objective Functions Evaluation
	
! Modification of VAR_MAIN.f90

    INTEGER,PARAMETER:: filhotstart_wb=127 !
    INTEGER,PARAMETER:: FILHOTSTART_Inert=128 !
	
	
! Create flood_WHRITEHOT.f90 : file with initials conditions 

		
! Create flood_READHOT.f90 : For futur simulation, read file with initials conditions (flood_WHRITEHOT)

		

		




