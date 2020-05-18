# MGB-Niger

Niger Source code with Write and/or Read initial conditions

! Processing  - SIMULA.f90

	- To read initial condition
	"if (1===1) then call flood_READHOT" 
	"if (0===1) then call flood_WHRITEHOT"
	
	- To write initial condition 	
	if (1===1) then call flood_WHRITEHOT" 
	if (0===1) then call flood_READHOT"
	
	- To write and read initial condition 
	"if (1===1) then call flood_WHRITEHOT" 
	"if (1===1) then call flood_READHOT
