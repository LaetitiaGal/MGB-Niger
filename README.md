# MGB-Niger

    ! To write initial condition (in SIMULA.f90) 
       - if (0==1) then  call flood_READHOT  (line 139)
       - if (1==1) then  call flood_WRITEHOT  (line 147)
    
    ! To read initial condition (in SIMULA.f90)
       - if (1==1) then  call flood_READHOT  (line 139)
       - if (0==1) then  call flood_WRITEHOT  (line 147)
    
    ! To read AND write initial condition (in SIMULA.f90)
       - if (1==1) then  call flood_READHOT  (line 139)
       - if (1==1) then  call flood_WRITEHOT  (line 147)
    
    ! To not read AND not write initial condition (in SIMULA.f90)
       - if (0==1) then  call flood_READHOT  (line 139)
       - if (0==1) then  call flood_WRITEHOT  (line 147)
       
    !!! WARNING !!!!
       initial condition files : StateVars_inertial.hot   AND   StateVars_WB.hot, are created in Output folder but Read in Input folder
      
		

