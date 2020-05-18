	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine is the main subroutine of Inertial model  (Inertial version).
    !
    !
    ! Usage: flood_timestep, flood_continuity and flood_discharge
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

    SUBROUTINE flood_inercial
    
    ! Variables and parameters
    USE VARS_MAIN
    USE VARS_INERC
    use IFPORT
    
    implicit none
    REAL(8) elapsed_time 
    
    
    !-------------------------------------------------------------------------------------
    tflood=0.0
    AFLTUDO(1,IT)=0.0 !inicializa a variável de areas inundadas no inland delta
    
    do while (tflood<dtflood0)
        call flood_timestep
        dtflood=min(dtflood,dtflood0-tflood)
        tflood=tflood+dtflood
        !write(*,*) 'Flood inundation iT=',iT, 100.*tflood/dtflood0,'%',', dt=',dtflood,' s'
        ! Inertial equation:
        call flood_discharge
        ! Continuity equation:
        call flood_continuity 
          
    enddo

    !Write results for each IC in IT timestep
	!if(IT==883) call flood_write 
 !   if(IT==975) call flood_write 
    !if(IT==1024) call flood_write
 !    if(IT==1097) call flood_write
    
	return
	endsubroutine
    
    
    