	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine writes the results for eatch catchment  (Inertial version).
    !
    !
    ! Usage:
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
    ! * output files
    !
    ! reads
    !
    ! * no files are created in this routine
    !
    ! creates
    !
    ! * output files
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
	
	subroutine flood_write
	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use VARS_INERC
	use VARS_MAIN
	implicit none
	character(5) trash
	!-------------------------------------------------------------------------------------

    
   ! Writing results for catchments
    if (int(real(iT)/1.0)==real(iT)/1.0.or.iT==1) then

	    ! Open output files:
	    if (iT<10) then
	        write(trash,'(A4,I1)') '0000', iT
	    elseif (iT<100) then
	        write(trash,'(A3,I2)') '000', iT
	    elseif (iT<1000) then
	        write(trash,'(A2,I3)') '00', iT
	    elseif (iT<10000) then
	        write(trash,'(A1,I4)') '0', iT
	    else
	        write(trash,'(I5)') iT
	    endif
    	
	    trash=trim(trash)
	    !open(13131,FILE=OUTPUT_DIRECTORY //'flood\Q2fl_'//trash//'.txt')
	    !open(13133,FILE=OUTPUT_DIRECTORY //'flood\Yfl_'//trash//'.txt')
	    open(13134,FILE=OUTPUT_DIRECTORY //'flood\HflHAND_'//trash//'.txt')
	    !open(13135,FILE=OUTPUT_DIRECTORY //'face\Q2face_'//trash//'.txt')
	    !open(13136,FILE=OUTPUT_DIRECTORY //'flood\Afl_'//trash//'.txt')

        ! Write results:
        do iC=1,nC
		 !   write(13131,'(F16.6)') (Q2fl(iC))
		    !write(13133,'(F16.6)') (Yfl(iC))
            write(13134,'(F16.6)')  (max(Hfl(iC)-HRIO(iC),0.0))
            
		    !write(13136,'(F16.6)') (Area2(iC))
	    enddo
	  !do iFACE=1,nFACE
	  !write(13135,'(F16.6)') (Q2face(iFACE))
	   !enddo	    
	  

        ! Close files:
        !close(13131)
        close(13133)
        !close(13134)
        !close(13135)
        !close(13136)
        
    endif

	endsubroutine
