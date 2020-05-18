!---------------------------------------------------------------------------------
!  Discussion:
! 
!    This sub-routine reads climate files.
!
!    uses modules and functions
!
!    * module VARS_MAIN in VARSMAIN.f90
!
!	 opens
!
!    * ARQCLI(KLI),  the climatic data file
!
!    reads
!
!    * ARQCLI(KLI),  the climatic data file
!
!    creates
!
!    * no files are created in this routine
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

SUBROUTINE ARQCLISUB

!  Variables and Parameters:
USE VARS_MAIN
IMPLICIT NONE
INTEGER KB,KLI,K
!*************************************************************
!Verify
!Checks the unit-catchment in each sub-catchment
KCB(1)=IEXUT(1)
DO KB=2,NB
    KCB(KB)=IEXUT(KB)-IEXUT(KB-1)
ENDDO

!IF flagclimed==1 uses CRU data.
if (flagaclimed==1) then
    TD(:,:)=-9999.0
    UD(:,:)=-9999.0
    VD(:,:)=-9999.0
    SOLD(:,:)=-9999.0
    PAD(:,:)=-9999.0
else

    !Opens daily climatic data files.
    DO KLI=1,NCLI
        NARQ(KLI)=KLI+30
        OPEN(NARQ(KLI),FILE=INPUT_DIRECTORY // ''//ARQCLI(KLI),STATUS='OLD')
        READ(NARQ(KLI),733)(CABE(K),K=1,4)
    ENDDO
    DO KLI=1,NCLI
        DO IT=1,NT
            TD(KLI,IT) = -1.0
            UD(KLI,IT) = -1.0
            VD(KLI,IT) = -1.0
            SOLD(KLI,IT) = -1.0
            PAD(KLI,IT) = -1.0
        ENDDO
    ENDDO
    DO KLI=1,NCLI
        NARQ(KLI)=KLI+30
        CLOSE (NARQ(KLI))
    ENDDO
endif
733	FORMAT(4A10)
734	FORMAT(5F10.2)
!__________________________________________________________
RETURN
END
