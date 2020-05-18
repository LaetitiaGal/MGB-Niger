    !*********************************************************************************
    !
    !  SUBROUTINE CALIBRA is the routine for hydrological parameters calibration    
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine calculates optimal parameters for MGB-IPH using the MOCOM-UA 
    !    optimization method. 
    !
    !    Calibration parameters are read before in LECALIB.F90
    !
    !    MOCOM-UA is applied for a vector of i=1,ns samples of parameters.
    !    
    !
    !    Parameters and objective-functions are stored in vectors PAR and FOBJ and  
    !    written in evolution.txt file.
    !
    !
    !  Usage:
    !
    !    CALL CALIBRA
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90
    !    * module     VARS_CALIB  in      VARS_CALIB.f90
    !    * subroutine SEMENTE     in      SEMENTE.f90
    !    * function   RAN1        in      RAN1.f90
    !    * subroutine CALIBPARAM  in      CALIBPARAM.F90
    !    * subroutine CONDINIC    in      CONDINIC.F90
    !    * subroutine MODELO      in      MODELO.F90
    !    * subroutine FOBJ        in      FOBJ.F90    
    !
    !	 opens
    !
    !    * Population0.TXT, contains initial parameters and objective-function for i=1,NS 
    !    population samples. attention: this is optional, manual check if(1==1)
    !
    !    reads
    !
    !    * PARHIG.HIG, the main simulation settings
    !
    !    creates
    !
    !    * EVOLUTION.TXT, the parameters and objective function file.
    !    
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the ...
    !
    !  Version/Modified: 
    !
    !    2014.17.11 - 17 November 2014 (By: Mino V. Sorribas)    
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
    
	
	SUBROUTINE CALIBRA
	USE VARS_MAIN
	USE VARS_CALIB
	IMPLICIT NONE

    ! Variables declaration

	INTEGER I,L,J,KB,K,IW,KS,IFUNC,IPAR,JF,IBOBO,IRMAX  !Variable indexes/counter
	REAL    XBOBO,RAN1                                  !Random number gen
	INTEGER ICONG,ITESTE                                !Not used
	INTEGER i1,i2,i3                                    !Not used
	INTEGER NPCOMPLEX                                   !Number of Parameter Complexes
	
	
	! Start Optimization parameters	
	NPCOMPLEX=5


    ! Opens Output Evolution.txt file
	OPEN(FILEVO,FILE=OUTPUT_DIRECTORY // 'EVOLUTION.TXT',STATUS='UNKNOWN')


    ! Generate random number seed	
	CALL SEMENTE(ISEED)
	IBOBO=MOD(ISEED,100)
	DO IW=1,IBOBO
		XBOBO=RAN1(ISEED)   !Generates random number between 0-1
	ENDDO

	! Saves original parameters values
	WMOLD=WM
	BOLD=B
	KIOLD=KINS
	KBOLD=KBAS
	PLAMOLD=PLAM
	CAPOLD=CAP
	WCOLD=WC
	CSOLD=CS
	CIOLD=CI
	CBOLD=CB

    ! Initializes parameters set for samples
	PAR=1.0
	if (1==1) then
		DO I=2,NS                                               !NS is the number os simulations to be done = number of samples
			DO L=1,NPAR
				PPAR(L,I)=RAN1(ISEED)
				PAR(L,I)=PMIN(L)+PPAR(L,I)*(PMAX(L)-PMIN(L))    !Generates i-th parameters set based on uniform distribution between min and max
			ENDDO
		ENDDO
	else
		! Read initial population from external file
		open(123456,FILE=INPUT_DIRECTORY // 'Population0.txt')

		DO I=1,NS
			read(123456,*) (PAR(L,I),L=1,NPAR) 	
			DO L=1,NPAR
				PPAR(L,I)=RAN1(ISEED)
				PAR(L,I)=PMIN(L)+PPAR(L,I)*(PMAX(L)-PMIN(L))

			ENDDO
		ENDDO
	endif
	
	
	KOUNTF=0
	write(*,*) 'Generates initial population:'	
	! Evaluate function for each initial point in sample
	DO I=1,NS
                
		PARX=PAR(:,I)
		call CalibParam		
		CALL CONDINIC
		CALL MODELO		
		CALL FOBJ
		KOUNTF=KOUNTF+1     !Counts how many times the function was evaluated (number of simulations)
		
		!Store function values
		DO JF=1,NF
			FO(I,JF)=VFO(JF)
		ENDDO
		write(*,'(A9,I5,A1,I3)') 'Individual',I,'/',NS
		WRITE(*,*)'F.O.=',(FO(I,JF),JF=1,NF)
	ENDDO	


    ! Starts optimization method
	FMIN=999999999999999999.9
	ISHUFFLE=0
	IRMAX=10
	
	DO WHILE(IRMAX.GT.1.and.ISHUFFLE<iMaxGen) ! Iterates while there is more than one hierarchical category     		
		                                      ! and less than the max number of generations
		
		ISHUFFLE=ISHUFFLE+1
		IRUIM=0

		! Initialize Solution Ranking
		IPARET=1
		ISOMA=0
		IRMAX=1
3000	DO I=1,NS
			IF(IPARET(I).EQ.IRMAX)THEN
				DO J=1,NS
					IF(J.NE.I.AND.IPARET(J).EQ.IRMAX)THEN	
						DO K=1,NF
							IF(FO(I,K).GT.FO(J,K))THEN
								IDOMIN=1
								IDD=J
							ELSE
								IDOMIN=0
								GOTO 1000
							ENDIF
						ENDDO
						IF(IDOMIN.EQ.1)THEN !This point is dominated
							GOTO 2000 
						ENDIF
1000				ENDIF
				ENDDO
2000			IF(IDOMIN.EQ.1)THEN
					!WRITE(2,*)I,X(I),' DOMINADO POR',IDD
					IPARET(I)=IPARET(I)+1
					ISOMA=1
				ELSE
					!WRITE(2,*)I,X(I),' NAO DOMINADO'
					IPARET(I)=IRMAX
				ENDIF
				IDOMIN=0
				IDD=0
			ENDIF
		ENDDO
		IF(ISOMA.EQ.1)THEN
			IRMAX=IRMAX+1
			ISOMA=0
			GOTO 3000
		ENDIF
		
		IF(IRMAX==1)EXIT ! When all sets are at level 1 calibration is ended
		
		!---------------------------------------------
		WRITE(FILEVO,*)
		WRITE(FILEVO,*)ISHUFFLE
		DO I=1,NS		   ! Save points in Evolution.txt
		  WRITE(FILEVO,172)(PAR(L,I),L=1,NPAR),(FO(I,JF),JF=1,NF),IPARET(I)
		ENDDO
172		FORMAT(<NPAR>F10.3,<NF>F12.4,I8)
		!---------------------------------------------

		! At this point there is NS sample points, which have been evaluated and ranked in IPARET
		! IPARET(I)=1 indicates the best points; IPARET(I)=IRMAX indicates the worst points.
		NPC=0
		DO I=1,NS
			IF(IPARET(I).EQ.IRMAX)THEN
				NPC=NPC+1               !Count the number of points that have ranking equal to the worst of them IRMAX
				IRUIM(NPC)=I            !Store the position of worst points
			ENDIF
		ENDDO

		! Sets probability of each one of NS points to be chosen constitue a simplex
		SPARET=0.0
		
		RMAX=IRMAX
		DO I=1,NS
			SPARET=SPARET+IPARET(I)
		ENDDO
		DO I=1,NS
			PROB(I)=(RMAX-IPARET(I)+1.)/(NS*(RMAX+1.)-SPARET)
		ENDDO
		DO I=2,NS
			PROB(I)=PROB(I-1)+PROB(I) !Acumulaed probability
		ENDDO
		IF(ABS(PROB(NS)-1.0).GT.0.001)STOP 'Something wrong with the probability'
		PROB(NS)=1.0 !Fix rounding errors if necessary

		! SPAR contais NPCOMPLEX points of SIMPLEX 1...NPC, defined by NPAR parameters
		!ALLOCATE (SPAR(NPC,NPAR+1,NPAR),FPLEX(NPC,NPAR+1,NF))
		ALLOCATE (SPAR(NPC,NPCOMPLEX,NPAR),FPLEX(NPC,NPCOMPLEX,NF))
        
        !Evolution of each complex
		DO K=1,NPC
			LRUIM=IRUIM(K)
			DO IPAR=1,NPAR
				SPAR(K,1,IPAR)=PAR(IPAR,LRUIM) !Store worst point
			ENDDO

			!DO J=2,NPAR+1  
			DO J=2,NPCOMPLEX        !Chooses the other points to each complex
				XRAND=RAN1(ISEED)
				PINI=0.0
				DO I=1,NS
					IF(XRAND.GE.PINI.AND.XRAND.LE.PROB(I))THEN
						! I is choosen
						DO IPAR=1,NPAR
							SPAR(K,J,IPAR)=PAR(IPAR,I)	     !Stores parameters values
							DO IFUNC=1,NF
								FPLEX(K,J,IFUNC)=FO(I,IFUNC) !Stores function values
							ENDDO
						ENDDO
					ENDIF
					PINI=PROB(I)
				ENDDO
			ENDDO
		ENDDO

107		FORMAT(I4,7F8.2)
        
        !At this point simplexes have been chosen!		
        
        ! The routine has a sample of NS points in space with dimension NPAR
        ! The objective function was evaluated and a  mutltiobjective ranking was established
    	! The sample was separated into NPC complexes, each one with NS+1 points
    	
    	
    	! Starts the  Competitive Complex Evolution		
		DO K=1,NPC              !Evolve each one of the simplexes
		
			LRUIM=IRUIM(K)
			SOMAPAR=0.0
			IREJECT=0
			! Calculates centroid of best points (average of point coordinates)
			!DO J=2,NPAR+1  
			DO J=2,NPCOMPLEX 
				DO L=1,NPAR
					SOMAPAR(L)=SOMAPAR(L)+SPAR(K,J,L)
				ENDDO
			ENDDO
			!SOMAPAR=SOMAPAR/(NPAR)
			SOMAPAR=SOMAPAR/(NPCOMPLEX-1)
			
			!write(*,*)' complexo:',k
			!write(*,*)' pior ponto: ',(spar(k,1,l),l=1,npar)
			!write(*,*)' centróide: ',(somapar(l),l=1,npar)
	
			DO L=1,NPAR
				DIFPAR=SOMAPAR(L)-SPAR(K,1,L)   !distance of the worst point to the centroid
				REFLEX(L)=SOMAPAR(L)+DIFPAR	    !reflection point coordinates
				CONTRA(L)=SOMAPAR(L)-DIFPAR/2.  !contraction point coordinates
			ENDDO
			!write(*,*)' reflex: ',(reflex(l),l=1,npar)
			
			! Check if point generated by reflex is inside valid region for parameters
			DO L=1,NPAR
				IF(REFLEX(L).LT.PMIN(L).OR.REFLEX(L).GE.PMAX(L))THEN
					write(*,*)' reject reflex on parameter',l
					!write(*,*)' contra: ',(contra(l),l=1,npar)
					IREJECT=1
					GOTO 4000
				ENDIF
			ENDDO


			PARX=REFLEX         ! Loads reflex parameters
			CALL CalibParam     ! Sets MGB-IPH parameters
			CALL CONDINIC       ! MGB-IPH Initial Condition
			CALL MODELO         ! Runs MGC-IPH			
			
			CALL FOBJ           ! Evaluate objective-function
			KOUNTF=KOUNTF+1     ! NUmber of simulations
			
			! Stores objective-function values
			DO JF=1,NF
				FO(LRUIM,JF)=VFO(JF)
			ENDDO
			
			! Checks if reflex point is dominated and accepts if it is not.
			!DO J=2,NPAR+1
			DO J=2,NPCOMPLEX !checks if j-th points dominates
				IDOMIN=0
				DO IFUNC=1,NF
					IF(FPLEX(K,J,IFUNC).LT.FO(LRUIM,IFUNC))THEN
						IDOMIN=IDOMIN+1	    !counts number of points where the function for the new point is worst
					ENDIF
				ENDDO
				IF(IDOMIN.EQ.NF)THEN        !if new point is the worst, rejects
					WRITE(*,*)'PONTO REJEITADO POR',(FPLEX(K,J,IFUNC),IFUNC=1,NF)
					IREJECT=1
					GOTO 4000
				ENDIF
			ENDDO
4000		IF(IREJECT.EQ.0)THEN
				! if it is here than reflex is acecepted
				WRITE(*,*)'PONTO ACEITO'
				DO IPAR=1,NPAR
					PAR(IPAR,LRUIM)=REFLEX(IPAR)
				ENDDO
	
			ELSE
				! If it is here, reflex was rejected, then uses contraction

				PARX=CONTRA          ! Loads contraction parameteres
				CALL CalibParam
				CALL CONDINIC
				CALL MODELO				
				CALL FOBJ 
				KOUNTF=KOUNTF+1				
				DO JF=1,NF
					FO(LRUIM,JF)=VFO(JF)
				ENDDO
	
				DO IPAR=1,NPAR
					PAR(IPAR,LRUIM)=CONTRA(IPAR)
				ENDDO
			ENDIF

		ENDDO !End of simplexes loop
	
		DEALLOCATE (SPAR,FPLEX)
	
		VMIN=99999999999990.0
		DO KS=1,NS
			DO JF=1,NF
				VMIN(JF)=MIN(VMIN(JF),FO(KS,JF))
			ENDDO
		ENDDO

		WRITE(*,79)ISHUFFLE,IRMAX,NPC,(VMIN(JF),JF=1,NF)
79		FORMAT(3I8,8F12.4)


	ENDDO !End of shuffle-evolution loop - only when irmax = 1
	!************************************************************************************************

	RETURN
	END
