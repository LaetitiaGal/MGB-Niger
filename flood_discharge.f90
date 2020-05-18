	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the flow for each catchment from Inertial equation   (Inertial version).
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

	subroutine flood_discharge

	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use VARS_INERC
	use VARS_MAIN
	implicit none
    INTEGER K,KHID 					!indexes and counters
    
    !$OMP PARALLEL NUM_THREADS(8)
    !$OMP DO PRIVATE(y1,y2,z1,z2,iCJus,hflow,dxflow,bflow,xMan,q0,Sflow,q)
	
    !-------------------------------------------------------------------------------------
    do iC=1,nC

        IB=IBAC(IC) 
		!IF((IB<SUBini).OR.(IB>SUBfim))CYCLE ! Manual controls for especialized calibration

        
        !SIMULATE ONLY RED FLOODS
        !IF(IBAC(IC)<5.AND.IBAC(IC)>7) cycle
        !IF(IBAC(IC).NE.4) cycle
        !IF(IBAC(IC)>4) cycle
        !IF((IBAC(IC)>4).AND.(IBAC(IC).NE.9)) cycle
        !IF(IBAC(IC)>9) cycle
        !IF(IBAC(IC).NE.9) cycle
       ! IF(IBAC(IC).NE.5) cycle
        !IF(IBAC(IC)<10.OR.IBAC(IC)>14) cycle   !modif cécile
!          IF((IBAC(IC)<10.OR.IBAC(IC)>14).AND.(IBAC(IC).NE.20)) cycle   !modif cécile

        
            ! Bottom level and water level of IC catchment:
            z1=ZTAB(1,iC)
            y1=Hfl(iC)+z1

            ! Bottom level and water level of downstream IC catchment:
            iCJus = CELJUS(iC)
            
            if(iCJus == -1)then
                
                !z2=z1-0.0005*dxart*1000.
                !y2=y1-0.0005*dxart*1000.
                
                z2=z1-0.0001*SRIO(IC)*1000. !assuming a 0.0001 m/m slope at the downstream boundary
                y2=y1-0.0001*SRIO(IC)*1000.
                
                !boundary condition for outflow in ocean
                !z2=z1
                !y2=0.0
            else
                z2=ZTAB(1,iCJus)
                y2=Hfl(iCJus)+z2
            endif
            
            ! Calculates the hflow variable:
            hflow=max(y2,y1)-max(z2,z1)
            hflow=max(hflow,0.0)
            
           
            if(iCJus /= -1)then
                dxflow=DBLE(SRIO(IC)*1000.) + DBLE(SRIO(iCJus)*1000.)       
                dxflow=dxflow/2.
            else
                !dxflow=DBLE(SRIO(IC)*1000.) + DBLE(dxart*1000.)
                dxflow=DBLE(SRIO(IC)*1000.)
                !dxflow=dxflow/2.
            endif

            !River width
            bflow=DBLE(BRIO(iC))
                      
            !River manning coefficient
            xMan=nMan(iC)
            
            ! Flow in the last time-step
            q0=Q2fl(iC)/bflow ! in m2/s
                    
            ! Water table slope:
            Sflow=-(y1-y2)/dxflow
            
            ! Calculates flow from Inertial equation (m2/s) for each IC:
            if (hflow>0.0) then
                q=(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*xMan*xMan*abs(q0)/(hflow**(10.0/3.0)))
                q=q*bflow !(m3/s)
            else
                q=0.0;
            endif
            
            ! Updates variable flow inertial:
            Q2fl(iC)=q ! in m3/s
            QJ2(iC)=Q2fl(iC) ! in m3/s
            
           
    enddo
    
    !$OMP END DO
    !$OMP END PARALLEL
    
        !Loop que calcula a vazão nas interconexões entre minibacias.
    if (1==1) then !liga conexões
    Q2viz=0.0
    !Q2face=0.0
    !$OMP PARALLEL DO NUM_THREADS(8) DEFAULT (SHARED)
    !PRIVATE ()
    do iFACE=1,nFACE
    
            KCAT=nFACECAT1(iFACE)
            KCAT2=nFACECAT2(iFACE)
            
            ! Nível de Fundo e Nível da Água da minibacia iC:
            z1=ZTAB(1,KCAT)
            !z1=nFACEY1(iFACE)
            !z1=z1-HRIO(KCAT) 
            y1=Hfl(KCAT)+z1
            z2=ZTAB(1,KCAT2)
            !z2=nFACEY2(iFACE)
            !z2=z2-HRIO(KCAT2) 
            y2=Hfl(KCAT2)+z2
            !if (iFace==122) then
            !write(*,*)iFACE,nFACEY1(iFACE),z1,HRIO(KCAT),nFACEY2(iFACE),z2,HRIO(KCAT2)
            !    write(*,*)z1,z2
                !write(*,*) KCAT,ZTAB(1,KCAT),ZTAB(1,KCAT2),z1,z2
            !pause
            !endif
            ! Cálculo da profundidade de escoamento:
            hflow=max(y2,y1)-max(z2,z1)
            !Limitador de fundo
!            hflow=hflow-1.0
            !Correção de valores negativos
            hflow=max(hflow,0.0)
            
            !A rotina DBLE transforma a variável de entrada em um real*8
            !Média dos dx de IC e ICJUS
            dxflow=DBLE(nFACEDX(iFACE))       !Verificar se precisa de um limitador do dx
            bflow=100.0
            !bflow=DBLE(BRIO(KCAT))
            !if (NORTHERN_DELTA(KCAT)==1.or.NORTHERN_DELTA(KCAT2)==1) bflow=200
            !if (SOUTHERN_DELTA(KCAT)==1.or.SOUTHERN_DELTA(KCAT2)==1) bflow=30
            if (NORTHERN_DELTA(KCAT)==1.or.NORTHERN_DELTA(KCAT2)==1) bflow=500
            if (SOUTHERN_DELTA(KCAT)==1.or.SOUTHERN_DELTA(KCAT2)==1) bflow=30
            
            !WIDTH FOR ESPECIFIC CONNECTIONS (E.G. RIVER DEFLUENCES)
            if(iFACE==230) bflow=600 !diaka distributary
            !if(iFACE==209) bflow=600 !diaka distributary
            xMan=nMan(KCAT)  
            !xMan=0.055 !modif Cécile - Manning pour le delta intérieur
               
            ! Vazão no tempo anterior:
            q0=Q2face(iFACE)/bflow ! em m2/s
                    
            ! Declividade da linha de água:
            Sflow=-(y1-y2)/dxflow
            
                
            ! Cálculo da vazão Inercial (por unidade de largura do rio) na face de jusante da minibacia iC:
            if (hflow>0.0) then
                q=(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*xMan*xMan*abs(q0)/(hflow**(10.0/3.0)))
                !q=(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*0*0*abs(q0)/(hflow**(10.0/3.0)))
                q=q*bflow
            else
                q=0.0;
            endif
            !write(*,*)(q0-(g*dtflood*hflow*Sflow))/(1.+g*dtflood*hflow*0.003*0.003*abs(q0)/(hflow**(10.0/3.0)))
            
            ! Calcula a nova vazão no próximo intervalo de tempo: EXCLUI FACES QUE NÃO ESTÃO LOCALIZADAS NO INLAND DELTA
            IF(NORTHERN_DELTA(KCAT)==1.or.NORTHERN_DELTA(KCAT2)==1.or.SOUTHERN_DELTA(KCAT)==1.or.SOUTHERN_DELTA(KCAT2)==1) THEN
            !IF(SOUTHERN_DELTA(KCAT)==1.or.SOUTHERN_DELTA(KCAT2)==1) THEN !DESABILITA NORTHERN_DELTA
            !IF(NORTHERN_DELTA(KCAT)==1.or.NORTHERN_DELTA(KCAT2)==1) THEN !DESABILITA SOUTHERN_DELTA
            !if(KCAT==364.OR.KCAT2==4155) THEN !iFACE==230 (DIAKA RIVER)
                Q2face(iFACE)=q ! em m3/s
                Q2viz(KCAT)=Q2viz(KCAT)- Q2face(iFACE)
                Q2viz(KCAT2)=Q2viz(KCAT2)+ Q2face(iFACE)
            ELSE
                Q2face(iFACE)=0
                Q2viz(KCAT)=0
                Q2viz(KCAT2)=0
            ENDIF
            
            
            DO K=1,NUMHIDG 				! store discharge in catchments defined in PARHIG.HIG
				KHID=IHIDG(K) 			! indexed catchment
				IF(KHID==KCAT) THEN
                     QRG_viz(K,IT)=Q2viz(KCAT)  
                ENDIF
            ENDDO
            

            
    enddo 
    !$OMP END PARALLEL DO
    endif

	endsubroutine
