	SUBROUTINE ALLOCA_VARS(IOP)

!---------------------------------------------------------------------------------
! ALLOCA_VARS.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This routine allocate or deallocate the MGB-IPH memory for the main variables
!
! Usage:
!
! CALL ALLOCA_VARS(IOP)
!
! uses modules, functions, and subroutines
!
! *Module VARS_MAIN !module with the model main variables
!
! opens
!
! * no files are opened in this routine
!
! reads
!
! * no files are read in this routine
!
! creates
!
! * no files are created in this routine
!
!---------------------------------------------------------------------------------
! Licensing:
!
! This code is distributed under the...
!
! Version/Modified:
!
! 2015.06.21 - 21 June 2016
! By: Fernando Mainardi Fan
!
! Authors:
!
! Original fortran version by Walter Collischonn
! Present fortran version by:
! * Walter Collischonn
! * Rodrigo Cauduro Dias de Paiva
! * Diogo da Costa Buarque
! * Paulo Pontes Rógenes
! * Mino Viana Sorribas
! * Fernando Mainardi Fan
! * Juan Martin Bravo
!
! Main References:
! COLLISCHONN, W. ; ALLASIA, D. G. ; SILVA, B. C. ; TUCCI, C. E. M. The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
! COLLISCHONN, W., TUCCI, C. E. M. Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n. 2, 2001.
! COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------

! Variables and Modules

	USE VARS_MAIN
	!use AUX_MOD !RP
	use vars_calib !RP
	IMPLICIT NONE	
	INTEGER IOP
    	SAVE

! Select the case how the subroutine is working: for allocation, or for deallocation.

	ALLOCV_CASE: SELECT CASE (IOP) !Checks for allocation or deallocation. IOP=0 means allocation, and IOP=1 means deallocation.
	CASE (0) ! Allocation
		ALLOCATE (QTUDO(NC,NT),HTUDO(NC,NT),YTUDO(NC,NT),QVIZTUDO(NC,NT)) !Streamflow at all minibasins flow at all time steps
		ALLOCATE (DIAH(NT),MESH(NT),ANOH(NT),HORAH(NT)) !Date and time corresponding to the time interval counter
		ALLOCATE (ICODMUSK(NC)) !Code that indicates linear or non-linear
		ALLOCATE (BPLAN(NC)) !Width of the floodplain (includes own river)
		ALLOCATE (HCALHA1(NC),HCALHA2(NC)) !Depth that starts the floodplain and that this floodplain is totally flooded
		ALLOCATE (QMUP(NC,NTMUP),AMUP(NC,NTMUP),BMUP(NC,NTMUP),CMUP(NC,NTMUP)) !Non-linear Muskingum table
		ALLOCATE (EVQ(NC)) !Direct evaporation of the liquid surface of the cell in m3 / s
		ALLOCATE (QRIOINI(NC,NUMUSK+1)) !Initial condition of muskingum Cunge rounting
		ALLOCATE (QCONTORM(NC+1,25),QCONTORJ(NC+1,25)) !Streamflow boundary condition for Muskingum Cunge
		ALLOCATE (SI(NC,NU))
		ALLOCATE (XYC(NCLI,2))!x and y coordinates of climatological stations
		ALLOCATE (TAMM(NCLI,12),URMM(NCLI,12),VVMM(NCLI,12))
		ALLOCATE (PAMM(NCLI,12),SOLMM(NCLI,12))
		ALLOCATE (PMB(NB)) !Average rainfall at each sub-basin
		ALLOCATE (KPM(NB)) !Auxiliar variable
		
		!if (flagaclimed/=0) then !#MV When using CRU file of South America (or not use time serie), do not allocate!
		 ALLOCATE (TD(NCLI,NT),UD(NCLI,NT),VD(NCLI,NT),SOLD(NCLI,NT))
		 ALLOCATE (PAD(NCLI,NT))
		!end if

		ALLOCATE (QREF(NC)) !Reference flow
		ALLOCATE (ICBOM(NC))!Number (id) of the nearest climatological station to the minibasin center
	 	ALLOCATE (PUSO(NC,NU))!Proportion of uses in minibasin
	 	ALLOCATE (ACEL(NC),ACUR(NC),SRIO(NC),DECL(NC))
		ALLOCATE (X(NC),Y(NC))!Coordinates of the cell center
		ALLOCATE (IBAC(NC),HCEL(NC),LCEL(NC),CELJUS(NC))!SubBasin,HMAX,HMIN,Downstream minibasin
	 	ALLOCATE (PM(NC)) !Average rainfall and rainfall at the timestep in the minibasin
		ALLOCATE (BRIO(NC))!River width
		ALLOCATE (IEXUT(NB)) !Indicate the outlets of the basin
		ALLOCATE (CI(NB),CB(NB),CS(NB)) !Pparameters of propagation in cell
		ALLOCATE (CIOLD(NB),CSOLD(NB)) !Pparameters of propagation in cell
		ALLOCATE (PLAM(NB,NU),CAP(NB,NU),WC(NB,NU))
		ALLOCATE (WMOLD(NB,NU),BOLD(NB,NU),KBOLD(NB,NU),KIOLD(NB,NU))
		ALLOCATE (WM(NB,NU),B(NB,NU),KINS(NB,NU),KBAS(NB,NU))
		ALLOCATE (NARQ(NCLI))
		ALLOCATE (NSUBT(NC)) !Number of stretchs for MUSKINGUM CUNGE routing (diferent at each minibasin)
		ALLOCATE (DT(NC)) !Time steps for MUSKINGUM CUNGE (diferent at each minibasin)
		ALLOCATE (CEL(NC),TIND(NC))
		ALLOCATE (ALB(NU,12),RIAF(NU,12),Z(NU,12),RS(NU,12))
		ALLOCATE (QESP(NB))	 !Specific baseflow (M3/S/KM2)
		ALLOCATE (QB(NT),QBI(NT),QBIS(NT))!Streamflow according to its origin
		ALLOCATE (DSUM(NB,NT),DINM(NB,NT),DBAM(NB,NT))!Averages over the SubBasins
		ALLOCATE (SIM(NB,NT),EM(NB,NT),WBM(NB,NT),PM2(NB,NT))!Averages over the SubBasins
		ALLOCATE (KCB(NB)) !Number of minibasins inside each subbasin
		ALLOCATE (QR(NOBS,NT),QOBS(NOBS,NT))!Calculated and Observed streamflow in the outlet of each microbasin
		ALLOCATE (QLIDO(NUMSUBST,NT)) !Flow that is read to replace calculated
		ALLOCATE (QMANTEIGA(NT),QSAOROM(NT),QSAOFRA(NT),QPMCRUZ(NT),QMANGA(NT),QCARI(NT))
		ALLOCATE (QBJLAPA(NT),QGAMEL(NT),QMORP(NT))
		ALLOCATE (QRG(NUMHIDG,NT)) !Hydrograph for recording
		ALLOCATE (QRB(NB,NT)) !Hydrographs of sub-basins
		ALLOCATE (QM1(NC+1),QJ1(NC),QM2(NC+1),QJ2(NC)) !Streamflows upstream and downstream in i
		ALLOCATE (QCEL1(NC),QCEL2(NC)) !Original streamflows at minibasins at the instants t E t+1 at minibasin i
		ALLOCATE (QCEL1PREV(NC),QCEL2PREV(NC)) !Originated streamflows at minibasins at the instants t E t+1 at minibasin i
		ALLOCATE (PMB2(NC+1),PMI2(NC+1),PMS2(NC+1),PJB2(NC+1),PJI2(NC+1),PJS2(NC+1)) !Proportions source of streamflow to the river
		ALLOCATE (VRB(NC),VRI(NC),VRS(NC)) !Volumes of proportions source of streamflow to the river
		ALLOCATE (ET(NC,NU)) 	!Total evapotranspiration
		ALLOCATE (CAF(NC,NU))		!Upward capillary flow
		ALLOCATE (W(NC,NU)) 	!Amount of water in the soil
		ALLOCATE (QBAS(NC),QINT(NC),QSUP(NC))		!Flows at the minibasins
		ALLOCATE (VBAS(NC),VINT(NC),VSUP(NC))		!Volume at the minibasin
		ALLOCATE (TONTEM(NC))		!Previous day temperature
		ALLOCATE (P(NC)) !Rainfall at the interval in the cell
		ALLOCATE (TA(NC),UR(NC),VV(NC),SOL(NC),PA(NC))	!Temperature, Umidity, Wind, Insolation, Pressure
		ALLOCATE (R2(NOBS),ERRV(NOBS),R2L(NOBS)) !COEF. NASH-SUTCLIFFE, VOLUME ERROR, COEF. NASH LOGARITHMS
		ALLOCATE (WPREV(NC,NU)) 	!Amoint of water in the soil at the beginning of the forecasting cycle
		ALLOCATE (VBASPREV(NC),VINTPREV(NC),VSUPPREV(NC)) !Volume of water in the minibasin at the beginning of the forecasting cycle
		ALLOCATE (TAPREV(NC)) !Temperature at the beginning of the forecasting cycle
		ALLOCATE (QM2PREV(NC+1),QJ2PREV(NC),QJ1PREV(NC)) !Streamflows at the beginning of the forecasting cycle
		ALLOCATE (QPREV(NUMHIDG,100)) !Forecasted streamfows at interest points with multiple lead-times
		ALLOCATE (SIPREV(NC,NU)) 	!Intercepeted amount of water (forecating cycle)
		ALLOCATE (QRIOINIPREV(NC,NUMUSK)) !Initial condition for MUSKINGUM CUNGE at forecast cycle
		ALLOCATE (QCONTORMPREV(NC+1,25),QCONTORJPREV(NC+1,25)) !Boundary condition for MUSKINGUM CUNGE at forecast cycle
		ALLOCATE (PMB2PREV(NC+1),PMI2PREV(NC+1),PMS2PREV(NC+1)) !Proportions of origin of the flows in the river
		ALLOCATE (PJB2PREV(NC+1),PJI2PREV(NC+1),PJS2PREV(NC+1)) !Proportions of origin of the flows in the river

		!allocate (OD(NC),hdFLAG(NC)) !RP
		allocate (OD(NC)) !Inertial modification
 
		ALLOCATE (CBOLD(NB),PLAMOLD(NB,NU),CAPOLD(NB,NU),WCOLD(NB,NU)) ! RP
		allocate (sbtFLAG(nC))
		
		ALLOCATE (PPREV(72,NC))

		! Variables to calculate indicators of maximum evapotranspiration:
		allocate (E0agua(nC),E0TOPO(nC),E0SUP(nC))

        !FMF 21/06/2015
        ALLOCATE (ATIVABACIA(NB)) !flag to activate sub-basins or not
        
        !FMF 21/06/2015        
        ALLOCATE (QITUDO(NC,NT)) !incremental discharges everytime and everywhere
	    ALLOCATE (EVAPTUDO(NC,NT)) !et everytime and everywhere
	    ALLOCATE (QBASTUDO(NC,NT)) !catchment baseflows everytime and everywhere
		ALLOCATE (WTUDO(NC,NT)) !soil moisture volme everytime and everywhere
		
		!FMF 09/09/2015 
		ALLOCATE (RUGMAN(NC)) !Manning
        
        !WHOLE DELTA
        ALLOCATE(EVAPTUDO_NIGER(NT)) !ET FOR INLAND DELTA NIGER
        ALLOCATE(EVAPTUDO_NIGER_INT(NT)) !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        ALLOCATE(EVAPTUDO_NIGER_SOLO(NT)) !ET FOR INLAND DELTA NIGER ON SOIL
        ALLOCATE(EVAPTUDO_NIGER_FLOOD(NT)) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(PRECTUDO_NIGER(NT)) !PREC FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(PRECTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(EVAPTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(PRECTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(EVAPTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(ETPOTENCIAL(NT)) !AVALIAÇÃO DE ET POTENCIAL
        
        !SOUTHERN DELTA
        ALLOCATE(S_EVAPTUDO_NIGER(NT)) !ET FOR INLAND DELTA NIGER
        ALLOCATE(S_EVAPTUDO_NIGER_INT(NT)) !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        ALLOCATE(S_EVAPTUDO_NIGER_SOLO(NT)) !ET FOR INLAND DELTA NIGER ON SOIL
        ALLOCATE(S_EVAPTUDO_NIGER_FLOOD(NT)) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(S_PRECTUDO_NIGER(NT)) !PREC FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(S_PRECTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(S_EVAPTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(S_PRECTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(S_EVAPTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        
        !NORTHERN DELTA
        ALLOCATE(N_EVAPTUDO_NIGER(NT)) !ET FOR INLAND DELTA NIGER
        ALLOCATE(N_EVAPTUDO_NIGER_INT(NT)) !ET FOR INLAND DELTA NIGER ON CANOPY INTERCEPTATION
        ALLOCATE(N_EVAPTUDO_NIGER_SOLO(NT)) !ET FOR INLAND DELTA NIGER ON SOIL
        ALLOCATE(N_EVAPTUDO_NIGER_FLOOD(NT)) !ET FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(N_PRECTUDO_NIGER(NT)) !PREC FOR INLAND DELTA NIGER ON FLOODPLAINS
        ALLOCATE(N_PRECTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(N_EVAPTUDO_NIGER_MES(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(N_PRECTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        ALLOCATE(N_EVAPTUDO_NIGER_ANO(NT/25))  !PREC mensal FOR INLAND DELTA NIGER
        
    CASE (1) ! Deallocation
	    DEALLOCATE (QTUDO,HTUDO) !Streamflow at all minibasins flow at all time steps
	    DEALLOCATE (DIAH,MESH,ANOH,HORAH) !Date and time corresponding to the time interval counter
		DEALLOCATE (ICODMUSK) !Code that indicates linear or non-linear
		DEALLOCATE (BPLAN) !Width of the floodplain (includes own river)
		DEALLOCATE (HCALHA1,HCALHA2) !Depth that starts the floodplain and that this floodplain is totally flooded
		DEALLOCATE (QMUP,AMUP,BMUP,CMUP) !Non-linear Muskingum table
		DEALLOCATE (QMANTEIGA,QSAOROM,QSAOFRA,QPMCRUZ,QMANGA,QCARI)
		DEALLOCATE (QBJLAPA,QGAMEL,QMORP)
		DEALLOCATE (EVQ) !Direct evaporation of the liquid surface of the cell in m3 / s
		DEALLOCATE (QCONTORM,QCONTORJ) !Boundary condition of muskingum Cunge rounting
		DEALLOCATE (WMOLD,BOLD,KBOLD,KIOLD)
		DEALLOCATE (XYC)!x and y coordinates of climatological stations
		DEALLOCATE (TAMM,URMM,VVMM)
		DEALLOCATE (PAMM,SOLMM)
		DEALLOCATE (PMB,KPM) !Average rainfall at each sub-basin
		
		!if (flagaclimed/=0) then !#MV When using CRU file of South America (or not use time serie), do not allocate!
		DEALLOCATE (TD,UD,VD,SOLD)
		!end if

		DEALLOCATE (PAD)
		DEALLOCATE (QREF) !Reference flow
		DEALLOCATE (ICBOM)!Number (id) of the nearest climatological station to the minibasin center
	 	DEALLOCATE (PUSO)!Proportion of uses in minibasin
	 	DEALLOCATE (ACEL,ACUR,SRIO,DECL)
		DEALLOCATE (X,Y)!Coordinates of the cell center
		DEALLOCATE (IBAC,HCEL,LCEL,CELJUS)!SubBasin,HMAX,HMIN,Downstream minibasin
	 	DEALLOCATE (PM) !Average rainfall and rainfall at the timestep in the minibasin
		DEALLOCATE (BRIO)!River width
		DEALLOCATE (IEXUT) !Indicate the outlets of the basin
		DEALLOCATE (CI,CB,CS) !Parameters of propagation in cell
		DEALLOCATE (CIOLD,CSOLD) !Parameters of propagation in cell
		DEALLOCATE (PLAM,CAP,WC)
		DEALLOCATE (WM,B,KINS,KBAS)
		DEALLOCATE (NARQ)
		DEALLOCATE (NSUBT,DT)
		DEALLOCATE (CEL,TIND)
		DEALLOCATE (ALB,RIAF,Z,RS)
		DEALLOCATE (QESP)	 !Specific baseflow (M3/S/KM2)
		DEALLOCATE (QB,QBI,QBIS)!Streamflow according to its origin
		DEALLOCATE (DSUM,DINM,DBAM)!Averages over the SubBasins
		DEALLOCATE (SIM,EM,WBM,PM2)!Averages over the SubBasins
		DEALLOCATE (KCB) !Number of minibasins inside each subbasin
		DEALLOCATE (QR,QOBS)!Calculated and Observed streamflow in the outlet of each microbasin
		DEALLOCATE (SI)
		DEALLOCATE (QRG) 	!Streamflow n the outlets of microbasins
		DEALLOCATE (QRB) !Hydrograph for recording
		DEALLOCATE (QM1,QJ1,QM2,QJ2) !Streamflows upstream and downstream in i
		DEALLOCATE (PMB2,PMI2,PMS2,PJB2,PJI2,PJS2) !Proportions source of streamflow to the river
		DEALLOCATE (VRB,VRI,VRS) !Volumes of proportions source of streamflow to the riverO
		DEALLOCATE (ET) 	!Total evapotranspiration
		DEALLOCATE (CAF)		!Upward capillary flow
		DEALLOCATE (W) 	!Amount of water in the soil
		DEALLOCATE (QBAS,QINT,QSUP)		!Flows at the minibasin
		DEALLOCATE (VBAS,VINT,VSUP)		!Volume at the minibasin
		DEALLOCATE (TONTEM)		!Previous day temperature
		DEALLOCATE (P) !Rainfall at the interval in the cell
		DEALLOCATE (TA,UR,VV,SOL,PA)	!Temperature, Umidity, Wind, Insolation, Pressure
		DEALLOCATE (R2,ERRV,R2L) !COEF. NASH-SUTCLIFFE, VOLUME ERROR, COEF. NASH LOGARITHMS
		DEALLOCATE (QCEL1,QCEL2) !Originated streamflows at minibasins at the instants t E t+1 at minibasin i
		DEALLOCATE (QLIDO) !Calculated streamflow for substitution
		DEALLOCATE (PMB2PREV,PMI2PREV,PMS2PREV) !Proportions of origin of the flows in the river
		DEALLOCATE (PJB2PREV,PJI2PREV,PJS2PREV) !Proportions of origin of the flows in the river

		deallocate (OD) !RP
		!deallocate (OD,hdFLAG) !Inertial modification
		
		DEALLOCATE (PPREV)
		
		!FMF 21/06/2015
        DEALLOCATE (ATIVABACIA) !flag to activate subbasins or not
        
        !FMF 21/06/2015
        DEALLOCATE (QITUDO) !incremental discharges everytime and everywhere
	    DEALLOCATE (EVAPTUDO) !et everytime and everywhere
	    DEALLOCATE (QBASTUDO) !catchment baseflows everytime and everywhere
		DEALLOCATE (WTUDO) !soil moisture volme everytime and everywhere
		
		!FMF 09/09/2015 
		DEALLOCATE (RUGMAN) !Manning
		
	CASE DEFAULT
		STOP ' ERROR: IOP UNKNOWN AT ROUTINE ALLOCA_CALIB!!!'
	END SELECT ALLOCV_CASE

	! End of the routine.
	END
