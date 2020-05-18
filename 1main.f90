    PROGRAM MGB_Inercial
!*********************************************************************************
!***************MODELO HIDROLÓGICO DE GRANDES BACIAS (MGB-IPH)********************
!*********************************************************************************
!*		Large Scale Hydrology Research Group          	      *
!*	      Instituto de Pesquisas Hidráulicas (IPH)      	  *
!*	 Universidade Federal do Rio Grande do Sul (UFRGS)        *
!*                  	   June 2015                          *
!*                  	  Version 2.0                         *
!*                                                            *
!*********************************************************************************
!*********************************************************************************
!
!  PROGRAM MGB_Inercial (1main.f90) is the main routine from MGB-IPH.
!  Modelo de Grandes Bacias - Instituto de Pesquisas Hidraulicas
!
!---------------------------------------------------------------------------------
!  Opening:
!
!    This model and routines were developed by Walter Collischonn
!    during his research thesis in early 2001, in Brazil.
!
!    From there on MGB-IPH was modified for research and project purposes.
!    The current package is the official version of MGB-IPH.
!
!   As references we list below some publications where the model was used:
!
!   MELLER, A. ; COLLISCHONN, W. ; FAN, F. M. ; BUARQUE, D. C. ; PAIVA, R. C. D. ; DIAS, P. ; MOREIRA, D. . Previsão de Cheias por Conjunto em Curto Prazo. Revista Brasileira de Recursos Hídricos, v. 19, p. 33-49, 2014.
!   FAN, F. M.; COLLISCHONN, W. ; MELLER, A. ; BOTELHO, L. C. M. Ensemble streamflow forecasting experiments in a tropical basin: The São Francisco river case study. Journal of Hydrology (Amsterdam), v. SI, p. 10.1016/j.jhydr, 2014.
!	FAN, F. M.; COLLISCHONN, W. .Integração do Modelo MGB-IPH com Sistema de Informação Geográfica. Revista Brasileira de Recursos Hídricos, v. 19, p. 243-254, 2014.
!	Paiva, R. C. D. ; Collischonn, W. ; Buarque, D. C. . Validation of a full hydrodynamic model for large scale hydrologic modelling in the Amazon. Hydrological Processes, 27, p. 333–346. DOI: 10.1002/hyp.8425, 2013.
!	Paiva, R. C. D. ; Paiva, R. C. D. ; Collischonn, W. ; Bonnet, M.-P. ; De Gonçalves, L. G. G. ; Calmant, S. ; Getirana, A. ; Santos Da Silva, J. . Assimilating in situ and radar altimetry data into a large-scale hydrologic-hydrodynamic model for streamflow forecast in the Amazon. Hydrology and Earth System Sciences Discussions (Online), v. 10, p. 2879-2925, 2013.
!	Paiva, R. C. D.; Buarque, D. C. ; Collischonn, W. ; Bonnet, M.-P. ; Frappart, F.; Calmant, S.; Bulhões Mendes, C. A.. Large-scale hydrologic and hydrodynamic modeling of the Amazon River basin. Water Resources Research, v. 49, p. 1226-1243, 2013.
!	Collischonn, W. ; Meller, A. ; Fan, F. ; Moreira, D.S. ; Silva Dias, P.L. ; Buarque, D.; Bravo, J. M. . Short-term Ensemble Flood Forecasting Experiments in Brazil. Geophysical Research Abstracts, v. 15, p. 11910, 2013.
!	Bravo, J. M. ; Allasia, D. ; Paz, A. R. ; Collischonn, W. ; Tucci, C. E. M. . Coupled Hydrologic-Hydraulic Modeling of the Upper Paraguay River Basin. Journal of Hydrologic Engineering, v. 17, p. 635, 2012.
!	Paiva, R. C. D.; Collischonn, W.; Bonnet, M. P.; De Gonçalves, L. G. G. . On the sources of hydrological prediction uncertainty in the Amazon. Hydrology and Earth System Sciences, v. 16, p. 3127-3137, 2012.
!	Paiva, R. C. D.; Collischonn, W.; Buarque, D. C.. Validation of a full hydrodynamic model for large-scale hydrologic modelling in the Amazon. Hydrological Processes (Print), v. -, p. n/a-n/a, 2012.
!	Pontes, P. R. M. ; Collischonn, W. . Conservação de Volume em Modelos Simplificados de Propagação de Vazão. Revista Brasileira de Recursos Hídricos, v. 17, p. 83-96, 2012.
!	Sorribas, M. V.; Collischonn, W.; Marques, D. M.; Fragoso Jr., C. R.; Castro, N. M. R.; Souza, R. S. . Modelagem Distribuída do Carbono em Bacias Hidrográficas. Revista Brasileira de Recursos Hídricos, v. 17, p. 225-240, 2012.
!	Meller, A. ; Bravo, J. M. ; Collischonn, W. . Assimilação de Dados de Vazão na Previsão de Cheias em Tempo-Real com o Modelo Hidrológico MGB-IPH. Revista Brasileira de Recursos Hídricos, v. 17, p. 209-224, 2012.
!	Buarque, D. C. ; Collischonn, W. ; Paiva, R. C. D. . Coupling a basin erosion and river sediment transport model into a large scale hydrological model: an application in the Amazon basin. Geophysical Research Abstracts, v. 14, p. 11935, 2012.
!	Nóbrega, M. T. ; Collischonn, W. ; Tucci, C. E. M. ; Paz, A. R. . Uncertainty in climate change impacts on water resources in the Rio Grande Basin, Brazil. Hydrology and Earth System Sciences, v. 15, p. 585-595, 2011.
!	Paz, Adriano Rolim da ; Collischonn, Walter ; Tucci, Carlos E. M. ; Padovani, Carlos R. . Large-scale modelling of channel flow and floodplain inundation dynamics and its application to the Pantanal (Brazil). Hydrological Processes (Print), v. 25, p. 1498-1516, 2011.
!	Paiva, Rodrigo C.D. ; Collischonn, Walter ; Tucci, Carlos E.M. . Large scale hydrologic and hydrodynamic modeling using limited data and a GIS based approach. Journal of Hydrology (Amsterdam), v. 406, p. 170-181, 2011.
!	Collischonn, B. ; Paiva, R. C. D. ; Meirelles, F. S. C. ; Collischonn, W. ; Fan, F. M. ; Camano, E. . Modelagem Hidrológica de Uma Bacia com Uso Intensivo de Água: Caso do Rio Quaraí-RS. Revista Brasileira de Recursos Hídricos, v. 16, p. 119-133, 2011.
!	Paiva, R. C. D. ; Buarque, D. C. ; Collischonn, W. ; Sorribas, M. ; Allasia, D. G. ; Mendes, C. A. B. ; Tucci, C. E. M. ; Bonnet, M. P. . Hydrologic and Hydrodynamic Modelling of the Amazon Basin using TRMM Rainfall Estimates. Geophysical Research Abstracts, v. 13, p. 12666, 2011.
!	Paiva, R. ; Buarque, Diogo ; Collischonn, W. ; Sorribas, M. ; Allasia, D. G. P. ; Mendes, C. A. B. ; Tucci, C. E. M. ; Tucci, C. E. M. ; Bonnet, M. . Using TRMM rainfall estimates in hydrological and hydrodynamic modelling of the Amazon Basin. IAHS-AISH Publication, v. 343, p. 72-77, 2011.
!	Paz, A. R. ; Bravo, J. M. ; Allasia, D. ; Collischonn, W. ; Tucci, C. E. M. . Large-Scale Hydrodynamic Modeling of a Complex River Network and Floodplains. Journal of Hydrologic Engineering, v. 15, p. 152-165, 2010.
!	Paz, A. R. ; Collischonn, W. ; Tucci, C. E. M. . Simulação hidrológica de rios com grandes planícies de inundação. Revista Brasileira de Recursos Hídricos, v. 15, p. 31-43, 2010.
!	Paz, A. R. ; Collischonn, W. . Derivação de rede de drenagem a partir de dados do SRTM. Revista Geográfica Acadêmica, v. 2, p. 84-95, 2008.
!	Larentis, D. G. ; Collischonn, W. ; Tucci, C. E. M. . Simulação da qualidade de água em grandes bacias: Rio Taquari-Antas, RS. Revista Brasileira de Recursos Hídricos, v. 13, p. 5-22, 2008.
!	Tucci, Carlos E. M. ; Collischonn, W. ; Clarke, R. T. ; Clarke, Robin T. ; Paz, Adriano R. ; Allasia, Daniel . Short- and long-term flow forecasting in the Rio Grande watershed (Brazil). Atmospheric Science Letters, v. 9, p. 53-56, 2008.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Chou, S. C. ; Guilhon, L. G. ; Cataldi, M. ; Allasia, D. G. . Medium-range reservoir inflow predictions based on quantitative precipitation forecasts. Journal of Hydrology, v. 344, p. 112-122, 2007.
!	Collischonn, W. ; Allasia, D. G. ; Silva, B. C. ; Tucci, C. E. M. . The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Delgado, M. C. ; Silva, B. C. ; Collischonn, B. ; Allasia, D. G. ; Paz, A. R. . Modelo hidrológico distribuído para previsão de vazões incrementais na bacia do rio Paranaíba entre Itumbiara e São Simão. Revista Brasileira de Recursos Hídricos, v. 12, p. 43-56, 2007.
!	Silva, B. C. ; Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Delgado, M. C. . Previsão hidroclimática de vazão de curto prazo na bacia do rio São Francisco. Revista Brasileira de Recursos Hídricos, v. 12, p. 31-42, 2007.
!	Collischonn, W. ; Silva, B. C. ; Tucci, C. E. M. ; Allasia, D. G. . Large basin simulation experience in South America. IAHS Pubblication n. 303, v. 303, p. 360-370, 2006.
!	Silva, B. C. ; Tucci, C. E. M. ; Collischonn, W. . Previsão de vazão com modelos hidroclimáticos. Revista Brasileira de Recursos Hídricos, v. 11, p. 15-30, 2006.
!	Andreolli, I. ; Collischonn, W. ; Tucci, C. E. M. ; Haas, R. ; Regina, J. V. M. . Previsão de vazão afluente a um reservatório utilizando previsão quantitativa de chuva. Revista Brasileira de Recursos Hídricos, v. 11, p. 55-70, 2006.
!	Ribeiro Neto, A. ; Collischonn, W. ; Silva, R. C. V. ; Tucci, C. E. . Hydrological modelling in Amazonia use of the MGB-IPH model and alternative databases. IAHS-AISH Publication, v. 303, p. 246-254, 2006.
!	Collischonn, W. ; Tucci, C. E. M. ; Haas, R. ; Andreolli, I. . Forecasting river Uruguay flow using rainfall forecasts from a regional weather-prediction model. Journal of Hydrology, v. 305, p. 87-98, 2005.
!	Collischonn, W. ; Tucci, C. E. M. . Previsão Sazonal de vazão na bacia do rio Uruguai 1: Ajuste e verificação do modelo hidrológico distribuído. Revista Brasileira de Recursos Hídricos, v. 10, n.4, p. 43-59, 2005.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. ; Dias, P. L. S. ; Sampaio, G. O. . Previsão sazonal de vazão na bacia do rio Uruguai 2: Previsão Climática-Hidrológica. Revista Brasileira de Recursos Hídricos, v. 10, n.4, p. 60-72, 2005.
!	Bravo, J. M. ; Collischonn, W. ; Pilar, J. V. ; Silva, B. C. ; Tucci, C. E. M. . Operação de um reservatório com múltiplos usos com base na previsão de curto prazo. Revista Brasileira de Energia, v. 11, p. 85-110, 2005.
!	Collischonn, W. ; Tucci, C. E. M. ; Clarke, R. T. . Variabilidade temporal no regime hidrológico da bacia do rio Paraguai. Revista Brasileira de Recursos Hídricos, Porto Alegre RS, v. 8, n.1, p. 201-211, 2003.
!	Tucci, C. E. M. ; Dias, P. L. S. ; Clarke, R. T. ; Sampaio, G. O. ; Collischonn, W. . Long-term flow forecasts based on climate and hydrologic modeling: Uruguay river basin. Water Resources Research, New York, v. 39, n.7, p. 1-2, 2003.
!	Collischonn, W. ; Tucci, C. E. M. . Ajuste multiobjetivo dos parâmetros de um modelo hidrológico. Revista Brasileira de Recursos Hídricos, Porto Alegre, v. 8, n.3, p. 27-39, 2003.
!	Collischonn, W. ; Tucci, C. E. M. . Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n.2, 2001.
!   COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001 
!
!    Older source codes have been distributed by some people in
!    past research, thus any other Large Scale Model with similar source
!    codes are non-official and doesn't respect copyrights and patents.
!
!    This distribution is a version of the model in English.
!
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!  1main.f90
!---------------------------------------------------------------------------------
!  Discussion:
! 
!    This is the main routine of the model.
!    This routines starts the model, open some input files and calls
!    model simulation modes (i.e. simulation, calibration) etc...
!
!  Usage:
!
!    PROGRAM MGB_Inercial
!
!    uses modules, functions, and subroutines
!
!	*Module	PORTLIB     !library to calculate the processing time
!	*Module	IFPORT      !library to calculate the processing time
!	*Module	VARS_MAIN   !module containing the main variables used in the model
!	*Module	VARS_CALIB  !module containing the variables used in the model calibration
!	*Module	VARS_INERC  !module containing the variables of the model Local Inertial flow routing model
!	*Function	TEMPO(0,ISEED)  ! INICIA CONTAGEM DO TEMPO
!	*Subroutine	LEFIX 	!Suboutine that reads the model main input files with the running main information 
!	*Subroutine	ALLOCA_VARS(0) !Allocate the model main variables
!	*Subroutine	ALLOCA_VARSINERC(0) !Allocate the variables of the model Local Inertial flow routing model
!	*Subroutine	LEVAR !Subroutine for reading the file with fixed parameters (ALBEDOS, IAF, RS, Z; mounthly) (File: ALBIAF.FIX)
!	*Subroutine	LEUSO !Subroutine for reading the file with calibrated parameters (File: PARUSO.CAL)
!	*Subroutine	LECELL !Subroutine for reading the file with the model MiniBasins information (File: MINI.GTP)
!	*Subroutine	LEFLOODPLN !Subroutine for reading the file with a table containing a table of water-depth vs. area at each MiniBasin obtained by DEM preprocessing (File: COTA_AREA.FLP)
!	*Subroutine	LECLIMED !Subroutine for reading the file with mensal climate information (File: MEDIAS.CLM)
!	*Subroutine	ARQCLISUB	  !Subroutine for reading the files with daily/hourly meteorological information (Files with extension .CLI)
!	*Subroutine	LEQOBS	!Subroutine for reading the file observed streamflow data (File with extension .QOB)
!	*Subroutine	LESUBST !Then calls the subroutine for the reading the file with to be substituded streamflow data (File: QSUBST.QSB) 
!	*Subroutine	PARCEL !Calculates parameters related to MicroBains and rivers
!	*Subroutine	PARCUNGE !Calcualtes Muskingum-Cunge parameters
!	*Subroutine	flood_TOPO !Creates a matrix with topology information that is used in the Local Inertial routing method
!	*Subroutine	flood_TAB !Creates a table with the volume of water in the floodplain from the table of water depth vs. area obtained by DEM preprocessing (File with extension .FLP)
!	*Subroutine	SIMULA !Subroutine that calls a standard model simulation run
!	*Subroutine	LeCalib !Subroutine that reads automatic calibration input files
!	*Subroutine	ALLOCA_CALIB(0) !Subroutine that allocates automatic calibration input files
!	*Subroutine	CALIBRA !Subroutine that alls a model automatic calibration run
!	*Subroutine	ALLOCA_CALIB(1) !Subroutine that cleans automatic calibration input files
!	*Subroutine	PREVISAO !Subroutine that calls a forecasting run
!	*Function	TEMPO(1,0)  !Finish time count
!	*Subroutine	ALLOCA_VARS(1) !Cleans main variables
!
!	 opens
!
!    * CHUVAbin.pbi, the precipitation binary file
!    * NOSOLO.tx, information about soil water at one minibasin
!    * AJUSTE.fob, model performance coeficcients
!
!    reads
!
!    * no files are read in this routine
!
!    creates
!
!    * NOSOLO.tx, information about soil water at one minibasin
!    * AJUSTE.fob, model performance coeficcients
!
!---------------------------------------------------------------------------------
!  Licensing:
!
!    This code is distributed under the...
!
!  Version/Modified: 
!
!    21 June 2015
!    By: Fernando Mainardi Fan
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
!  Main References:
!  COLLISCHONN, W. ; ALLASIA, D. G. ; SILVA, B. C. ; TUCCI, C. E. M. The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
!  COLLISCHONN, W., TUCCI, C. E. M. Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n. 2, 2001.
!  COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001 
!
!---------------------------------------------------------------------------------
!  Variables and Parameters:
!  *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------

!  Calls used libraries and modules 

	USE PORTLIB     !library to calculate the processing time
	USE IFPORT      !library to calculate the processing time
	USE VARS_MAIN   !module containing the main variables used in the model
	USE VARS_CALIB  !module containing the variables used in the model calibration
	USE VARS_INERC  !module containing the variables of the model Local Inertial flow routing model

!  Variables declaration 

	IMPLICIT NONE !States that non declared variables cannot be used
	INTEGER JULDAY !Variable used for day of the year calculation on a Julian calendar
	REAL TIME_SECONDS !Variable to calculate the processing time

!  Statements used to calculate the processing time

	CALL TEMPO(0,ISEED)  ! INICIA CONTAGEM DO TEMPO
	TIME_SECONDS = SECNDS(0.0)	

!   Input and Output Directories and log file
    !INPUT_DIRECTORY = TRIM(adjustl(('C:\MGB\Input\')))
    !OUTPUT_DIRECTORY = TRIM(adjustl(('C:\MGB\Output\')))
    INPUT_DIRECTORY = TRIM(adjustl(('.\Input\')))
    OUTPUT_DIRECTORY = TRIM(adjustl(('.\Output\')))
	OPEN(FILLOG,FILE=OUTPUT_DIRECTORY // 'LOG_MGB.txt',STATUS='UNKNOWN') !Output file with simulation LOG

!  Write the model name in the black box
    WRITE(*,*)'Modelo Hidrologico de Grandes Bacias: MGB-IPH'
    WRITE(FILLOG,*)'Modelo Hidrologico de Grandes Bacias: MGB-IPH'

!  Call subroutines related to input files that are readen only once reading and MGB-IPH variables allocation
	WRITE(*,*)'Reading infoMGB.sim file...'
	WRITE(FILLOG,*)'Reading infoMGB.sim file...'
	CALL LEFIX 	!Suboutine that reads the model main input files with the running main information 
	WRITE(*,*)'Allocating variables...'
	WRITE(FILLOG,*)'Allocating variables...'
	CALL ALLOCA_VARS(0) !Allocate the model main variables
	CALL ALLOCA_VARSINERC(0) !Allocate the variables of the model Local Inertial flow routing model
    	! WRITE(*,*)'ALOCOU VARIAVEIS PRINCIPAIS' !Warning that the model main variables were allocated
	IDINI=JULDAY(IMES,IDIA,IANO) !Calculate the first day of the simulation in the Julian calendar
	WRITE(*,*)'Reading fixed parameters file...'
	WRITE(FILLOG,*)'Reading fixed parameters file...'
	CALL LEVAR !Subroutine for reading the file with fixed parameters (ALBEDOS, IAF, RS, Z; mounthly) (File: ALBIAF.FIX)
	WRITE(*,*)'Reading calibrated parameters file...'
	WRITE(FILLOG,*)'Reading calibrated parameters file...'
	CALL LEUSO !Subroutine for reading the file with calibrated parameters (File: PARUSO.CAL)
	WRITE(*,*)'Reading model topology (MINI.GTP) file...'
	WRITE(FILLOG,*)'Reading calibrated parameters file...'
	CALL LECELL !Subroutine for reading the file with the model catchments (minibasins) information (File: MINI.GTP)
	
	WRITE(*,*)'Reading observed flows, substituted flows, and climate files...'
	WRITE(FILLOG,*)'Reading observed flows, substituted flows, and climate files...'
	CALL LECLIMED !Subroutine for reading the file with mensal climate information (File: MEDIAS.CLM)
	WRITE(FILLOG,*)'CLIMED OK'
	CALL ARQCLISUB	  !Subroutine for reading the files with daily/hourly meteorological information (Files with extension .CLI)
	WRITE(FILLOG,*)'CLISUB OK'
	CALL LEQOBS	!Subroutine for reading the file observed streamflow data (File with extension .QOB)
	WRITE(FILLOG,*)'QOBS OK'
	IF(NUMSUBST>0)THEN !If it was identified in the LEFIX routine that streamflow substitution is being done
		CALL LESUBST !Then calls the subroutine for the reading the file with to be substituded streamflow data (File: QSUBST.QSB) 
	ENDIF 
	WRITE(FILLOG,*)'QSUBS OK'

!*****************************************************************
	! Hiydrodynamic modeling:
	hdFLAG=1 !Use 0 to deactivate Hydrodynamic modeling; Use 1 to activate.
	hdFLAG0=sum(hdFLAG) ! Test
	If(hdFLAG(1)==1)THEN
	CALL LEFLOODPLN !Subroutine for reading the file with a table containing a table of water-depth vs. area at each MiniBasin obtained by DEM preprocessing (File: COTA_AREA.FLP)
	ENDIF
!*****************************************************************

! Delimiters
    SUBini=0
    SUBfim=99999999999999
! Enf of delimiters

!  End of subroutines related to input files that are readen only once reading and MGB-IPH variables allocation

!  Call subroutines that calculates MicroBasins, rivers channels, and floodplains characteristics
    WRITE(*,*)'Calculating river geometry informations...'
    WRITE(FILLOG,*)'Calculating river geometry informations...'
	CALL PARCEL !Calculates parameters related to MicroBains and rivers
	CALL PARCUNGE !Calcualtes Muskingum-Cunge parameters
	CALL flood_TOPO !Creates a matrix with topology information that is used in the Local Inertial routing method
	CALL flood_TAB !Creates a table with the volume of water in the floodplain from the table of water depth vs. area obtained by DEM preprocessing (File with extension .FLP)

!  End of subroutines that calculates MicroBasins, rivers channels, and floodplains characteristics

!  Openning of MGB-IPH Input and Output files that are read and written during the simulation
    WRITE(*,*)'Openning input/output files...'
    WRITE(FILLOG,*)'Openning input/output files...'
	OPEN(FILPLU,FILE=INPUT_DIRECTORY // 'CHUVAbin.pbi',STATUS='old',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') !Input binary file with interpolated rainfall data
	OPEN(FILSOL,FILE=OUTPUT_DIRECTORY // 'NOSOLO.txt',STATUS='UNKNOWN') !Output file with soil data from one single definied MicroBasin
	OPEN(FILAJU,FILE=OUTPUT_DIRECTORY // 'AJUSTE.fob',STATUS='UNKNOWN') !Output file with objective functions valuess for the locations with observed data

!  End of openning MGB-IPH main Input and Output files that are read and written over the simulation

!  Ask for verification
!  Informs the model general setup
    WRITE(FILLOG,*)'Input data reading process ended'
    WRITE(FILLOG,*)''
    WRITE(FILLOG,'(A30,I8)')'Number of catchments:',NC
    WRITE(FILLOG,'(A30,I8)')'Number of subbasins:',NB
    WRITE(FILLOG,'(A30,I8)')'Number of HRUs:',NU
    WRITE(FILLOG,'(A30,4I5)')'Start date (dd mm yyyy hh):',IDIA,IMES,IANO,HORAINI 
    WRITE(FILLOG,'(A30,F16.2)')'Time step (s):',DTP
    WRITE(FILLOG,'(A30,I8)')'Number of time steps:',NT
    WRITE(FILLOG,'(A30,I8)')'Number of observation gauges:',NOBS
    WRITE(FILLOG,'(A30,I8)')'Number of substitution gauges:',NUMSUBST
    WRITE(FILLOG,'(A30,I8)')'Number of climate gauges:',NCLI
    WRITE(FILLOG,*)'ICALIB (simulation=0; autocalibration=1; forecasting=2)'
    WRITE(FILLOG,'(A30,I8)')'ICALIB :',ICALIB
    WRITE(FILLOG,*)'PLEASE VERIFY INPUT INFORMATION AND PRESS ENTER TO START RUN...'
    WRITE(FILLOG,*)''

    WRITE(*,*)'Input data reading process ended'
    WRITE(*,*)''
    WRITE(*,'(A30,I8)')'Number of catchments:',NC
    WRITE(*,'(A30,I8)')'Number of subbasins:',NB
    WRITE(*,'(A30,I8)')'Number of HRUs:',NU
    WRITE(*,'(A30,4I5)')'Start date (dd mm yyyy hh):',IDIA,IMES,IANO,HORAINI 
    WRITE(*,'(A30,F16.2)')'Time step (s):',DTP
    WRITE(*,'(A30,I8)')'Number of time steps:',NT
    WRITE(*,'(A30,I8)')'Number of observation gauges:',NOBS
    WRITE(*,'(A30,I8)')'Number of substitution gauges:',NUMSUBST
    WRITE(*,'(A30,I8)')'Number of climate gauges:',NCLI
    WRITE(*,*)'ICALIB (simulation=0; autocalibration=1; forecasting=2)'
    WRITE(*,'(A30,I8)')'ICALIB :',ICALIB
    WRITE(*,*)'PLEASE VERIFY INPUT INFORMATION AND PRESS ENTER TO START RUN...'
    WRITE(*,*)''
    !PAUSE

!  Select the mode of MGB-IPH run: Simulation (ICALIB=0), Automatic Calibration (ICALIB=1), Forecasting (ICALIB=2)

	CALIBRA_CASE: SELECT CASE (ICALIB) !Verify the mode
	CASE (0) !Standard model simulation run
		WRITE(*,*)'SIMULATION MODE' !Warns about the chosen mode
		WRITE(FILLOG,*)'SIMULATION MODE' !Warns about the chosen mode
		CALL SIMULA !Subroutine that calls a standard model simulation run
	CASE (1) !Automatic Calibration run
		WRITE(*,*)'AUTOMATIC CALIBRATION MODE' !Warns about the chosen mode
		WRITE(FILLOG,*)'AUTOMATIC CALIBRATION MODE' !Warns about the chosen mode
		CALL LeCalib !Subroutine that reads automatic calibration input files
		CALL ALLOCA_CALIB(0) !Subroutine that allocates automatic calibration input files
		CALL CALIBRA !Subroutine that alls a model automatic calibration run
		CALL ALLOCA_CALIB(1) !Subroutine that cleans automatic calibration input files
	CASE (2) !Forecasting run
		WRITE(*,*)'FORECASTING MODE' !Warns about the chosen mode
		WRITE(FILLOG,*)'FORECASTING MODE' !Warns about the chosen mode
		CALL PREVISAO !Subroutine that calls a forecasting run
	CASE DEFAULT
		PAUSE ' ERROR: UNKNOWM ICALIB!!!' !Warns about unknown ICALIB's
	END SELECT CALIBRA_CASE

!  End of MGB-IPH

!  After running the model, finish the time count, close all files, and clean all variables 

	CALL TEMPO(1,0)  !Finish time count
	CLOSE (FILPLU) !Close input rainfall binary file
	CLOSE (FILSOL) !Close output file with soil data from one single definied MicroBasin
	CLOSE (FILAJU) !Close output file with objective functions valuess for the locations with observed data
	CALL ALLOCA_VARS(1) !Cleans main variables
	
! End of closings and cleanings

! Inform the user that the model run is over and informations about total processing time

	WRITE(*,*)'MGB-IPH ENDED SUCCESSFULLY'
	WRITE(FILLOG,*)'MGB-IPH ENDED SUCCESSFULLY'
	print *, char(7)
	TIME_SECONDS = SECNDS(TIME_SECONDS)
	WRITE(*,*)'TOTAL TIME:',TIME_SECONDS/60.,'MINUTES'
	WRITE(*,*)'PRESS ENTER TO CLOSE...'
	WRITE(FILLOG,*)'TOTAL TIME:',TIME_SECONDS/60.,'MINUTES'
	WRITE(FILLOG,*)'PRESS ENTER TO CLOSE...'
	CLOSE(FILLOG)
	!PAUSE
	!READ(*,*)
	
	END         
