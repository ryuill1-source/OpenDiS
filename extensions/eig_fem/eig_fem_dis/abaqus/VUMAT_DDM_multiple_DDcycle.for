!DIR$ FREEFORM
!====================================================================================
!   Developed by Taejoon Park 
!	July 01, 2019.
!
!   Modified by Cuong Nguyen (02/02/2020)
!   Modification: added subroutine VUAMP to control the amplitude of applied stress
!                 modified VEXTERNALDB to check if the Hold Stress flag exist
!   Modified by Cuong Nguyen (10/20/2022)
!   Modification: allow multiple DD cycles for each DDM step
!====================================================================================
! PROPERTIES
! PROPS(1)= EMOD : Young's modulus
! PROPS(2)= ENU : Poisson's ratio
!------------------------------------------------------------------------------------
! Solution dependent variables
! STATEOLD(1)= Crystal_COORDS(1) : Coord(1) in the crystal coordinate system.
! STATEOLD(2)= Crystal_COORDS(1) : Coord(2) in the crystal coordinate system.
! STATEOLD(3)= Crystal_COORDS(2) : Coord(3) in the crystal coordinate system.
! STATEOLD(4~12)= P(1,1), P(1,2), ... , P(3,2), P(3,3) : Rotation tensor which 
! 		rotates the components in the material coord system into crystal coord system.
! STATEOLD(13) : Parameter to check whether the integration point is inside the 
!		activation volume during the current step or not.
! STATEOLD(14) : Parameter to check how many dislocation line segments generated
!		 plastic strain to the integration point.
! STATEOLD(15~20) : (total) Plastic strain increment in the crystal coordinate system.
!		11, 22, 33, 12, 23, 31
! STATEOLD(21~26) : Stress component in the crystal coordinate system.
!		11, 22, 33, 12, 23, 31
! STATEOLD(27~29) : Plastic spin increment in the crystal coordinate system.
!		27: w12, 28: w13, 29: w23
! STATEOLD(30) : Accumulated equivalent plastic strain
! STATEOLD(31) : Dislocation density
! STATEOLD(32) : Accumulated pe in z-direction
!
!
!	COMMON /sharedf/ de_p, t_step, swept_normal, swept_center, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8, EffSlip, dw_p!,rn1_old, rn2_old, rn1_new, rn2_new
!	COMMON /sharedi/ n_segment, k_iteration
! rn_id:    COL(1) node ID starts from 0, COL(2) flag, COL(3) num. neighbors	
! rn_nodes: COL(1-3): x, y, z coordinates
! node_id : rn id of all segments (large id followed by smaller id)
!====================================================================================
SUBROUTINE VUMAT(&
!  READ ONLY (UNMODIFIABLE) VARIABLES -
            jBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,&
            STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,&
            PROPS, DENSITY, STRAININC, RELSPININC,&
            TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,&
            STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,&
            TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,&
!  WRITE ONLY (MODIFIABLE) VARIABLES -
            STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
!
!    INCLUDE 'VABA_PARAM.INC'
	IMPLICIT DOUBLE PRECISION (a-h,o-z)
	INTEGER, PARAMETER :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length
!	The above two lines are needed only for ABAQUS 2018! 
	DIMENSION JBLOCK(*), PROPS(NPROPS),DENSITY(*), COORDMP(*),&
		CHARLENGTH(*), STRAININC(*),&
		RELSPININC(*), TEMPOLD(*),&
		STRETCHOLD(*),&
		DEFGRADOLD(*),&
		FIELDOLD(*), STRESSOLD(*),&
		STATEOLD(*), ENERINTERNOLD(*),&
		ENERINELASOLD(*), TEMPNEW(*),&
		STRETCHNEW(*),&
		DEFGRADNEW(*),&
		FIELDNEW(*),&
		STRESSNEW(*), STATENEW(*),&
		ENERINTERNNEW(*), ENERINELASNEW(*)
    CHARACTER*80 CMNAME
	INTEGER, PARAMETER :: I_UMT_NBLOCK=1, I_UMT_NPT=2, I_UMT_LAYER=3,I_UMT_KSPT=4, I_UMT_NOEL=5
!
	INTEGER, PARAMETER :: n_segment_max=80000
	INTEGER:: n_segment, k_iteration
	DOUBLE PRECISION :: EffSlip
	DOUBLE PRECISION:: t_step, t_nextstep
	! t_step: total time increment in the current step
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: rn1_old, rn2_old, rn1_new, rn2_new
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: swept_normal, swept_center
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: dw_p
	DOUBLE PRECISION, DIMENSION(n_segment_max,6):: de_p
	COMMON /sharedf/ de_p, dw_p, t_step, t_nextstep, EffSlip, rn1_old, rn2_old, rn1_new, rn2_new!, swept_normal, swept_center, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	COMMON /sharedi/ n_segment, k_iteration
!
	CALL VUMATXTRARG (JBLOCK(I_UMT_NBLOCK),&
		NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,&
		STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,&
		PROPS, DENSITY, STRAININC, RELSPININC,&
		TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,&
		STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,&
		TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,&
		STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW,&
		JBLOCK(I_UMT_NOEL), JBLOCK(I_UMT_NPT),&
		JBLOCK(I_UMT_LAYER), JBLOCK(I_UMT_KSPT))
	RETURN
END SUBROUTINE VUMAT
!
! Subroutine VUAMP added 2/2/2020
SUBROUTINE VUAMP(&
! passed in for information and state variables
			ampName, time, ampValueOld, dt, nprops, props, nSvars,& 
			svars, lFlagsInfo, nSensor, sensorValues, sensorNames,& 
			jSensorLookUpTable,& 
! to be defined
			ampValueNew,& 
			lFlagsDefine,&
			AmpDerivative, AmpSecDerivative, AmpIncIntegral)
!      
	IMPLICIT DOUBLE PRECISION (a-h,o-z)
	INTEGER, PARAMETER :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length
! Above two lines needed for ABAQUS 2018
! svars - additional state variables, similar to (V)UEL
    DIMENSION sensorValues(*), svars(*)
    DIMENSION props(*)
	DIMENSION jSensorLookUpTable(*)
	CHARACTER*80 sensorNames(nSensor)
    CHARACTER*80 ampName
!	CHARACTER(*), PARAMETER :: filename_stress_flag = "HoldStress"
! time indices
    INTEGER, PARAMETER :: iStepTime = 1  
	INTEGER, PARAMETER :: iTotalTime = 2
	INTEGER, PARAMETER :: nTime = 2
! flags passed in for information
    INTEGER, PARAMETER :: iInitialization = 1
    INTEGER, PARAMETER :: iRegularInc = 2
    INTEGER, PARAMETER :: ikStep = 3
    INTEGER, PARAMETER :: nFlagsInfo = 3
! optional flags to be defined
    INTEGER, PARAMETER :: iComputeDeriv = 1
    INTEGER, PARAMETER :: iComputeSecDeriv = 2
    INTEGER, PARAMETER :: iComputeInteg = 3
    INTEGER, PARAMETER :: iStopAnalysis = 4
    INTEGER, PARAMETER :: iConcludeStep = 5
    INTEGER, PARAMETER :: nFlagsDefine = 5
!	
	INTEGER :: LoadType ! 1: dynamic (synchronization) 2: quasi-static
	DIMENSION time(nTime), lFlagsInfo(nFlagsInfo), lFlagsDefine(nFlagsDefine)
	DOUBLE PRECISION :: tim
	DOUBLE PRECISION :: ampValueNew
	DOUBLE PRECISION :: amp_increment
	DOUBLE PRECISION :: paradis_nextdt, FEM_Du
!	DOUBLE PRECISION :: time_1, time_2, time_3, time_4, time_5
!	LOGICAL HOLD_ONE, HOLD_TWO
! HOLD_STRESS determined by VEXTERNALDB	
	COMMON /sharednextdt/ LoadType, FEM_Du
!	
	IF (lFlagsInfo(iInitialization) .EQ. 1) THEN
		ampValueNew = ampValueOld
	ELSE
		IF (LoadType.eq.1) THEN
			amp_increment=FEM_Du*1.0D0
		ELSEIF (LoadType.eq.2) THEN
			amp_increment=FEM_Du*0.01D0
		END IF
		tim = dble(time(iStepTime)) 
		IF (tim .le. 0.0) THEN
			ampValueNew = ampValueOld+0.004D0
        ELSE
            ampValueNew = ampValueOld+amp_increment
		END IF		

		! write(*,*) tim, ampValueOld, ampValueNew

		if (ampValueNew > 70.d0) then
			ampValueNew = ampValueOld
		endif


		
	END IF
!
	RETURN
END SUBROUTINE VUAMP
!
SUBROUTINE VUMATXTRARG(&
!	READ ONLY (UNMODIFIABLE) VARIABLES -
			NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,&
			STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,&
			PROPS, DENSITY, STRAININC, RELSPININC,&
			TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,&
			STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,&
			TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,&
!	WRITE ONLY (MODIFIABLE) VARIABLES -
			STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW,&
!	READ ONLY EXTRA ARGUMENTS -
			NELEMENT, NMATPOINT, NLAYER, NSECPOINT)
!
!    INCLUDE 'VABA_PARAM.INC'
	IMPLICIT DOUBLE PRECISION (a-h,o-z)
	INTEGER, PARAMETER :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length
!	The above two lines are needed only for ABAQUS 2018! 		
    DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),&
            CHARLENGTH(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR),&
            RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),&
            STRETCHOLD(NBLOCK,NDIR+NSHR),&
            DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),&
            FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),&
            STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),&
            ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),&
            STRETCHNEW(NBLOCK,NDIR+NSHR),&
            DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),&
            FIELDNEW(NBLOCK,NFIELDV),&
            STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),&
            ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK),&
			NELEMENT(NBLOCK)
!
!  
	INTEGER, DIMENSION(200000) :: Element_id
    DOUBLE PRECISION, PARAMETER :: PI=acos(-1.0_16)
    DOUBLE PRECISION, DIMENSION(3) :: Crystal_COORDS, domega
    DOUBLE PRECISION, DIMENSION(6) :: dfyds, dstranp
    DOUBLE PRECISION, DIMENSION(3,3) :: dStrainGlobal, dStrainLocal, StressGlobal, StressLocal, Pnew, P, dstranp33
    DOUBLE PRECISION, DIMENSION(3,3,3,3) :: CE
    DOUBLE PRECISION :: EMOD, ENU, c11, c12, c44
    INTEGER:: NBLOCK, iBlock, k1, k2, k3, k4, i, j, k, l, LENOUTDIR
	CHARACTER(*), PARAMETER :: filename_stress_head = "Output_stress"
	CHARACTER(*), PARAMETER :: filename_surfaceEl_head = "Surface_stress"
	CHARACTER(LEN=1024) :: foldername, filename_output, filename_stress, filename_surface
	INTEGER :: KPROCESSNUM
!
	INTEGER, PARAMETER :: n_segment_max=80000
	INTEGER:: n_segment, k_iteration, num_elements
	DOUBLE PRECISION :: EffSlip
	DOUBLE PRECISION:: t_step, t_nextstep
	! t_step: total time increment in the current step
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: rn1_old, rn2_old, rn1_new, rn2_new
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: swept_normal, swept_center
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	INTEGER, PARAMETER :: n_elements_max=600000
	INTEGER, DIMENSION (n_elements_max,2):: surface_elements
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: dw_p
	DOUBLE PRECISION, DIMENSION(n_segment_max,6):: de_p
	COMMON /sharedf/ de_p, dw_p, t_step, t_nextstep, EffSlip, rn1_old, rn2_old, rn1_new, rn2_new!, swept_normal, swept_center, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	COMMON /sharedi/ n_segment, k_iteration
	COMMON /sharedsurface/ surface_elements, num_elements
	COMMON /elem_id/ Element_id
!    
    SAVE
!====================================================================================
!   Input Material properties
    EMOD=PROPS(1)
    ENU=PROPS(2)
    IF(ENU.GT.0.4999.AND.ENU.LT.0.5001) ENU=0.499
!
!====================================================================================
!
	c11=EMOD*(1.0-ENU)/(1.0-2.0*ENU)/(1.0+ENU)
	c12=EMOD*ENU/(1.0-2.0*ENU)/(1.0+ENU)
	c44=EMOD/2.0/(1.0+ENU)
!	Elastic stiffness tensor in the crystal coordinate system.
    Call E_modulus(c11,c12,c44,CE)
!==================================== ====================================
!                             DO  Nblock loop
!==================================== ====================================
    DO iblock=1,NBLOCK
!-------------------------------------------------------------------------
!	Define Initial state variables
		IF (TOTALTIME.EQ.0.0) THEN
			! initial positions at the crystal coordinate system.
			! (Initially, the crystal coord system is aligned to mat coord system.)
			!STATEOLD(iblock,1)=COORDMP(iblock,1)
			!STATEOLD(iblock,2)=COORDMP(iblock,2)
			!STATEOLD(iblock,3)=COORDMP(iblock,3)
			!
			! Initial rotation of the crystal coord system (aligned to mat coord sysrtem)
			!
			!
			!STATEOLD(iblock,4)=1.0 !P(1,1)
			!STATEOLD(iblock,5)=0.0 !P(1,2)
			!STATEOLD(iblock,6)=0.0 !P(1,3)
			!STATEOLD(iblock,7)=0.0 !P(2,1)
			!STATEOLD(iblock,8)=1.0 !P(2,2)
			!STATEOLD(iblock,9)=0.0 !P(2,3)
			!STATEOLD(iblock,10)=0.0 !P(3,1)
			!STATEOLD(iblock,11)=0.0 !P(3,2)
			!STATEOLD(iblock,12)=1.0 !P(3,3)
		END IF
			!
			!
!	End initial state variables
!-------------------------------------------------------------------------
!	Reset parameter when new step starts
!
!
!
!-------------------------------------------------------------------------
		
!	Obtain the rotation matrix P
		P(1,1)=STATEOLD(iblock,4)
		P(1,2)=STATEOLD(iblock,5)
		P(1,3)=STATEOLD(iblock,6)
		P(2,1)=STATEOLD(iblock,7)
		P(2,2)=STATEOLD(iblock,8)
		P(2,3)=STATEOLD(iblock,9)
		P(3,1)=STATEOLD(iblock,10)
		P(3,2)=STATEOLD(iblock,11)
		P(3,3)=STATEOLD(iblock,12)
!	Calculate stress and strain increment components in the crystal coordinate system
		!------------------------------------ ------------------------------------
		dStrainGlobal(1,1)=STRAININC(iblock,1)
		dStrainGlobal(2,2)=STRAININC(iblock,2)
		dStrainGlobal(3,3)=STRAININC(iblock,3)
		dStrainGlobal(1,2)=STRAININC(iblock,4)
		dStrainGlobal(2,3)=STRAININC(iblock,5)
		dStrainGlobal(3,1)=STRAININC(iblock,6)
		dStrainGlobal(2,1)=dStrainGlobal(1,2)
		dStrainGlobal(3,2)=dStrainGlobal(2,3)
		dStrainGlobal(1,3)=dStrainGlobal(3,1)

		StressGlobal(1,1)=STRESSOLD(iblock,1)
		StressGlobal(2,2)=STRESSOLD(iblock,2)
		StressGlobal(3,3)=STRESSOLD(iblock,3)
		StressGlobal(1,2)=STRESSOLD(iblock,4)
		StressGlobal(2,3)=STRESSOLD(iblock,5)
		StressGlobal(3,1)=STRESSOLD(iblock,6)
		StressGlobal(2,1)=StressGlobal(1,2)
		StressGlobal(3,2)=StressGlobal(2,3)
		StressGlobal(1,3)=StressGlobal(3,1)
		!------------------------------------ ------------------------------------
		!   dStrain in mat coord to crystal coord (dstrain_cr=P*dstrain_mat*P')
		!   Stress  in mat coord to crystal coord (Stress_cr=P*Stress_mat*P')
		DO j=1,3
			DO i=1,3
				dStrainLocal(i,j)=0.0
				StressLocal(i,j)=0.0
				DO k=1,3
					DO l=1,3
						dStrainLocal(i,j)=dStrainLocal(i,j)+P(i,k)*dStrainGlobal(k,l)*P(j,l)
						StressLocal(i,j)=StressLocal(i,j)+P(i,k)*StressGlobal(k,l)*P(j,l)
					END DO
				END DO
			END DO
		END DO
		!------------------------------------ ------------------------------------
		IF (TOTALTIME.EQ.0.0) THEN
			dstranp33(1,1)=0.0
			dstranp33(2,1)=0.0
			dstranp33(3,1)=0.0
			dstranp33(1,2)=0.0
			dstranp33(2,2)=0.0
			dstranp33(3,2)=0.0
			dstranp33(1,3)=0.0
			dstranp33(2,3)=0.0
			dstranp33(3,3)=0.0
			domega(1)=0.0
			domega(2)=0.0
			domega(3)=0.0
		ELSE
			dstranp33(1,1)=STATEOLD(iBlock,15)*DT
			dstranp33(2,2)=STATEOLD(iBlock,16)*DT
			dstranp33(3,3)=STATEOLD(iBlock,17)*DT
			!
			dstranp33(1,2)=STATEOLD(iBlock,18)*DT
			dstranp33(2,1)=dstranp33(1,2)
			!
			dstranp33(2,3)=STATEOLD(iBlock,19)*DT
			dstranp33(3,2)=dstranp33(2,3)
			!
			dstranp33(1,3)=STATEOLD(iBlock,20)*DT
			dstranp33(3,1)=dstranp33(1,3)
			
			domega(1)=STATEOLD(iBlock,27)*DT
			domega(2)=STATEOLD(iBlock,28)*DT
			domega(3)=STATEOLD(iBlock,29)*DT
		END IF
		!------------------------------------ ------------------------------------
		Do j=1,3
			Do i=1,3
				DO k=1,3
					DO l=1,3
						StressLocal(i,j)=StressLocal(i,j)+CE(i,j,k,l)*(dStrainLocal(k,l)-dstranp33(k,l))
					END DO
				END DO
			End Do
		End Do
		!
		!------------------------------------ ------------------------------------
		!
		!------------------------------------ --------------------------------
		!      Update Rotation matrix
		!------------------------------------ -------------------------------- 
		CALL Update_Crystal_Orientation(domega,P,Pnew)
		!
		!------------------------------------ --------------------------------
		!       END: Rotation of grain
		!------------------------------------ -------------------------------- 
		!   Stress in crystal coord to mat coord (Stress_mat=P'*Stress_cr*P)       
		DO j=1,3
			DO i=1,3
				StressGlobal(i,j)=0.0
				DO k=1,3
					DO l=1,3
						StressGlobal(i,j)=StressGlobal(i,j)+Pnew(k,i)*StressLocal(k,l)*Pnew(l,j) 
					END DO
				END DO
			END DO
		END DO
		!
		!   Solution dependent variables update
		!
		!
		STATENEW(iblock,1)=STATEOLD(iblock,1)
		STATENEW(iblock,2)=STATEOLD(iblock,2)
		STATENEW(iblock,3)=STATEOLD(iblock,3)
        !
		STATENEW(iblock,4)=Pnew(1,1)
		STATENEW(iblock,5)=Pnew(1,2)
		STATENEW(iblock,6)=Pnew(1,3)
		STATENEW(iblock,7)=Pnew(2,1)
		STATENEW(iblock,8)=Pnew(2,2)
		STATENEW(iblock,9)=Pnew(2,3)
		STATENEW(iblock,10)=Pnew(3,1)
		STATENEW(iblock,11)=Pnew(3,2)
		STATENEW(iblock,12)=Pnew(3,3)
		! STATENEW(iblock,4)=P(1,1)
		! STATENEW(iblock,5)=P(1,2)
		! STATENEW(iblock,6)=P(1,3)
		! STATENEW(iblock,7)=P(2,1)
		! STATENEW(iblock,8)=P(2,2)
		! STATENEW(iblock,9)=P(2,3)
		! STATENEW(iblock,10)=P(3,1)
		! STATENEW(iblock,11)=P(3,2)
		! STATENEW(iblock,12)=P(3,3)
		!
		!
		STRESSNEW(iblock,1) = StressGlobal(1,1)
		STRESSNEW(iblock,2) = StressGlobal(2,2)
		STRESSNEW(iblock,3) = StressGlobal(3,3)
		STRESSNEW(iblock,4) = StressGlobal(1,2)
		STRESSNEW(iblock,5) = StressGlobal(2,3)
		STRESSNEW(iblock,6) = StressGlobal(3,1)
		!
		!
		STATENEW(iblock,13)=STATEOLD(iblock,13)
		STATENEW(iblock,14)=STATEOLD(iblock,14)
		!
		STATENEW(iblock,15)=STATEOLD(iblock,15)
		STATENEW(iblock,16)=STATEOLD(iblock,16)
		STATENEW(iBlock,17)=STATEOLD(iblock,17)
		STATENEW(iBlock,18)=STATEOLD(iblock,18)
		STATENEW(iBlock,19)=STATEOLD(iblock,19)
		STATENEW(iBlock,20)=STATEOLD(iblock,20)
        !
		! Update State variable for stress components in the crystal coordinate system
		STATENEW(iblock,21)=StressLocal(1,1)
		STATENEW(iblock,22)=StressLocal(2,2)
		STATENEW(iBlock,23)=StressLocal(3,3)
		STATENEW(iBlock,24)=StressLocal(1,2)
		STATENEW(iBlock,25)=StressLocal(2,3)
		STATENEW(iBlock,26)=StressLocal(3,1)
		!
		STATENEW(iBlock,27)=STATEOLD(iBlock,27)
		STATENEW(iBlock,28)=STATEOLD(iBlock,28)
		STATENEW(iBlock,29)=STATEOLD(iBlock,29)
		STATENEW(iBlock,30)=STATEOLD(iBlock,30)
		STATENEW(iBlock,31)=STATEOLD(iBlock,31)
		STATENEW(iBlock,32)=STATEOLD(iBlock,32)
		!
    END DO
	!
	IF((TOTALTIME>=t_nextstep).AND.(t_nextstep>0.0)) THEN
		CALL VGETOUTDIR(foldername, LENOUTDIR)
		CALL VGETRANK(KPROCESSNUM)
		foldername=foldername(1:LENOUTDIR)
		write (filename_output, "(I4.4)") KPROCESSNUM+1
		!
		filename_stress=TRIM(foldername)//"/"//TRIM(filename_stress_head)//TRIM(filename_output)//".txt"
		filename_surface=TRIM(foldername)//"/"//TRIM(filename_surfaceEl_head)//TRIM(filename_output)//".txt"
		OPEN(16,file=TRIM(filename_stress),POSITION='append')
		OPEN(17,file=TRIM(filename_surface),POSITION='append')
		DO iBlock=1,NBLOCK
			WRITE(16,"(I8.8, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)")&
			Element_id(iBlock), STATENEW(iblock,21), STATENEW(iblock,22), STATENEW(iblock,23),&
			STATENEW(iblock,25), STATENEW(iblock,26), STATENEW(iblock,24), TEMPNEW(iBlock)
!			
			IF (surface_elements(Element_id(iBlock),2) .ge. 1) THEN
				IF (STATEOLD(iblock,14).EQ.1) THEN ! Activated Element: Not allowed nucleation
					WRITE(17,"(I8.8, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)")&
					Element_id(iBlock), 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
				ELSE
					WRITE(17,"(I8.8, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)")&
					Element_id(iBlock), STATENEW(iblock,21), STATENEW(iblock,22), STATENEW(iblock,23),&
					STATENEW(iblock,25), STATENEW(iblock,26), STATENEW(iblock,24)
				END IF
			END IF
		END DO
		CLOSE(17)
		CLOSE(16)
	END IF
	!
!==================================== ====================================
!                             END DO Nblock loop
!==================================== ====================================
RETURN
END SUBROUTINE VUMATXTRARG
!   END VUMATXTRARG


!==================================== ====================================
!   Subroutine for Elastic stiffness tensor
!   for tetragonal crystal, note that c66 is for c1212
!==================================== ====================================
Subroutine E_modulus(c11,c12,c44,C3333)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(3,3,3,3):: C3333
    DOUBLE PRECISION:: c11, c33, c12, c13, c44, c66
    INTEGER:: i, j, k, l
    !
    DO i=1,3
        DO j=1,3
            DO k=1,3
                DO l=1,3
                    C3333(i,j,k,l)=0.0
                END DO
            END DO
        END DO
    END DO
	!
	c33=c11
	c13=c12
	c66=c44
    !
    C3333(1,1,1,1)=c11
    C3333(2,2,2,2)=c11
    C3333(3,3,3,3)=c33
    !
    C3333(1,1,2,2)=c12
    C3333(2,2,1,1)=c12
    C3333(2,2,3,3)=c13
    C3333(3,3,2,2)=c13
    C3333(3,3,1,1)=c13
    C3333(1,1,3,3)=c13
    !
    C3333(1,2,1,2)=c66
    C3333(1,2,2,1)=c66
    C3333(2,1,1,2)=c66
    C3333(2,1,2,1)=c66
    C3333(2,3,2,3)=c44
    C3333(2,3,3,2)=c44
    C3333(3,2,2,3)=c44
    C3333(3,2,3,2)=c44
    C3333(3,1,3,1)=c44
    C3333(3,1,1,3)=c44
    C3333(1,3,3,1)=c44
    C3333(1,3,1,3)=c44
!
RETURN
END SUBROUTINE E_modulus
!====================================================================
!
!
!
!====================================================================
!     Input : xx (eight nodal points for the activation volume)
!     Output: invxv (inverse matrix to calculate the shape function)
!====================================================================
SUBROUTINE CAL_shapefunc8(xx,invxv)
	IMPLICIT NONE
	DOUBLE PRECISION:: detinv
	DOUBLE PRECISION, DIMENSION(3,3):: xv, invxv
	DOUBLE PRECISION, DIMENSION(8,3):: xx
	INTEGER:: K1
!
!	Three vectors for eta calculation (basis for the shape function)
	DO K1=1,3
		xv(K1,1)=0.125*(xx(2,K1)+xx(3,K1)+xx(6,K1)+xx(7,K1)-xx(1,K1)-xx(4,K1)-xx(5,K1)-xx(8,K1))
		xv(K1,2)=0.125*(xx(3,K1)+xx(4,K1)+xx(7,K1)+xx(8,K1)-xx(1,K1)-xx(2,K1)-xx(5,K1)-xx(6,K1))
		xv(K1,3)=0.125*(xx(5,K1)+xx(6,K1)+xx(7,K1)+xx(8,K1)-xx(1,K1)-xx(2,K1)-xx(3,K1)-xx(4,K1))
	END DO
!
	detinv = 1/(xv(1,1)*xv(2,2)*xv(3,3) - xv(1,1)*xv(2,3)*xv(3,2)&
	-xv(1,2)*xv(2,1)*xv(3,3)+xv(1,2)*xv(2,3)*xv(3,1)&
	+xv(1,3)*xv(2,1)*xv(3,2)-xv(1,3)*xv(2,2)*xv(3,1))
!
! Calculate the inverse of the matrix xv
    invxv(1,1) = +detinv * (xv(2,2)*xv(3,3) - xv(2,3)*xv(3,2))
    invxv(2,1) = -detinv * (xv(2,1)*xv(3,3) - xv(2,3)*xv(3,1))
    invxv(3,1) = +detinv * (xv(2,1)*xv(3,2) - xv(2,2)*xv(3,1))
    invxv(1,2) = -detinv * (xv(1,2)*xv(3,3) - xv(1,3)*xv(3,2))
    invxv(2,2) = +detinv * (xv(1,1)*xv(3,3) - xv(1,3)*xv(3,1))
    invxv(3,2) = -detinv * (xv(1,1)*xv(3,2) - xv(1,2)*xv(3,1))
    invxv(1,3) = +detinv * (xv(1,2)*xv(2,3) - xv(1,3)*xv(2,2))
    invxv(2,3) = -detinv * (xv(1,1)*xv(2,3) - xv(1,3)*xv(2,1))
    invxv(3,3) = +detinv * (xv(1,1)*xv(2,2) - xv(1,2)*xv(2,1))
!	
	RETURN
END SUBROUTINE CAL_shapefunc8

SUBROUTINE CAL_shapefunc4(xx,invxv)
	IMPLICIT NONE
	DOUBLE PRECISION:: detinv
	DOUBLE PRECISION, DIMENSION(3,3):: xv, invxv
	DOUBLE PRECISION, DIMENSION(4,3):: xx
	INTEGER:: K1
!
!	Three vectors for eta calculation (basis for the shape function)
	DO K1=1,3
		xv(K1,1)=xx(2,K1)-xx(1,K1)
		xv(K1,2)=xx(3,K1)-xx(1,K1)
		xv(K1,3)=xx(4,K1)-xx(1,K1)
	END DO
!
	detinv = 1/(xv(1,1)*xv(2,2)*xv(3,3) - xv(1,1)*xv(2,3)*xv(3,2)&
	-xv(1,2)*xv(2,1)*xv(3,3)+xv(1,2)*xv(2,3)*xv(3,1)&
	+xv(1,3)*xv(2,1)*xv(3,2)-xv(1,3)*xv(2,2)*xv(3,1))
!
! Calculate the inverse of the matrix xv
    invxv(1,1) = +detinv * (xv(2,2)*xv(3,3) - xv(2,3)*xv(3,2))
    invxv(2,1) = -detinv * (xv(2,1)*xv(3,3) - xv(2,3)*xv(3,1))
    invxv(3,1) = +detinv * (xv(2,1)*xv(3,2) - xv(2,2)*xv(3,1))
    invxv(1,2) = -detinv * (xv(1,2)*xv(3,3) - xv(1,3)*xv(3,2))
    invxv(2,2) = +detinv * (xv(1,1)*xv(3,3) - xv(1,3)*xv(3,1))
    invxv(3,2) = -detinv * (xv(1,1)*xv(3,2) - xv(1,2)*xv(3,1))
    invxv(1,3) = +detinv * (xv(1,2)*xv(2,3) - xv(1,3)*xv(2,2))
    invxv(2,3) = -detinv * (xv(1,1)*xv(2,3) - xv(1,3)*xv(2,1))
    invxv(3,3) = +detinv * (xv(1,1)*xv(2,2) - xv(1,2)*xv(2,1))
!	
	RETURN
END SUBROUTINE CAL_shapefunc4
!==================================== ====================================
!
!
!==================================== ====================================
!   Update the rotation matrix P
!==================================== ====================================
SUBROUTINE Update_Crystal_Orientation(dOmega,P,Pnew)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: PI=acos(-1.0_16)
    DOUBLE PRECISION :: w12, w13, w23, determinant
	DOUBLE PRECISION, DIMENSION(3) :: dOmega
    DOUBLE PRECISION, DIMENSION(3,3) :: Qr, P, Pnew
    INTEGER :: i, j, k
	!
    w12=dOmega(1)
	w13=dOmega(2)
	w23=dOmega(3)
	determinant=1.0+w12*w12+w13*w13+w23*w23
	!
	Qr(1,1)=(1.0+w23*w23-w12*w12-w13*w13)/determinant
    Qr(1,2)=(2.0*w12-2.0*w13*w23)/determinant
    Qr(1,3)=(2.0*w13+2.0*w12*w23)/determinant
    Qr(2,1)=(-2.0*w12-2.0*w13*w23)/determinant
    Qr(2,2)=(1.0+w13*w13-w12*w12-w23*w23)/determinant
    Qr(2,3)=(2.0*w23-2.0*w12*w13)/determinant
    Qr(3,1)=(-2.0*w13+2.0*w12*w23)/determinant
    Qr(3,2)=(-2.0*w23-2.0*w12*w13)/determinant
    Qr(3,3)=(1.0+w12*w12-w13*w13-w23*w23)/determinant
    !
    Pnew(1,1)=Qr(1,1)*P(1,1)+Qr(1,2)*P(2,1)+Qr(1,3)*P(3,1)
    Pnew(1,2)=Qr(1,1)*P(1,2)+Qr(1,2)*P(2,2)+Qr(1,3)*P(3,2)
    Pnew(1,3)=Qr(1,1)*P(1,3)+Qr(1,2)*P(2,3)+Qr(1,3)*P(3,3)
    Pnew(2,1)=Qr(2,1)*P(1,1)+Qr(2,2)*P(2,1)+Qr(2,3)*P(3,1)
    Pnew(2,2)=Qr(2,1)*P(1,2)+Qr(2,2)*P(2,2)+Qr(2,3)*P(3,2)
    Pnew(2,3)=Qr(2,1)*P(1,3)+Qr(2,2)*P(2,3)+Qr(2,3)*P(3,3)
    Pnew(3,1)=Qr(3,1)*P(1,1)+Qr(3,2)*P(2,1)+Qr(3,3)*P(3,1)
    Pnew(3,2)=Qr(3,1)*P(1,2)+Qr(3,2)*P(2,2)+Qr(3,3)*P(3,2)
    Pnew(3,3)=Qr(3,1)*P(1,3)+Qr(3,2)*P(2,3)+Qr(3,3)*P(3,3)
    !
    RETURN
END SUBROUTINE Update_Crystal_Orientation

!==================================== ====================================
!
! Find cross product of two vector
SUBROUTINE CAL_cross_product(cross, a, b)
	DOUBLE PRECISION, DIMENSION(3) :: cross, a, b
	
	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
	
	RETURN
END SUBROUTINE CAL_cross_product
!==================================== ====================================
!
! Find the volume of the element (in Burger's scale)
SUBROUTINE CAL_element_volume(elem_volume, xx)
	IMPLICIT NONE
	DOUBLE PRECISION:: elem_volume
	DOUBLE PRECISION, DIMENSION(3) :: A, B, C, B_cross_C
	DOUBLE PRECISION, DIMENSION(8,3):: xx, nodal_xx
	INTEGER:: i, j	
	DO i=1,3
		nodal_xx(1,i)=xx(1,i)
		nodal_xx(2,i)=xx(2,i)
		nodal_xx(3,i)=xx(4,i)
		nodal_xx(4,i)=xx(3,i)
		nodal_xx(5,i)=xx(5,i)
		nodal_xx(6,i)=xx(6,i)
		nodal_xx(7,i)=xx(8,i)
		nodal_xx(8,i)=xx(7,i)
	END DO
	! First term of the sum
	DO i=1,3
		A(i)=nodal_xx(8,i)-nodal_xx(1,i)
		B(i)=nodal_xx(2,i)-nodal_xx(1,i)
		C(i)=nodal_xx(4,i)-nodal_xx(6,i)
	END DO
	CALL CAL_cross_product(B_cross_C, B, C)
	elem_volume=abs(A(1)*B_cross_C(1) + A(2)*B_cross_C(2) + A(3)*B_cross_C(3))
	! Second term of the sum
	DO i=1,3
		B(i)=nodal_xx(5,i)-nodal_xx(1,i)
		C(i)=nodal_xx(6,i)-nodal_xx(7,i)
	END DO
	CALL CAL_cross_product(B_cross_C, B, C)
	elem_volume=elem_volume + abs(A(1)*B_cross_C(1) + A(2)*B_cross_C(2) + A(3)*B_cross_C(3))
	! Third term of the sum
	DO i=1,3
		B(i)=nodal_xx(3,i)-nodal_xx(1,i)
		C(i)=nodal_xx(7,i)-nodal_xx(4,i)
	END DO
	CALL CAL_cross_product(B_cross_C, B, C)
	elem_volume=elem_volume + abs(A(1)*B_cross_C(1) + A(2)*B_cross_C(2) + A(3)*B_cross_C(3))	
	
	elem_volume=elem_volume/6.0
	RETURN
END SUBROUTINE CAL_element_volume
!==================================== ====================================
!
! Find projection of a surface node on the surface plane
SUBROUTINE Find_Projection(projection, xx, rn1_xx, rn2_xx)
	IMPLICIT NONE
	DOUBLE PRECISION:: min_distance, i_distance
	DOUBLE PRECISION, DIMENSION(3):: projection, rn1_xx, rn2_xx, node1, node2, node3, node4, i_intersect
	DOUBLE PRECISION, DIMENSION(8,3):: xx
	INTEGER:: i, j
	
	! Surface from element node 1, 2, 3, 4
	DO i=1,3 
		node1(i)=xx(1,i)
		node2(i)=xx(2,i)
		node3(i)=xx(3,i)
		node4(i)=xx(4,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	min_distance = i_distance
	projection(1)= i_intersect(1)
	projection(2)= i_intersect(2)
	projection(3)= i_intersect(3)
	! Surface from element node 1, 4, 8, 5
	DO i=1,3 
		node1(i)=xx(1,i)
		node2(i)=xx(4,i)
		node3(i)=xx(8,i)
		node4(i)=xx(5,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IF (min_distance.GT.i_distance) THEN
		projection(1)= i_intersect(1)
		projection(2)= i_intersect(2)
		projection(3)= i_intersect(3)
		min_distance=i_distance
	END IF
	! Surface from element node 1, 2, 6, 5
	DO i=1,3 
		node1(i)=xx(1,i)
		node2(i)=xx(2,i)
		node3(i)=xx(6,i)
		node4(i)=xx(5,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IF (min_distance.GT.i_distance) THEN
		projection(1)= i_intersect(1)
		projection(2)= i_intersect(2)
		projection(3)= i_intersect(3)
		min_distance=i_distance
	END IF	
	! Surface from element node 3, 4, 8, 7
	DO i=1,3 
		node1(i)=xx(3,i)
		node2(i)=xx(4,i)
		node3(i)=xx(8,i)
		node4(i)=xx(7,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IF (min_distance.GT.i_distance) THEN
		projection(1)= i_intersect(1)
		projection(2)= i_intersect(2)
		projection(3)= i_intersect(3)
		min_distance=i_distance
	END IF	
	! Surface from element node 2, 3, 7, 6
	DO i=1,3 
		node1(i)=xx(2,i)
		node2(i)=xx(3,i)
		node3(i)=xx(7,i)
		node4(i)=xx(6,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IF (min_distance.GT.i_distance) THEN
		projection(1)= i_intersect(1)
		projection(2)= i_intersect(2)
		projection(3)= i_intersect(3)
		min_distance=i_distance
	END IF	
	! Surface from element node 5, 6, 7, 8
	DO i=1,3 
		node1(i)=xx(5,i)
		node2(i)=xx(6,i)
		node3(i)=xx(7,i)
		node4(i)=xx(8,i)
		i_intersect(i)=0.0
	END DO
	i_distance=0.0
	CALL Find_Intersection_Surface(i_intersect, i_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IF (min_distance.GT.i_distance) THEN
		projection(1)= i_intersect(1)
		projection(2)= i_intersect(2)
		projection(3)= i_intersect(3)
		min_distance=i_distance
	END IF		
	
	RETURN
END SUBROUTINE	Find_Projection


!==================================== ====================================
!
! Find intersection of a segment to the surface
SUBROUTINE Find_Intersection_Surface(intersection, min_distance, node1, node2, node3, node4, rn1_xx, rn2_xx)
	IMPLICIT NONE
	DOUBLE PRECISION:: min_distance, i_distance
	DOUBLE PRECISION, DIMENSION(3):: intersection, rn1_xx, rn2_xx, node1, node2, node3, node4, point1, point2, point3, i_intersect
	INTEGER:: i, j
	! Plane from point 1, 2, 3
	DO i=1,3 
		point1(i)=node1(i)
		point2(i)=node2(i)
		point3(i)=node3(i)
		i_intersect(i)=0.0
	END DO	
	i_distance=0.0
	CALL Find_Intersection_Plane(i_intersect, i_distance, point1, point2, point3, rn1_xx, rn2_xx)
	min_distance = i_distance
	intersection(1)= i_intersect(1)
	intersection(2)= i_intersect(2)
	intersection(3)= i_intersect(3)	
	! Plane from point 1, 2, 4
	DO i=1,3 
		point1(i)=node1(i)
		point2(i)=node2(i)
		point3(i)=node4(i)
		i_intersect(i)=0.0
	END DO	
	i_distance=0.0
	CALL Find_Intersection_Plane(i_intersect, i_distance, point1, point2, point3, rn1_xx, rn2_xx)	
	IF (min_distance.GT.i_distance) THEN
		intersection(1)= i_intersect(1)
		intersection(2)= i_intersect(2)
		intersection(3)= i_intersect(3)	
		min_distance=i_distance
	END IF
	! Plane from point 1, 3, 4
	DO i=1,3 
		point1(i)=node1(i)
		point2(i)=node3(i)
		point3(i)=node4(i)
		i_intersect(i)=0.0
	END DO	
	i_distance=0.0
	CALL Find_Intersection_Plane(i_intersect, i_distance, point1, point2, point3, rn1_xx, rn2_xx)	
	IF (min_distance.GT.i_distance) THEN
		intersection(1)= i_intersect(1)
		intersection(2)= i_intersect(2)
		intersection(3)= i_intersect(3)	
		min_distance=i_distance
	END IF	
	! Plane from point 2, 3, 4
	DO i=1,3 
		point1(i)=node2(i)
		point2(i)=node3(i)
		point3(i)=node4(i)
		i_intersect(i)=0.0
	END DO	
	i_distance=0.0
	CALL Find_Intersection_Plane(i_intersect, i_distance, point1, point2, point3, rn1_xx, rn2_xx)	
	IF (min_distance.GT.i_distance) THEN
		intersection(1)= i_intersect(1)
		intersection(2)= i_intersect(2)
		intersection(3)= i_intersect(3)	
		min_distance=i_distance
	END IF	
	
	RETURN
END SUBROUTINE Find_Intersection_Surface

!==================================== ====================================
!
! Find intersection of a segment to a plane
SUBROUTINE Find_Intersection_Plane(intersection, distance, point1, point2, point3, rn1_xx, rn2_xx)
	IMPLICIT NONE
	DOUBLE PRECISION:: distance, normal_length, rn_dot_n, d
	DOUBLE PRECISION, DIMENSION(3):: intersection, rn1_xx, rn2_xx, point1, point2, point3
	DOUBLE PRECISION, DIMENSION(3):: vector_a, vector_b, vector_rn, rn_center, plane_normal
	INTEGER:: LENOUTDIR, KPROCESSNUM
	CHARACTER(LEN=1024) :: foldername, filename_output, filename_intersect
	INTEGER:: i, j
	
	DO i=1,3
		vector_a(i)=point2(i)-point1(i)
		vector_b(i)=point3(i)-point1(i)
		vector_rn(i)=rn2_xx(i)-rn1_xx(i)
		rn_center(i)=0.5*(rn1_xx(i)+rn2_xx(i))
	END DO
	CALL CAL_cross_product(plane_normal, vector_a, vector_b)
	!
	rn_dot_n=vector_rn(1)*plane_normal(1)+vector_rn(2)*plane_normal(2)+vector_rn(3)*plane_normal(3)
	IF (abs(rn_dot_n).GT.0.0) THEN! There is one intersection
		d=(point1(1)-rn1_xx(1))*plane_normal(1)+(point1(2)-rn1_xx(2))*plane_normal(2)+(point1(3)-rn1_xx(3))*plane_normal(3)
		d=d/rn_dot_n
		intersection(1)=rn1_xx(1)+d*vector_rn(1)
		intersection(2)=rn1_xx(2)+d*vector_rn(2)
		intersection(3)=rn1_xx(3)+d*vector_rn(3)
		distance=sqrt((intersection(1)-rn_center(1))*(intersection(1)-rn_center(1))+(intersection(2)-rn_center(2))*(intersection(2)-rn_center(2))+(intersection(3)-rn_center(3))*(intersection(3)-rn_center(3)))
	ELSE ! Segment parallel or coincide with plane
		distance=1.0e9
		intersection(1)=0.0
		intersection(2)=0.0
		intersection(3)=0.0
	END IF
	CALL VGETOUTDIR(foldername, LENOUTDIR)
	CALL VGETRANK(KPROCESSNUM)
	foldername=foldername(1:LENOUTDIR)
	write (filename_output, "(I4.4)") KPROCESSNUM+1
	filename_intersect=TRIM(foldername)//"/plane_intersection"//TRIM(filename_output)//".txt"
	!OPEN(17,file=TRIM(filename_intersect),POSITION='append')
	!WRITE(17,"(E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)")  rn1_xx(3), rn2_xx(3), distance, intersection(1), intersection(2), intersection(3)
	!CLOSE(17)
	
	RETURN
END SUBROUTINE Find_Intersection_Plane
!==================================== ====================================
!
!====================================================================
!     vucharlength: used for assigning plastinc strains for each int point
!====================================================================
SUBROUTINE vucharlength(&
!	 Read only variables-
		nblock, nfieldv, nprops, ncomp, ndim, nnode, nstatev,&
		kSecPt, kLayer, kIntPt, jElType, jElem,&
		totalTime, stepTime, dt,&
		cmname, coordMp, coordNode, direct, T, props,&
		field, stateOld,&
!	Write only variables-
		charLength)
!
	IMPLICIT DOUBLE PRECISION (a-h,o-z)
!	include 'vaba_param.inc'
!
	DIMENSION jElType(3), jElem(nblock), coordMp(nblock,ndim),&
		coordNode(nblock,nnode,ndim),&
		direct(nblock,3,3), T(nblock,3,3), props(nprops),&
		stateOld(nblock,nstatev), charLength(nblock,ncomp),&
		field(nblock, nfieldv)
!
	CHARACTER*80 cmname
!
!
	DOUBLE PRECISION, PARAMETER :: tolerance=1.001, burgers=1.0, shape_fn_tol=1.1, shape_fn_tol_rn=1.15, b_mag=2.864e-10
	DOUBLE PRECISION, PARAMETER :: PI=acos(-1.0_16)
	! Here, we assume that the length unit for FEM is the same as the size of burgers vector.
	DOUBLE PRECISION:: normalizer
	DOUBLE PRECISION:: temp_rn_length, temp_rn_r, temp_rn_z, common_omega, temp_omega, determinant, rn_old_dot_new, seg_length, seg_inside_length, new_seg_length, elem_volume
	DOUBLE PRECISION, DIMENSION(3):: Crystal_COORDS, domega, Ralra, eta!xx_min, xx_max, temp_swept_center, Ralra, eta
	DOUBLE PRECISION, DIMENSION(3):: temp_rn1, temp_rn2, temp_rn_v, temp_rn_h, temp_rn_m, temp_rn_center, rn_min, rn_max, segment_center, temp_link_center
	DOUBLE PRECISION, DIMENSION(3):: rn_surface_old, projection, rn1_xx, rn2_xx
	DOUBLE PRECISION, DIMENSION(6):: dstranp
	DOUBLE PRECISION, DIMENSION(3,3):: invxv
	DOUBLE PRECISION, DIMENSION(8,3):: xx
	DOUBLE PRECISION, DIMENSION(4,3):: xx4by3
	!
	INTEGER:: LENOUTDIR
	CHARACTER(*), PARAMETER :: file_element_head = "Element_center"
	CHARACTER(*), PARAMETER :: file_element_output_head = "Output_element_center"
	CHARACTER(*), PARAMETER :: file_link_head = "Output_link"
	CHARACTER(*), PARAMETER :: file_output_head = "Output_nodes"
!	CHARACTER(*), PARAMETER :: file_CoordNode_head = "OriginalCoordNodes"
	CHARACTER(*), PARAMETER :: file_rn_head = "rn_check"
	CHARACTER(*), PARAMETER :: file_FEM_Dt_head = "FEM_Dt"
	CHARACTER(LEN=1024) :: foldername, filename_output, filename_link, filename_element, filename_output_element, filename_node, filename_CoordNode, filename_rn_check, filename_FEM_Dt
	INTEGER :: KPROCESSNUM
!	DOUBLE PRECISION, DIMENSION(3,3):: invxv
!	DOUBLE PRECISION, DIMENSION(8,3):: xx
!
	INTEGER:: iBlock, num_Act_vol, K1, K2, K3, i_segment, i_rn, wt_Act_Vol, num_elements, rn_surface_check, segment_elem, nbr_elem, nbr_jElem, j, k, m, n, rn1_id, rn2_id
!
	INTEGER, DIMENSION(200000) :: Element_id
	INTEGER, PARAMETER :: n_segment_max=80000
	INTEGER, PARAMETER :: n_rn_max=80000
	INTEGER:: n_segment, n_segment_final, k_iteration, link_element, num_rn
	DOUBLE PRECISION :: EffSlip
	DOUBLE PRECISION:: t_step, t_nextstep
	! t_step: total time increment in the current step
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: rn1_old, rn2_old, rn1_new, rn2_new, rn1_old_final, rn2_old_final, rn1_new_final, rn2_new_final
	INTEGER, DIMENSION (n_segment_max, 2):: link_elem_id
	DOUBLE PRECISION, DIMENSION(n_rn_max,3):: rn_nodes
	INTEGER, PARAMETER :: n_elements_max=600000
	INTEGER, DIMENSION (n_elements_max,2):: surface_elements
	DOUBLE PRECISION, DIMENSION(n_elements_max,8,3):: original_coordNodes
	DOUBLE PRECISION, DIMENSION(nblock,8,ndim):: Origin_CoordNode
!	INTEGER, DIMENSION(n_segment_max):: link_flag
	INTEGER, DIMENSION (n_segment_max,2):: node_id_final
	INTEGER, DIMENSION (n_rn_max,3):: rn_id
	INTEGER, DIMENSION (n_rn_max,2):: rn_element
	DOUBLE PRECISION, DIMENSION(nblock,3):: element_center, element_min, element_max 
	DOUBLE PRECISION, DIMENSION(nblock):: element_volume, element_density
	INTEGER, DIMENSION (nblock):: element_within
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: swept_normal, swept_center
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: dw_p
	DOUBLE PRECISION, DIMENSION(n_segment_max,6):: de_p
	COMMON /sharedf/ de_p, dw_p, t_step, t_nextstep, EffSlip, rn1_old, rn2_old, rn1_new, rn2_new, rn1_old_final, rn2_old_final, rn1_new_final, rn2_new_final, link_elem_id!
	COMMON /sharedi/ n_segment, n_segment_final, k_iteration
	COMMON /sharenode/ node_id_final, rn_id, rn_nodes, num_rn, rn_element
	COMMON /sharedsurface/ surface_elements, num_elements
	COMMON /elem_id/ Element_id
	COMMON /sharedoriginalnodes/ original_coordNodes	
!
	IF (TOTALTIME.EQ.0.0) THEN
		DO iBlock = 1, nblock
			! initial positions at the crystal coordinate system.
			! (Initially, the crystal coord system is aligned to mat coord system.)
			STATEOLD(iblock,1)=COORDMP(iblock,1)
			STATEOLD(iblock,2)=COORDMP(iblock,2)
			STATEOLD(iblock,3)=COORDMP(iblock,3)
			!
			! Initial rotation of the crystal coord system (aligned to mat coord sysrtem)
			!
			STATEOLD(iblock,4)=1.0 !P(1,1)
			STATEOLD(iblock,5)=0.0 !P(1,2)
			STATEOLD(iblock,6)=0.0 !P(1,3)
			STATEOLD(iblock,7)=0.0 !P(2,1)
			STATEOLD(iblock,8)=1.0 !P(2,2)
			STATEOLD(iblock,9)=0.0 !P(2,3)
			STATEOLD(iblock,10)=0.0 !P(3,1)
			STATEOLD(iblock,11)=0.0 !P(3,2)
			STATEOLD(iblock,12)=1.0 !P(3,3)
			STATEOLD(iBlock,13)=0
			STATEOLD(iBlock,31)=0.0
		END DO
	END IF
	IF (stepTime.EQ.0.0) THEN
		DO iBlock = 1, nblock
			STATEOLD(iBlock,14)=0
		END DO
	END IF
	!
	common_omega=t_step*PI*EffSlip*EffSlip/burgers
	!
	DO iBlock = 1, nblock! initialization of plastic strain increment and rotation
		STATEOLD(iBlock,15)=0! dstranp(1)=dstranp(11)
		STATEOLD(iBlock,16)=0! dstranp(2)=dstranp(22)
		STATEOLD(iBlock,17)=0! dstranp(3)=dstranp(33)
		STATEOLD(iBlock,18)=0! dstranp(4)=dstranp(12)
		STATEOLD(iBlock,19)=0! dstranp(5)=dstranp(23)
		STATEOLD(iBlock,20)=0! dstranp(6)=dstranp(31)
		!
		STATEOLD(iBlock,27)=0! domega(1)=domega(12) half
		STATEOLD(iBlock,28)=0! domega(2)=domega(13) half
		STATEOLD(iBlock,29)=0! domega(3)=domega(23) half
		element_within(iBlock)=0
		Element_id(iBlock)=jElem(iBlock)
	END DO
	!
	DO i_segment=1, n_segment
		! Find a temporary position of rn1 and rn2 for each time increment
		temp_rn1(1)=((t_nextstep-totalTime)*rn1_old(i_segment,1)+(totalTime+t_step-t_nextstep)*rn1_new(i_segment,1))/t_step
		temp_rn1(2)=((t_nextstep-totalTime)*rn1_old(i_segment,2)+(totalTime+t_step-t_nextstep)*rn1_new(i_segment,2))/t_step
		temp_rn1(3)=((t_nextstep-totalTime)*rn1_old(i_segment,3)+(totalTime+t_step-t_nextstep)*rn1_new(i_segment,3))/t_step
		temp_rn2(1)=((t_nextstep-totalTime)*rn2_old(i_segment,1)+(totalTime+t_step-t_nextstep)*rn2_new(i_segment,1))/t_step
		temp_rn2(2)=((t_nextstep-totalTime)*rn2_old(i_segment,2)+(totalTime+t_step-t_nextstep)*rn2_new(i_segment,2))/t_step
		temp_rn2(3)=((t_nextstep-totalTime)*rn2_old(i_segment,3)+(totalTime+t_step-t_nextstep)*rn2_new(i_segment,3))/t_step
		! Calculate center of the temporary rn1 and rn2
		temp_rn_center(1)=0.5*(temp_rn1(1)+temp_rn2(1))
		temp_rn_center(2)=0.5*(temp_rn1(2)+temp_rn2(2))
		temp_rn_center(3)=0.5*(temp_rn1(3)+temp_rn2(3))
		! Calculate segment vector of the temporary rn1 and rn2: v=(rn2-rn1)
		temp_rn_v(1)=temp_rn2(1)-temp_rn1(1)
		temp_rn_v(2)=temp_rn2(2)-temp_rn1(2)
		temp_rn_v(3)=temp_rn2(3)-temp_rn1(3)
		temp_rn_length=sqrt(temp_rn_v(1)*temp_rn_v(1)+temp_rn_v(2)*temp_rn_v(2)+temp_rn_v(3)*temp_rn_v(3))
		IF(temp_rn_length.GT.0) THEN!normalization
			temp_rn_v(1)=temp_rn_v(1)/temp_rn_length
			temp_rn_v(2)=temp_rn_v(2)/temp_rn_length
			temp_rn_v(3)=temp_rn_v(3)/temp_rn_length
		END IF
		!
		! Find max and min of rn1 and rn2 
		rn_max(1)=max(temp_rn1(1),temp_rn2(1))+EffSlip
		rn_max(2)=max(temp_rn1(2),temp_rn2(2))+EffSlip
		rn_max(3)=max(temp_rn1(3),temp_rn2(3))+EffSlip
		rn_min(1)=min(temp_rn1(1),temp_rn2(1))-EffSlip
		rn_min(2)=min(temp_rn1(2),temp_rn2(2))-EffSlip
		rn_min(3)=min(temp_rn1(3),temp_rn2(3))-EffSlip
		!
		! Calculate effective moving distance of the dislocation segment during t_step
		temp_omega=(rn1_new(i_segment,1)+rn2_new(i_segment,1)-rn1_old(i_segment,1)-rn2_old(i_segment,1))**2.0
		temp_omega=temp_omega+(rn1_new(i_segment,2)+rn2_new(i_segment,2)-rn1_old(i_segment,2)-rn2_old(i_segment,2))**2.0
		temp_omega=temp_omega+(rn1_new(i_segment,3)+rn2_new(i_segment,3)-rn1_old(i_segment,3)-rn2_old(i_segment,3))**2.0
		temp_omega=0.5*sqrt(temp_omega)
		IF(temp_omega/=0) THEN
			temp_omega=common_omega/temp_omega
		END IF
		
		!
		!! Element loop start
		DO iBlock = 1, nblock
			wt_Act_Vol=0
			Crystal_COORDS(1)=STATEOLD(iBlock,1)
			Crystal_COORDS(2)=STATEOLD(iBlock,2)
			Crystal_COORDS(3)=STATEOLD(iBlock,3)
			! Determine inside or outside of the current act vol by calculating eta.
			IF ((Crystal_COORDS(1).LE.rn_max(1)).AND.(Crystal_COORDS(1).GE.rn_min(1))) THEN
				IF ((Crystal_COORDS(2).LE.rn_max(2)).AND.(Crystal_COORDS(2).GE.rn_min(2))) THEN
					IF ((Crystal_COORDS(3).LE.rn_max(3)).AND.(Crystal_COORDS(3).GE.rn_min(3))) THEN
						!
						temp_rn_m(1)=Crystal_COORDS(1)-temp_rn_center(1)
						temp_rn_m(2)=Crystal_COORDS(2)-temp_rn_center(2)
						temp_rn_m(3)=Crystal_COORDS(3)-temp_rn_center(3)
						!
						temp_rn_z=temp_rn_m(1)*temp_rn_v(1)+temp_rn_m(2)*temp_rn_v(2)+temp_rn_m(3)*temp_rn_v(3)
						!
						IF(abs(temp_rn_z).GT.(0.5*temp_rn_length)) THEN
							IF(temp_rn_z.GT.0.0) THEN!temp_rn2 side
								!temp_rn_r=sqrt((Crystal_COORDS(1)-temp_rn2(1))**2+(Crystal_COORDS(2)-temp_rn2(2))**2+(Crystal_COORDS(3)-temp_rn2(3))**2)
							ELSE!temp_rn1 side
								!temp_rn_r=sqrt((Crystal_COORDS(1)-temp_rn1(1))**2+(Crystal_COORDS(2)-temp_rn1(2))**2+(Crystal_COORDS(3)-temp_rn1(3))**2)
							END IF
							temp_rn_r=2*EffSlip
						ELSE
							temp_rn_r=sqrt((temp_rn_m(1)-temp_rn_z*temp_rn_v(1))**2+(temp_rn_m(2)-temp_rn_z*temp_rn_v(2))**2+(temp_rn_m(3)-temp_rn_z*temp_rn_v(3))**2)
						END IF
						!
						IF(temp_rn_r.LE.EffSlip) THEN
							wt_Act_Vol=1
						END IF
						!
					END IF
				END IF
			END IF
			! If the current element is inside the act vol, accumulate the plastic strain and spin.
			IF (wt_Act_Vol.EQ.1) THEN
				IF(STATEOLD(iblock,14).EQ.0) THEN
					STATEOLD(iblock,13)=STATEOLD(iblock,13)+1
				END IF
				STATEOLD(iblock,14)=1
				element_within(iBlock)=1
				!
				IF(temp_omega/=0) THEN
					STATEOLD(iBlock,15)=STATEOLD(iBlock,15)+de_p(i_segment,1)/temp_omega! dstranp(1)=dstranp(11)
					STATEOLD(iBlock,16)=STATEOLD(iBlock,16)+de_p(i_segment,2)/temp_omega! dstranp(2)=dstranp(22)
					STATEOLD(iBlock,17)=STATEOLD(iBlock,17)+de_p(i_segment,3)/temp_omega! dstranp(3)=dstranp(33)
					STATEOLD(iBlock,18)=STATEOLD(iBlock,18)+de_p(i_segment,4)/temp_omega! dstranp(4)=dstranp(12)
					STATEOLD(iBlock,19)=STATEOLD(iBlock,19)+de_p(i_segment,5)/temp_omega! dstranp(5)=dstranp(23)
					STATEOLD(iBlock,20)=STATEOLD(iBlock,20)+de_p(i_segment,6)/temp_omega! dstranp(6)=dstranp(31)
					!
					STATEOLD(iBlock,27)=STATEOLD(iBlock,27)+dw_p(i_segment,1)/temp_omega! domega(1)=domega(12) half
					STATEOLD(iBlock,28)=STATEOLD(iBlock,28)+dw_p(i_segment,2)/temp_omega! domega(2)=domega(13) half
					STATEOLD(iBlock,29)=STATEOLD(iBlock,29)+dw_p(i_segment,3)/temp_omega! domega(3)=domega(23) half
				END IF
				!
			END IF
		END DO
		!!
	END DO
	DO iBlock = 1, nblock
		IF (element_within(iBlock) .EQ. 0) THEN
			STATEOLD(iblock,14)=0
		END IF
!       SDV32 stores the accumulated plastic strain e_zz        
		STATEOLD(iBlock,32)=STATEOLD(iBlock,32)+DT*STATEOLD(iBlock,17)
		IF (STATEOLD(iblock,14).EQ.1) THEN
			! Accumulation of equivalent plastic strain (based on von Mises yield function)
			STATEOLD(iBlock,30)=STATEOLD(iBlock,30)+DT*sqrt(2.0*(STATEOLD(iBlock,15)*STATEOLD(iBlock,15)&
			+STATEOLD(iBlock,16)*STATEOLD(iBlock,16)+STATEOLD(iBlock,17)*STATEOLD(iBlock,17)&
			+2.0*STATEOLD(iBlock,18)*STATEOLD(iBlock,18)+2.0*STATEOLD(iBlock,19)*STATEOLD(iBlock,19)&
			+2.0*STATEOLD(iBlock,20)*STATEOLD(iBlock,20))/3.0)
		END IF
	END DO
	DO iBlock = 1, nblock
        charLength(iBlock,ncomp)=1.0
	END DO		
!	Print out the data for next ParaDiS run
	IF((TOTALTIME>=t_nextstep).AND.(t_nextstep>0.0)) THEN
		CALL VGETOUTDIR(foldername, LENOUTDIR)
		CALL VGETRANK(KPROCESSNUM)
		foldername=foldername(1:LENOUTDIR)
		write (filename_output, "(I4.4)") KPROCESSNUM+1 
! 	Print out the stable time increment	
		IF (KPROCESSNUM.eq.0) THEN
			filename_FEM_Dt=TRIM(foldername)//"/"//TRIM(file_FEM_Dt_head)//TRIM(filename_output)//".txt"
			OPEN(16,file=TRIM(filename_FEM_Dt),POSITION='append')
				WRITE(16,"(E16.8E2)") DT
			CLOSE (16)
		END IF 			
!		Import original coordinates
		DO iBlock=1,NBLOCK
			DO j=1,8
				DO k=1,3
					Origin_CoordNode(iBlock,j,k)=original_coordNodes(jElem(iBlock),j,k)
				END DO
			END DO		
		END DO
! 		Write nodal coordinates (no longer needed)
		filename_node=TRIM(foldername)//"/"//TRIM(file_output_head)//TRIM(filename_output)//".txt"
		OPEN(16,file=TRIM(filename_node),POSITION='append')
		IF(nnode.EQ.4) THEN
			DO iBlock=1,NBLOCK
				WRITE(16,"(I8.8, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2,&
				E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)") jElem(iBlock),&
				Origin_CoordNode(iBlock,1,1), Origin_CoordNode(iBlock,1,2), Origin_CoordNode(iBlock,1,3),&
				Origin_CoordNode(iBlock,2,1), Origin_CoordNode(iBlock,2,2), Origin_CoordNode(iBlock,2,3),&
				Origin_CoordNode(iBlock,3,1), Origin_CoordNode(iBlock,3,2), Origin_CoordNode(iBlock,3,3),&
				Origin_CoordNode(iBlock,4,1), Origin_CoordNode(iBlock,4,2), Origin_CoordNode(iBlock,4,3)
			END DO
		ELSEIF(nnode.EQ.8) THEN
			DO iBlock=1,NBLOCK
				WRITE(16,"(I8.8, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2,&
				E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2,&
				E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2,&
				E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2, E16.8E2)") jElem(iBlock),&
				Origin_CoordNode(iBlock,1,1), Origin_CoordNode(iBlock,1,2), Origin_CoordNode(iBlock,1,3),& 
				Origin_CoordNode(iBlock,2,1), Origin_CoordNode(iBlock,2,2), Origin_CoordNode(iBlock,2,3),&
				Origin_CoordNode(iBlock,3,1), Origin_CoordNode(iBlock,3,2), Origin_CoordNode(iBlock,3,3),&
				Origin_CoordNode(iBlock,4,1), Origin_CoordNode(iBlock,4,2), Origin_CoordNode(iBlock,4,3),&
				Origin_CoordNode(iBlock,5,1), Origin_CoordNode(iBlock,5,2), Origin_CoordNode(iBlock,5,3),&
				Origin_CoordNode(iBlock,6,1), Origin_CoordNode(iBlock,6,2), Origin_CoordNode(iBlock,6,3),&
				Origin_CoordNode(iBlock,7,1), Origin_CoordNode(iBlock,7,2), Origin_CoordNode(iBlock,7,3),&
				Origin_CoordNode(iBlock,8,1), Origin_CoordNode(iBlock,8,2), Origin_CoordNode(iBlock,8,3)
			END DO
		END IF
		CLOSE(16)		
! Check the links and print out the element in which lies the segment center 
		IF(nnode.EQ.4) THEN
			DO iBlock=1,nblock
				element_center(iBlock,1)=0.25D0*(Origin_CoordNode(iBlock,1,1)+Origin_CoordNode(iBlock,2,1)+Origin_CoordNode(iBlock,3,1)+Origin_CoordNode(iBlock,4,1))
				element_center(iBlock,2)=0.25D0*(Origin_CoordNode(iBlock,1,2)+Origin_CoordNode(iBlock,2,2)+Origin_CoordNode(iBlock,3,2)+Origin_CoordNode(iBlock,4,2))
				element_center(iBlock,3)=0.25D0*(Origin_CoordNode(iBlock,1,3)+Origin_CoordNode(iBlock,2,3)+Origin_CoordNode(iBlock,3,3)+Origin_CoordNode(iBlock,4,3))
				element_max(iBlock,1)=max(Origin_CoordNode(iBlock,1,1),Origin_CoordNode(iBlock,2,1),Origin_CoordNode(iBlock,3,1),Origin_CoordNode(iBlock,4,1))
				element_max(iBlock,2)=max(Origin_CoordNode(iBlock,1,2),Origin_CoordNode(iBlock,2,2),Origin_CoordNode(iBlock,3,2),Origin_CoordNode(iBlock,4,2))
				element_max(iBlock,3)=max(Origin_CoordNode(iBlock,1,3),Origin_CoordNode(iBlock,2,3),Origin_CoordNode(iBlock,3,3),Origin_CoordNode(iBlock,4,3))
				element_min(iBlock,1)=min(Origin_CoordNode(iBlock,1,1),Origin_CoordNode(iBlock,2,1),Origin_CoordNode(iBlock,3,1),Origin_CoordNode(iBlock,4,1))
				element_min(iBlock,2)=min(Origin_CoordNode(iBlock,1,2),Origin_CoordNode(iBlock,2,2),Origin_CoordNode(iBlock,3,2),Origin_CoordNode(iBlock,4,2))
				element_min(iBlock,3)=min(Origin_CoordNode(iBlock,1,3),Origin_CoordNode(iBlock,2,3),Origin_CoordNode(iBlock,3,3),Origin_CoordNode(iBlock,4,3))
			END DO
		ELSEIF(nnode.EQ.8) THEN
			DO iBlock=1,nblock
				element_center(iBlock,1)=0.125D0*(Origin_CoordNode(iBlock,1,1)+Origin_CoordNode(iBlock,2,1)+Origin_CoordNode(iBlock,3,1)+Origin_CoordNode(iBlock,4,1)+Origin_CoordNode(iBlock,5,1)+Origin_CoordNode(iBlock,6,1)+Origin_CoordNode(iBlock,7,1)+Origin_CoordNode(iBlock,8,1))
				element_center(iBlock,2)=0.125D0*(Origin_CoordNode(iBlock,1,2)+Origin_CoordNode(iBlock,2,2)+Origin_CoordNode(iBlock,3,2)+Origin_CoordNode(iBlock,4,2)+Origin_CoordNode(iBlock,5,2)+Origin_CoordNode(iBlock,6,2)+Origin_CoordNode(iBlock,7,2)+Origin_CoordNode(iBlock,8,2))
				element_center(iBlock,3)=0.125D0*(Origin_CoordNode(iBlock,1,3)+Origin_CoordNode(iBlock,2,3)+Origin_CoordNode(iBlock,3,3)+Origin_CoordNode(iBlock,4,3)+Origin_CoordNode(iBlock,5,3)+Origin_CoordNode(iBlock,6,3)+Origin_CoordNode(iBlock,7,3)+Origin_CoordNode(iBlock,8,3))
				element_max(iBlock,1)=max(Origin_CoordNode(iBlock,1,1),Origin_CoordNode(iBlock,2,1),Origin_CoordNode(iBlock,3,1),Origin_CoordNode(iBlock,4,1),Origin_CoordNode(iBlock,5,1),Origin_CoordNode(iBlock,6,1),Origin_CoordNode(iBlock,7,1),Origin_CoordNode(iBlock,8,1))
				element_max(iBlock,2)=max(Origin_CoordNode(iBlock,1,2),Origin_CoordNode(iBlock,2,2),Origin_CoordNode(iBlock,3,2),Origin_CoordNode(iBlock,4,2),Origin_CoordNode(iBlock,5,2),Origin_CoordNode(iBlock,6,2),Origin_CoordNode(iBlock,7,2),Origin_CoordNode(iBlock,8,2))
				element_max(iBlock,3)=max(Origin_CoordNode(iBlock,1,3),Origin_CoordNode(iBlock,2,3),Origin_CoordNode(iBlock,3,3),Origin_CoordNode(iBlock,4,3),Origin_CoordNode(iBlock,5,3),Origin_CoordNode(iBlock,6,3),Origin_CoordNode(iBlock,7,3),Origin_CoordNode(iBlock,8,3))
				element_min(iBlock,1)=min(Origin_CoordNode(iBlock,1,1),Origin_CoordNode(iBlock,2,1),Origin_CoordNode(iBlock,3,1),Origin_CoordNode(iBlock,4,1),Origin_CoordNode(iBlock,5,1),Origin_CoordNode(iBlock,6,1),Origin_CoordNode(iBlock,7,1),Origin_CoordNode(iBlock,8,1))
				element_min(iBlock,2)=min(Origin_CoordNode(iBlock,1,2),Origin_CoordNode(iBlock,2,2),Origin_CoordNode(iBlock,3,2),Origin_CoordNode(iBlock,4,2),Origin_CoordNode(iBlock,5,2),Origin_CoordNode(iBlock,6,2),Origin_CoordNode(iBlock,7,2),Origin_CoordNode(iBlock,8,2))
				element_min(iBlock,3)=min(Origin_CoordNode(iBlock,1,3),Origin_CoordNode(iBlock,2,3),Origin_CoordNode(iBlock,3,3),Origin_CoordNode(iBlock,4,3),Origin_CoordNode(iBlock,5,3),Origin_CoordNode(iBlock,6,3),Origin_CoordNode(iBlock,7,3),Origin_CoordNode(iBlock,8,3))
			END DO
		END IF
		filename_link=TRIM(foldername)//"/"//TRIM(file_link_head)//TRIM(filename_output)//".txt"
		OPEN(16,file=TRIM(filename_link),POSITION='append')
		DO i_segment=1, n_segment_final
			temp_link_center(1)=0.5D0*(rn1_new_final(i_segment,1)+rn2_new_final(i_segment,1))
			temp_link_center(2)=0.5D0*(rn1_new_final(i_segment,2)+rn2_new_final(i_segment,2))
			temp_link_center(3)=0.5D0*(rn1_new_final(i_segment,3)+rn2_new_final(i_segment,3))
			link_element=0
			iBlock=1
			DO WHILE ((link_element.eq.0).and.(iBlock.le.nblock))
				IF ((element_max(iBlock,1).ge.temp_link_center(1)).and.(element_min(iBlock,1).le.temp_link_center(1))) THEN
					IF ((element_max(iBlock,2).ge.temp_link_center(2)).and.(element_min(iBlock,2).le.temp_link_center(2))) THEN
						IF ((element_max(iBlock,3).ge.temp_link_center(3)).and.(element_min(iBlock,3).le.temp_link_center(3))) THEN
							IF(nnode.EQ.4) THEN
								DO j=1,4
									DO k=1,3
										xx4by3(j,k)=Origin_CoordNode(iBlock,j,k)
									END DO
								END DO
								! Reset inverse vx
								DO j=1,3
									DO k=1,3
										invxv(j,k)=0.0D0
									END DO
								END DO
								!Call function to get inverse of vx
								CALL CAL_shapefunc4(xx4by3,invxv)
								DO m=1,3
									Ralra(m)=temp_link_center(m)-xx4by3(1,m)
								END DO
								DO n=1,3
									eta(n)=invxv(n,1)*Ralra(1)+invxv(n,2)*Ralra(2)+invxv(n,3)*Ralra(3)
								END DO
								! Check if center within shape function
								IF ((eta(1).ge.0.0).and.(eta(2).ge.0.0).and.(eta(3).ge.0.0)) THEN
									IF ((eta(1)+eta(2)+eta(3)).le.shape_fn_tol) THEN
										link_element=jElem(iBlock)
										link_elem_id(i_segment,1)=iBlock
										link_elem_id(i_segment,2)=jElem(iBlock)										
									END IF
								END IF								
							ELSEIF(nnode.EQ.8) THEN
								DO j=1,8
									DO k=1,3
										xx(j,k)=Origin_CoordNode(iBlock,j,k)
									END DO
								END DO
								! Reset inverse vx
								DO j=1,3
									DO k=1,3
										invxv(j,k)=0.0D0
									END DO
								END DO
								!Call function to get inverse of vx
								CALL CAL_shapefunc8(xx,invxv)
								DO m=1,3
									Ralra(m)=temp_link_center(m)-element_center(iBlock,m)
								END DO
								DO n=1,3
									eta(n)=invxv(n,1)*Ralra(1)+invxv(n,2)*Ralra(2)+invxv(n,3)*Ralra(3)
								END DO
								! Check if center within shape function
								IF ((abs(eta(1)).le.shape_fn_tol).and.(abs(eta(2)).le.shape_fn_tol).and.(abs(eta(3)).le.shape_fn_tol)) THEN
									link_element=jElem(iBlock)
									link_elem_id(i_segment,1)=iBlock
									link_elem_id(i_segment,2)=jElem(iBlock)
								END IF
							END IF
						END IF			
					END IF
				END IF
				iBlock=iBlock+1
			END DO
			IF (link_element.gt.0) THEN
!				link_flag(i_segment)=link_element
				WRITE(16,"(I6.6, I10.6, I12.8)") node_id_final(i_segment,1), node_id_final(i_segment,2), link_element
			END IF	
		END DO
		CLOSE(16)
! Check the rn and fill rn_element matrix; COL(1): iBlock, COL(2): FEM element ID
		DO i_rn=1, num_rn
			iBlock=1
			DO WHILE ((rn_element(rn_id(i_rn,1)+1,1).eq.0).and.(iBlock.le.nblock))
				IF ((element_max(iBlock,1).ge.rn_nodes(i_rn,1)).and.(element_min(iBlock,1).le.rn_nodes(i_rn,1))) THEN
					IF ((element_max(iBlock,2).ge.rn_nodes(i_rn,2)).and.(element_min(iBlock,2).le.rn_nodes(i_rn,2))) THEN
						IF ((element_max(iBlock,3).ge.rn_nodes(i_rn,3)).and.(element_min(iBlock,3).le.rn_nodes(i_rn,3))) THEN
							IF(nnode.EQ.4) THEN
								DO j=1,4
									DO k=1,3
										xx4by3(j,k)=Origin_CoordNode(iBlock,j,k)
									END DO
								END DO
								! Reset inverse vx
								DO j=1,3
									DO k=1,3
										invxv(j,k)=0.0D0
									END DO
								END DO
								!Call function to get inverse of vx
								CALL CAL_shapefunc4(xx4by3,invxv)
								DO m=1,3
									Ralra(m)=rn_nodes(i_rn,m)-xx4by3(1,m)
								END DO
								DO n=1,3
									eta(n)=invxv(n,1)*Ralra(1)+invxv(n,2)*Ralra(2)+invxv(n,3)*Ralra(3)
								END DO
								! Check if center within shape function
								IF ((eta(1).ge.0.0).and.(eta(2).ge.0.0).and.(eta(3).ge.0.0)) THEN
									IF ((eta(1)+eta(2)+eta(3)).le.shape_fn_tol_rn) THEN
										rn_element(rn_id(i_rn,1)+1,1)=iBlock
										rn_element(rn_id(i_rn,1)+1,2)=jElem(iBlock)
									END IF
								END IF						
							ELSEIF(nnode.EQ.8) THEN
								DO j=1,8
									DO k=1,3
										xx(j,k)=Origin_CoordNode(iBlock,j,k)
									END DO
								END DO
								! Reset inverse vx
								DO j=1,3
									DO k=1,3
										invxv(j,k)=0.0D0
									END DO
								END DO
								!Call function to get inverse of vx
								CALL CAL_shapefunc8(xx,invxv)
								DO m=1,3
									Ralra(m)=rn_nodes(i_rn,m)-element_center(iBlock,m)
								END DO
								DO n=1,3
									eta(n)=invxv(n,1)*Ralra(1)+invxv(n,2)*Ralra(2)+invxv(n,3)*Ralra(3)
								END DO
								! Check if center within shape function
								IF ((abs(eta(1)).le.shape_fn_tol_rn).and.(abs(eta(2)).le.shape_fn_tol_rn).and.(abs(eta(3)).le.shape_fn_tol_rn)) THEN
									rn_element(rn_id(i_rn,1)+1,1)=iBlock
									rn_element(rn_id(i_rn,1)+1,2)=jElem(iBlock)
!									WRITE (*,"(I6.2, I8.6, I8.6)") rn_id(i_rn,1), rn_element(rn_id(i_rn,1)+1,1), rn_element(rn_id(i_rn,1)+1,2)
								END IF
							END IF
						END IF			
					END IF
				END IF
				iBlock=iBlock+1
			END DO
		END DO
! Update the density of each element
		IF(nnode.EQ.8) THEN
			DO iBlock=1, nblock
				STATEOLD(iBlock,31)=0.0
				DO j=1,8
					DO k=1,3
						xx(j,k)=Origin_CoordNode(iBlock,j,k)
					END DO
				END DO	
				CALL CAL_element_volume(elem_volume, xx)
				element_volume(iBlock)=elem_volume
			END DO
			DO i_segment=1, n_segment_final
				segment_elem=link_elem_id(i_segment,1)
				IF ((segment_elem.gt.0).and.(segment_elem.le.nblock).and.(jElem(segment_elem).eq.link_elem_id(i_segment,2))) THEN				
					segment_center(1)=0.5D0*(rn1_new_final(i_segment,1)+rn2_new_final(i_segment,1))
					segment_center(2)=0.5D0*(rn1_new_final(i_segment,2)+rn2_new_final(i_segment,2))
					segment_center(3)=0.5D0*(rn1_new_final(i_segment,3)+rn2_new_final(i_segment,3))	
					DO j=1,3
						rn1_xx(j)=rn1_new_final(i_segment,j)
						rn2_xx(j)=rn2_new_final(i_segment,j)
					END DO			
					! Determne the density contribution from segment center to rn1
					rn1_id=node_id_final(i_segment,1)
					nbr_elem=rn_element(rn1_id+1,1)
					nbr_jElem=rn_element(rn1_id+1,2)!							
					seg_length=sqrt((segment_center(1)-rn1_xx(1))*(segment_center(1)-rn1_xx(1)))
					seg_length=seg_length+sqrt((segment_center(2)-rn1_xx(2))*(segment_center(2)-rn1_xx(2)))
					seg_length=seg_length+sqrt((segment_center(3)-rn1_xx(3))*(segment_center(3)-rn1_xx(3)))
					IF (jElem(segment_elem).eq.nbr_jElem) THEN ! rn1 lies in the same element as the segment center
						STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_length/(element_volume(segment_elem)*b_mag*b_mag)
					ELSE
						DO j=1,8
							DO k=1,3
								xx(j,k)=Origin_CoordNode(segment_elem,j,k)
							END DO
						END DO	
						projection(1)=0.0
						projection(2)=0.0
						projection(3)=0.0					
						CALL Find_Projection(projection, xx, segment_center, rn1_xx)
						seg_inside_length=sqrt((segment_center(1)-projection(1))*(segment_center(1)-projection(1)))
						seg_inside_length=seg_inside_length+sqrt((segment_center(2)-projection(2))*(segment_center(2)-projection(2)))
						seg_inside_length=seg_inside_length+sqrt((segment_center(3)-projection(3))*(segment_center(3)-projection(3)))
						IF (seg_inside_length.LT.seg_length) THEN
							STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_inside_length/(element_volume(segment_elem)*b_mag*b_mag)
							IF ((nbr_elem.gt.0).AND.(nbr_elem.le.nblock).and.(jElem(nbr_elem).eq.nbr_jElem)) THEN
								STATEOLD(nbr_elem,31)=STATEOLD(nbr_elem,31)+(seg_length-seg_inside_length)/(element_volume(nbr_elem)*b_mag*b_mag)
							END IF						
						ELSE
							STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_length/(element_volume(segment_elem)*b_mag*b_mag)
						END IF
					END IF
					! Determne the density contribution from segment center to rn2
					rn2_id=node_id_final(i_segment,2)
					nbr_elem=rn_element(rn2_id+1,1)
					nbr_jElem=rn_element(rn2_id+1,2)
					IF (jElem(segment_elem).eq.nbr_jElem) THEN ! rn2 lies in the same element as the segment center
						STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_length/(element_volume(segment_elem)*b_mag*b_mag)
					ELSE		
						projection(1)=0.0
						projection(2)=0.0
						projection(3)=0.0				
						CALL Find_Projection(projection, xx, segment_center, rn2_xx)
						seg_inside_length=sqrt((segment_center(1)-projection(1))*(segment_center(1)-projection(1)))
						seg_inside_length=seg_inside_length+sqrt((segment_center(2)-projection(2))*(segment_center(2)-projection(2)))
						seg_inside_length=seg_inside_length+sqrt((segment_center(3)-projection(3))*(segment_center(3)-projection(3)))
						IF (seg_inside_length.LT.seg_length) THEN
							STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_inside_length/(element_volume(segment_elem)*b_mag*b_mag)
							IF ((nbr_elem.gt.0).AND.(nbr_elem.le.nblock).and.(jElem(nbr_elem).eq.nbr_jElem)) THEN
								STATEOLD(nbr_elem,31)=STATEOLD(nbr_elem,31)+(seg_length-seg_inside_length)/(element_volume(nbr_elem)*b_mag*b_mag)
							END IF						
						ELSE
							STATEOLD(segment_elem,31)=STATEOLD(segment_elem,31)+seg_length/(element_volume(segment_elem)*b_mag*b_mag)
						END IF
					END IF	
				END IF			
			END DO
		END IF
! Perform surface node projection on rn flagged 6 and print out the new coordinates of the surface node
! If a node is found in side FEM domain, the coordinates is print out as three 0's
		filename_rn_check=TRIM(foldername)//"/"//TRIM(file_rn_head)//TRIM(filename_output)//".txt"
		OPEN(16,file=TRIM(filename_rn_check),POSITION='append')		
		DO i_rn=1, num_rn
			rn_surface_old(1)=0.0
			rn_surface_old(2)=0.0
			rn_surface_old(3)=0.0		
			IF (rn_element(rn_id(i_rn,1)+1,1).gt.0) THEN ! node found inside FEM domain
				IF (rn_id(i_rn,2).EQ.6) THEN
					i_segment=1
					rn_surface_check=0
					DO WHILE ((i_segment.LE.n_segment_final).AND.(rn_surface_check.EQ.0))
						IF (rn_id(i_rn,1).EQ.node_id_final(i_segment,1)) THEN
							rn_surface_old(1)=rn1_old_final(i_segment,1)
							rn_surface_old(2)=rn1_old_final(i_segment,2)
							rn_surface_old(3)=rn1_old_final(i_segment,3)
							rn_surface_check=1
						ELSEIF (rn_id(i_rn,1).EQ.node_id_final(i_segment,2)) THEN
							rn_surface_old(1)=rn2_old_final(i_segment,1)
							rn_surface_old(2)=rn2_old_final(i_segment,2)
							rn_surface_old(3)=rn2_old_final(i_segment,3)
							rn_surface_check=1
						END IF
						i_segment=i_segment+1
					END DO				
				END IF
!				rn_surface_old(1)=0.0
!				rn_surface_old(2)=0.0
!				rn_surface_old(3)=0.0
				if (totalTime > 0.19347) then 
					rn_id(i_rn,2)=7
				endif
				WRITE(16,"(I6.6, I6.2, E16.8E2, E16.8E2, E16.8E2)") rn_id(i_rn,1), rn_id(i_rn,2), rn_surface_old(1), rn_surface_old(2), rn_surface_old(3)
			ELSEIF (rn_element(rn_id(i_rn,1)+1,1).EQ.0) THEN ! node not found in FEM domain
				IF ((rn_id(i_rn,2).EQ.6)) THEN
					projection(1)=0.0
					projection(2)=0.0
					projection(3)=0.0
					! Find previous coordinates of the surface node
					i_segment=1
					rn_surface_check=0
					DO WHILE ((i_segment.LE.n_segment_final).AND.(rn_surface_check.EQ.0))
						IF (rn_id(i_rn,1).EQ.node_id_final(i_segment,1)) THEN
							rn_surface_old(1)=rn1_old_final(i_segment,1)
							rn_surface_old(2)=rn1_old_final(i_segment,2)
							rn_surface_old(3)=rn1_old_final(i_segment,3)
							rn_surface_check=1
						ELSEIF (rn_id(i_rn,1).EQ.node_id_final(i_segment,2)) THEN
							rn_surface_old(1)=rn2_old_final(i_segment,1)
							rn_surface_old(2)=rn2_old_final(i_segment,2)
							rn_surface_old(3)=rn2_old_final(i_segment,3)
							rn_surface_check=1
						END IF					
						i_segment=i_segment+1
					END DO					
					IF ((rn_id(i_rn,3).LE.1).AND.(nnode.EQ.8)) THEN ! Find projection for this surface node
						nbr_elem=0
						nbr_jElem=0
						i_segment=1
						rn_surface_check=0
						DO WHILE ((i_segment.LE.n_segment_final).AND.(rn_surface_check.EQ.0))
							IF (rn_id(i_rn,1).EQ.node_id_final(i_segment,1)) THEN
								nbr_elem = rn_element(node_id_final(i_segment,2)+1,1)
								nbr_jElem = rn_element(node_id_final(i_segment,2)+1,2)
								DO j=1,3
									rn1_xx(j)=rn2_new_final(i_segment,j)
									rn2_xx(j)=rn1_new_final(i_segment,j)
								END DO
								rn_surface_check=1			
							ELSEIF (rn_id(i_rn,1).EQ.node_id_final(i_segment,2)) THEN
								nbr_elem = rn_element(node_id_final(i_segment,1)+1,1)
								nbr_jElem = rn_element(node_id_final(i_segment,1)+1,2)
								DO j=1,3
									rn1_xx(j)=rn1_new_final(i_segment,j)
									rn2_xx(j)=rn2_new_final(i_segment,j)
								END DO							
								rn_surface_check=1
							END IF					
							i_segment=i_segment+1
						END DO
!						WRITE (16,"(I6.2, I8.6, I8.6)") rn_id(i_rn,1), nbr_elem, nbr_jElem
						IF ((nbr_elem.GT.0).AND.(nbr_elem.LE.nblock)) THEN
							IF ((jElem(nbr_elem).EQ.nbr_jElem).AND.(surface_elements(nbr_jElem,2) .GE. 1)) THEN
								DO j=1,8
									DO k=1,3
										xx(j,k)=Origin_CoordNode(nbr_elem,j,k)
									END DO
								END DO
								!Call function to find projection of surface node on surface plane (rn2_xx is coordinates of the surface node)
								CALL Find_Projection(projection, xx, rn1_xx, rn2_xx)
								rn_old_dot_new=(projection(1)-rn1_xx(1))*(rn2_xx(1)-rn1_xx(1))+(projection(2)-rn1_xx(2))*(rn2_xx(2)-rn1_xx(2))+(projection(3)-rn1_xx(3))*(rn2_xx(3)-rn1_xx(3))
								!
								seg_length=sqrt((rn2_xx(1)-rn1_xx(1))*(rn2_xx(1)-rn1_xx(1))+(rn2_xx(2)-rn1_xx(2))*(rn2_xx(2)-rn1_xx(2))+(rn2_xx(3)-rn1_xx(3))*(rn2_xx(3)-rn1_xx(3)))
								new_seg_length=sqrt((projection(1)-rn1_xx(1))*(projection(1)-rn1_xx(1))+(projection(2)-rn1_xx(2))*(projection(2)-rn1_xx(2))+(projection(3)-rn1_xx(3))*(projection(3)-rn1_xx(3)))
								IF ((rn_old_dot_new.LE.0.0).OR.(new_seg_length.GT.seg_length).OR.(new_seg_length.LT.80.0)) THEN
									projection(1)=0.0
									projection(2)=0.0
									projection(3)=0.0
								END IF
							END IF
						END IF
					ELSEIF (rn_id(i_rn,3).GT.1) THEN
						projection(1)=rn_surface_old(1)
						projection(2)=rn_surface_old(2)
						projection(3)=rn_surface_old(3)		
					END IF
					IF ((abs(projection(1)).GT.0.0).OR.(abs(projection(2)).GT.0.0).OR.(abs(projection(3)).GT.0.0)) THEN
						rn_surface_old(1)=projection(1)
						rn_surface_old(2)=projection(2)
						rn_surface_old(3)=projection(3)
						WRITE(16,"(I6.6, I6.2, E16.8E2, E16.8E2, E16.8E2)") rn_id(i_rn,1), rn_id(i_rn,2), projection(1), projection(2), projection(3)
					ELSE
!						WRITE(16,"(I6.6, I6.2, E16.8E2, E16.8E2, E16.8E2)") rn_id(i_rn,1), rn_id(i_rn,2), rn_surface_old(1), rn_surface_old(2), rn_surface_old(3)
					END IF					
				END IF
			END IF
		END DO
		CLOSE(16)
!Output element center
		filename_element=TRIM(foldername)//"/"//TRIM(file_element_head)//TRIM(filename_output)//".txt"
		filename_output_element=TRIM(foldername)//"/"//TRIM(file_element_output_head)//TRIM(filename_output)//".txt"
		OPEN(16,file=TRIM(filename_element),POSITION='append')
		OPEN(17,file=TRIM(filename_output_element),POSITION='append')
		DO iBlock=1,nblock
			WRITE(17,"(I8.8, E16.8E2, E16.8E2, E16.8E2)") jElem(iBlock), element_center(iBlock,1), element_center(iBlock,2), element_center(iBlock,3)
			IF (surface_elements(jElem(iBlock),2) .ge. 1) THEN
				WRITE(16,"(I8.8, E16.8E2, E16.8E2, E16.8E2)") jElem(iBlock), element_center(iBlock,1), element_center(iBlock,2), element_center(iBlock,3)
			END IF
		END DO
		CLOSE(16)
		CLOSE(17)
	ELSE
!		DO i_segment=1, n_segment
!			link_flag(i_segment)=0
!		END DO
	END IF
!	
	RETURN
END SUBROUTINE vucharlength
!==================================== ====================================
!
!
!====================================================================
!     Input dislocation segment information from external db 
!====================================================================
SUBROUTINE vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
!    INCLUDE 'VABA_PARAM.INC'
	USE IFPORT
	Implicit double precision (a-h,o-z)
	INTEGER, PARAMETER :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length
!	The above two lines are needed only for ABAQUS 2018!
!
!	Contents of i_Array
	parameter( i_int_nTotalNodes     = 1,&
				i_int_nTotalElements  = 2,&
				i_int_kStep           = 3,&
				i_int_kInc            = 4,&
				i_int_iStatus         = 5,&
				i_int_lWriteRestart   = 6  )
!
!	Possible values for the lOp argument
	parameter( j_int_StartAnalysis    = 0,&
				j_int_StartStep        = 1,&
				j_int_SetupIncrement   = 2,&
				j_int_StartIncrement   = 3,&
				j_int_EndIncrement     = 4,&
				j_int_EndStep          = 5,&
				j_int_EndAnalysis      = 6 )     
!
!	Possible values for i_Array(i_int_iStatus)
	parameter( j_int_Continue          = 0,&
				j_int_TerminateStep     = 1,&
				j_int_TerminateAnalysis = 2)      
!
!	Contents of r_Array
	parameter( i_flt_TotalTime   = 1,&
				i_flt_StepTime    = 2,&
				i_flt_dTime       = 3 )
!
	dimension i_Array(niArray),r_Array(nrArray)
!
	DOUBLE PRECISION:: normalizer
	DOUBLE PRECISION, DIMENSION(3,3):: tensor
	DOUBLE PRECISION, DIMENSION(3,3):: tensor_final
	CHARACTER(*), PARAMETER :: filename_segs_num = "Data_segs_num.txt"
	CHARACTER(*), PARAMETER :: filename_segs = "Data_segs.txt"
	CHARACTER(*), PARAMETER :: filename_segs_final = "Data_segs_final.txt"
	CHARACTER(*), PARAMETER :: filename_tensor = "Data_tensor.txt"
	CHARACTER(*), PARAMETER :: filename_ABAQUS_flag = "ABAQUSworking"
	CHARACTER(*), PARAMETER :: filename_Matlab_flag = "Matlabworking"
	CHARACTER(*), PARAMETER :: filename_hold_stress_flag = "HoldStress"
	CHARACTER(*), PARAMETER :: filename_rn_num = "nodes_num.txt"
	CHARACTER(*), PARAMETER :: filename_FEM_deltaU = "FEM_deltaU.txt"
	CHARACTER(*), PARAMETER :: filename_rn_nodes = "rn_matrix.txt"
!	CHARACTER(*), PARAMETER :: filename_elements_num = "elements_num.txt"
	CHARACTER(*), PARAMETER :: filename_surface_elements = "surface_elements.txt"
	CHARACTER(*), PARAMETER :: file_CoordNode = "OriginalCoordNodes.txt"
	CHARACTER(LEN=1024) :: foldername, filename_flag
	INTEGER :: KPROCESSNUM
	INTEGER :: K1, K2, K3, LENOUTDIR, MATLABWORKING
	INTEGER :: LoadType 
	LOGICAL FILEEXISTANCE,FILECHANGERESULT, HOLD_ONE, HOLD_TWO
!
	INTEGER, PARAMETER :: n_segment_max=80000
	INTEGER, PARAMETER :: n_rn_max=80000
	INTEGER, DIMENSION (n_segment_max,2):: node_id, node_id_final
	INTEGER, DIMENSION (n_rn_max,3):: rn_id
	INTEGER, PARAMETER :: n_elements_max=600000
	INTEGER, DIMENSION (n_elements_max,2):: surface_elements
	DOUBLE PRECISION, DIMENSION (n_elements_max,8,3):: original_coordNodes
	INTEGER, DIMENSION (n_rn_max, 2):: rn_element
	INTEGER:: n_segment, n_segment_final, k_iteration, num_rn, num_elements
	DOUBLE PRECISION :: EffSlip
	DOUBLE PRECISION:: time_now, t_step, t_nextstep, paradis_nextdt, FEM_Du
	! t_step: total time increment in the current step	
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: rn1_old, rn2_old, rn1_new, rn2_new, rn1_old_final, rn2_old_final, rn1_new_final, rn2_new_final
	INTEGER, DIMENSION (n_segment_max, 2):: link_elem_id
	DOUBLE PRECISION, DIMENSION(n_rn_max,3):: rn_nodes
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: swept_normal, swept_center
	! DOUBLE PRECISION, DIMENSION(n_segment_max,3):: xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8
	DOUBLE PRECISION, DIMENSION(n_segment_max,3):: dw_p
	DOUBLE PRECISION, DIMENSION(n_segment_max,6):: de_p
	COMMON /sharedf/ de_p, dw_p, t_step, t_nextstep, EffSlip, rn1_old, rn2_old, rn1_new, rn2_new, rn1_old_final, rn2_old_final, rn1_new_final, rn2_new_final, link_elem_id!
	COMMON /sharedi/ n_segment, n_segment_final, k_iteration
!	COMMON /sharedamp/ HOLD_ONE, HOLD_TWO
	COMMON /sharenode/ node_id_final, rn_id, rn_nodes, num_rn, rn_element
	COMMON /sharedsurface/ surface_elements, num_elements
	COMMON /sharednextdt/ LoadType, FEM_Du
	COMMON /sharedoriginalnodes/ original_coordNodes
	kStep = i_Array(i_int_kStep)
	kInc  = i_Array(i_int_kInc)
	num_elements=i_Array(i_int_nTotalElements)
!
!--------------------------------------------------------------------
!	Start of the analysis
	IF (lOp .eq. j_int_StartAnalysis) THEN
		t_nextstep=0.0
		CALL VGETOUTDIR(foldername, LENOUTDIR)
		foldername=foldername(1:LENOUTDIR)
!		Import surface elements flags
		OPEN(16,file=TRIM(foldername)//"/"//filename_surface_elements)
		DO K1=1,num_elements
			READ(16,*) (surface_elements(K1,K2), K2=1,2)
		END DO
		CLOSE(16)
! 		Import original nodal coordinates
		OPEN(16,file=TRIM(foldername)//"/"//file_CoordNode)
		DO K1=1,num_elements
			READ(16,*) K3, (original_coordNodes(K1,1,K2), K2=1,3), (original_coordNodes(K1,2,K2), K2=1,3), (original_coordNodes(K1,3,K2), K2=1,3), (original_coordNodes(K1,4,K2), K2=1,3), (original_coordNodes(K1,5,K2), K2=1,3), (original_coordNodes(K1,6,K2), K2=1,3), (original_coordNodes(K1,7,K2), K2=1,3), (original_coordNodes(K1,8,K2), K2=1,3)
		END DO
		CLOSE(16)
	END IF
!
	IF (lOp .eq. j_int_StartIncrement) THEN
		k_iteration=kInc
	END IF
!
	time_now = r_Array(1) ! the current total simulation time
	IF(time_now.GT.t_nextstep) THEN!Run ParaDis through Matlab and wait until ABAQUSworking.flag exists
		CALL VGETOUTDIR(foldername, LENOUTDIR)
		CALL VGETRANK(KPROCESSNUM)
		foldername=foldername(1:LENOUTDIR)
		write (filename_flag, "(I4.4)") KPROCESSNUM+1
		! change ABAQUSworking.flag to Matlabworking.flag
!		FILECHANGERESULT=SYSTEMQQ("move "//TRIM(foldername)//"\"//filename_ABAQUS_flag//TRIM(filename_flag)//".flag "//TRIM(foldername)//"\"//filename_Matlab_flag//TRIM(filename_flag)//".flag")
		FILECHANGERESULT=SYSTEMQQ("mv "//TRIM(foldername)//"/"//filename_ABAQUS_flag//TRIM(filename_flag)//".flag "//TRIM(foldername)//"/"//filename_Matlab_flag//TRIM(filename_flag)//".flag")
		!
		MATLABWORKING=1
		DO WHILE (MATLABWORKING>0)
			CALL SLEEPQQ(250)
			INQUIRE (FILE=TRIM(foldername)//"/"//filename_ABAQUS_flag//TRIM(filename_flag)//".flag", EXIST=FILEEXISTANCE)
			IF(FILEEXISTANCE) THEN
				MATLABWORKING=0
			END IF
		END DO
		!
		OPEN(16,file=TRIM(foldername)//"/"//filename_segs_num)
		READ(16,*) n_segment, n_segment_final, t_step, EffSlip
		CLOSE(16)
		OPEN(16,file=TRIM(foldername)//"/"//filename_segs)
		DO K1=1,n_segment
			READ(16,*) (node_id(K1,K2), K2=1,2), (rn1_old(K1,K2), K2=1,3), (rn2_old(K1,K2), K2=1,3), (rn1_new(K1,K2), K2=1,3),&
                       (rn2_new(K1,K2), K2=1,3),(tensor(1,K2), K2=1,3),(tensor(2,K2), K2=1,3),(tensor(3,K2), K2=1,3)
			IF (ISNAN(tensor(1,1))) THEN
				de_p(K1,1)=0.0
				de_p(K1,2)=0.0
				de_p(K1,3)=0.0
				de_p(K1,4)=0.0
				de_p(K1,5)=0.0
				de_p(K1,6)=0.0
				!
				dw_p(K1,1)=0.0
				dw_p(K1,2)=0.0
				dw_p(K1,3)=0.0
			ELSE
				de_p(K1,1)=tensor(1,1)
				de_p(K1,2)=tensor(2,2)
				de_p(K1,3)=tensor(3,3)
				de_p(K1,4)=0.5*(tensor(1,2)+tensor(2,1))!ep12 Not Engineering strain quantities for VUMAT
				de_p(K1,5)=0.5*(tensor(2,3)+tensor(3,2))!ep23 Not Engineering strain quantities 
				de_p(K1,6)=0.5*(tensor(3,1)+tensor(1,3))!ep31 Not Engineering strain quantities
				!
				dw_p(K1,1)=0.25*(tensor(1,2)-tensor(2,1))!w12 !Here half of the spin tensor is used for numerical efficiency
				dw_p(K1,2)=0.25*(tensor(1,3)-tensor(3,1))!w13
				dw_p(K1,3)=0.25*(tensor(2,3)-tensor(3,2))!w23
				!
			END IF
			!
		END DO
		CLOSE(16)
! import the segments from the latest cycle in ParaDiS 		
		OPEN(16,file=TRIM(foldername)//"/"//filename_segs_final)
		DO K1=1,n_segment_final
!			READ(16,*) (node_id_final(K1,K2), K2=1,2), (rn1_old_final(K1,K2), K2=1,3), (rn2_old_final(K1,K2), K2=1,3), (rn1_new_final(K1,K2), K2=1,3), (rn2_new_final(K1,K2), K2=1,3)
			READ(16,*) (node_id_final(K1,K2), K2=1,2), (rn1_old_final(K1,K2), K2=1,3), (rn2_old_final(K1,K2), K2=1,3), &
                       (rn1_new_final(K1,K2), K2=1,3),(rn2_new_final(K1,K2), K2=1,3),(tensor_final(1,K2), K2=1,3),&
                       (tensor_final(2,K2),K2=1,3),(tensor_final(3,K2), K2=1,3)
		END DO
		DO K1=1, n_segment_final
			link_elem_id(K1,1)=0
			link_elem_id(K1,2)=0
		END DO
		! Get rn
		OPEN(16,file=TRIM(foldername)//"/"//filename_rn_num)
		READ(16,*) num_rn
		CLOSE(16)
		! Read loading condition. 1: dynamic (synchronization). 2: quasi-static
		OPEN(16,file=TRIM(foldername)//"/"//filename_FEM_deltaU)
		READ(16,*) LoadType, FEM_Du
		CLOSE(16)
		
		OPEN(17,file=TRIM(foldername)//"/"//filename_rn_nodes)
! rn_id:    COL(1) node ID starts from 0, COL(2) flag, COL(3) num. neighbors	
! rn_nodes: COL(1-3): x, y, z coordinates		
		DO K1=1,num_rn
			READ(17,*) (rn_id(K1,K2), K2=1,3), (rn_nodes(K1,K2), K2=1,3)
		END DO
		DO K1=1, num_rn
			rn_element(rn_id(K1,1)+1,1)=0
			rn_element(rn_id(K1,1)+1,2)=0
		END DO
		CLOSE(17)
! VUAMP input			
		INQUIRE (FILE=TRIM(foldername)//"/"//filename_hold_stress_flag//"01.flag", EXIST=HOLD_ONE)
		INQUIRE (FILE=TRIM(foldername)//"/"//filename_hold_stress_flag//"02.flag", EXIST=HOLD_TWO)
		t_nextstep=time_now+t_step
	END IF
!
	RETURN
END SUBROUTINE vexternaldb
!====================================================================================
!
!====================================================================================
