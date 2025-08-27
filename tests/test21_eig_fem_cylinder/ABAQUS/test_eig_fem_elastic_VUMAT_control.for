!DIR$ FREEFORM
subroutine vumat(&
    nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
    stepTime, totalTime, dt, cmname, coordMp, charLength, &
    props, density, strainInc, relSpinInc, &
    tempOld, stretchOld, defgradOld, fieldOld, &
    stressOld, stateOld, enerInternOld, enerInelasOld, &
    tempNew, stretchNew, defgradNew, fieldNew, &
    stressNew, stateNew, enerInternNew, enerInelasNew)

    include 'vaba_param.inc'
    !IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
    !Elastic material model only

    dimension props(nprops), density(nblock), &
        coordMp(nblock,*), charLength(*), strainInc(nblock,ndir+nshr), &
        relSpinInc(*), tempOld(*), stretchOld(*), defgradOld(*), &
        fieldOld(*), stressOld(nblock,ndir+nshr), &
        stateOld(nblock,nstatev), enerInternOld(nblock), &
        enerInelasOld(nblock), tempNew(*), stretchNew(*), &
        defgradNew(*), fieldNew(*), stressNew(nblock,ndir+nshr), &
        stateNew(nblock,nstatev), enerInternNew(nblock), &
        enerInelasNew(nblock)

    character*80 cmname
    parameter (zero=0.0d0, one=1.0d0)

    !Read properties: Young's modulus and Poisson's ratio
    E = props(1)
    NU = props(2)

    !Compute Lame constants
    G = E / (2.0d0 * (1.0d0 + NU))
    LAMBDA = E * NU / ((1.0d0 + NU)*(1.0d0 - 2.0d0 * NU))

    do i = 1, nblock

        !Compute strain trace (volumetric strain)
        trace_eps = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)

        !Compute stress update: sigma = sigma_old + lambda*trace*I + 2G*eps
        stressNew(i,1) = stressOld(i,1) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,1)
        stressNew(i,2) = stressOld(i,2) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,2)
        stressNew(i,3) = stressOld(i,3) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,3)
        stressNew(i,4) = stressOld(i,4) + 2.0d0* G * strainInc(i,4)
        if (nshr .ge. 2) stressNew(i,5) = stressOld(i,5) + 2.0d0 * G * strainInc(i,5)
        if (nshr .ge. 3) stressNew(i,6) = stressOld(i,6) + 2.0d0 * G * strainInc(i,6)

        !Copy state variables unchanged
        do j = 1, nstatev
            stateNew(i,j) = stateOld(i,j)
        end do
    end do

    return
end subroutine vumat

subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
    USE IFPORT
    IMPLICIT NONE
    INTEGER            :: lOp, niArray, nrArray
    INTEGER            :: i_Array(niArray)
    DOUBLE PRECISION   :: r_Array(nrArray)
    LOGICAL            :: FILEEXISTANCE
    INTEGER            :: kInc
!  i_Array index map
    INTEGER, PARAMETER :: i_int_nTotalNodes    = 1
    INTEGER, PARAMETER :: i_int_nTotalElements = 2
    INTEGER, PARAMETER :: i_int_kStep          = 3
    INTEGER, PARAMETER :: i_int_kInc           = 4
    INTEGER, PARAMETER :: i_int_iStatus        = 5
    INTEGER, PARAMETER :: i_int_lWriteRestart  = 6
!  lOp codes
    INTEGER, PARAMETER :: j_int_StartAnalysis  = 0
    INTEGER, PARAMETER :: j_int_StartStep      = 1
    INTEGER, PARAMETER :: j_int_SetupIncrement = 2
    INTEGER, PARAMETER :: j_int_StartIncrement = 3
    INTEGER, PARAMETER :: j_int_EndIncrement   = 4
    INTEGER, PARAMETER :: j_int_EndStep        = 5
    INTEGER, PARAMETER :: j_int_EndAnalysis    = 6
!  r_Array indices
    INTEGER, PARAMETER :: i_flt_TotalTime=1, i_flt_StepTime=2, i_flt_dTime=3
!  flag
	CHARACTER(*), PARAMETER :: f_stress_ready_flag = "ABAQUS_stress_ready"
	CHARACTER(*), PARAMETER :: f_pause_flag = "ABAQUS_pause"
    CHARACTER(*), PARAMETER :: f_running_flag = "ABAQUS_running"
!   [CAUTION] foldername should be ABSOLUTE PATH !!!!
    CHARACTER(*), PARAMETER :: foldername = "/home/wnswjswk1/Codes/OpenDiS/tests/test21_eig_fem_cylinder/ABAQUS"
    CHARACTER(*), PARAMETER :: f_flag = ""
    CHARACTER(LEN=256) :: full_f_pause
    CHARACTER(LEN=256) :: full_f_stress_ready
    CHARACTER(LEN=256) :: full_f_running

    full_f_pause = TRIM(foldername)//"/"//TRIM(f_pause_flag)//TRIM(f_flag)//".flag"
    full_f_stress_ready = TRIM(foldername)//"/"//TRIM(f_stress_ready_flag)//TRIM(f_flag)//".flag"
    full_f_running = TRIM(foldername)//"/"//TRIM(f_running_flag)//TRIM(f_flag)//".flag"
    
    SELECT CASE (lOp)
    
    CASE (j_int_EndIncrement)
        PRINT *, "EndInc: kInc=", i_Array(i_int_kInc), &
           " StepTime=", r_Array(i_flt_StepTime), " dt=", r_Array(i_flt_dTime)
        ! Check the step number
        kInc = i_Array(i_int_kInc)
        IF (kInc == 1) THEN
        ! remove ABAQUS_running.flag & write ABAQUS_pause.flag only
            PRINT *, "[EndStep] kInc = 1. write ABAQUS_PAUSE.flag only"
            
            OPEN(UNIT=97, FILE=TRIM(full_f_running), STATUS="REPLACE", ACTION="WRITE")
            CLOSE(97, STATUS="DELETE")
            PRINT *, "VEXTXERNALDB deleted empty file: ", TRIM(f_running_flag)
            
            OPEN(UNIT=99, FILE=TRIM(full_f_pause), STATUS="REPLACE", ACTION="WRITE")
            CLOSE(99)
            PRINT *, "VEXTERNALDB made empty file: ", TRIM(f_pause_flag)
            
            DO WHILE (.TRUE.)
                PRINT *, "Checking for file: ", TRIM(f_pause_flag)

                INQUIRE (FILE=full_f_pause, EXIST=FILEEXISTANCE)
            
                IF(FILEEXISTANCE) THEN
                    PRINT *, "ABAQUS_pause.flag exists"
                    CALL SLEEPQQ(1000)
                ELSE
                    PRINT *, "ABAQUS_pause.flag NOT exists, keep calculation"
                    EXIT
                END IF
            END DO
        ELSE
        ! remove ABAQUS_running.flag & write ABAQUS_pause.flag, ABAQUS_stress_ready.flag
            PRINT *, "[EndStep] kInc = ", kInc, ". write ABAQUS_PAUSE.flag and ABAQUS_stress_ready.flag"
            
            OPEN(UNIT=97, FILE=TRIM(full_f_running), STATUS="REPLACE", ACTION="WRITE")
            CLOSE(97, STATUS="DELETE")
            PRINT *, "VEXTERNALDB deleted empty file: ", TRIM(f_running_flag)
            
            OPEN(UNIT=99, FILE=TRIM(full_f_pause), STATUS="REPLACE", ACTION="WRITE")
            CLOSE(99)
            PRINT *, "VEXTERNALDB made empty file: ", TRIM(f_pause_flag)
     
            OPEN(UNIT=99, FILE=TRIM(full_f_stress_ready), STATUS="REPLACE", ACTION="WRITE")
            CLOSE(99)
            PRINT *, "VEXTERNALDB made empty file: ", TRIM(f_stress_ready_flag)    
            
            DO WHILE (.TRUE.)
                PRINT *, "Checking for file: ", TRIM(f_pause_flag)

                INQUIRE (FILE=full_f_pause, EXIST=FILEEXISTANCE)
            
                IF(FILEEXISTANCE) THEN
                    PRINT *, "ABAQUS_pause.flag exists"
                    CALL SLEEPQQ(1000)
                ELSE
                    PRINT *, "ABAQUS_pause.flag NOT exists, keep calculation"
                    EXIT
                END IF
            END DO
        ENDIF

    CASE (j_int_EndAnalysis)
    ! Clean up the flags
        PRINT *, "[EndAnalysis] Clean up flags ..."
        
        INQUIRE(FILE=TRIM(full_f_pause), EXIST=FILEEXISTANCE)
        IF (FILEEXISTANCE) THEN
            OPEN(99, FILE=TRIM(full_f_pause), STATUS="OLD")
            CLOSE(99, STATUS="DELETE")
            PRINT *, "Deleted: ", TRIM(full_f_pause)
        END IF

        INQUIRE(FILE=TRIM(full_f_stress_ready), EXIST=FILEEXISTANCE)
        IF (FILEEXISTANCE) THEN
            OPEN(98, FILE=TRIM(full_f_stress_ready), STATUS="OLD")
            CLOSE(98, STATUS="DELETE")
            PRINT *, "Deleted: ", TRIM(full_f_stress_ready)
        END IF

        INQUIRE(FILE=TRIM(full_f_running), EXIST=FILEEXISTANCE)
        IF (FILEEXISTANCE) THEN
            OPEN(97, FILE=TRIM(full_f_running), STATUS="OLD")
            CLOSE(97, STATUS="DELETE")
            PRINT *, "Deleted: ", TRIM(full_f_running)
        END IF

        PRINT *, "TotalTime=", r_Array(i_flt_TotalTime), &
                 ", LastStepTime=", r_Array(i_flt_StepTime)

        CALL FLUSH(6)
    END SELECT
    RETURN
END subroutine vexternaldb