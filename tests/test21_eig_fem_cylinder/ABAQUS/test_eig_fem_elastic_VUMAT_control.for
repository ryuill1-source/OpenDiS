      subroutine vumat(
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'
	!IMPLICIT DOUBLE PRECISION (a-h,o-z)
	CHARACTER(*), PARAMETER :: filename_ABAQUS_stress_ready_flag = "ABAQUS_stress_ready"
	CHARACTER(*), PARAMETER :: filename_ABAQUS_pause_flag = "ABAQUS_pause"
      CHARACTER(*), PARAMETER :: filename_ABAQUS_running_flag = "ABAQUS_running"
      
C     [CAUTION] foldername should be ABSOLUTE PATH !!!!
      CHARACTER(*), PARAMETER :: foldername = "/home/kyeongmi/Codes/OpenDiS.eig_fem.git/tests/test21_eig_fem_cylinder/ABAQUS"
      
      CHARACTER(*), PARAMETER :: filename_flag = ""
      CHARACTER(LEN=256) :: full_filename_pause
      CHARACTER(LEN=256) :: full_filename_stress_ready
      CHARACTER(LEN=256) :: full_filename_running
      LOGICAL :: FILEEXISTANCE

C     Elastic material model only

      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*), charLength(*), strainInc(nblock,ndir+nshr),
     2  relSpinInc(*), tempOld(*), stretchOld(*), defgradOld(*),
     3  fieldOld(*), stressOld(nblock,ndir+nshr),
     4  stateOld(nblock,nstatev), enerInternOld(nblock),
     5  enerInelasOld(nblock), tempNew(*), stretchNew(*),
     6  defgradNew(*), fieldNew(*), stressNew(nblock,ndir+nshr),
     7  stateNew(nblock,nstatev), enerInternNew(nblock),
     8  enerInelasNew(nblock)

      character*80 cmname

      parameter (zero=0.0d0, one=1.0d0)

C     Read properties: Young's modulus and Poisson's ratio
      E = props(1)
      NU = props(2)

C     Compute Lame constants
      G = E / (2.0d0 * (1.0d0 + NU))
      LAMBDA = E * NU / ((1.0d0 + NU)*(1.0d0 - 2.0d0 * NU))

      do i = 1, nblock

C       Compute strain trace (volumetric strain)
        trace_eps = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)

C       Compute stress update: sigma = sigma_old + lambda*trace*I + 2G*eps
        stressNew(i,1) = stressOld(i,1) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,1)
        stressNew(i,2) = stressOld(i,2) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,2)
        stressNew(i,3) = stressOld(i,3) + LAMBDA*trace_eps + 2.0d0*G*strainInc(i,3)
        stressNew(i,4) = stressOld(i,4) + 2.0d0* G * strainInc(i,4)
        if (nshr .ge. 2) stressNew(i,5) = stressOld(i,5) + 2.0d0 * G * strainInc(i,5)
        if (nshr .ge. 3) stressNew(i,6) = stressOld(i,6) + 2.0d0 * G * strainInc(i,6)

C       Copy state variables unchanged
        do j = 1, nstatev
          stateNew(i,j) = stateOld(i,j)
        end do

      end do


C  The following needs to be placed in vexternalDB
C  Remove ABAQUS_running.flag & Making ABAQUS_pause.flag, ABAQUS_stress_ready.flag *********************
      full_filename_running = TRIM(foldername)//"/"//filename_ABAQUS_running_flag//TRIM(filename_flag)//".flag"
      OPEN(UNIT=97, FILE=TRIM(full_filename_running), STATUS="REPLACE", ACTION="WRITE")
      CLOSE(97, STATUS="DELETE")
      PRINT *, "VUMAT deleted empty file: ", TRIM(filename_ABAQUS_running_flag)

      full_filename_stress_ready = TRIM(foldername)//"/"//filename_ABAQUS_stress_ready_flag//TRIM(filename_flag)//".flag"
      OPEN(UNIT=98, FILE=TRIM(full_filename_stress_ready), STATUS="REPLACE", ACTION="WRITE")
      CLOSE(98)
      PRINT *, "VUMAT made empty file: ", TRIM(filename_ABAQUS_stress_ready_flag)
      
      full_filename_pause = TRIM(foldername)//"/"//filename_ABAQUS_pause_flag//TRIM(filename_flag)//".flag"
      OPEN(UNIT=99, FILE=TRIM(full_filename_pause), STATUS="REPLACE", ACTION="WRITE")
      CLOSE(99)
      PRINT *, "VUMAT made empty file: ", TRIM(filename_ABAQUS_pause_flag)
C  *************************************************************************

C  Control ABAQUS run (pause&go)  ***********************************  
	DO WHILE (.TRUE.)
            PRINT *, "Checking for file: ", TRIM(filename_ABAQUS_pause_flag)

            INQUIRE (FILE=full_filename_pause, EXIST=FILEEXISTANCE)
		
            IF(FILEEXISTANCE) THEN
                  PRINT *, "ABAQUS_pause.flag exists"
		      CALL SLEEPQQ(1000)
            ELSE
                  PRINT *, "ABAQUS_pause.flag NOT exists, keep calculation"
			EXIT
            END IF
	END DO
C  End of Control ABAQUS run (pause&go)  *********************************** 

      return
      end
