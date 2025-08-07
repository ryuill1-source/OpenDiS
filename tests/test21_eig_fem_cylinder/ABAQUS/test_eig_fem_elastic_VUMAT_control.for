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
	CHARACTER(*), PARAMETER :: filename_ABAQUS_flag = "ABAQUSrunning"
      
C     [CAUTION] foldername should be ABSOLUTE PATH !!!!
      CHARACTER(*), PARAMETER :: foldername = "/home/kyeongmi/Codes/OpenDiS.eig_fem.git/tests/test21_eig_fem_cylinder/ABAQUS"
	
      CHARACTER(*), PARAMETER :: filename_flag = ""
      CHARACTER(LEN=200) :: full_filename	
      INTEGER :: MATLABWORKING = 10
	
      LOGICAL FILEEXISTANCE

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

C  Control ABAQUS run (pause&go)  ***********************************V   
C  The following needs to be placed in vexternalDB
	DO WHILE (MATLABWORKING>0)
		full_filename = TRIM(foldername)//"/"//filename_ABAQUS_flag//TRIM(filename_flag)//".flag"
            PRINT *, "Checking for file: ", TRIM(full_filename)

            INQUIRE (FILE=full_filename, EXIST=FILEEXISTANCE)
		
            IF(FILEEXISTANCE) THEN
                  PRINT *, "FILEEXISTANCE"
		      !CALL SLEEPQQ(10000)
		      MATLABWORKING=0
            ELSE
                  PRINT *, "Not FILEEXISTANCE"
			MATLABWORKING=0
            END IF
	END DO
C  End of Control ABAQUS run (pause&go)  ***********************************V   
 
      return
      end
