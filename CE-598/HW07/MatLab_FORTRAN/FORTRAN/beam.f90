      program beam
!     peridynamic plane stress beam example
      implicit none
      integer maxNumPtcls, numTimeSteps, numberOfMovieFrames
      parameter (maxNumPtcls = 300) 
      parameter (numTimeSteps = 2000)
      parameter (numberOfMovieFrames = 50)
      real :: L = 0.01 ! lattice spacing
      real :: density = 2323.0 ! material density of concrete, kg/m^3
      real :: E = 24.86E9 ! Young's modulus of concrete
      real :: loadFactor = -5.E4 ! (times beam self-weight)
      real :: dampRatio = 0.5 ! external damping coefficient
      real :: beamSpan = 0.30 ! span of beam between supports, m
      real :: beamHeight = 0.06  ! width of beam in Y-direction, m
      real :: beamThickness = 0.03 ! thickness of beam in Z-direction, m
      real :: LatRotAngle = 0.   !lattice rotation angle about Z-axis, in degrees
      integer :: writeInterval = numTimeSteps/numberOfMovieFrames ! interval for writing restart file
      integer :: PtclToMonitor = 1 ! particle to monitor for time history output
      real :: Xref(2, maxNumPtcls) = 0. ! reference positions of particles
      real  dT ! time step
      real beamLength ! length of beam in X-direction, m

      real mPtcl, volPtcl  ! mass and volume of each particle
      real c ! stiffness (force per unit stretch) of bond between particles.
      real xCur(2, maxNumPtcls)   ! current posiitons of particles
      integer :: bondList(6, maxNumPtcls) = 0 !list of neighboring particles for each particle
      integer numPtcls, iStep, iBond, iDir, iPtcl
      integer :: BCcodes(2, maxNumPtcls) = 0 ! boundary condition codes: 
                                    !  0 = free to move
                                    !  1 = applied displacement
      real :: BCvals(2, maxNumPtcls) = 0.    ! applied force or displacement values.
      real :: velocity(2, maxNumPtcls) = 0.  ! initial velocities
      real internForce(2, maxNumPtcls)   ! internal forces acting upon particles.
      real dx(2), curBondLength, dc(2), stretch, externalForce
      real domLimits(2,2) !minimum and maximum domain limits in X and Y 
      real dampForce, forcePtcl, accel, xI(2), xF(2)
      real histPtclMonitored(4, numTimeSteps)
      real curTime, massBeam, Inertia, dampCoeff
      real beamOverhang, g, loadPerPtcl
      real deltaAnalytical, fundAngularFreqBeam, fundPeriod
      real latticeRadius    !spherical radius centered at orginal containing the lattice domain
      real latticeOrigin(2) !origin in physical space of lattice coordinate system
      real latBmatrix(2, 2)  !lattice basis matrix
      integer i, j, nX, nY, leftSupport, rightSupport
      real, parameter :: pi   = 3.141592653589793
      real :: tol = 0.0001   ! a small tolerance
      
      beamLength = beamSpan + 4.*L
      volPtcl = (sqrt(3.)/2.)*beamThickness*(L**2) !Eq. 7.5 of Practical Peridynamics
      mPtcl = density*volPtcl
      c = sqrt((3.)/2.)*E*L*beamThickness  !see pg. 210 of Practical Peridynamics
      dT = sqrt(mPtcl/c)/20.
      latticeOrigin(1) = beamSpan/2.
      latticeOrigin(2) = beamHeight/2.
      LatRotAngle = LatRotAngle*pi/180.  !rotation about X-axis, in radians
      call setupLattice(L, LatRotAngle, latBmatrix)
      latticeRadius = 2.0*sqrt((beamSpan/2.)**2 + (beamHeight/2.)**2)
      nX = floor(latticeRadius/L)
      nY = floor(latticeRadius/L)
      domLimits(1,1) = +1E10 ! minimum X
      domLimits(2,1) = -1E10  ! maximum X
      domLimits(1,2) = +1E10 ! minimum Y
      domLimits(2,2) = -1E10  ! maximum Y
      !
      ! Insert particles
      !
      beamOverhang = (beamLength - beamSpan)/2
      numPtcls = 0
      do i = -nX, nX
        do j = -nY, nY
          xI(1) = i
          xI(2) = j
          xF = matmul(latBmatrix, xI) + latticeOrigin
          if(((xF(1) .gt. -beamOverhang) .and. &
                      (xF(1) .lt. beamSpan + beamOverhang) ) .and.  &
                   ((xF(2) .gt. tol) .and.  &
                       (xF(2) .lt. beamHeight - tol) ) ) then
            numPtcls = numPtcls + 1
            if(numPtcls .gt. maxNumPtcls) then
              write(*,*)'number of particles exceeds maximum allowed'
              return
            endif
            Xref(:, numPtcls) = xF(:)
            domLimits(1,1) = min(Xref(1, numPtcls), domLimits(1,1))
            domLimits(2,1) = max(Xref(1, numPtcls), domLimits(2,1))
            domLimits(1,2) = min(Xref(2, numPtcls), domLimits(1,2))
            domLimits(2,2) = max(Xref(2, numPtcls), domLimits(2,2))
            if((abs(xF(1)) .lt. 0.5*L) .and. (xF(2) .lt. L)) then
              BCcodes(1:2, numPtcls) = 1 !pin the lower left-most particle
              leftSupport = numPtcls
            elseif((abs(xF(1) - beamSpan) .lt. L/2.) &
                                       .and. (xF(2) .lt. L)) then
              BCcodes(2, numPtcls) = 1 !roller at the lower right-most particle
              rightSupport = numPtcls
            else
              BCvals(2, numPtcls) = 1.
            endif 
            if((i .eq. 0) .and. (j .eq. 0)) then
              PtclToMonitor = numPtcls
            endif
          endif
        enddo
      enddo
      g = 9.81 ! m/s^2
      loadPerPtcl = loadFactor*mPtcl*g
      BCvals = BCvals*loadPerPtcl
      xCur = Xref      !initial positions of the particles
      velocity = 0
      massBeam = mPtcl*numPtcls
      Inertia = beamThickness*(beamHeight**3)/12
      deltaAnalytical = loadFactor*5.*massBeam*g*(beamSpan**3)/  &
                             (384.*E*Inertia)
      fundAngularFreqBeam = sqrt((48*E*Inertia/(beamSpan**3))/(massBeam/2))
      fundPeriod = 2*pi/fundAngularFreqBeam
      dampCoeff = -2*mPtcl*dampRatio*fundAngularFreqBeam
      call setupBondList(Xref, latBmatrix, numPtcls, L, bondList) !determine the bond list for each particle
      call writeData(iStep, writeInterval, curTime, numPtcls, &
                           numTimeSteps, bcCodes, xCur, Xref)
      do iStep = 1, numTimeSteps
        curTime = iStep*dT
        ! Find internal forces acting upon each particle:
        internForce = 0.
        do iPtcl = 1, numPtcls
          do iBond = 1, 6
            if(bondList(iBond, iPtcl) .ne. 0) then
              dx = xCur(:, bondList(iBond, iPtcl)) - xCur(:, iPtcl)
              curBondLength = sqrt(dx(1)**2 + dx(2)**2)
              dc = dx/curBondLength
              stretch = (curBondLength - L)/L
              internForce(:, iPtcl) = internForce(:, iPtcl) + dc*(c*stretch)
            endif
          enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        ! Update the particle positions x:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do iPtcl = 1, numPtcls
          do iDir = 1, 2
            if(BCcodes(iDir,iPtcl) .eq. 0.) then
                externalForce = BCvals(iDir,iPtcl)
                dampForce = dampCoeff*velocity(iDir,iPtcl)
                forcePtcl = internForce(iDir,iPtcl) + &
                    externalForce + dampForce
                accel = forcePtcl/mPtcl
                velocity(iDir,iPtcl) = velocity(iDir,iPtcl) + dT*accel
                xCur(iDir,iPtcl)= xCur(iDir,iPtcl)+ &
                    dT*velocity(iDir,iPtcl)
            else
                xCur(iDir,iPtcl)= Xref(iDir,iPtcl) + BCvals(iDir,iPtcl)
                velocity(iDir,iPtcl) = 0.
            endif                           
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
          enddo
        enddo
        call writeData(iStep, writeInterval, curTime, numPtcls, &
                           numTimeSteps, bcCodes, xCur, Xref)
        histPtclMonitored(1:2, iStep) = xCur(:, PtclToMonitor)
        histPtclMonitored(3:4, iStep) = Xref(:, PtclToMonitor)
      enddo
      
      open(14, file = 'timehist.pdq', status = 'replace')
      write(14, *)dT, PtclToMonitor
      do iStep = 1, numTimeSteps
        write(14,'(4ES23.14E3)') histPtclMonitored(1:4, iStep)
      enddo
      close(14)      
      open(13, file = 'postProcParams.pdq', status = 'replace')
      domLimits(1,1) = domLimits(1,1)
      domLimits(2,1) = domLimits(2,1)
      domLimits(1,2) = domLimits(1,2)
      domLimits(2,2) = domLimits(2,2)
      write(13,*)numTimeSteps, writeInterval, &
                 domLimits(1,1), domLimits(2,1),domLimits(1,2), &
                 domLimits(2,2), PtclToMonitor
      close(13)
      end program
      
!----------------------------------------------------------------------------------------
      subroutine setupLattice(L, LatRotAngle, latB)
!----------------------------------------------------------------------------------------
!
! This subroutine sets up the rotated lattice matrix latB for a hexagonal lattice
!
      implicit none
      real L            ! lattice spacing, input
      real LatRotAngle  !lattice rotation about the Z axis, in radians, input
      real latB(2,2)    !lattice matrix, returned
      
      real latRotZ(2,2)   !lattice basis matrix
      
      latB(1,1) = L*1.
      latB(1,2) = L/2.
      latB(2,1) = 0.
      latB(2,2) = L*sqrt(3.)/2.
      latRotZ(1,1) = cos(LatRotAngle)
      latRotZ(1,2) = -sin(LatRotAngle)
      latRotZ(2,1) = sin(LatRotAngle)
      latRotZ(2,2) = cos(LatRotAngle)
      latB = matmul(latRotZ, LatB)
      return
      end

!----------------------------------------------------------------------------------------
      subroutine setupBondList(Xref, latBmatrix, numPtcls, L, bondList)
!----------------------------------------------------------------------------------------
!
! This subroutine sets up the lattice bond list for a hexagonal lattice body
!
      implicit none
      real Xref(2, numPtcls) ! current positions of particles
      real latBmatrix(2, 2)  !lattice basis matrix
      integer numPtcls       ! number of particles
      real L ! lattice spacing
      integer bondList(6, numPtcls) !list of neighboring particles for each particle (returned)

      integer iPtcl, jPtcl, iBond
      real Xi(2), Xj(2), dX(2), bondLength
      real bond(2, 6), dBonds(2), dist
      real tol
      
      tol = 0.1*L
      Xi(1) = 1.
      Xi(2) = 0.
      bond(:, 1) =  matmul(latBmatrix, Xi)
      Xi(1) = -1.
      Xi(2) = 0.
      bond(:, 2) =  matmul(latBmatrix, Xi)
      Xi(1) = 0.
      Xi(2) = 1.
      bond(:, 3) =  matmul(latBmatrix, Xi)
      Xi(1) = 0.
      Xi(2) = -1.
      bond(:, 4) =  matmul(latBmatrix, Xi)
      Xi(1) = -1.
      Xi(2) = 1.
      bond(:, 5) =  matmul(latBmatrix, Xi)
      Xi(1) = 1.
      Xi(2) = -1.
      bond(:, 6) =  matmul(latBmatrix, Xi)
      bondList = 0
      do iPtcl = 1, numPtcls
        Xi = Xref(:, iPtcl)
        do jPtcl = 1, numPtcls
          Xj = Xref(:, jPtcl)
          dX = Xj - Xi
          bondLength = sqrt(dX(1)**2 + dX(2)**2)
          if(abs(bondLength - L) .lt. tol*L) then
            do iBond = 1, 6
              dBonds = dX - bond(:, iBond)
              dist = dBonds(1)**2 + dBonds(2)**2
              if(dist .lt. tol*L) then
                bondList(iBond, iPtcl) = jPtcl
              endif
            enddo
          endif
        enddo
      enddo
      return
      end
      
 !----------------------------------------------------------------------------------------
      subroutine writeData(iStep, writeInterval, curTime, numPtcls, &
                           numTimeSteps, bcCodes, xCur, Xref)
!----------------------------------------------------------------------------------------
!
! This subroutine writes out restart files and the time history file
!
      implicit none
      integer iStep
      integer writeInterval
      real curTime
      integer numPtcls
      integer numTimeSteps
      integer bcCodes(2, numPtcls)
      real xCur(2, numPtcls)
      real Xref(2, numPtcls)
      
      integer iPtcl
      character*60 Resfile

      if((mod(iStep, writeInterval) .eq. 0) .or. (iStep .eq. 1)) then
        write(*,*) 'at time step: ', istep, ' of ', numTimeSteps
        if (istep.lt.10) then
          write(Resfile,'(A,I1, A)') 'restart.',istep,'.pdq'
        elseif (istep.lt.100) then
          write(Resfile,'(A,I2, A)') 'restart.',iStep,'.pdq'
        elseif (istep.lt.1000) then
          write(Resfile,'(A,I3, A)') 'restart.',iStep,'.pdq'
        elseif (istep.lt.10000) then
          write(Resfile,'(A,I4, A)') 'restart.',iStep,'.pdq'
        endif
        open(18, file = Resfile, status = 'replace')
        write(18, *)curTime
        do iPtcl = 1, numPtcls
          write(18,'(1I6, 2I3, 4ES23.14E3)') &
                iPtcl, bcCodes(:, iPtcl), xCur(:, iPtcl), Xref(:, iPtcl)
        enddo
        close(18)
      endif
      return
      end
