      program prob3
!     peridynamic cable example in 2D space

      implicit none
      integer maxNumPtcls, maxNumSteps
      parameter (maxNumPtcls = 20) 
      parameter (maxNumSteps = 2000) 
      real :: L = 1. ! lattice spacing
      real :: M = 1. ! mass of each particle
      real :: K = 100. ! stiffness (force per unit stretch) of bond between particles.
      real :: dampCoeff = 0.05 ! external damping coefficient
      real  dT ! time step
      integer :: numTimeSteps = 2000 ! number of time steps
      integer :: numPtcls = 12 ! number of lattice particles
      real :: Xref(maxNumPtcls, 2) = 0. ! reference positions of particles
      real xCur(maxNumPtcls, 2)   ! current posiitons of particles
      integer iPtcl, iStep, iBond, iDir
      integer :: BCcodes(maxNumPtcls, 2) = 0 ! boundary condition codes: 
                                    !  0 = free to move
                                    !  1 = applied displacement
      real :: BCvals(maxNumPtcls, 2) = 0.    ! applied force or displacement values.
      real :: velocity(maxNumPtcls, 2) = 0.  ! initial velocities
      real internForce(maxNumPtcls, 2)   ! internal forces acting upon particles.
      real dx(2), curBondLength, dc(2), stretch, externalForce
      real dampForce, forcePtcl, accel
      real history(maxNumSteps, 4)
      real curTime
      integer :: writeInterval = maxNumSteps/100 ! interval for writing restart file
      character*60 Resfile

      dT = 8*sqrt(M/K)/6.
      do iPtcl = 1, numPtcls
        Xref(iPtcl, 1) = L*(iPtcl-1)! !initial horizontal positions of particles
      enddo
      BCcodes(1, 1:2) = 1 !fix the left-most particle
      
      BCcodes(numPtcls, 1) = 1 !fix the left-most particle
      BCcodes(numPtcls, 2) = 0 !fix the left-most particle
      
      BCvals(2:numPtcls, 2) = -1.! !apply a vertical force free particles.
      xCur = Xref      !initial positions of the particles
      do iStep = 1, numTimeSteps
        curTime = iStep*dT
        ! Find internal forces acting upon each particle:
        internForce = 0. 
        do iBond = 1, numPtcls - 1
          dx = xCur(iBond + 1, :) - xCur(iBond, :)!
          curBondLength = sqrt(dx(1)**2 + dx(2)**2)!
          dc = dx/curBondLength!
          stretch = (curBondLength - L)/L!
          internForce(iBond, :) = internForce(iBond, :) + dc*(K*stretch)!
          internForce(iBond+1, :) = internForce(iBond+1, :) - dc*(K*stretch)!
        enddo
        ! Update the particle positions x:
        do iPtcl = 1, numPtcls
          do iDir = 1, 2
            if(BCcodes(iPtcl, iDir) .eq. 0) then
              externalForce = BCvals(iPtcl, iDir)!
              dampForce = -dampCoeff*velocity(iPtcl, iDir)!
              forcePtcl = internForce(iPtcl, iDir) + externalForce &
                                                          + dampForce
              accel = (forcePtcl)/M
              velocity(iPtcl, iDir) = velocity(iPtcl, iDir) + dT*accel!
              xCur(iPtcl, iDir)= xCur(iPtcl, iDir) &
                                   + dT*velocity(iPtcl, iDir)!
            else
              xCur(iPtcl, iDir)= xRef(iPtcl, iDir) + BCvals(iPtcl, iDir)!
              velocity(iPtcl, iDir) = 0.
            endif
          enddo
        enddo
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
            write(18,'(1I3, 2I3, 4ES23.14E3)') &
                 iPtcl, bcCodes(iPtcl,:), xCur(iPtcl,:), Xref(iPtcl, :)
          enddo
          close(18)
        endif
        history(iStep, 1:2) = xCur(numPtcls, :)
        history(iStep, 3:4) = Xref(numPtcls, :)
      enddo
      
      open(14, file = 'timehist.pdq', status = 'replace')
      write(14, *)dT, numPtcls
      do iStep = 1, numTimeSteps
        write(14,'(4ES23.14E3)') history(iStep, 1:4)
      enddo
      close(14)      
      
      open(13, file = 'postProcParams.pdq', status = 'replace')
      write(13,*)numTimeSteps, writeInterval
      close(13)
      end program prob3
