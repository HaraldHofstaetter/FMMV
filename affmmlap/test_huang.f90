program test_huang
    use, intrinsic :: iso_c_binding
    use fmmv3d_module
    use test_utilities
    implicit none
    integer, parameter :: prec = 8

    real(prec), pointer :: particles(:,:)
    real(prec), pointer :: charges(:)
    real(prec), pointer :: pot(:)
    real(prec), pointer :: targets(:,:)
    real(prec), pointer :: dipoleMoments(:,:)
    real(prec), pointer :: grad(:,:)
    real(prec), pointer :: exPot(:)
    real(prec), pointer :: exGrad(:,:)
    integer :: NParticles
    integer :: NTargets
    integer :: NDirect
    logical :: printResult
    logical :: directEvalExtendedPrecision
    type( c_ptr ) :: errorMessage
    real(prec) :: errL2, errL2grad
    real(8) :: exTime, time0, time1
    integer :: exAcc
    integer :: ier
    integer :: i
    real(8) ::second

    printResult = .false.
    directEvalExtendedPrecision = .false.
    NParticles = 10000
    NDirect = 100
    call set_from_command_line("NParticles", NParticles)
    call set_from_command_line("NDirect", NDirect)
    call set_from_command_line("printResult", printResult)
    call set_from_command_line("directEvalExtendedPrecision", directEvalExtendedPrecision)

    allocate( particles(3, NParticles) )
    call random_number(particles)
    particles = particles -0.5d0 ! Huang expects particles in [-.5,.5]^3 ...
    allocate( charges(NParticles) )
    charges = 1.0d0

    NTargets = NParticles
    allocate( pot(NTargets) )
    allocate( grad(3, NTargets) )

    write(*,*)
    write(*,*) "FMM STARTED"
!
!-----3. main fmm call, call fmm to calculate the potential
!
      TIME0 = SECOND()
      CALL FMMLAP_A(NParticles,particles,charges,pot,grad,IER)
      TIME1=SECOND()


      write(*,*) "FMM FINISHED"
      print *, "ier = ", ier
      WRITE(*,554)TIME1-TIME0
554   FORMAT(' time = ',F8.2)


   write(*,*)

    if (NTargets < NDirect) then 
        NDirect = NTargets;
    end if 
    allocate( exPot(NDirect) )
    allocate( exGrad(3, NDirect) )
    
    write(*,*) "DIRECT METHOD STARTED"
    if (directEvalExtendedPrecision) then
       call directMethodCoulomb3d_xp( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       & 
           NDirect ,            & 
           targets,             & 
           exPot,               & 
           exGrad,              & 
           exTime,              &
           errorMessage         &
           )
    else
       call directMethodCoulomb3d( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       & 
           NDirect ,            & 
           targets,             & 
           exPot,               & 
           exGrad,              & 
           exAcc,               &
           exTime,              &
           errorMessage         &
           )
    end if   
    write(*,*) "DIRECT METHOD FINISHED"

        call relL2Err(0)  ! initialize
        call relL2Err3(0) ! initialize

        if (printResult) then 
             write(*,*)
             write(*,*) "\n        Pot(FMM)   Pot(exact)     rel.err.  rel.err.(grad)"
             write(*,*) "============================================================"
        end if
      
        do i=1,NDirect 
             call relL2Err(1, pot(i), exPot(i))    ! accumulate 
             call relL2Err3(1, grad(:,i), -exGrad(:,i)) ! accumulate 
             if (printResult) then
                   write(*,'(I3,1X,E12.4,1X,E12.4,1X,E12.4,1X,E12.4)') &
                              i, pot(i), exPot(i), relErr(pot(i), exPot(i)),  &
                              relErr3(grad(:,i), -exGrad(:,i))
             end if
        end do
        call relL2Err(2, res=errL2)      ! finalize
        call relL2Err3(2, res=errL2grad) ! finalize

    write(*,*)
       100 format ('err_L2(pot) = ',E12.4,'  err_L2(grad) = ',E12.4) 
       write (*,100) errL2, errL2grad
    

    deallocate( particles, charges, pot, exPot )
    deallocate(grad, exGrad)

end program test_huang

