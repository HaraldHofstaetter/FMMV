#if (FMM_DIM==2)
program test_fmmv2d
    use fmmv2d_module
#else
program test_fmmv3d
    use fmmv3d_module
#endif
    use, intrinsic :: iso_c_binding
#ifdef USE_SINGLE_PRECISION
    use test_utilities_f
#else		
    use test_utilities
#endif	
    implicit none
#ifdef USE_SINGLE_PRECISION
    integer, parameter :: prec = 4
#else		
    integer, parameter :: prec = 8
#endif	

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
    type( FmmvOptions ) :: options
    type( FmmvStatistics ) :: statistics
    logical :: with_dipoleMoments
    logical :: with_gradients
    logical :: with_targets
    logical :: printResult
    logical :: splitCall
    type( c_ptr ) :: errorMessage
    type( c_ptr ) :: fmmvHandle
    integer( kind = c_int ) :: typeSources
    integer( kind = c_int ) :: typeTargets
    real(prec) :: errL2, errL2grad
    real(8) :: exTime
    integer :: exAcc
    integer :: i

    NParticles = 10000
    NTargets = 0
    NDirect = 100
    with_dipoleMoments = .false.
    with_gradients = .false.
    with_targets = .false.
    printResult = .false.
    splitCall = .false.
    call set_from_command_line("dipoleMoments", with_dipoleMoments)
    call set_from_command_line("gradients", with_gradients)
    call set_from_command_line("printResult", printResult)
    call set_from_command_line("splitCall", splitCall)

    options = fmmvGetDefaultOptions()
    options%splitThreshold=190
    options%levels=51
    options%directEvalThreshold=-1
    options%periodicBoundaryConditions = 0
    options%extrinsicCorrection = 0
    options%pM = 6
    options%pL = 6
    options%s = 8
    options%ws = 1
    options%reducedScheme = 0
    options%useHilbertOrder = 0
#ifdef USE_SINGLE_PRECISION
    options%directEvalAccuracy = 1
    exAcc = 1
#else
    options%directEvalAccuracy = 2
    exAcc = 2
#endif
    options%scale = 1.0
    options%useFarfieldNearfieldThreads = 0
    options%beta = 0.0


    call set_from_command_line("NParticles", NParticles)
    call set_from_command_line("NTargets", NTargets)
    call set_from_command_line("NDirect", NDirect)

    call set_from_command_line("pM", options%pM)
    call set_from_command_line("pL", options%pL)
    call set_from_command_line("s", options%s)
    
    call set_from_command_line("ws", options%ws)
    call set_from_command_line("reducedScheme", options%reducedScheme)
   
    call set_from_command_line("scale", options%scale)
    call set_from_command_line("beta", options%beta)
    call set_from_command_line("splitThreshold", options%splitThreshold)
    call set_from_command_line("splitTargetThreshold", options%splitTargetThreshold)
    call set_from_command_line("levels", options%levels)
    call set_from_command_line("directEvalThreshold", options%directEvalThreshold)
    call set_from_command_line("periodicBoundaryConditions", options%periodicBoundaryConditions)
    call set_from_command_line("extrinsicCorrection", options%extrinsicCorrection)
    call set_from_command_line("useHilbertOrder", options%useHilbertOrder)
    call set_from_command_line("directEvalAccuracy", options%directEvalAccuracy)
    call set_from_command_line("useFarfieldNearfieldThreads", options%useFarfieldNearfieldThreads)


    allocate( particles(FMM_DIM, NParticles) )
    call random_number(particles)
    allocate( charges(NParticles) )
    charges = 1.0d0

    with_targets = (NTargets>0)
    if (with_targets) then
        allocate( targets(FMM_DIM, NTargets) )
        call random_number(targets)
    else
        targets => null()
        NTargets = NParticles
    end if

    if(with_dipoleMoments) then
        allocate( dipoleMoments(FMM_DIM, NParticles) )
        call random_number(dipoleMoments)
    else
        dipoleMoments => null()
    end if

    allocate( pot(NTargets) )
    if(with_gradients) then
        allocate( grad(FMM_DIM, NTargets) )
    else
       grad => null()
    end if


    write(*,*)
    if (splitCall) then
      typeSources = FMMV_CHARGES
      if(with_dipoleMoments) then
         typeSources = ior(typeSources, FMMV_DIPOLEMOMENTS)
      end if
      typeTargets = FMMV_POTENTIALS
      if(with_gradients) then
         typeTargets = ior(typeTargets, FMMV_GRADIENTS)
      end if
      write(*,*) "FMM INITIALIZE STARTED"
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
      call fmmv2df_initialize( &
#else
      call fmmv2d_initialize( &
#endif
#else
#ifdef USE_SINGLE_PRECISION
      call fmmv3df_initialize( &
#else
      call fmmv3d_initialize( &
#endif
#endif
           fmmvHandle,          &
           NParticles,          & 
           particles,           &
           typeSources,         &
           NTargets,            &
           targets,             &
           typeTargets,         &
           options,             &
           statistics,          &
           errorMessage         &
           )    
      write(*,*) "FMM EVALUATE STARTED"
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
      call fmmv2df_evaluate( &
#else
      call fmmv2d_evaluate( &
#endif
#else
#ifdef USE_SINGLE_PRECISION
      call fmmv3df_evaluate( &
#else
      call fmmv3d_evaluate( &
#endif
#endif
           fmmvHandle,          & 
           charges,             &
           dipoleMoments,       &
           pot,                 &
           grad,                &
           statistics,          &
           errorMessage         &
           ) 
      write(*,*) "FMM FINALIZE STARTED"
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
      call fmmv2df_finalize( &
#else
      call fmmv2d_finalize( &
#endif
#else
#ifdef USE_SINGLE_PRECISION
      call fmmv3df_finalize( &
#else
      call fmmv3d_finalize( &
#endif
#endif
           fmmvHandle,          & 
           statistics,          &
           errorMessage         &
           )
      write(*,*) "FMM FINISHED"
    else
      write(*,*) "FMM STARTED"
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
      call fmmv2df(      &
#else
      call fmmv2d(       &
#endif
#else
#ifdef USE_SINGLE_PRECISION
      call fmmv3df(      &
#else
      call fmmv3d(       &
#endif
#endif
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       & 
           NTargets,            & 
           targets,             & 
           pot,                 & 
           grad,                & 
           options,             &
           statistics,          &
           errorMessage         &
           )
      write(*,*) "FMM FINISHED"
    end if       
    write(*,*)

    call printFmmvStatistics(statistics)
    write(*,*)

    if (NTargets < NDirect) then 
        NDirect = NTargets;
    end if 
    allocate( exPot(NDirect) )
    if(with_gradients) then
        allocate( exGrad(3, NDirect) )
    else
       exGrad => null()
    end if
    
    write(*,*) "DIRECT METHOD STARTED"
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
    call fmmv2df_direct(&
#else
    call fmmv2d_direct( &
#endif
#else
#ifdef USE_SINGLE_PRECISION
       call fmmv3df_direct(&
#else
       call fmmv3d_direct( &
#endif
#endif
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       & 
           NDirect ,            & 
           targets,             & 
           exPot,               & 
           exGrad,              & 
           exAcc,               &
#if (FMM_DIM==3)
           options%beta,        &
#endif
           exTime,              &
           errorMessage         &
           )
    write(*,*) "DIRECT METHOD FINISHED"

    if (with_gradients) then
        call relL2Err(0)  ! initialize
#if (FMM_DIM==2)
        call relL2Err2(0) ! initialize
#else
        call relL2Err3(0) ! initialize
#endif

        if (printResult) then 
             write(*,*)
             write(*,*) "\n        Pot(FMM)   Pot(exact)     rel.err.  rel.err.(grad)"
             write(*,*) "============================================================"
        end if
      
        do i=1,NDirect 
             call relL2Err(1, pot(i), exPot(i))    ! accumulate 
#if (FMM_DIM==2)
             call relL2Err2(1, grad(:,i), exGrad(:,i)) ! accumulate 
#else
             call relL2Err3(1, grad(:,i), exGrad(:,i)) ! accumulate 
#endif
             if (printResult) then
                   write(*,'(I3,1X,E12.4,1X,E12.4,1X,E12.4,1X,E12.4)') &
                              i, pot(i), exPot(i), relErr(pot(i), exPot(i)),  &
#if (FMM_DIM==2)
                              relErr2(grad(:,i), &
#else
                              relErr3(grad(:,i), &
#endif
                              exGrad(:,i))
             end if
        end do
        call relL2Err(2, res=errL2)      ! finalize
#if (FMM_DIM==2)
        call relL2Err2(2, res=errL2grad) ! finalize
#else
        call relL2Err3(2, res=errL2grad) ! finalize
#endif
    else 
        call relL2Err(0); ! initialize 
        if (printResult) then 
             write(*,*)
             write(*,*) "        Pot(FMM)   Pot(exact)     rel.err."
             write(*,*) "============================================"
        end if
        
        do i=1,NDirect
             call relL2Err(1, pot(i), exPot(i)) ! accumulate 
             if (printResult) then 
                   write(*,'(I3,1X,E24.16,1X,E24.16,1X,E12.4)') &
                              i, pot(i), exPot(i), relErr(pot(i), exPot(i))
             end if
        end do
        call relL2Err(2, res=errL2) ! finalize
    end if

    write(*,*)
    if (with_gradients) then
       100 format ('err_L2(pot) = ',E12.4,'  err_L2(grad) = ',E12.4) 
       write (*,100) errL2, errL2grad
    else 
       101 format ('err_L2(pot) = ',E12.4) 
       write (*,101) errL2
    end if
    

    deallocate( particles, charges, pot, exPot )
    if (with_targets) deallocate(targets)
    if (with_dipoleMoments) deallocate(dipoleMoments)
    if (with_gradients) deallocate(grad, exGrad)

#if (FMM_DIM==2)
end program test_fmmv2d
#else
end program test_fmmv3d
#endif
