module fmmv3d_types
   use, intrinsic :: iso_c_binding
   implicit none

   integer, parameter :: FMMV_STAT_MAX = 20
   integer, parameter :: FMMV_MAX_NUM_PAPI_EVENTS = 10
  
   ! Fortran 2003 feature enum not implemented ?
!   enum, bind(x)  
!     enumerator :: &
!        FMMV_STAT_TOTAL=1,      &
!        FMMV_STAT_BUILD_TREE,   &
!        FMMV_STAT_GEN_M,        &
!        FMMV_STAT_M2M,  &
!        FMMV_STAT_M2L,  &
!        FMMV_STAT_L2L,  &
!        FMMV_STAT_EVAL_L,       &
!        FMMV_STAT_LIST1,        &
!        FMMV_STAT_LIST3,        &
!        FMMV_STAT_LIST4,        &
!        FMMV_STAT_LIST34,       &
!        FMMV_STAT_FARFIELD,     &
!        FMMV_STAT_NEARFIELD,    &
!        FMMV_STAT_INITIALIZE,   &
!        FMMV_STAT_EVALUATE,     &
!        FMMV_STAT_FINALIZE,     &
!        FFMV_STAT_LAST_ ! mark last entry, do not remove!
!   end enum
 
   type, bind( c ) :: FmmvStatistics 
        real(kind=c_double) :: time(FMMV_STAT_MAX)
        real(kind=c_double) :: etime(FMMV_STAT_MAX)

        integer(kind=c_int) :: PAPIeventSet
        integer(kind=c_long_long) :: PAPIvalues(FMMV_MAX_NUM_PAPI_EVENTS,FMMV_STAT_MAX)

        integer(kind=c_int) :: pM
        integer(kind=c_int) :: pL
        integer(kind=c_int) :: s_eps
        integer(kind=c_int) :: s_exp
        
        integer(kind=c_int) :: noOfParticles
        integer(kind=c_int) :: noOfTargets
        integer(kind=c_int) :: noOfSourceLevels
        integer(kind=c_int) :: noOfTargetLevels
        
        integer(kind=c_int) :: maxNoOfStoredXin
        integer(kind=c_int) :: maxAllocatedMemory
        integer(kind=c_long_long) :: noOfDirectInteractions
        
        integer(kind=c_int) :: noOfSourceBoxes
        integer(kind=c_int) :: noOfTargetBoxes
        integer(kind=c_int) :: noOfSourceLeafBoxes
        integer(kind=c_int) :: noOfTargetLeafBoxes
        real(kind=c_float) :: averageNoOfParticlesPerLeafBox
        real(kind=c_float) :: averageNoOfTargetsPerLeafBox
        
        integer(kind=c_int) :: noOfParticlesInLevel(52)
        integer(kind=c_int) :: noOfTargetsInLevel(52)
        integer(kind=c_int) :: noOfSourceBoxesInLevel(52)
        integer(kind=c_int) :: noOfTargetBoxesInLevel(52)
        integer(kind=c_int) :: noOfSourceLeafBoxesInLevel(52)
        integer(kind=c_int) :: noOfTargetLeafBoxesInLevel(52)
        real(kind=c_float) :: averageNoOfParticlesPerLeafBoxInLevel(52)
        real(kind=c_float) :: averageNoOfTargetsPerLeafBoxInLevel(52)
   end type FmmvStatistics 

   type, bind( c ) :: FmmvOptions
       real(kind=c_double) :: beta

       integer(kind=c_int) :: pM
       integer(kind=c_int) :: pL
       integer(kind=c_int) :: s
       
       integer(kind=c_int) :: ws
       integer(kind=c_int) :: reducedScheme
       
       real(kind=c_double) :: scale
       integer(kind=c_int) :: splitThreshold
       integer(kind=c_int) :: splitTargetThreshold
       integer(kind=c_int) :: levels
       integer(kind=c_int) :: directEvalThreshold
       integer(kind=c_int) :: periodicBoundaryConditions
       integer(kind=c_int) :: extrinsicCorrection
       integer(kind=c_int) :: useHilbertOrder
       integer(kind=c_int) :: directEvalAccuracy
       integer(kind=c_int) :: useFarfieldNearfieldThreads;
       
       integer(kind=c_int) :: PAPIeventSet;


       type(c_ptr) :: x
       type(c_ptr) :: w
       type(c_ptr) :: M
   end type FmmvOptions
end module fmmv3d_types

module fmmv3d_module
   use, intrinsic :: iso_c_binding
   use fmmv3d_types
   implicit none

  integer, parameter :: FMMV_DIPOLEMOMENTS = 1
  integer, parameter :: FMMV_CHARGES = 2
  integer, parameter :: FMMV_POTENTIALS = 1
  integer, parameter :: FMMV_GRADIENTS = 2

   interface
       function fmmvGetDefaultOptions() bind( c, name="fmmvGetDefaultOptions" )
           use, intrinsic :: iso_c_binding
           use fmmv3d_types 
           implicit none
           type( FmmvOptions ) :: fmmvGetDefaultOptions
       end function fmmvGetDefaultOptions

       subroutine printFmmvStatistics(statistics) bind( c, name="printFmmvStatistics" )
           use, intrinsic :: iso_c_binding
           use fmmv3d_types 
           implicit none
           type( FmmvStatistics ) :: statistics
       end subroutine printFmmvStatistics

       subroutine fmmv3d( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       &           
           NTargets,            &
           targets,             &
           potentials,          &
           gradients,           &
           options,             &
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3d" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          integer( kind = c_int ), value :: NParticles
          real( kind = c_double ) :: particles(3,*) 
          real( kind = c_double ) :: charges(*) 
          real( kind = c_double ) :: dipoleMoments(3,*) 
          integer( kind = c_int ), value :: NTargets
          real( kind = c_double ) :: targets(3,*) 
          real( kind = c_double ) :: potentials(*) 
          real( kind = c_double ) :: gradients(3,*)
          type( FmmvOptions ) :: options
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3d

       subroutine fmmv3df( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       &           
           NTargets,            &
           targets,             &
           potentials,          &
           gradients,           &
           options,             &
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3df" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types
          implicit none
          integer( kind = c_int ), value :: NParticles
          real( kind = c_float ) :: particles(3,*) 
          real( kind = c_float ) :: charges(*) 
          real( kind = c_float ) :: dipoleMoments(3,*) 
          integer( kind = c_int ), value :: NTargets
          real( kind = c_float ) :: targets(3,*) 
          real( kind = c_float ) :: potentials(*) 
          real( kind = c_float ) :: gradients(3,*)
          type( FmmvOptions ) :: options
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3df

       subroutine fmmv3d_initialize( &
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
           ) bind( c, name="fmmv3d_initialize" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr) :: fmmvHandle
          integer( kind = c_int ), value :: NParticles
          real( kind = c_double ) :: particles(3,*) 
          integer( kind = c_int ), value :: typeSources
          integer( kind = c_int ), value :: NTargets
          real( kind = c_double ) :: targets(3,*) 
          integer( kind = c_int ), value :: typeTargets
          type( FmmvOptions ) :: options
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3d_initialize
       
       subroutine fmmv3df_initialize( &
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
           ) bind( c, name="fmmv3df_initialize" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr) :: fmmvHandle
          integer( kind = c_int ), value :: NParticles
          real( kind = c_float ) :: particles(3,*) 
          integer( kind = c_int ), value :: typeSources
          integer( kind = c_int ), value :: NTargets
          real( kind = c_float ) :: targets(3,*) 
          integer( kind = c_int ), value :: typeTargets
          type( FmmvOptions ) :: options
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3df_initialize
       
       subroutine fmmv3d_evaluate( &
           fmmvHandle,          & 
           charges,             &
           dipoleMoments,       &
           potentials,          &
           gradients,           &
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3d_evaluate" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr), value :: fmmvHandle
          real( kind = c_double ) :: charges(*) 
          real( kind = c_double ) :: dipoleMoments(3,*) 
          real( kind = c_double ) :: potentials(*) 
          real( kind = c_double ) :: gradients(3,*)
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3d_evaluate

       subroutine fmmv3df_evaluate( &
           fmmvHandle,          & 
           charges,             &
           dipoleMoments,       &
           potentials,          &
           gradients,           &
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3df_evaluate" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr), value :: fmmvHandle
          real( kind = c_float ) :: charges(*) 
          real( kind = c_float ) :: dipoleMoments(3,*) 
          real( kind = c_float ) :: potentials(*) 
          real( kind = c_float ) :: gradients(3,*)
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3df_evaluate

       subroutine fmmv3d_finalize( &
           fmmvHandle,          & 
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3d_finalize" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr), value :: fmmvHandle
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3d_finalize
       
       subroutine fmmv3df_finalize( &
           fmmvHandle,          & 
           statistics,          &
           errorMessage         &
           ) bind( c, name="fmmv3df_finalize" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types 
          implicit none
          type (c_ptr), value :: fmmvHandle
          type( FmmvStatistics ) :: statistics
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3df_finalize


       subroutine fmmv3d_direct( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       &           
           NTargets,            &
           targets,             &
           potentials,          &
           gradients,           &
           accuracy,            &
           beta,                &
           time,                &
           errorMessage         &
           ) bind( c, name="fmmv3d_direct" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types !, only: FmmvOptions, FmmvStatistics
          implicit none
          integer( kind = c_int), value :: NParticles
          real( kind = c_double) :: particles(3,*)
          real( kind = c_double) :: charges(*)
          real( kind = c_double) :: dipoleMoments(3,*)
          integer( kind = c_int), value :: NTargets
          real( kind = c_double) :: targets(3,*)
          real( kind = c_double) :: potentials(*)
          real( kind = c_double) :: gradients(3,*)
          integer( kind = c_int), value :: accuracy
          real( kind = c_double), value :: beta
          real( kind = c_double) :: time
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3d_direct

       subroutine fmmv3df_direct( &
           NParticles,          & 
           particles,           &
           charges,             &
           dipoleMoments,       &           
           NTargets,            &
           targets,             &
           potentials,          &
           gradients,           &
           accuracy,            &
           beta,                &
           time,                &
           errorMessage         &
           ) bind( c, name="fmmv3df_direct" )
          use, intrinsic :: iso_c_binding
          use fmmv3d_types !, only: FmmvOptions, FmmvStatistics
          implicit none
          integer( kind = c_int), value :: NParticles
          real( kind = c_float) :: particles(3,*)
          real( kind = c_float) :: charges(*)
          real( kind = c_float) :: dipoleMoments(3,*)
          integer( kind = c_int), value :: NTargets
          real( kind = c_float) :: targets(3,*)
          real( kind = c_float) :: potentials(*)
          real( kind = c_float) :: gradients(3,*)
          integer( kind = c_int), value :: accuracy
          real( kind = c_double), value :: beta
          real( kind = c_double) :: time
          type( c_ptr ) :: errorMessage
       end subroutine fmmv3df_direct

   end interface

end module fmmv3d_module



