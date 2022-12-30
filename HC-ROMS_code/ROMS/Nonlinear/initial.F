#include "cppdefs.h"
      SUBROUTINE initial
!
!svn $Id: initial.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all model variables.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
      USE close_io_mod,      ONLY : close_inp
      USE dateclock_mod,     ONLY : time_string
#if defined FOUR_DVAR && !defined CORRELATION
      USE def_ini_mod,       ONLY : def_ini
#endif
#ifdef DISTRIBUTE
      USE distribute_mod,    ONLY : mp_bcasti, mp_bcastf, mp_bcasts
#endif
      USE get_state_mod,     ONLY : get_state
#ifdef WET_DRY
      USE get_wetdry_mod,    ONLY : get_wetdry
#endif
#ifdef TLM_CHECK
      USE ini_adjust_mod,    ONLY : ini_perturb
#endif
      USE ini_hmixcoef_mod,  ONLY : ini_hmixcoef
#ifdef NESTING
      USE nesting_mod,       ONLY : nesting
#endif
#ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth0, set_depth
      USE omega_mod,         ONLY : omega
      USE rho_eos_mod,       ONLY : rho_eos
      USE set_massflux_mod,  ONLY : set_massflux
#endif
#if defined OBSERVATIONS && !defined RBL4DVAR_FCT_SENSITIVITY
      USE obs_initial_mod,   ONLY : obs_initial
#endif
#ifdef MASKING
      USE set_masks_mod,     ONLY : set_masks
#endif
      USE stiffness_mod,     ONLY : stiffness
#if defined WAV_COUPLING && defined MCT_LIB
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
#endif
      USE strings_mod,       ONLY : FoundError
#ifdef WET_DRY
      USE wetdry_mod,        ONLY : wetdry
#endif
#if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod,       ONLY : wpoints
#endif
#ifdef HCGCM
      USE netcdf
      USE mpi
      USE hcgcm_filter_tools
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical :: update = .FALSE.
!
      integer :: Fcount
      integer :: ng, thread, tile
#ifdef NESTING
      integer :: ig, nl
      integer :: cr, i, m
#endif

      integer, dimension(Ngrids) :: IniRec, Tindex

#if defined ADJUST_BOUNDARY || \
    defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
      integer :: irec
#endif
!
#ifdef SP4DVAR
      real(dp) :: my_dstart
!
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__

#ifdef HCGCM
      ! integer,              parameter :: ocn_iunit     = 10
      ! character( len = * ), parameter :: ocn_drv_file  = 'HCGCM.in' ! by yman
      ! integer,       parameter :: ocn_max_grids = 1
      ! integer,  parameter :: filename_len = 1024
      ! integer, parameter  :: tamu_filter_opt_none      = 0, &
      !                        tamu_filter_opt_sst_sma   = 1, &
      !                        tamu_filter_opt_sst_loess = 2
      ! integer, parameter  :: tamu_filter_opt_miss_zonal_avg = 0, &
      !                        tamu_filter_opt_miss_local_avg = 1, &
      !                        tamu_filter_opt_miss_original  = 99
      integer :: ierr, io_stat, ip, jp,ni, rc
      integer :: ncid, lat_id, lon_id, lat_ndims, lon_ndims, nx, ny
      character(len = filename_len), dimension(ocn_max_grids) :: &
            filter_name, grid_input, lon_name, lat_name
      integer , dimension(ocn_max_grids) :: filter_nx, &
                                            filter_ny, &
                                            filter_miss_opt
      integer, dimension(nf90_max_var_dims) :: lat_dimids, lon_dimids
      real(dp), dimension(ocn_max_grids) :: filter_dx, &
                                            filter_dy, &
                                            filter_tolerance
      real(r8), dimension(:,:), allocatable :: lat, lon
      logical :: exists

      ! type t_hcgcm_filter
      !    character( len = filename_len ) :: name
      !    integer  :: option
      !    integer  :: nx, ny
      !    real(r8) :: dx, dy
      !    real(r8) :: tol
      !    integer  :: miss
      ! end type t_hcgcm_filter

      ! type (t_hcgcm_filter), dimension(Ngrids) :: hcgcm_filter

      namelist /ocn_filter/     &
            grid_input,  &
            lon_name,    &
            lat_name,    &
            filter_name, &
            filter_nx,   &
            filter_ny,   &
            filter_dx,   &
            filter_dy,   &
            filter_tolerance, &
            filter_miss_opt
#endif
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
!$OMP MASTER
      IF (Master) THEN
#if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
#elif defined I4DVAR      || defined RBL4DVAR          || \
      defined R4DVAR      || defined SENSITIVITY_4DVAR || \
      defined TL_RBL4DVAR || defined TL_R4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
#endif
        WRITE (stdout,20) 'INITIAL: Configuring and initializing ',     &
     &                    'forward nonlinear model ...'
 20     FORMAT (/,1x,a,a,/,1x,'*******')
      END IF
!$OMP END MASTER
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        iif(ng)=1
        indx1(ng)=1
        kstp(ng)=1
        krhs(ng)=1
        knew(ng)=1
        PREDICTOR_2D_STEP(ng)=.FALSE.
!
        iic(ng)=0
        nstp(ng)=1
        nrhs(ng)=1
        nnew(ng)=1
#ifdef FLOATS
        nf(ng)=0
        nfp1(ng)=1
        nfm1(ng)=4
        nfm2(ng)=3
        nfm3(ng)=2
#endif
!
        IniRec(ng)=nrrec(ng)
        Tindex(ng)=1
!
        synchro_flag(ng)=.TRUE.
        first_time(ng)=0
#ifdef SP4DVAR
        IF (INItime(ng).lt.0.0_dp) THEN
          my_dstart=dstart
        ELSE
          my_dstart=INItime(ng)/86400.0_dp
        END IF
#else
        tdays(ng)=dstart
#endif
        time(ng)=tdays(ng)*day2sec
!$OMP MASTER
        ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
        ntend(ng)=ntstart(ng)+ntimes(ng)-1
        ntfirst(ng)=ntstart(ng)
!$OMP END MASTER
!$OMP BARRIER
        step_counter(ng)=0
      END DO
!
!  Initialize global diagnostics variables.
!
      avgke=0.0_dp
      avgpe=0.0_dp
      avgkp=0.0_dp
      volume=0.0_dp

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 2, __LINE__, MyFile)
        END DO
      END DO
!$OMP BARRIER
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------
!
!  Set initial conditions record to process. If applicable open existing
!  nonlinear model initial conditions NetCDF file and, if needed, define
!  new variables. Then, inquire about available variables.
!
      DO ng=1,Ngrids
        IF (LdefINI(ng)) THEN
          LdefINI(ng)=.FALSE.          ! needed to inquire variables IDs
# ifndef CORRELATION
          CALL def_ini (ng)
# endif
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IniRec(ng)=nrrec(ng)
          INI(ng)%Rindex=IniRec(ng)
        ELSE
          IniRec(ng)=INI(ng)%Rindex
        END IF
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      DO ng=1,Ngrids
        OBCcount(ng)=0
      END DO
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      DO ng=1,Ngrids
        SFcount(ng)=0
      END DO
# endif
!
!  Reset nonlinear history time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
      DO ng=1,Ngrids
        HIS(ng)%Rindex=0
        Fcount=HIS(ng)%Fcount
        HIS(ng)%Nrec(Fcount)=0
      END DO

# ifdef I4DVAR
!
!  Activate switches to writting data into average, history and
!  restart files.
!
      DO ng=1,Ngrids
        LwrtAVG(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        LwrtRST(ng)=.TRUE.
      END DO
# endif
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  Set application grid, metrics, and associated variables and
!  parameters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (SetGridConfig(ng)) THEN
          CALL set_grid (ng, iNLM)
          SetGridConfig(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
#ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
#endif
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_hmixcoef (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO

#ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (Lsponge(ng)) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_sponge (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================

#ifdef TLM_CHECK
!
!  Clear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_depth (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------

#ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_initial (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_passive (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_biology (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_SEDIMENT && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!

      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_sediment (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined INI_FILE && !defined RBL4DVAR_FCT_SENSITIVITY
!
!  Read in initial conditions from initial NetCDF file.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL get_state (ng, iNLM, 1, INI(ng), IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        time(ng)=io_time                     ! needed for shared-memory
      END DO
#else
!
!  If restart, read in initial conditions restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_state (ng, 0, 1, INI(ng), IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          time(ng)=io_time                   ! needed for shared-memory
        END IF
      END DO
#endif

#ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!
!  If restart, read in wet/dry masks.
!
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
# ifdef DISTRIBUTE
          CALL get_wetdry (ng, MyRank, iNLM, IniRec(ng))
# else
          CALL get_wetdry (ng, -1, iNLM, IniRec(ng))
# endif
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined OBSERVATIONS && !defined RBL4DVAR_FCT_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations. Need to be done after processing initial conditions since
!  the correct initial time is needed to determine the first "ObsTime"
!  to process.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL obs_initial (ng, iNLM, .FALSE.)
!$OMP END MASTER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!$OMP BARRIER
#endif

#if (defined ADJUST_BOUNDARY          || \
     defined ADJUST_STFLUX            || \
     defined ADJUST_WSTRESS)          && \
    (defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR)
!
!-----------------------------------------------------------------------
!  Read in the surface forcing and or open boundary conditions
!  increments for RBL4D-Var from record IniRec of the NLM initial
!  NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_state (ng, 5, 5, INI(ng), IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#endif

#ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variable according to the outer
!  loop iteration with the steepest descent direction of the gradient
!  (adjoint state).
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_perturb (ng, tile, Lnew(ng), Tindex(ng))
          END DO
!$OMP BARRIER
        END DO
      END IF
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time independent (Zt_avg1=0) anf initial time dependent
!  depths and level thicknesses.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_depth0 (ng, tile, iNLM)
          CALL set_depth  (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_massflux (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL omega (ng, tile, iNLM)
          CALL rho_eos (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#ifdef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_psource (ng, tile, iNLM)
          END DO
        END IF
!$OMP BARRIER
      END DO
#endif

#if defined FOUR_DVAR || !defined TANGENT || !defined ADJOINT
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------

# ifdef ADJUST_BOUNDARY
!
!  If first pass of iteration loop, set time of open boundary
!  adjustment.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          OBC_time(1,ng)=time(ng)
          DO irec=2,Nbrec(ng)
            OBC_time(irec,ng)=OBC_time(irec-1,ng)+nOBC(ng)*dt(ng)
          END DO
        END DO
      END IF
!$OMP END MASTER
!$OMP BARRIER
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  If first pass of iteration loop, set time of surface forcing
!  adjustment.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          SF_time(1,ng)=time(ng)
          DO irec=2,Nfrec(ng)
            SF_time(irec,ng)=SF_time(irec-1,ng)+nSFF(ng)*dt(ng)
          END DO
        END DO
      END IF
!$OMP END MASTER
!$OMP BARRIER
# endif
# if !defined CORRELATION
!
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL close_inp (ng, iNLM)
        CALL check_multifile (ng, iNLM)
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  If applicable, read in input data.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL get_idata (ng)
        CALL get_data (ng)
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
#endif

#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_masks (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#if !defined CORRELATION
# ifdef NESTING
#  if defined MASKING || defined WET_DRY
!
!-----------------------------------------------------------------------
!  If nesting and Land/Sea masking, scale horizontal interpolation
!  weights to account for land contact points.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL nesting (ng, iNLM, nmask)
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  If nesting, process state fields initial conditions in the contact
!  regions.
!-----------------------------------------------------------------------
!
!  Free-surface and 2D-momentum.
!
      DO nl=1,NestLayers
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          IF (ANY(CompositeGrid(:,ng))) THEN
            CALL nesting (ng, iNLM, nFSIC)        ! free-surface
#  ifndef SOLVE3D
            CALL nesting (ng, iNLM, n2dIC)        ! 2d momentum
#  endif
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
!
!  Determine vertical indices and vertical interpolation weights in
!  the contact zone using initial unperturbed depth arrays.
!
      DO ng=1,Ngrids
        CALL nesting (ng, iNLM, nzwgt)
      END DO
!
!  3D-momentum and tracers.
!
      DO nl=1,NestLayers
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          IF (ANY(CompositeGrid(:,ng))) THEN
            CALL nesting (ng, iNLM, n3dIC)        ! 3D momentum
            CALL nesting (ng, iNLM, nTVIC)        ! Tracer variables
          END IF
        END DO
      END DO
#  endif
# endif
#endif

#if defined PROPAGATOR || \
   (defined MASKING    && (defined READ_WATER || defined WRITE_WATER ))
!
!-----------------------------------------------------------------------
!  Set variables associated with the processing water points and/or
!  size of packed state arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL wpoints (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#if defined NLM_OUTER                || \
    defined RBL4DVAR                 || \
    defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined TL_RBL4DVAR
!
!-----------------------------------------------------------------------
!  Read in convolved adjoint impulse forcing (first record) and its
!  application time.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (SporadicImpulse(ng)) THEN
          FrcRec(ng)=1
!$OMP MASTER
          CALL get_state (ng, 7, 7, TLF(ng), FrcRec(ng), 1)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#endif

#if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_drag (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL stiffness (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
      END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL grid_coords (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
      END DO
#endif

# if defined WAV_COUPLING && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ocn2wav_coupling (ng, tile)
        END DO
!$OMP BARRIER
        IF (Master) WRITE (stdout,'(/)')
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and date/time string. Save NLM
!  initial conditions time.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        INItime(ng)=time(ng)
        iic(ng)=ntstart(ng)
        CALL time_string (time(ng), time_code(ng))
      END DO

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 2, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
!

#ifdef HCGCM
!
!-----------------------------------------------------------------------
!  inital filter
!  by yman
!-----------------------------------------------------------------------
!
      it_ou = 0
      ierr = 0
      IF (Master) WRITE (stdout,*) 'HCGCM inital (by yman):'
      if (Ngrids > ocn_max_grids)  WRITE (stdout,*) 'error!'
      if (Master) then
        inquire(file = trim(ocn_drv_file), exist = exists)
        WRITE (stdout,*) 'exists:', exists
        if (exists) then
          open(unit = ocn_iunit, file = trim(ocn_drv_file),     &
          form = 'formatted', action = 'read', iostat = io_stat)
          if (io_stat .ne. 0) then
            ierr = -1
            write(stdout, 10) trim(ocn_drv_file)
          endif
        endif
        if (ierr .eq. 0) then
          rewind(unit = ocn_iunit)
          read(unit = ocn_iunit, nml = ocn_filter , iostat = io_stat)
          if (io_stat .ne. 0) then
            ierr = -1
            write(stdout, 40) '&ocn_filter', trim(ocn_drv_file)
          end if
        end if
      endif

      ! call MPI_Bcast(ierr, 1, MPI_INTEGER, MyMaster, mpi_comm, io_stat)
      call mp_bcasti(ng, iNLM, ierr)

      if (ierr .eq. 0) then
            do ng = 1, Ngrids
                  call mp_bcasts(ng, iNLM, grid_input(ng))
                  call mp_bcasts(ng, iNLM, lon_name(ng))
                  call mp_bcasts(ng, iNLM, lat_name(ng))
                  call mp_bcasts(ng, iNLM, filter_name(ng))
                  call mp_bcasti(ng, iNLM, filter_nx(ng))
                  call mp_bcasti(ng, iNLM, filter_ny(ng))
                  call mp_bcastf(ng, iNLM, filter_dx(ng))
                  call mp_bcastf(ng, iNLM, filter_dy(ng))
                  call mp_bcastf(ng, iNLM, filter_tolerance(ng))
                  call mp_bcasti(ng, iNLM, filter_miss_opt(ng))
                  ! call MPI_Bcast(filter_name     , ng*filename_len,
                  ! MPI_CHARACTER       , MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_nx       , ng,
                  ! MPI_INTEGER         , MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_ny       , ng,
                  ! MPI_INTEGER         , MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_dx       , ng,
                  ! MPI_DOUBLE_PRECISION, MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_dy       , ng,
                  ! MPI_DOUBLE_PRECISION, MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_tolerance, ng,
                  ! MPI_DOUBLE_PRECISION, MyMaster, mpi_comm, io_stat)
                  ! call MPI_Bcast(filter_miss_opt , ng,
                  ! MPI_INTEGER         , MyMaster, mpi_comm, io_stat)
            end do

      endif

      if (MASTER) then
            WRITE (stdout,*) 'grid_input=', grid_input
            WRITE (stdout,*) 'lon_rho=', lon_name
            WRITE (stdout,*) 'lat_rho=', lat_name
            WRITE (stdout,*) 'filter_name=', filter_name
            WRITE (stdout,*) 'filter_nx=', filter_nx
            WRITE (stdout,*) 'filter_ny=', filter_ny
            WRITE (stdout,*) 'filter_dx=', filter_dx
            WRITE (stdout,*) 'filter_dy=', filter_dy
            WRITE (stdout,*) 'filter_tolerance=', filter_tolerance
            WRITE (stdout,*) 'filter_miss_opt=', filter_miss_opt
      endif

      allocate(hcgcm_filter(Ngrids))
      allocate(hcgcm_grid(Ngrids))

      do ng = 1, Ngrids
            if (MASTER) then

                  call check( nf90_open(trim(grid_input(ng)), nf90_nowrite,ncid) )
                  call check( nf90_inq_varid(ncid, trim(lon_name(ng)), lon_id) )
                  call check( nf90_inquire_variable(ncid, lon_id, ndims = lon_ndims, dimids = lon_dimids) )
                  call check( nf90_inquire_dimension(ncid, lon_dimids(1), len = nx) )
                  call check( nf90_inquire_dimension(ncid, lon_dimids(2), len = ny) )

                  call check( nf90_inq_varid(ncid, trim(lat_name(ng)), lat_id) )
                  call check( nf90_inquire_variable(ncid, lat_id, ndims = lat_ndims, dimids = lat_dimids) )

                  allocate(lon(nx, ny))
                  allocate(lat(nx, ny))

                  call check( nf90_get_var(ncid, lon_id, lon) )
                  call check( nf90_get_var(ncid, lat_id, lat) )

                  do jp = 1, ny
                        do ip = 2, nx
                              if (lon(ip,jp)-lon(ip-1,jp) .lt. 0._r8) &
                                    lon(ip,jp) = lon(ip,jp) + 360._r8
                        end do
                  end do
                  call check( nf90_close(ncid) )
                  WRITE (stdout,*) 'i am master'
                  WRITE (stdout,*) 'myrank=', myrank, 'nx=', nx, 'ny=',ny, 'shape=', shape(lon), 'shape2=', shape(lat)
                  WRITE (stdout,*) 'i am master'

            endif
            call mp_bcasti(ng, iNLM, nx)
            call mp_bcasti(ng, iNLM, ny)

            hcgcm_filter(ng) % name   = ""
            hcgcm_filter(ng) % option = 0
            hcgcm_filter(ng) % nx     = 0
            hcgcm_filter(ng) % ny     = 0
            hcgcm_filter(ng) % dx     = 0._r8
            hcgcm_filter(ng) % dy     = 0._r8
            hcgcm_filter(ng) % tol    = 0._r8
            hcgcm_filter(ng) % miss   = 0

            hcgcm_grid(ng) % LBi = BOUNDS(ng)%LBi(MyRank)
            hcgcm_grid(ng) % UBi = BOUNDS(ng)%UBi(MyRank)
            hcgcm_grid(ng) % LBj = BOUNDS(ng)%LBj(MyRank)
            hcgcm_grid(ng) % UBj = BOUNDS(ng)%UBj(MyRank)

            hcgcm_grid(ng) % west_east = nx
            hcgcm_grid(ng) % south_north = ny

            allocate(hcgcm_grid(ng) % Long &
                  (hcgcm_grid(ng) % west_east,   &
                   hcgcm_grid(ng) % south_north))
            allocate(hcgcm_grid(ng) % Lat &
                  (hcgcm_grid(ng) % west_east,   &
                   hcgcm_grid(ng) % south_north))

            if (MASTER) then
              hcgcm_grid(ng) % Long = lon
              hcgcm_grid(ng) % Lat = lat
              WRITE (stdout,*) 'i am master'
              WRITE (stdout,*) 'myrank=', myrank, 'nx=', nx, 'ny=',ny, &
                'shape=', shape(hcgcm_grid(ng) % Long), 'shape2=', shape(hcgcm_grid(ng) % Lat)
              WRITE (stdout,*) 'i am master'
            endif

            call mp_bcastf(ng, iNLM, hcgcm_grid(ng) % Long)
            call mp_bcastf(ng, iNLM, hcgcm_grid(ng) % Lat)

            ! hcgcm_grid(ng) % Long = lon
            ! hcgcm_grid(ng) % Lat = lat
            ! WRITE (stdout,*) 'myrank=', myrank, 'nx=', nx, 'ny=',ny, &
            !   'shape=', shape(hcgcm_grid(ng) % Long), 'shape2=',  shape(hcgcm_grid(ng) % Lat)

            select case (trim(filter_name(ng)))
                  case ('','none','None','NONE')
                        hcgcm_filter(ng) % option = tamu_filter_opt_none
                  case ('sst_sma','SST_SMA')
                        hcgcm_filter(ng) % option = tamu_filter_opt_sst_sma
                  case('sst_loess','SST_LOESS')
                        hcgcm_filter(ng) % option = tamu_filter_opt_sst_loess
                  case default
                        hcgcm_filter(ng) % option = tamu_filter_opt_none
            end select

            if (hcgcm_filter(ng) % option .ne. tamu_filter_opt_none) then
                  if ((filter_nx(ng).gt.0) .and. (filter_ny(ng).gt.0)) then
                        hcgcm_filter(ng) % nx   = filter_nx(ng)
                        hcgcm_filter(ng) % ny   = filter_ny(ng)
                  else if ((filter_dx(ng).gt.0) .and. (filter_dy(ng).gt.0)) then
                        hcgcm_filter(ng) % name = trim(filter_name(ng))
                        hcgcm_filter(ng) % dx   = filter_dx(ng)
                        hcgcm_filter(ng) % dy   = filter_dy(ng)
                  else
                        hcgcm_filter(ng) % option = tamu_filter_opt_none
                  end if
            end if
            hcgcm_filter(ng) % miss = filter_miss_opt(ng)
            ! if (hcgcm_filter(ng) % option .ne. tamu_filter_opt_none) then
            !       hcgcm_filter(ng) % name = trim(filter_name(ng))
            !       if (hcgcm_filter(ng) % option .eq.
            !       tamu_filter_opt_sst_loess) then
            !             hcgcm_filter(ng) % tol  = filter_tolerance(ng)
            !             if (      (filter_miss_opt(ng) .ne.
            !             tamu_filter_opt_miss_zonal_avg) &
            !                   .and. (filter_miss_opt(ng) .ne.
            !                   tamu_filter_opt_miss_local_avg) &
            !                   .and. (filter_miss_opt(ng) .ne.
            !                   tamu_filter_opt_miss_original )) then
            !                   ierr = -1
            !                   if (OutThread) then
            !                         write(stdout, '("ocn_init: Unknown
            !                         filter_miss_opt.")')
            !                         write(stdout, '("ocn_init: Error
            !                         initializing ROMS. Aborting...")')
            !                   end if
            !                   return
            !             end if
            !       end if
            ! end if
      end do

      ! deallocate(lat)
      ! deallocate(lon)

      ! WRITE (stdout,*) 'myrank=', MyRank, &
      ! 'west=', hcgcm_grid%west_east, &
      ! 'south=',  hcgcm_grid%south_north

10     format(1x,"Error opening file ",a)
40     format(1x,"Error reading namelist ",a," in file ",a)
#endif


      RETURN
      END SUBROUTINE initial
