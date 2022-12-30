#include "cppdefs.h"

   module mod_tamu_mp

#ifdef HCGCM
      implicit none

      public

   contains

      subroutine tamu_gather_2d_field(ng, model, gtype, &
                                      LBi, UBi, LBj, UBj, &
                                      local_field, global_field)

          use mod_kinds
          use mod_scalars, only : exit_flag, NoError
          use mod_param,   only : IOBOUNDS
          use mod_ncparam, only : p2dvar, r2dvar, u2dvar, v2dvar
          use mod_grid,    only : GRID

#if defined DISTRIBUTE
          use distribute_mod, ONLY : mp_bcasti, mp_gather2d
#endif

          implicit none

          integer,  intent(in)  :: ng, model, gtype
          integer,  intent(in)  :: LBi, UBi, LBj, UBj
          real(r8), intent(in)  :: local_field(LBi:UBi,LBj:UBj)
          real(r8), intent(out) :: global_field(:)

          integer :: Imin, Imax, Jmin, Jmax, Npts

# ifdef MASKING
          ! Masking while retrieving?
          real(r8), dimension(:,:), pointer :: mask
# endif

!
!        Set first and last global grid point according to
!        staggered C-grid classification.
!
          SELECT CASE (ABS(gtype))
            CASE (p2dvar)
              Imin=IOBOUNDS(ng)%ILB_psi
              Imax=IOBOUNDS(ng)%IUB_psi
              Jmin=IOBOUNDS(ng)%JLB_psi
              Jmax=IOBOUNDS(ng)%JUB_psi
              mask => GRID(ng) % pmask_full
            CASE (r2dvar)
              ! No change for atmosphere
              Imin=IOBOUNDS(ng)%ILB_rho
              Imax=IOBOUNDS(ng)%IUB_rho
              Jmin=IOBOUNDS(ng)%JLB_rho
              Jmax=IOBOUNDS(ng)%JUB_rho
              mask => GRID(ng) % rmask_full
            CASE (u2dvar)
              Imin=IOBOUNDS(ng)%ILB_u
              Imax=IOBOUNDS(ng)%IUB_u
              Jmin=IOBOUNDS(ng)%JLB_u
              Jmax=IOBOUNDS(ng)%JUB_u
              mask => GRID(ng) % umask_full
            CASE (v2dvar)
              Imin=IOBOUNDS(ng)%ILB_v
              Imax=IOBOUNDS(ng)%IUB_v
              Jmin=IOBOUNDS(ng)%JLB_v
              Jmax=IOBOUNDS(ng)%JUB_v
              mask => GRID(ng) % vmask_full
            CASE DEFAULT
              Imin=IOBOUNDS(ng)%ILB_rho
              Imax=IOBOUNDS(ng)%IUB_rho
              Jmin=IOBOUNDS(ng)%JLB_rho
              Jmax=IOBOUNDS(ng)%JUB_rho
              mask => GRID(ng) % rmask_full
          END SELECT


# ifdef DISTRIBUTE
!
!       If distributed memory, collects data from all tasks and store it
!       into a global scratch 1D array, packed in column-major order
!
          CALL mp_gather2d(ng, model, LBi, UBi, LBj, UBj,        &
                           0, gtype, 1.0_r8,                     &
#  ifdef MASKING
                           mask,                                 &
#  endif
                           local_field,    &
                           Npts, global_field, SetFillVal = .false.)

          NULLIFY(mask)
          IF (Npts.ne.(Jmax-Jmin+1)*(Imax-Imin+1)) THEN
            exit_flag=6
            CALL mp_bcasti(ng, model, exit_flag)
          END IF

# else
!
!       If serial or shared-memory mode and serial output, store data
!       into a global scratch 1D array, packed in column-major order
!

          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              global_field(ic)=local_field(i,j)
            END DO
          END DO
          Npts=ic
# endif

      end subroutine tamu_gather_2d_field


      subroutine tamu_scatter_2d_field(ng, model, gtype, &
                                      LBi, UBi, LBj, UBj, &
                                      global_field, local_field)

          use mod_kinds
          use mod_param,   only : IOBOUNDS, iADM, NghostPoints
          use mod_ncparam, only : p2dvar, r2dvar, u2dvar, v2dvar
          use mod_iounits,     only : stdout

#if defined DISTRIBUTE
          use distribute_mod, ONLY : mp_bcasti, mp_scatter2d
#endif

          implicit none

          integer,  intent(in)  :: ng, model, gtype
          integer,  intent(in)  :: LBi, UBi, LBj, UBj
          real(r8), intent(out)  :: local_field(LBi:UBi,LBj:UBj)
          real(r8), intent(inout) :: global_field(:)

          integer  :: Imin, Imax, Ilen, Jmin, Jmax, Jlen, Nghost, Nopts
          integer  :: i, ic, ibuffer(1), j
          real(r8) :: f_min, f_max

!
!        Set first and last global grid point according to
!        staggered C-grid classification.
!
          SELECT CASE (ABS(gtype))
            CASE (p2dvar)
              Imin=IOBOUNDS(ng)%ILB_psi
              Imax=IOBOUNDS(ng)%IUB_psi
              Jmin=IOBOUNDS(ng)%JLB_psi
              Jmax=IOBOUNDS(ng)%JUB_psi
            CASE (r2dvar)
              ! No change for atmosphere
              Imin=IOBOUNDS(ng)%ILB_rho
              Imax=IOBOUNDS(ng)%IUB_rho
              Jmin=IOBOUNDS(ng)%JLB_rho
              Jmax=IOBOUNDS(ng)%JUB_rho
            CASE (u2dvar)
              Imin=IOBOUNDS(ng)%ILB_u
              Imax=IOBOUNDS(ng)%IUB_u
              Jmin=IOBOUNDS(ng)%JLB_u
              Jmax=IOBOUNDS(ng)%JUB_u
            CASE (v2dvar)
              Imin=IOBOUNDS(ng)%ILB_v
              Imax=IOBOUNDS(ng)%IUB_v
              Jmin=IOBOUNDS(ng)%JLB_v
              Jmax=IOBOUNDS(ng)%JUB_v
            CASE DEFAULT
              Imin=IOBOUNDS(ng)%ILB_rho
              Imax=IOBOUNDS(ng)%IUB_rho
              Jmin=IOBOUNDS(ng)%JLB_rho
              Jmax=IOBOUNDS(ng)%JUB_rho
          END SELECT

          Ilen=Imax-Imin+1
          Jlen=Jmax-Jmin+1

          Nopts = Ilen * Jlen
          global_field(Nopts + 1) = 0._r8
          global_field(Nopts + 2) = 0._r8


# ifdef DISTRIBUTE
!
!  Set the number of tile ghost points, Nghost, to scatter in
!  distributed-memory applications. If Nghost=0, the ghost points
!  are not processed.  They will be processed elsewhere by the
!  appropriate call to any of the routines in "mp_exchange.F".
!
#  ifdef NO_READ_GHOST
          Nghost=0
#  else
          IF (model.eq.iADM) THEN
            Nghost=0
          ELSE
            Nghost = NghostPoints
          END IF
#  endif

          CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,   &
                             Nghost, gtype, f_min, f_max,     &
#  if defined READ_WATER && defined MASKING
                             0, ibuffer,                      &
#  endif
                             Nopts, global_field, local_field)

# else
!
!       If serial or shared-memory mode and serial output, store data
!       into a global scratch 1D array, packed in column-major order
!
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              local_field(i,j) = global_field(ic)
            END DO
          END DO

# endif

      end subroutine tamu_scatter_2d_field


#endif
   end module mod_tamu_mp
