#include "cppdefs.h"

module mod_tamu_filter

#ifdef HCGCM1

      use mod_kinds,       only : r8
      use mod_param,       only : NtileI, NtileJ
      use mod_scalars,     only : exit_flag, NoError
      use mod_iounits,     only : stdout
# if defined DISTRIBUTE
      use mod_parallel,    only : MyRank, MyMaster, OCN_COMM_WORLD, &
                                  MPI_SUCCESS, MPI_SUM,  &
                                  MP_FLOAT, MPI_INTEGER, &
                                  OutThread
      use distribute_mod,  only : mp_barrier, mp_bcastf
      use hcgcm_filter_type, only : hcgcm_filter, log_error, hcgcm_filter
# endif

      implicit none

      integer, parameter  :: tamu_filter_opt_none      = 0, &
                             tamu_filter_opt_sst_sma   = 1, &
                             tamu_filter_opt_sst_loess = 2

      integer, parameter  :: tamu_filter_opt_miss_zonal_avg = 0, &
                             tamu_filter_opt_miss_local_avg = 1, &
                             tamu_filter_opt_miss_original  = 99


      private

      public :: tamu_filter_opt_none,    &
                tamu_filter_opt_sst_sma, &
                tamu_filter_opt_sst_loess

      public :: tamu_filter_opt_miss_zonal_avg, &
                tamu_filter_opt_miss_local_avg, &
                tamu_filter_opt_miss_original

      public :: tamu_run_filter

      public :: tamu_filter_simple_moving_average, &
                tamu_filter_loess_2d

      interface tamu_filter_simple_moving_average
          module procedure tamu_filter_sma_ij
          module procedure tamu_filter_sma_xy
      end interface

      interface tamu_filter_loess_2d
          module procedure tamu_filter_loess_ij
          module procedure tamu_filter_loess_xy
      end interface

  contains

      subroutine tamu_run_filter(ng, model, v, rc)

         implicit none

         ! --- Arguments ---
         integer,                  intent(in)    :: ng, model
         real(r8), dimension(:,:), intent(inout) :: v
         integer,                  intent(out)   :: rc

         ! --- Local variables
         integer  :: is, ie, ip, js, je, jp, nx, ny
         integer  :: ierr, nmiss, nmissed
         real(r8) :: wtime_beg, wtime_end
         real(r8), dimension(2) :: wtime
         real(r8), external     :: my_wtime
         real(r8), dimension(:,:), pointer :: x, y
         character( len = * ), parameter :: subname = 'tamu_run_filter'

         ! --- Begin

         rc = 0
         if (hcgcm_filter(ng) % option &
             .eq. tamu_filter_opt_none) return

#if defined DISTRIBUTE
         call mp_barrier(ng)
#endif
         if (OutThread) wtime_beg = my_wtime(wtime)

         nmiss   = 0
         nmissed = 0
         nullify(x,y)

         nx = TAMU_CPL(ng) % Atmo % west_east
         ny = TAMU_CPL(ng) % Atmo % south_north

         if (nx .ne. size(v, dim = 1)) then
            call log_error(subname, &
                           'Horizontal size of field /= grid size')
            rc = -1
         end if

         if (rc .eq. 0) then
            if (ny .ne. size(v, dim = 2)) then
               call log_error(subname, &
                              'Vertical size of field /= grid size')
               rc = -1
            end if
         end if

         if (rc .eq. 0) then
            ! --- Compute tiling
# ifdef DISTRIBUTE
            jp = MyRank / NtileI(ng)
            ip = MyRank - jp * NtileI(ng)

            is = ip * nx / NtileI(ng) + 1
            ie = (ip+1) * nx / NtileI(ng)
            if (ip .eq. NtileI(ng)-1) ie = nx

            js = jp * ny / NtileJ(ng) + 1
            je = (jp+1) * ny / NtileJ(ng)
            if (jp .eq. NtileJ(ng)-1) je = ny
# else
            is = 1
            ie = nx
            js = 1
            je = ny
# endif
            allocate(x(nx, ny), stat=ierr)
            if (ierr .ne. 0) then
                call log_error(subname, 'Insufficient memory')
                rc = -1
            end if
         end if

         if (rc .eq. 0) then
            ! --- Make sure longitude increases monotonically from W to E
            x = TAMU_CPL(ng) % Atmo % Long
            do jp = 1, ny
               do ip = 2, nx
                  if (x(ip,jp)-x(ip-1,jp) .lt. 0._r8) &
                      x(ip,jp) = x(ip,jp) + 360._r8
               end do
            end do

            y => TAMU_CPL(ng) % Atmo % Lat

# ifdef DISTRIBUTE
            call mp_bcastf (ng, model, v)
            if (exit_flag .ne. NoError) then
                call log_error(subname, 'Error broadcasting field values')
                deallocate(x)
                nullify(x, y)
                rc = -1
            end if
         end if
# endif

         if (rc .eq. 0) then
           select case (hcgcm_filter(ng) % option)

           case (tamu_filter_opt_sst_sma  )

             if ((hcgcm_filter(ng) % dx .gt. 0._r8) .and. &
                 (hcgcm_filter(ng) % dy .gt. 0._r8)) then
                call tamu_filter_simple_moving_average(v, x, y, &
                             hcgcm_filter(ng) % dx, &
                             hcgcm_filter(ng) % dy, &
                             mask=0._r8, its=is, ite=ie, jts=js, jte=je)
             else if ((hcgcm_filter(ng) % nx .gt. 0) .and. &
                      (hcgcm_filter(ng) % ny .gt. 0)) then
                call tamu_filter_simple_moving_average(v,       &
                             hcgcm_filter(ng) % nx, &
                             hcgcm_filter(ng) % ny, &
                             mask=0._r8, its=is, ite=ie, jts=js, jte=je)
             end if

           case (tamu_filter_opt_sst_loess)

             if ((hcgcm_filter(ng) % dx .gt. 0._r8) .and. &
                 (hcgcm_filter(ng) % dy .gt. 0._r8)) then
                call tamu_filter_loess_2d(v, x, y, &
                             hcgcm_filter(ng) % tol,  &
                             hcgcm_filter(ng) % miss, &
                             nmiss, &
                             hcgcm_filter(ng) % dx,   &
                             hcgcm_filter(ng) % dy,   &
                             ierr, &
                             mask=0._r8, its=is, ite=ie, jts=js, jte=je)
             else if ((hcgcm_filter(ng) % nx .gt. 0) .and. &
                      (hcgcm_filter(ng) % ny .gt. 0)) then
                call tamu_filter_loess_2d(v, &
                             hcgcm_filter(ng) % tol,  &
                             hcgcm_filter(ng) % miss, &
                             nmiss, &
                             hcgcm_filter(ng) % nx,   &
                             hcgcm_filter(ng) % ny,   &
                             ierr, &
                             mask=0._r8, its=is, ite=ie, jts=js, jte=je)
             end if
             if (ierr .ne. 0) rc = -1

           end select
         end if

#if defined DISTRIBUTE
         if (rc .eq. 0) then
           x = 0._r8
           call mpi_reduce(v, x, nx*ny, MP_FLOAT, MPI_SUM, &
                           MyMaster, OCN_COMM_WORLD, ierr)
           if (ierr .ne. MPI_SUCCESS) then
              call log_error(subname, 'Error collecting filter output')
              deallocate(x)
              nullify(x, y)
              rc = -1
           end if
         end if
         if (rc .eq. 0)  v = x
#endif

         if (rc .eq. 0) then
            deallocate(x, stat=ierr)
            if (ierr .ne. 0) then
                call log_error(subname, 'Error freeing memory')
                rc = -1
            end if
         end if

         nullify(x, y)

#if defined DISTRIBUTE

         if (rc .eq. 0) then
            call mpi_reduce(nmiss, nmissed, 1, MPI_INTEGER, MPI_SUM, &
                            MyMaster, OCN_COMM_WORLD, ierr)
            if (ierr .ne. MPI_SUCCESS) then
                call log_error(subname, &
                                'Error collecting filter missed values')
                rc = -1
            end if
         end if

#else
         nmissed = nmiss
#endif

#if defined DISTRIBUTE
         ip = 0
         call mpi_allreduce(rc, ip, 1, MPI_INTEGER, &
                            MPI_SUM, OCN_COMM_WORLD, ierr)
         if (ierr .ne. MPI_SUCCESS) then
            rc = -1
         else if (ip .ne. 0) then
            rc = -1
         end if

         call mp_barrier(ng)
#endif

         if (rc .ne. 0) exit_flag = 8

         if (OutThread) then
            wtime_end = my_wtime(wtime)
            write(stdout, &
             '(31x,".. Timing for filter on domain ",i0,&
             &": ",f12.5," elapsed seconds")') ng, wtime_end-wtime_beg
            if (nmissed .gt. 0) write(stdout, &
                      '(34x,"* replaced ",i0," missed values")') nmissed
         end if

      end subroutine tamu_run_filter


      subroutine tamu_filter_sma_ij(v, idx, idy, mask, its, ite, jts, jte)

         implicit none

         real(r8), dimension(:,:), intent(inout) :: v
         integer,                  intent(in)    :: idx, idy
         real(r8), optional,       intent(in)    :: mask
         integer,  optional,       intent(in)    :: its, ite, jts, jte

         integer  :: i, ia, i_beg, i_end, i_size, &
                     j, ja, j_beg, j_end, j_size, n
         integer  :: is, ie, js, je
         real(r8) :: s, v_mask
         real(r8), dimension(:,:), allocatable :: vp


         if (present(mask)) then
            v_mask = mask
         else
            v_mask = 0._r8
         end if

         i_size = size(v, dim = 1)
         j_size = size(v, dim = 2)

         allocate(vp(i_size, j_size))
         vp = 0._r8

         is = 1
         ie = i_size
         js = 1
         je = j_size

         if (present(its)) is = its
         if (present(ite)) ie = ite
         if (present(jts)) js = jts
         if (present(jte)) je = jte

         do j = js, je
            j_beg = max(1,     j-idy)
            j_end = min(j_size,j+idy)
            do i = is, ie
               i_beg = max(1,      i-idx)
               i_end = min(i_size, i+idx)
               vp(i,j) = v(i,j)
               if (v(i,j) .ne. v_mask) then
                  n = 0
                  s = 0._r8
                  do ja = j_beg, j_end
                     do ia = i_beg, i_end
                        if (v(ia,ja) .ne. v_mask) then
                           s = s + v(ia,ja)
                           n = n + 1
                        end if
                     end do
                  end do
                  if (n .gt. 0) vp(i,j) = s / n
               end if
            end do
         end do

         v = vp

         deallocate(vp)

      end subroutine tamu_filter_sma_ij

      subroutine tamu_filter_sma_xy(v, x, y, dx, dy, mask, &
                                    its, ite, jts, jte)

         implicit none

         real(r8), dimension(:,:), intent(inout) :: v
         real(r8), dimension(:,:), intent(in)    :: x, y
         real(r8),                 intent(in)    :: dx, dy
         real(r8), optional,       intent(in)    :: mask
         integer,  optional,       intent(in)    :: its, ite, jts, jte

         integer  :: i, ia, i_beg, i_end, i_size, idx, idy, &
                     j, ja, j_beg, j_end, j_size, n
         integer  :: is, ie, js, je
         real(r8) :: dmin, s, sn, sx, sxn, sy, syn, v_mask
         real(r8), dimension(:,:), allocatable :: vp


         if (present(mask)) then
            v_mask = mask
         else
            v_mask = 0._r8
         end if

         i_size = size(v, dim = 1)
         j_size = size(v, dim = 2)

         if (i_size .ne. size(x, dim = 1)) return
         if (i_size .ne. size(y, dim = 1)) return
         if (j_size .ne. size(x, dim = 2)) return
         if (j_size .ne. size(y, dim = 2)) return

         is = 1
         ie = i_size
         js = 1
         je = j_size

         if (present(its)) is = its
         if (present(ite)) ie = ite
         if (present(jts)) js = jts
         if (present(jte)) je = jte

         ! --- Check max size for regression

         dmin = minval(abs(x(2:,:)-x(1:i_size-1,:)))
         if (dmin .eq. 0._r8) return
         idx = int(dx/dmin)+1

         dmin = minval(abs(y(:,2:)-y(:,1:j_size-1)))
         if (dmin .eq. 0._r8) return
         idy = int(dy/dmin)+1

         allocate(vp(i_size, j_size))
         vp = 0._r8

         do j = js, je
            j_beg = max(1,     j-idy)
            j_end = min(j_size,j+idy)
            do i = is, ie
               i_beg = max(1,      i-idx)
               i_end = min(i_size, i+idx)
               vp(i,j) = v(i,j)
               if (v(i,j) .ne. v_mask) then
                  n = 0
                  s = 0._r8
                  do ja = j_beg, j_end
                     do ia = i_beg, i_end
                        if (v(ia,ja) .ne. v_mask) then
                           sx = x(ia,ja)-x(i,j)
                           sy = y(ia,ja)-y(i,j)
                           sxn = sx / dx
                           syn = sy / dy
                           sn = sqrt(sxn * sxn + syn * syn)
                           if (sn .le. 1._r8) then
                              s = s + v(ia,ja)
                              n = n + 1
                           end if
                        end if
                     end do
                  end do
                  if (n .gt. 0) vp(i,j) = s / n
               end if
            end do
         end do

         v = vp

         deallocate(vp)

      end subroutine tamu_filter_sma_xy

      subroutine tamu_filter_loess_ij(v, tol, opt, nmiss, idx, idy, rc,&
                                      mask, its, ite, jts, jte)

         implicit none

         real(r8), dimension(:,:), intent(inout) :: v
         real(r8),                 intent(in)    :: tol
         integer,                  intent(in)    :: opt
         integer,                  intent(in)    :: idx, idy
         real(r8), optional,       intent(in)    :: mask
         integer,  optional,       intent(in)    :: its, ite, jts, jte
         integer,                  intent(out)   :: nmiss
         integer,                  intent(out)   :: rc

         integer  :: i, ia, i_beg, i_end, i_size, &
                     j, ja, j_beg, j_end, j_size, &
                     l_size, n
         integer  :: is, ie, js, je
         real(r8) :: sx, sy, sn, v_mask, wn
         real(r8), dimension(:),   allocatable :: z, w
         real(r8), dimension(:,:), allocatable :: d, vp

         integer,  parameter :: reg_size_min = 10
         real(r8), parameter :: v_miss = -1.e+35_r8

         rc = -1

         i_size = size(v, dim = 1)
         j_size = size(v, dim = 2)

         l_size = (2 * idx + 1) * (2 * idy + 1)
         if (l_size .lt. reg_size_min) return

         if (present(mask)) then
            v_mask = mask
         else
            v_mask = 0._r8
         end if

         is = 1
         ie = i_size
         js = 1
         je = j_size

         if (present(its)) is = its
         if (present(ite)) ie = ite
         if (present(jts)) js = jts
         if (present(jte)) je = jte

         allocate(d(l_size,2))
         allocate(z(l_size))
         allocate(w(l_size))
         allocate(vp(i_size, j_size))
         d  = 0._r8
         z  = 0._r8
         w  = 0._r8
         vp = 0._r8

         nmiss = 0

         do j = js, je
            j_beg = max(1,     j-idy)
            j_end = min(j_size,j+idy)
            do i = is, ie
               i_beg = max(1,      i-idx)
               i_end = min(i_size, i+idx)
               vp(i,j) = v(i,j)
               if (v(i,j) .ne. v_mask) then
                  n = 0
                  do ja = j_beg, j_end
                     do ia = i_beg, i_end
                        if (v(ia,ja) .ne. v_mask) then
                           n = n + 1
                           sx = real(ia - i, kind = r8)
                           sy = real(ja - j, kind = r8)
                           d(n, 1) = sx
                           d(n, 2) = sy
                           z(n) = v(ia,ja)
                           sx = sx / idx
                           sy = sy / idy
                           sn = sqrt(sx * sx + sy * sy)
                           wn = 1._r8 - sn * sn * sn
                           w(n) = wn * wn * wn
                        end if
                     end do
                  end do
                  if (n .ge. reg_size_min) then
                     call regsm2(n,z,d,w,vp(i,j),v_miss,tol,l_size)
                  else
                     vp(i,j) = v_miss
                  end if
                  if (vp(i,j) .eq. v_miss) then
                     nmiss = nmiss + 1
                     select case (opt)
                     case (tamu_filter_opt_miss_zonal_avg)
                         ! --- Replace missing value with zonal average ---
                         n = count(v(:,j) .ne. v_mask)
                         if (n .gt. 0) then
                            vp(i,j) = &
                            sum(v(:,j), mask = (v(:,j).ne.v_mask))/n
                         end if
                     case (tamu_filter_opt_miss_local_avg)
                         ! --- Replace missing value with local average ---
                         if (n .gt. 0) vp(i,j) = sum(z(1:n))/n
                     case (tamu_filter_opt_miss_original)
                         ! --- Replace missing value with original value ---
                         vp(i,j) = v(i,j)
                     case default
                         ! --- Replace missing value with original value ---
                         vp(i,j) = v(i,j)
                     end select
                  end if
               else
                  vp(i,j) = v_mask
               end if
            end do
         end do

         v = vp

         rc = 0

         deallocate(d)
         deallocate(vp)
         deallocate(z)
         deallocate(w)

      end subroutine tamu_filter_loess_ij

      subroutine tamu_filter_loess_xy(v, x, y, tol, opt, nmiss, dx, dy,&
                                      rc, mask, its, ite, jts, jte)

         implicit none

         real(r8), dimension(:,:), intent(inout) :: v
         real(r8), dimension(:,:), intent(in)    :: x, y
         real(r8),                 intent(in)    :: tol
         integer,                  intent(in)    :: opt
         real(r8),                 intent(in)    :: dx, dy
         real(r8), optional,       intent(in)    :: mask
         integer,  optional,       intent(in)    :: its, ite, jts, jte
         integer,                  intent(out)   :: nmiss
         integer,                  intent(out)   :: rc

         integer  :: i, ia, i_beg, i_end, i_size, idx, &
                     j, ja, j_beg, j_end, j_size, idy, &
                     l_size, n
         integer  :: is, ie, js, je
         real(r8) :: dmin, sx, sxn, sy, syn, sn, v_mask, wn
         real(r8), dimension(:),   allocatable :: z, w
         real(r8), dimension(:,:), allocatable :: d, vp

         integer,  parameter :: reg_size_min = 10
         real(r8), parameter :: v_miss = -1.e+35_r8

         rc = -1

         i_size = size(v, dim = 1)
         j_size = size(v, dim = 2)

         if (i_size .ne. size(x, dim = 1)) return
         if (i_size .ne. size(y, dim = 1)) return
         if (j_size .ne. size(x, dim = 2)) return
         if (j_size .ne. size(y, dim = 2)) return

         ! --- Check max size for regression
         dmin = minval(abs(x(2:,:)-x(1:i_size-1,:)))
         if (dmin .eq. 0._r8) return
         idx = int(dx/dmin)+1

         dmin = minval(abs(y(:,2:)-y(:,1:j_size-1)))
         if (dmin .eq. 0._r8) return
         idy = int(dy/dmin)+1

         l_size = (2 * idx + 1) * (2 * idy + 1)
         if (l_size .lt. reg_size_min) return

         if (present(mask)) then
            v_mask = mask
         else
            v_mask = 0._r8
         end if

         is = 1
         ie = i_size
         js = 1
         je = j_size

         if (present(its)) is = its
         if (present(ite)) ie = ite
         if (present(jts)) js = jts
         if (present(jte)) je = jte

         allocate(d(l_size,2))
         allocate(z(l_size))
         allocate(w(l_size))
         allocate(vp(i_size, j_size))
         d  = 0._r8
         z  = 0._r8
         w  = 0._r8
         vp = 0._r8

         nmiss = 0

         do j = js, je
            j_beg = max(1,     j-idy)
            j_end = min(j_size,j+idy)
            do i = is, ie
               i_beg = max(1,      i-idx)
               i_end = min(i_size, i+idx)
               vp(i,j) = v(i,j)
               if (v(i,j) .ne. v_mask) then
                  n = 0
                  do ja = j_beg, j_end
                     do ia = i_beg, i_end
                        if (v(ia,ja) .ne. v_mask) then
                           sx = x(ia,ja)-x(i,j)
                           sy = y(ia,ja)-y(i,j)
                           sxn = sx / dx
                           syn = sy / dy
                           sn = sqrt(sxn * sxn + syn * syn)
                           if (sn .le. 1._r8) then
                              n = n + 1
                              d(n, 1) = sx
                              d(n, 2) = sy
                              z(n) = v(ia,ja)
                              wn = 1._r8 - sn * sn * sn
                              w(n) = wn * wn * wn
                           end if
                        end if
                     end do
                  end do
                  if (n .ge. reg_size_min) then
                     call regsm2(n,z,d,w,vp(i,j),v_miss,tol,l_size)
                  else
                     vp(i,j) = v_miss
                  end if
                  if (vp(i,j) .eq. v_miss) then
                     nmiss = nmiss + 1
                     select case (opt)
                     case (tamu_filter_opt_miss_zonal_avg)
                         ! --- Replace missing value with zonal average ---
                         n = count(v(:,j) .ne. v_mask)
                         if (n .gt. 0) then
                            vp(i,j) = &
                            sum(v(:,j), mask = (v(:,j).ne.v_mask))/n
                         end if
                     case (tamu_filter_opt_miss_local_avg)
                         ! --- Replace missing value with local average ---
                         if (n .gt. 0) vp(i,j) = sum(z(1:n))/n
                     case (tamu_filter_opt_miss_original)
                         ! --- Replace missing value with original value ---
                         vp(i,j) = v(i,j)
                     case default
                         ! --- Replace missing value with original value ---
                         vp(i,j) = v(i,j)
                     end select
                  end if
               else
                  vp(i,j) = v_mask
               end if
            end do
         end do

         v = vp

         rc = 0

         deallocate(d)
         deallocate(vp)
         deallocate(z)
         deallocate(w)

      end subroutine tamu_filter_loess_xy

      ! === Loess auxiliary subroutines: begin ===
      subroutine regsm2(n, z, des, w, zhat, zmiss, ztol, lddes)

         integer,                      intent(in)    :: n
         integer,                      intent(in)    :: lddes
         real(r8), dimension(lddes),   intent(in)    :: z
         real(r8), dimension(lddes,2), intent(in)    :: des
         real(r8), dimension(lddes),   intent(inout) :: w
         real(r8),                     intent(in)    :: zmiss
         real(r8),                     intent(in)    :: ztol
         real(r8),                     intent(out)   :: zhat

         integer,  parameter :: np = 6
         real(r8), parameter :: z_large = 1.e+35_r8

         real(r8) :: d
         integer,  dimension(np)    :: indx
         real(r8), dimension(np)    :: a, t
         real(r8), dimension(np,np) :: xtx

         integer  :: i, ierr, j, k
         real(r8) :: zi, zmax, zmin

         a    = 0._r8
         xtx  = 0._r8
         indx = 0

         zmax = -z_large
         zmin =  z_large

         do i = 1, n
            zmax = max(zmax, z(i))
            zmin = min(zmin, z(i))

            zi   = w(i) * z(i)
            t(1) = w(i)
            t(2) = des(i,1) * t(1)
            t(3) = des(i,1) * t(2)
            t(4) = des(i,2) * t(1)
            t(5) = des(i,2) * t(4)
            t(6) = des(i,2) * t(2)

            do k = 1, np
               a(k) = a(k) + t(k) * zi
               do j = 1, np
                  xtx(j,k) = xtx(j,k) + t(j) * t(k)
               end do
            end do
         end do


         ! --- solve

         call ludcmp(xtx, np, np, indx, d, w, ierr)
         if (ierr .ne. 0) then
            zhat = zmiss
            return
         end if

         call lubksb(xtx, np, np, indx, a)
         zhat = a(1)
         if ((zhat + ztol .gt. zmax) .or. (zhat - ztol .lt. zmin)) then
            zhat = zmiss
         end if

      end subroutine regsm2

      SUBROUTINE LUDCMP(A,N,NP,INDX,D,VV,IER)

         IMPLICIT NONE

         INTEGER,                    INTENT(IN)    :: N, NP
         INTEGER,  DIMENSION(N),     INTENT(OUT)   :: INDX
         REAL(r8), DIMENSION(NP,NP), INTENT(INOUT) :: A
         REAL(r8),                   INTENT(OUT)   :: D
         REAL(r8), DIMENSION(N),     INTENT(INOUT) :: VV
         INTEGER,                    INTENT(OUT)   :: IER

         REAL(r8), PARAMETER :: TINY = 1.0E-20_r8

         INTEGER  :: I, IMAX, J, K
         REAL(r8) :: AAMAX, DUM, SUM

         VV = 0._r8

         IER = 0
         D = 1._r8
         DO I = 1, N
            AAMAX = 0._r8
            DO J = 1, N
               IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
            END DO
            IF (AAMAX.EQ.0._r8) THEN
               IER = 1
               RETURN
            END IF
            VV(I) = 1._r8 / AAMAX
         END DO

         DO J = 1, N
            IF (J.GT.1) THEN
               DO I = 1, J-1
                  SUM = A(I,J)
                  IF (I.GT.1)THEN
                     DO K = 1, I-1
                        SUM = SUM - A(I,K) * A(K,J)
                     END DO
                     A(I,J) = SUM
                  END IF
               END DO
            ENDIF

            AAMAX = 0._r8

            DO I = J, N
               SUM = A(I,J)
               IF (J.GT.1)THEN
                  DO K = 1, J-1
                     SUM = SUM - A(I,K) * A(K,J)
                  END DO
                  A(I,J) = SUM
               END IF
               DUM = VV(I) * ABS(SUM)
               IF (DUM.GE.AAMAX) THEN
                  IMAX  = I
                  AAMAX = DUM
               END IF
            END DO

            IF (J.NE.IMAX)THEN
               DO K = 1, N
                  DUM = A(IMAX,K)
                  A(IMAX,K) = A(J,K)
                  A(J,K) = DUM
               END DO
               D = -D
               VV(IMAX) = VV(J)
            ENDIF

            INDX(J) = IMAX
            IF (J.NE.N) THEN
               IF (A(J,J).EQ.0._r8) A(J,J)=TINY
               DUM = 1._r8 / A(J,J)
               DO I = J+1, N
                  A(I,J) = A(I,J) * DUM
               END DO
            END IF
         END DO

         IF (A(N,N).EQ.0._r8) A(N,N) = TINY

      END SUBROUTINE LUDCMP

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

         IMPLICIT NONE

         INTEGER,                    INTENT(IN)    :: N, NP
         INTEGER,  DIMENSION(N),     INTENT(IN)    :: INDX
         REAL(r8), DIMENSION(NP,NP), INTENT(IN)    :: A
         REAL(r8), DIMENSION(N),     INTENT(INOUT) :: B

         INTEGER  :: I, II, J, LL
         REAL(r8) :: SUM

         II=0

         DO I = 1, N
            LL    = INDX(I)
            SUM   = B(LL)
            B(LL) = B(I)
            IF (II.NE.0) THEN
               DO J = II, I-1
                  SUM = SUM - A(I,J) * B(J)
               END DO
            ELSE IF (SUM.NE.0._r8) THEN
               II = I
            END IF
            B(I) = SUM
         END DO

         DO I = N, 1, -1
            SUM = B(I)
            IF (I.LT.N) THEN
               DO J = I+1, N
                  SUM = SUM - A(I,J) * B(J)
               END DO
            END IF
            B(I) = SUM / A(I,I)
         END DO

      END SUBROUTINE LUBKSB
      ! === Loess auxiliary subroutines: end ===

#endif
end module mod_tamu_filter
