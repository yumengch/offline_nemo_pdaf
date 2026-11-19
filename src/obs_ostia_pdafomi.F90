!> PDAF-OMI observations from daily OSTIA foundation SST analyses.
!!
!! Valid ocean pixels are binned into the nearest AMM7 surface cell.
!! Each cell contributes one area-weighted SST and one error estimate.
!! Pixel errors are treated as fully correlated within a model cell:
!! averaging must not reduce the error by sqrt(number of pixels).
module obs_ostia_pdafomi
   use mod_kind_pdaf
   use PDAF, only: obs_f, obs_l
   use obs_grid_pdaf
   use, intrinsic :: ieee_arithmetic
   implicit none
   private
   public :: read_config_ostia, print_ostia_configuration
   public :: init_dim_obs_ostia, obs_op_ostia, init_dim_obs_l_ostia

   ! *** Observation configuration ***
   logical :: assim_ostia = .false.
   integer :: id_temp = 1 ! Index of thetao in sfields.
   integer :: min_pixels = 1 ! Minimum valid pixels in a super-observation.
   character(lc) :: path_ostia = 'obsproc'
   character(lc) :: suffix_ostia = &
      '120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB_REP-v02.0-fv02.0.nc'
   real(pwp) :: error_floor = 0.5_pwp ! degC; tune for representation error.
   real(pwp) :: max_ice_fraction = 0.0_pwp ! Exclude ice by default.
   real(pwp) :: lradius = 20000.0_pwp ! Horizontal cut-off in metres.
   real(pwp) :: sradius = 20000.0_pwp ! Horizontal support in metres.
   namelist /ostia_nml/ assim_ostia, id_temp, path_ostia, suffix_ostia, &
      min_pixels, error_floor, max_ice_fraction, lradius, sradius

   type(obs_f), target :: thisobs
   type(obs_l), target :: thisobs_l
!$OMP THREADPRIVATE(thisobs_l)

contains
   !> Read an optional namelist; absence retains disabled defaults.
   subroutine read_config_ostia(unit)
      integer, intent(in) :: unit
      integer :: status
      if (.not. namelist_present(unit, 'ostia_nml')) return
      read (unit, nml=ostia_nml, iostat=status)
      call obs_require(status == 0, 'Cannot read ostia_nml')
      call obs_require(min_pixels >= 1, 'OSTIA min_pixels must be positive')
      call obs_require(ieee_is_finite(error_floor) .and. error_floor > 0, &
                       'Invalid OSTIA error floor')
      call obs_require(max_ice_fraction >= 0 .and. max_ice_fraction <= 1, &
                       'Invalid OSTIA maximum ice fraction')
      call obs_require(lradius > 0 .and. sradius > 0, &
                       'OSTIA localization radii must be positive')
   end subroutine read_config_ostia

   !> Print super-observation and uncertainty settings for the cycle log.
   subroutine print_ostia_configuration()
      write (*, '(a,l3)') 'NEMO-PDAF: OSTIA SST ', assim_ostia
      if (.not. assim_ostia) return
      write (*, '(a)') 'NEMO-PDAF: OSTIA path '//trim(path_ostia)
      write (*, '(a)') 'NEMO-PDAF: OSTIA suffix '//trim(suffix_ostia)
      write (*, '(a,2i5)') 'NEMO-PDAF: OSTIA field/min_pixels ', &
                           id_temp, min_pixels
      write (*, '(a,2f10.4)') 'NEMO-PDAF: OSTIA error floor/ice limit ', &
                              error_floor, max_ice_fraction
      write (*, '(a,2f12.1)') 'NEMO-PDAF: OSTIA radii (m) ', lradius, sradius
   end subroutine print_ostia_configuration

   !> Read daily SST and pack one super-observation per wet surface cell.
   subroutine init_dim_obs_ostia(step, dim_obs)
      use netcdf
      use nemo_pdaf, only: year, month, day, nwet, ni_p, nj_p, &
                          nav_lon, nav_lat, deg2rad
      integer, intent(in) :: step
      integer, intent(out) :: dim_obs
      character(lc) :: filename
      character(8) :: date
      real(pwp), allocatable :: lon(:), lat(:), sst(:,:), error(:,:)
      real(pwp), allocatable :: mask(:,:), ice(:,:)
      real(pwp), allocatable :: sums(:,:), errors(:,:), weights(:,:)
      real(pwp), allocatable :: values(:), ivar(:), coords(:,:)
      integer, allocatable :: counts(:,:), indices(:)
      integer :: ncid, first(2), extent(2), i, j, im, jm, id, n, capacity
      real(pwp) :: area, value, sd
      dim_obs = 0
      thisobs%doassim = 0
      if (.not. assim_ostia) return
      call validate_field(id_temp, .true.)
      write (date, '(i4.4,i2.2,i2.2)') year, month, day
      filename = trim(path_ostia)//'/'//date//trim(suffix_ostia)
      write (*, '(a)') 'NEMO-PDAF: OSTIA '//trim(filename)
      capacity = max(1, nwet)
      allocate(values(capacity), ivar(capacity), coords(2,capacity))
      allocate(indices(capacity))
      values = 0
      ivar = 0
      coords = 0
      indices = 0
      n = 0
      if (nwet > 0) then
         call obs_check(nf90_open(filename, nf90_nowrite, ncid), filename)
         call check_date(ncid, .false.)
         call check_units(ncid, 'analysed_sst', 'kelvin')
         call check_units(ncid, 'analysis_error', 'kelvin')
         call check_units(ncid, 'sea_ice_fraction', '1')
         call read_window(ncid, lon, lat, first, extent)
         if (all(extent > 0)) then
            call read_plane(ncid, 'analysed_sst', first, extent, sst)
            call read_plane(ncid, 'analysis_error', first, extent, error)
            call read_plane(ncid, 'mask', first, extent, mask)
            call read_plane(ncid, 'sea_ice_fraction', first, extent, ice)
         end if
         call obs_check(nf90_close(ncid), 'close OSTIA')
      end if
      if (allocated(sst)) then
         allocate(sums(ni_p,nj_p), errors(ni_p,nj_p))
         allocate(weights(ni_p,nj_p), counts(ni_p,nj_p))
         sums = 0
         errors = 0
         weights = 0
         counts = 0
         do j = 1, size(lat)
            jm = nearest_axis(nav_lat, lat(j))
            if (jm == 0) cycle
            ! OSTIA has uniform angular spacing; area is proportional to cos.
            area = cos(lat(j)*deg2rad)
            do i = 1, size(lon)
               im = nearest_axis(nav_lon, lon(i))
               if (im == 0) cycle
               id = state_index(id_temp, im, jm, 1)
               if (id == 0) cycle
               if (.not. ieee_is_finite(mask(i,j))) cycle
               ! GHRSST flags: accept water only, excluding land/lake/ice.
               if (nint(mask(i,j)) /= 1) cycle
               if (.not. ieee_is_finite(ice(i,j))) cycle
               if (ice(i,j) < 0 .or. ice(i,j) > max_ice_fraction) cycle
               value = sst(i,j) - 273.15_pwp
               sd = error(i,j)
               if (.not. ieee_is_finite(value)) cycle
               if (.not. ieee_is_finite(sd)) cycle
               if (value < -3 .or. value > 45 .or. sd <= 0) cycle
               sums(im,jm) = sums(im,jm) + area*value
               errors(im,jm) = errors(im,jm) + area*sd
               weights(im,jm) = weights(im,jm) + area
               counts(im,jm) = counts(im,jm) + 1
            end do
         end do
         do j = 1, nj_p
            do i = 1, ni_p
               if (counts(i,j) < min_pixels) cycle
               if (weights(i,j) <= 0) cycle
               n = n + 1
               values(n) = sums(i,j)/weights(i,j)
               sd = max(error_floor, errors(i,j)/weights(i,j))
               ivar(n) = 1.0_pwp/(sd*sd)
               coords(:,n) = [real(nav_lon(i),pwp), &
                              real(nav_lat(j),pwp)] * deg2rad
               indices(n) = state_index(id_temp, i, j, 1)
            end do
         end do
      end if
      call finish_observations(thisobs, n, values, ivar, coords, &
                               indices, lradius, dim_obs)
   end subroutine init_dim_obs_ostia

   !> Compare super-observations with model surface potential temperature.
   subroutine obs_op_ostia(dim_p, dim_obs, state_p, ostate)
      use PDAF, only: PDAFomi_obs_op_gridpoint
      integer, intent(in) :: dim_p, dim_obs
      real(pwp), intent(in) :: state_p(dim_p)
      real(pwp), intent(inout) :: ostate(dim_obs)
      if (thisobs%doassim /= 1) return
      call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
   end subroutine obs_op_ostia

   !> Initialize horizontal localization around the analysis column.
   subroutine init_dim_obs_l_ostia(domain_p, step, dim_obs, dim_obs_l)
      use PDAF, only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, only: domain_coords
      use config_pdaf, only: locweight
      integer, intent(in) :: domain_p, step, dim_obs
      integer, intent(inout) :: dim_obs_l
      if (thisobs%doassim /= 1) return
      call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
           locweight, lradius, sradius, dim_obs_l)
   end subroutine init_dim_obs_l_ostia
end module obs_ostia_pdafomi
