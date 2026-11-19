!> PDAF-OMI observations from EN4.2.2 monthly gridded analyses.
!!
!! These are objective analyses, not the EN4 profile product. Temperature
!! and salinity have separate OMI types, in that order. Native EN4 points
!! are mapped to the nearest wet AMM7 cell and level, without extrapolation
!! below the deepest wet model centre. Localization is horizontal only.
module obs_en4_pdafomi
   use mod_kind_pdaf
   use PDAF, only: obs_f, obs_l
   use obs_grid_pdaf
   use, intrinsic :: ieee_arithmetic
   implicit none
   private
   public :: read_config_en4, print_en4_configuration
   public :: init_dim_obs_en4, obs_op_en4, init_dim_obs_l_en4

   ! *** Observation configuration ***
   logical :: assim_en4_temp = .false.
   logical :: assim_en4_salt = .false.
   integer :: assim_day = 0 ! Zero: month end; otherwise calendar day.
   integer :: id_temp = 1, id_salt = 2 ! Indices in sfields.
   character(lc) :: path_en4 = 'obsproc'
   character(lc) :: prefix_en4 = 'EN.4.2.2.f.analysis.g10.'
   real(pwp) :: min_weight = 0.5_pwp ! Minimum observation influence.
   real(pwp) :: temp_error_floor = 0.5_pwp ! degC; tune for experiment.
   real(pwp) :: salt_error_floor = 0.1_pwp ! psu; tune for experiment.
   real(pwp) :: lradius = 20000.0_pwp ! Horizontal cut-off in metres.
   real(pwp) :: sradius = 20000.0_pwp ! Horizontal support in metres.
   namelist /en4_nml/ assim_en4_temp, assim_en4_salt, id_temp, id_salt, &
      path_en4, prefix_en4, assim_day, min_weight, temp_error_floor, &
      salt_error_floor, lradius, sradius

   type(obs_f), target :: thisobs(2)
   type(obs_l), target :: thisobs_l(2)
!$OMP THREADPRIVATE(thisobs_l)

contains
   !> Read an optional namelist; absence retains disabled defaults.
   subroutine read_config_en4(unit)
      integer, intent(in) :: unit
      integer :: status
      if (.not. namelist_present(unit, 'en4_nml')) return
      read (unit, nml=en4_nml, iostat=status)
      call obs_require(status == 0, 'Cannot read en4_nml')
      call obs_require(assim_day >= 0 .and. assim_day <= 28, &
                       'EN4 assim_day must be zero or between 1 and 28')
      call obs_require(min_weight >= 0 .and. min_weight <= 1, &
                       'EN4 min_weight must be between zero and one')
      call obs_require(ieee_is_finite(temp_error_floor) .and. &
                       temp_error_floor > 0, 'Invalid EN4 temperature floor')
      call obs_require(ieee_is_finite(salt_error_floor) .and. &
                       salt_error_floor > 0, 'Invalid EN4 salinity floor')
      call obs_require(lradius > 0 .and. sradius > 0, &
                       'EN4 localization radii must be positive')
   end subroutine read_config_en4

   !> Print product selection and error assumptions for reproducibility.
   subroutine print_en4_configuration()
      write (*, '(a,2l3)') 'NEMO-PDAF: EN4 temperature/salinity ', &
                           assim_en4_temp, assim_en4_salt
      if (.not. (assim_en4_temp .or. assim_en4_salt)) return
      write (*, '(a,i3)') 'NEMO-PDAF: EN4 day (0=month end) ', assim_day
      write (*, '(a)') 'NEMO-PDAF: EN4 path '//trim(path_en4)
      write (*, '(a)') 'NEMO-PDAF: EN4 prefix '//trim(prefix_en4)
      write (*, '(a,2i5)') 'NEMO-PDAF: EN4 field indices ', id_temp, id_salt
      write (*, '(a,3f10.4)') 'NEMO-PDAF: EN4 weight and error floors ', &
                              min_weight, temp_error_floor, salt_error_floor
      write (*, '(a,2f12.1)') 'NEMO-PDAF: EN4 radii (m) ', lradius, sradius
   end subroutine print_en4_configuration

   !> Initialize enabled types collectively on every filter rank.
   subroutine init_dim_obs_en4(step, dim_obs)
      use nemo_pdaf, only: year, month, day
      integer, intent(in) :: step
      integer, intent(out) :: dim_obs
      integer :: dimension, selected_day
      character(lc) :: filename
      character(6) :: date
      thisobs(:)%doassim = 0
      dim_obs = 0
      if (.not. (assim_en4_temp .or. assim_en4_salt)) return
      selected_day = assim_day
      if (selected_day == 0) selected_day = days_in_month(year, month)
      if (day /= selected_day) return
      write (date, '(i4.4,i2.2)') year, month
      filename = trim(path_en4)//'/'//trim(prefix_en4)//date//'.nc'
      write (*, '(a)') 'NEMO-PDAF: EN4 '//trim(filename)
      if (assim_en4_temp) then
         call init_en4_field(filename, 'temperature', id_temp, &
              temp_error_floor, thisobs(1), dimension)
         dim_obs = dim_obs + dimension
      end if
      if (assim_en4_salt) then
         call init_en4_field(filename, 'salinity', id_salt, &
              salt_error_floor, thisobs(2), dimension)
         dim_obs = dim_obs + dimension
      end if
   end subroutine init_dim_obs_en4

   !> Read a regional slab, apply QC, and pack one EN4 observation type.
   subroutine init_en4_field(filename, variable, id_field, floor_error, &
                             observation, dim_obs)
      use netcdf
      use nemo_pdaf, only: nwet, nav_lon, nav_lat, nav_lev, deg2rad, tmask
      character(*), intent(in) :: filename, variable
      integer, intent(in) :: id_field
      real(pwp), intent(in) :: floor_error
      type(obs_f), intent(inout) :: observation
      integer, intent(out) :: dim_obs
      real(pwp), allocatable :: lon(:), lat(:), depth(:)
      real(pwp), allocatable :: data(:,:,:), error(:,:,:), weight(:,:,:)
      real(pwp), allocatable :: values(:), ivar(:), coords(:,:)
      integer, allocatable :: indices(:)
      integer :: ncid, first(2), extent(2), capacity, n
      integer :: i, j, k, im, jm, km, id, deepest
      real(pwp) :: value, sd, influence
      logical :: temperature
      temperature = variable == 'temperature'
      call validate_field(id_field, temperature)
      n = 0
      capacity = 1
      if (nwet > 0) then
         call obs_check(nf90_open(filename, nf90_nowrite, ncid), filename)
         call check_date(ncid, .true.)
         if (temperature) then
            call check_units(ncid, variable, 'kelvin')
            call check_units(ncid, variable//'_uncertainty', 'kelvin')
         else
            call check_units(ncid, variable, '1')
            call check_units(ncid, variable//'_uncertainty', '1')
         end if
         call read_window(ncid, lon, lat, first, extent)
         call read_axis(ncid, 'depth', depth)
         call check_units(ncid, 'depth', 'metres')
         call obs_require(all(depth >= 0) .and. &
                          all(depth(2:) > depth(:size(depth)-1)), &
                          'EN4 depths must be increasing metres')
         if (all(extent > 0)) then
            call read_cube(ncid, variable, first, extent, size(depth), data)
            call read_cube(ncid, variable//'_uncertainty', &
                           first, extent, size(depth), error)
            call read_cube(ncid, variable//'_observation_weights', &
                           first, extent, size(depth), weight)
            capacity = max(1, size(data))
         end if
         call obs_check(nf90_close(ncid), 'close EN4')
      end if
      allocate(values(capacity), ivar(capacity), coords(2,capacity))
      allocate(indices(capacity))
      values = 0
      ivar = 0
      coords = 0
      indices = 0
      if (allocated(data)) then
         do j = 1, size(lat)
            jm = nearest_axis(nav_lat, lat(j))
            if (jm == 0) cycle
            do i = 1, size(lon)
               im = nearest_axis(nav_lon, lon(i))
               if (im == 0) cycle
               deepest = count(tmask(im,jm,:) == 1.0_pwp)
               if (deepest == 0) cycle
               do k = 1, size(depth)
                  if (depth(k) > real(nav_lev(deepest), pwp)) cycle
                  km = minloc(abs(real(nav_lev, pwp)-depth(k)), dim=1)
                  id = state_index(id_field, im, jm, km)
                  if (id == 0) cycle
                  value = data(i,j,k)
                  sd = error(i,j,k)
                  influence = weight(i,j,k)
                  if (.not. ieee_is_finite(value)) cycle
                  if (.not. ieee_is_finite(sd)) cycle
                  if (.not. ieee_is_finite(influence)) cycle
                  if (influence < min_weight .or. influence > 1.01_pwp) &
                     cycle
                  if (sd <= 0) cycle
                  if (temperature) value = value - 273.15_pwp
                  ! Reject implausible values; never clip observations.
                  if (temperature) then
                     if (value < -5 .or. value > 45) cycle
                  else
                     if (value < 0 .or. value > 48) cycle
                  end if
                  n = n + 1
                  values(n) = value
                  ivar(n) = 1.0_pwp / max(sd, floor_error)**2
                  coords(:,n) = [lon(i), lat(j)] * deg2rad
                  indices(n) = id
               end do
            end do
         end do
      end if
      call finish_observations(observation, n, values, ivar, coords, &
                               indices, lradius, dim_obs)
   end subroutine init_en4_field

   !> Apply nearest-model-grid-point operators in the initialization order.
   subroutine obs_op_en4(dim_p, dim_obs, state_p, ostate)
      use PDAF, only: PDAFomi_obs_op_gridpoint
      integer, intent(in) :: dim_p, dim_obs
      real(pwp), intent(in) :: state_p(dim_p)
      real(pwp), intent(inout) :: ostate(dim_obs)
      integer :: i
      do i = 1, 2
         if (thisobs(i)%doassim /= 1) cycle
         call PDAFomi_obs_op_gridpoint(thisobs(i), state_p, ostate)
      end do
   end subroutine obs_op_en4

   !> Initialize horizontal localization around the model analysis column.
   subroutine init_dim_obs_l_en4(domain_p, step, dim_obs, dim_obs_l)
      use PDAF, only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, only: domain_coords
      use config_pdaf, only: locweight
      integer, intent(in) :: domain_p, step, dim_obs
      integer, intent(inout) :: dim_obs_l
      integer :: i
      do i = 1, 2
         if (thisobs(i)%doassim /= 1) cycle
         call PDAFomi_init_dim_obs_l(thisobs_l(i), thisobs(i), &
              domain_coords, locweight, lradius, sradius, dim_obs_l)
      end do
   end subroutine init_dim_obs_l_en4
end module obs_en4_pdafomi
