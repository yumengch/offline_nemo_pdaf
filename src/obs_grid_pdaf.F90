!> Shared NetCDF and grid utilities for EN4 and OSTIA observations.
!!
!! The offline model exposes separable longitude/latitude axes. These
!! routines deliberately reject unsupported grids and transformed fields.
!! NetCDF Fortran reads packed values: decode before converting units.
module obs_grid_pdaf
   use mod_kind_pdaf
   use netcdf
   use, intrinsic :: ieee_arithmetic
   use parallel_pdaf, only: abort_parallel
   implicit none
   private
   public :: obs_check, obs_require, read_axis, read_cube, read_plane
   public :: read_window, check_units, check_date, state_index
   public :: nearest_axis, validate_field, finish_observations
   public :: days_in_month, namelist_present

contains
   !> Report NetCDF failures and stop all ranks, avoiding MPI deadlocks.
   subroutine obs_check(status, context)
      integer, intent(in) :: status
      character(*), intent(in) :: context
      if (status == nf90_noerr) return
      write (*, '(a)') 'NEMO-PDAF: '//trim(context)//': '// &
                      trim(nf90_strerror(status))
      call abort_parallel()
   end subroutine obs_check

   !> Fail clearly when an input does not satisfy the file contract.
   subroutine obs_require(condition, message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (condition) return
      write (*, '(a)') 'NEMO-PDAF: '//trim(message)
      call abort_parallel()
   end subroutine obs_require

   !> Distinguish a missing optional group from a malformed namelist.
   logical function namelist_present(unit, name) result(found)
      integer, intent(in) :: unit
      character(*), intent(in) :: name
      character(lc) :: line
      integer :: status, i, letter, length
      found = .false.
      rewind(unit)
      do
         read (unit, '(a)', iostat=status) line
         if (status < 0) exit
         call obs_require(status == 0, 'Cannot scan observation namelists')
         line = adjustl(line)
         do i = 1, len_trim(line)
            letter = iachar(line(i:i))
            if (letter >= iachar('A') .and. letter <= iachar('Z')) &
               line(i:i) = achar(letter + iachar('a') - iachar('A'))
         end do
         length = len_trim(name) + 1
         if (line(:length) /= '&'//trim(name)) cycle
         if (index(' ,!'//achar(9), line(length+1:length+1)) == 0) cycle
         found = .true.
         exit
      end do
      rewind(unit)
   end function namelist_present

   !> Read a one-dimensional coordinate axis.
   subroutine read_axis(ncid, name, axis)
      integer, intent(in) :: ncid
      character(*), intent(in) :: name
      real(pwp), allocatable, intent(out) :: axis(:)
      integer :: varid, ndims, dims(nf90_max_var_dims), n
      call obs_check(nf90_inq_varid(ncid, name, varid), name)
      call obs_check(nf90_inquire_variable(ncid, varid, ndims=ndims, &
                     dimids=dims), name)
      call obs_require(ndims == 1, name//' must be one-dimensional')
      call obs_check(nf90_inquire_dimension(ncid, dims(1), len=n), name)
      allocate(axis(n))
      call obs_check(nf90_get_var(ncid, varid, axis), name)
      call obs_require(all(ieee_is_finite(axis)), name//' is not finite')
   end subroutine read_axis

   !> Read only the local regional slab, including half-cell boundaries.
   !! Longitude wrapping can make the EN4 slab span its full 360 columns.
   subroutine read_window(ncid, lon, lat, first, extent)
      use nemo_pdaf, only: nav_lon, nav_lat
      integer, intent(in) :: ncid
      real(pwp), allocatable, intent(out) :: lon(:), lat(:)
      integer, intent(out) :: first(2), extent(2)
      real(pwp), allocatable :: x(:), y(:)
      real(pwp) :: dx, dy
      integer :: i, xfirst, xlast, yfirst, ylast
      call read_axis(ncid, 'lon', x)
      call read_axis(ncid, 'lat', y)
      call obs_require(size(x) > 1 .and. size(y) > 1, &
                       'Observation grid requires at least 2x2 cells')
      call obs_require(all(x(2:) > x(:size(x)-1)) .and. &
                       all(y(2:) > y(:size(y)-1)), &
                       'Observation coordinates must increase')
      x = modulo(x + 180.0_pwp, 360.0_pwp) - 180.0_pwp
      dx = abs(real(nav_lon(2) - nav_lon(1), pwp)) / 2
      dy = abs(real(nav_lat(2) - nav_lat(1), pwp)) / 2
      xfirst = size(x) + 1
      xlast = 0
      yfirst = size(y) + 1
      ylast = 0
      do i = 1, size(x)
         if (x(i) < minval(nav_lon)-dx) cycle
         if (x(i) >= maxval(nav_lon)+dx) cycle
         xfirst = min(xfirst, i)
         xlast = max(xlast, i)
      end do
      do i = 1, size(y)
         if (y(i) < minval(nav_lat)-dy) cycle
         if (y(i) >= maxval(nav_lat)+dy) cycle
         yfirst = min(yfirst, i)
         ylast = max(ylast, i)
      end do
      extent = [max(0, xlast-xfirst+1), max(0, ylast-yfirst+1)]
      first = [xfirst, yfirst]
      allocate(lon(extent(1)), lat(extent(2)))
      if (all(extent > 0)) then
         lon = x(xfirst:xlast)
         lat = y(yfirst:ylast)
      end if
   end subroutine read_window

   !> Read a (longitude, latitude, depth, time) variable at time index 1.
   subroutine read_cube(ncid, name, first, extent, nz, field)
      integer, intent(in) :: ncid, first(2), extent(2), nz
      character(*), intent(in) :: name
      real(pwp), allocatable, intent(out) :: field(:,:,:)
      real(pwp), allocatable :: raw(:,:,:,:)
      integer :: varid
      call field_id(ncid, name, ['lon  ', 'lat  ', 'depth', 'time '], &
                    varid)
      allocate(raw(extent(1), extent(2), nz, 1))
      call obs_check(nf90_get_var(ncid, varid, raw, &
                     start=[first, 1, 1], count=[extent, nz, 1]), name)
      allocate(field(extent(1), extent(2), nz))
      field = raw(:,:,:,1)
      call decode_field(ncid, varid, field)
   end subroutine read_cube

   !> Read a (longitude, latitude, time) variable at time index 1.
   subroutine read_plane(ncid, name, first, extent, field)
      integer, intent(in) :: ncid, first(2), extent(2)
      character(*), intent(in) :: name
      real(pwp), allocatable, intent(out) :: field(:,:)
      real(pwp), allocatable :: raw(:,:,:)
      integer :: varid
      call field_id(ncid, name, ['lon ', 'lat ', 'time'], varid)
      allocate(raw(extent(1), extent(2), 1))
      call obs_check(nf90_get_var(ncid, varid, raw, &
                     start=[first, 1], count=[extent, 1]), name)
      call decode_field(ncid, varid, raw)
      allocate(field(extent(1), extent(2)))
      field = raw(:,:,1)
   end subroutine read_plane

   !> Require the product's expected dimension order and one time record.
   subroutine field_id(ncid, name, expected, varid)
      integer, intent(in) :: ncid
      character(*), intent(in) :: name, expected(:)
      integer, intent(out) :: varid
      integer :: ndims, dims(nf90_max_var_dims), i, n
      character(nf90_max_name) :: dimname
      call obs_check(nf90_inq_varid(ncid, name, varid), name)
      call obs_check(nf90_inquire_variable(ncid, varid, ndims=ndims, &
                     dimids=dims), name)
      call obs_require(ndims == size(expected), name//' rank mismatch')
      do i = 1, ndims
         call obs_check(nf90_inquire_dimension(ncid, dims(i), &
                        name=dimname, len=n), name)
         call obs_require(trim(dimname) == trim(expected(i)), &
                          name//' dimension order mismatch')
         if (i == ndims) call obs_require(n == 1, &
                                          name//' requires one time record')
      end do
   end subroutine field_id

   !> Read optional scalar attributes without hiding other NetCDF errors.
   subroutine read_attribute(ncid, varid, name, value, present)
      integer, intent(in) :: ncid, varid
      character(*), intent(in) :: name
      real(pwp), intent(inout) :: value
      logical, intent(out) :: present
      integer :: status
      real(pwp) :: attribute
      status = nf90_get_att(ncid, varid, name, attribute)
      present = status /= nf90_enotatt
      if (present) then
         call obs_check(status, 'attribute '//name)
         value = attribute
      end if
   end subroutine read_attribute

   !> Mask packed fill/range values before applying scale and offset.
   subroutine decode_field(ncid, varid, field)
      integer, intent(in) :: ncid, varid
      real(pwp), intent(inout) :: field(:,:,:)
      real(pwp) :: scale, offset, fill, missing, lower, upper, nan
      real(pwp) :: range(2), value
      logical :: has_fill, has_missing, present
      integer :: i, j, k, status
      nan = ieee_value(0.0_pwp, ieee_quiet_nan)
      fill = 0
      missing = 0
      lower = -huge(0.0_pwp)
      upper = huge(0.0_pwp)
      call read_attribute(ncid, varid, '_FillValue', fill, has_fill)
      call read_attribute(ncid, varid, 'missing_value', missing, &
                          has_missing)
      call read_attribute(ncid, varid, 'valid_min', lower, present)
      call read_attribute(ncid, varid, 'valid_max', upper, present)
      status = nf90_get_att(ncid, varid, 'valid_range', range)
      if (status /= nf90_enotatt) then
         call obs_check(status, 'valid_range')
         lower = range(1)
         upper = range(2)
      end if
      scale = 1
      offset = 0
      call read_attribute(ncid, varid, 'scale_factor', scale, present)
      call read_attribute(ncid, varid, 'add_offset', offset, present)
      call obs_require(ieee_is_finite(scale) .and. &
                       ieee_is_finite(offset), 'Invalid packing attributes')
      do k = 1, size(field,3)
         do j = 1, size(field,2)
            do i = 1, size(field,1)
               value = field(i,j,k)
               field(i,j,k) = nan
               if (.not. ieee_is_finite(value)) cycle
               if (has_fill .and. value == fill) cycle
               if (has_missing .and. value == missing) cycle
               if (value < lower .or. value > upper) cycle
               field(i,j,k) = value*scale + offset
            end do
         end do
      end do
   end subroutine decode_field

   !> Confirm units before applying the product-specific conversion.
   subroutine check_units(ncid, name, expected)
      integer, intent(in) :: ncid
      character(*), intent(in) :: name, expected
      integer :: varid
      character(lc) :: units
      call obs_check(nf90_inq_varid(ncid, name, varid), name)
      call obs_check(nf90_get_att(ncid, varid, 'units', units), name)
      call obs_require(trim(units) == expected, name//' units mismatch')
   end subroutine check_units

   !> Convert a Gregorian date to an integer day number.
   integer function day_number(year, month, day) result(number)
      integer, intent(in) :: year, month, day
      integer :: a, y, m
      a = (14-month)/12
      y = year + 4800 - a
      m = month + 12*a - 3
      number = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400
   end function day_number

   !> Number of days in a Gregorian calendar month.
   integer function days_in_month(year, month) result(days)
      integer, intent(in) :: year, month
      integer :: next_year, next_month
      next_year = year
      next_month = month + 1
      if (month == 12) then
         next_year = year + 1
         next_month = 1
      end if
      days = day_number(next_year, next_month, 1) - &
             day_number(year, month, 1)
   end function days_in_month

   !> Check the file time against the selected monthly or daily cycle.
   !! Only the Gregorian CF epochs used by these two products are accepted.
   subroutine check_date(ncid, monthly)
      use nemo_pdaf, only: year, month, day
      integer, intent(in) :: ncid
      logical, intent(in) :: monthly
      integer :: varid, iy, im, id, status, origin, target, ending
      real(pwp) :: time(1), factor, lower, upper
      character(lc) :: units, calendar
      call obs_check(nf90_inq_varid(ncid, 'time', varid), 'time')
      call obs_check(nf90_get_att(ncid, varid, 'units', units), 'time')
      calendar = 'standard'
      status = nf90_get_att(ncid, varid, 'calendar', calendar)
      if (status == nf90_enotatt) then
         calendar = 'standard'
      else
         call obs_check(status, 'calendar')
      end if
      call obs_require(trim(calendar) == 'standard' .or. &
                       trim(calendar) == 'gregorian' .or. &
                       trim(calendar) == 'proleptic_gregorian', &
                       'Unsupported observation calendar')
      if (index(units, 'days since ') == 1) then
         factor = 1
         units = units(12:)
      else if (index(units, 'seconds since ') == 1) then
         factor = 86400
         units = units(15:)
      else
         call obs_require(.false., 'Unsupported observation time units')
         return
      end if
      read (units(1:4), *, iostat=status) iy
      call obs_require(status == 0, 'Invalid observation time epoch')
      read (units(6:7), *, iostat=status) im
      call obs_require(status == 0, 'Invalid observation time epoch')
      read (units(9:10), *, iostat=status) id
      call obs_require(status == 0, 'Invalid observation time epoch')
      origin = day_number(iy, im, id)
      target = day_number(year, month, day)
      ending = target + 1
      if (monthly) then
         target = day_number(year, month, 1)
         im = mod(month, 12) + 1
         iy = year
         if (month == 12) iy = year + 1
         ending = day_number(iy, im, 1)
      end if
      lower = real(target-origin, pwp)
      upper = real(ending-origin, pwp)
      call obs_check(nf90_get_var(ncid, varid, time), 'time')
      call obs_require(ieee_is_finite(time(1)), 'Invalid observation time')
      call obs_require(time(1)/factor >= lower .and. &
                       time(1)/factor < upper, 'Observation date mismatch')
   end subroutine check_date

   !> Locate a point within a regular model axis; outer edges are half-open.
   integer function nearest_axis(axis, point) result(index_point)
      real(4), intent(in) :: axis(:)
      real(pwp), intent(in) :: point
      real(pwp) :: delta, position
      index_point = 0
      if (.not. ieee_is_finite(point)) return
      delta = real(axis(2)-axis(1), pwp)
      position = (point-real(axis(1), pwp))/delta + 1
      if (position < 0.5_pwp .or. &
          position >= real(size(axis), pwp)+0.5_pwp) return
      index_point = floor(position+0.5_pwp)
   end function nearest_axis

   !> Validate state units, packing and the regular AMM7 model axes.
   subroutine validate_field(id_field, temperature)
      use statevector_pdaf, only: sfields
      use nemo_pdaf, only: nav_lon, nav_lat, nav_lev, ni_p, nj_p
      integer, intent(in) :: id_field
      logical, intent(in) :: temperature
      real(pwp) :: dx, dy
      call obs_require(id_field >= 1 .and. id_field <= size(sfields), &
                       'Observation field index outside sfields')
      call obs_require(sfields(id_field)%ndims == 3, &
                       'EN4/OSTIA requires a three-dimensional field')
      call obs_require(sfields(id_field)%transform == 0, &
                       'EN4/OSTIA requires untransformed physical fields')
      if (temperature) then
         call obs_require(trim(sfields(id_field)%variable) == 'thetao', &
                          'Temperature observations require thetao')
         call obs_require(trim(sfields(id_field)%unit) == 'degC', &
                          'Temperature state must use degC')
      else
         call obs_require(trim(sfields(id_field)%variable) == 'so', &
                          'Salinity observations require so')
         call obs_require(trim(sfields(id_field)%unit) == 'psu', &
                          'Salinity state must use practical salinity')
      end if
      call obs_require(ni_p >= 2 .and. nj_p >= 2, &
                       'Observation mapping requires at least 2x2 cells')
      dx = real(nav_lon(2)-nav_lon(1), pwp)
      dy = real(nav_lat(2)-nav_lat(1), pwp)
      call obs_require(dx > 0 .and. dy > 0, 'Model axes must increase')
      call obs_require(all(abs(real(nav_lon(2:)-nav_lon(:ni_p-1), &
                       pwp)-dx) < 1.e-4_pwp), 'Irregular model longitude')
      call obs_require(all(abs(real(nav_lat(2:)-nav_lat(:nj_p-1), &
                       pwp)-dy) < 1.e-4_pwp), 'Irregular model latitude')
      call obs_require(all(nav_lon >= -180) .and. &
                       all(nav_lon <= 180), 'Model longitude must be signed')
      call obs_require(all(nav_lev >= 0), 'Model depth must be in metres')
      call obs_require(all(nav_lev(2:) > nav_lev(:size(nav_lev)-1)), &
                       'Model depths must increase')
   end subroutine validate_field

   !> Map a wet grid cell into all three supported state-vector layouts.
   integer function state_index(id_field, i, j, k) result(index_state)
      use statevector_pdaf, only: sfields, dim_state_p
      use nemo_pdaf, only: tmask, ni_p, nj_p, nk_p, use_wet_state, &
                          idx_nwet, nwet, halo0, halo1
      integer, intent(in) :: id_field, i, j, k
      integer :: column
      index_state = 0
      if (i < 1 .or. i > ni_p .or. j < 1 .or. j > nj_p) return
      ! Halo cells are not owned by this rank: avoid duplicate observations.
      if (i <= halo0(1) .or. i > ni_p-halo1(1)) return
      if (j <= halo0(2) .or. j > nj_p-halo1(2)) return
      if (k < 1 .or. k > nk_p) return
      if (tmask(i,j,k) /= 1.0_pwp) return
      if (tmask(i,j,1) /= 1.0_pwp) return
      select case (use_wet_state)
      case (1)
         index_state = idx_nwet(i,j) + nwet*(k-1)
      case (2)
         ! The current packing assumes contiguous wet surface columns.
         column = count(tmask(i,j,:) == 1.0_pwp)
         call obs_require(all(tmask(i,j,:column) == 1.0_pwp), &
                          'Wet-only packing requires contiguous columns')
         index_state = idx_nwet(i,j) + k - 1
      case default
         index_state = i + ni_p*(j-1) + ni_p*nj_p*(k-1)
      end select
      index_state = index_state + sfields(id_field)%off
      call obs_require(index_state <= dim_state_p, &
                       'Observation index exceeds local state dimension')
   end function state_index

   !> Gather all ranks, including ranks with zero valid observations.
   subroutine finish_observations(thisobs, n, values, ivar, coords, &
                                  indices, radius, dim_obs)
      use PDAF, only: obs_f, PDAFomi_gather_obs
      use parallel_pdaf, only: mype_filter
      type(obs_f), intent(inout) :: thisobs
      integer, intent(in) :: n, indices(:)
      real(pwp), intent(in) :: values(:), ivar(:), coords(:,:), radius
      integer, intent(out) :: dim_obs
      integer :: allocated_size
      allocated_size = max(1, n)
      if (allocated(thisobs%id_obs_p)) deallocate(thisobs%id_obs_p)
      allocate(thisobs%id_obs_p(1, allocated_size))
      thisobs%id_obs_p = 0
      if (n > 0) thisobs%id_obs_p(1,:n) = indices(:n)
      thisobs%doassim = 1
      thisobs%disttype = 3 ! Haversine; horizontal radians, distance metres.
      thisobs%ncoord = 2
      thisobs%use_global_obs = 1
      call PDAFomi_gather_obs(thisobs, n, values(:allocated_size), &
           ivar(:allocated_size), coords(:,:allocated_size), 2, &
           radius, dim_obs)
      write (*, '(a,i0,a,i0,a,i0)') 'NEMO-PDAF: rank ', mype_filter, &
           ', local observations ', n, ', full observations ', dim_obs
   end subroutine finish_observations
end module obs_grid_pdaf
