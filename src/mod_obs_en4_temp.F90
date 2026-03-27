!> PDAF-OMI observation module for Monthly global Phytoplankton Carbon,
!! between 1998-2020 at 9 km resolution
!! (derived from the Ocean Colour Climate Change Initiative v5.0 dataset)
!!
!! The subroutines in this module are for the particular handling of
!! Phytoplankton Carbon observations available on the lat-lon grid.
!!
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (`init_dim_obs`)
!! and for the observation operator (`obs_op`).
!!
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initialized by the generic routines from `PDAFomi`.
!!
MODULE mod_obs_en4_temp
   USE mod_kind_pdaf
   USE pdafomi, ONLY: obs_f, obs_l
   IMPLICIT NONE

   real(pwp), parameter :: fill_value = -9999.0
   ! Variables which are inputs to the module (usually set in init_pdaf)
   integer   :: screen = 2
   logical   :: assim_en4_temp          !< Whether to assimilate this phytoplankton chlorophyll data
   real(pwp) :: lradius_en4_temp = 1.0  !< Localization cut-off radius
   real(pwp) :: sradius_en4_temp = 1.0  !< Support radius for weight function
   character(lc)    :: path_en4_temp = '.'    !< Path to phytoplankton carbon data
   !> Instance of full observation data type - see `PDAFomi` for details.
   TYPE(obs_f), TARGET, PUBLIC :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   TYPE(obs_l), TARGET, PUBLIC :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)
CONTAINS
   SUBROUTINE init_dim_obs_en4_temp(step, dim_obs)
      USE pdafomi, ONLY: PDAFomi_gather_obs
      USE mod_parallel_pdaf,  ONLY: mype_filter, npes_filter


      INTEGER, INTENT(in)    :: step    !< Current time step
      INTEGER, INTENT(inout) :: dim_obs !< Dimension of full observation vector
      ! count number of observations
      integer                :: i, j, cnt
      integer                :: dim_obs_p

      real(pwp), allocatable :: model_obs(:, :)
      real(pwp), allocatable :: model_obs_unc(:, :)

      real(pwp), allocatable :: obs_p(:)
      real(pwp), allocatable :: ivar_obs_p(:)
      real(pwp), allocatable :: ocoord_p(:, :)

      ! *****************************
      ! *** Global setting config ***
      ! *****************************
      ! Store whether to assimilate this observation type (used in routines below)
      thisobs%doassim = 0
      IF (assim_en4_temp) thisobs%doassim = 1
      if (thisobs%doassim == 0) then
         ! This is for the case that we do not assimilate this data at the current step
         dim_obs = 0
         dim_obs_p = 0
         return
      end if
      IF (mype_filter == 0) &
         WRITE (*, '(a,4x,a)') 'NEMO-PDAF', 'Assimilate observations - EN4 temperature profiles'
      ! Specify type of distance computation
      thisobs%disttype = 2   ! 2=Geographic
      ! Number of coordinates used for distance computation.
      ! The distance compution starts from the first row
      thisobs%ncoord = 3
      ! In case of MPI parallelization restrict observations to sub-domains
      if (npes_filter>1) thisobs%use_global_obs = 0
      ! **********************************
      ! *** Read PE-local observations ***
      ! **********************************
      IF (mype_filter == 0) &
         WRITE (*, '(a,4x,a)') 'NEMO-PDAF', 'Reading - EN4 temperature profiles'

      ! Set observation dimension
      cnt = 0
      do j = 1, nlats
         do i = 1, nlons
            if ((idx_nwet(i, j) <= 0) .or. &
                (model_obs(i, j) > 1e14)) cycle
            cnt = cnt + 1
         end do
      end do
      dim_obs_p = cnt
      dim_obs = dim_obs_p

      if (npes_filter==1) then
         write (6,'(a, 4x, a, i7)') 'NEMO-PDAF', '--- number of observations from phytoplankton_chlorophyll: ', dim_obs
      else
         if (screen>2) then
            write (6,'(a, 4x, a, i4, 2x, a, i7)') 'NEMO-PDAF', 'PE', mype_filter, &
               '--- number of observations from phytoplankton_chlorophyll: ', dim_obs
         end if
      end if

      ! *** Initialize vector of observations on the process sub-domain ***
      ! *** Initialize coordinate array of observations                 ***
      ! *** Initialize process local index array                        ***
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p))
         allocate(ivar_obs_p(dim_obs_p))
         allocate(ocoord_p(thisobs%ncoord, dim_obs_p))
         ! *** Super-Obbing ***
         ! Allocate process-local index array
         allocate(thisobs%id_obs_p(1, dim_obs_p))
         call get_obs_p_superob(dim_obs_p, model_obs, model_obs_unc, ocoord_p, obs_p, ivar_obs_p)
         deallocate(model_obs, model_obs_unc)
      else
         ! for DIM_OBS_P=0
         allocate(obs_p(1))
         allocate(ivar_obs_p(1))
         allocate(ocoord_p(thisobs%ncoord, 1))
         allocate(thisobs%id_obs_p(1, dim_obs_p))
      end if

      ! This routine is generic for the case that only the observations,
      ! inverse variances and observation coordinates are gathered
      call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, lradius_bgc_chlo, dim_obs)
      ! Deallocate all local arrays
      deallocate(obs_p, ocoord_p, ivar_obs_p)

   END SUBROUTINE init_dim_obs_bgc_chlo

   subroutine read_en4_qc()
      call read_en4_profiles(time, depth, lat, lon, potm, pos_qc, prof_potm_qc, potm_qc)
   end subroutine read_en4_qc


   subroutine read_en4_profiles(n_prof, n_level, time, depth, lat, lon, potm, pos_qc, prof_potm_qc, potm_qc)
      use netcdf
      use io_utils, only: check
      use mod_nemo_pdaf, only: year, month
      implicit none
      ! output variables
      integer,      intent(out)              :: n_prof
      integer,      intent(out)              :: n_level
      real(kind=8), intent(out), allocatable :: time(:)
      real(kind=4), intent(out), allocatable :: depth(:,:)
      real(kind=4), intent(out), allocatable :: lat(:)
      real(kind=4), intent(out), allocatable :: lon(:)
      real(kind=4), intent(out), allocatable :: potm(:,:)
      integer,      intent(out), allocatable :: pos_qc(:),
      integer,      intent(out), allocatable :: prof_potm_qc(:)
      integer,      intent(out), allocatable :: potm_qc(:,:)
      ! local variables
      character(256) :: ncfile
      integer :: ncid, retval
      integer :: dimid_prof, dimid_level
      integer :: var_time, var_lat, var_lon
      integer :: var_pos_qc, var_prof_potm_qc
      integer :: var_depth, var_potm, var_potm_qc

      ncfile = trim(path_en4_temp) // '/EN.4.2.2.f.profiles.g10.' // &
               trim(adjustl(char(year))) // &
               trim(adjustl(char(month))) // trim(adjustl(char(day))) // '.nc'
      call check(nf90_open(ncfile, NF90_NOWRITE, ncid), 'open')

      call check(nf90_inq_dimid(ncid, 'N_PROF', dimid_prof))
      call check(nf90_inq_dimlen(ncid, dimid_prof, n_prof))
      call check(nf90_inq_dimid(ncid, 'N_LEVELS', dimid_level))
      call check(nf90_inq_dimlen(ncid, dimid_level, n_level))

      call check(nf90_inq_varid(ncid, 'TIME', var_time))
      call check(nf90_inq_varid(ncid, 'LATITUDE', var_lat))
      call check(nf90_inq_varid(ncid, 'LONGITUDE', var_lon))
      call check(nf90_inq_varid(ncid, 'DEPTH', var_depth))

      call check(nf90_inq_varid(ncid, 'POSITION_QC', var_pos_qc))
      call check(nf90_inq_varid(ncid, 'PROFILE_POTM_QC', var_prof_potm_qc))
      call check(nf90_inq_varid(ncid, 'POTM_QC', var_potm_qc))

      call check(nf90_inq_varid(ncid, 'POTM_CORRECTED', var_potm))

      allocate(time(n_prof))
      allocate(lat(n_prof))
      allocate(lon(n_prof))
      allocate(depth(n_level, n_prof))

      allocate(pos_qc(n_prof))
      allocate(prof_potm_qc(n_prof))
      allocate(potm(n_level, n_prof))
      allocate(potm_qc(n_level, n_prof))
      ! read coordinates
      call check(nf90_get_var(ncid, var_time, time))
      call check(nf90_get_var(ncid, var_lat, lat))
      call check(nf90_get_var(ncid, var_lon, lon))
      call check(nf90_get_var(ncid, var_depth, depth))
      ! read QC variables and profiles
      call check(nf90_get_var(ncid, var_pos_qc, pos_qc))
      call check(nf90_get_var(ncid, var_prof_potm_qc, prof_potm_qc))
      call check(nf90_get_var(ncid, var_potm_qc, potm_qc))
      ! read profiles
      call check(nf90_get_var(ncid, var_potm, potm))
      call check(nf90_close(ncid))
   end subroutine read_en4_profiles

   subroutine perform_qc_en4_temp(n_prof, n_level, prof_potm_qc)
      use nemo_pdaf, only: nav_lon, nav_lat, time_days, ni_p, nj_p, dlon, dlat
      use parallel_pdaf, only: abort_parallel
      integer,    intent(in)  :: n_prof
      integer,    intent(in)  :: n_level
      real(pwp),  intent(in)  :: prof_potm_qc(:)
      integer,    intent(out) :: n_prof_qc

      integer :: i
      real(kind=4), parameter :: fill_value = -9999.0
      ! Count number of observations that pass QC
      dim_obs_p = 0
      do i = 1, n_prof
         if (.not. profile_pass_qc(time(i), lat(i), lon(i), pos_qc(i), prof_potm_qc(i))) cycle
         do k = 1, n_level
            if (potm_qc(k, i) == 1 .and. potm(k, i) /= fill_value) dim_obs_p = dim_obs_p + 1
         end do
      end do

      ! allocate arrays for QC-passed observations
      allocate(obs_p(dim_obs_p))
      allocate(ivar_obs_p(dim_obs_p))
      allocate(ocoord_p(thisobs%ncoord, dim_obs_p))
      allocate(thisobs%id_obs_p(8, dim_obs_p))
      cnt = 0
      do i = 1, n_prof
         if (.not. profile_pass_qc(time(i), lat(i), lon(i), pos_qc(i), prof_potm_qc(i))) cycle
         do k = 1, n_level
            if (potm_qc(k, i) == 1 .and. potm(k, i) /= fill_value) then
               cnt = cnt + 1
               obs_p(cnt) = potm(k, i)
               ivar_obs_p(cnt) = potm_unc(k, i)
               ocoord_p(1, cnt) = lat(i)
               ocoord_p(2, cnt) = lon(i)
               ocoord_p(3, cnt) = depth(k, i)
               ! get the grid cell corners corresponding to the observation location
               n_j = int( ( lat(i) - nav_lat(1) ) / dlat ) + 1
               n_i = int( ( lon(i) - nav_lon(1) ) / dlon ) + 1
               do while (nav_lon(n_i) > lon(i))
                  n_i = n_i + 1
               end do
               do while (nav_lon(n_i + 1) < lon(i))
                  n_i = n_i - 1
               end do
               do while (nav_lat(n_j) > lat(i))
                  n_j = n_j + 1
               end do
               do while (nav_lat(n_j + 1) < lat(i))
                  n_j = n_j - 1
               end do
               n_k = 1
               do while (nav_lev(n_k) > depth(k, i))
                  n_k = n_k + 1
               end do
               thisobs%id_obs_p(:, cnt) = []
               if ((idx_nwet(i, j)<=0.0)) cycle
               cnt = cnt + 1
               ! Set index of grid point
               if (use_wet_state==1 .or. use_wet_state==2) then
                  thisobs%id_obs_p(1, cnt) = &
                     idx_nwet(i, j) + sfields(ibgc_chlo)%off
               else
                  thisobs%id_obs_p(1, cnt) = &
                     i + nlons*(j-1) + sfields(ibgc_chlo)%off
               end if
            end if
         end do
      end do
   end subroutine perform_qc_en4_temp

   logical function profile_pass_qc(time, lat, lon, pos_qc, prof_potm_qc)
      use nemo_pdaf, only: nav_lon, nav_lat, time_days, ni_p, nj_p, dlon, dlat
      real, intent(in) :: time, lat, lon
      integer, intent(in) :: pos_qc, prof_potm_qc

      profile_pass_qc = .true.
      if (time < time_days .or. time >= time_days + 1) profile_pass_qc = .false.
      if (lat < nav_lat(1) .or. lat > nav_lat(nj_p) + 0.5*dlat) profile_pass_qc = .false.
      if (lon < nav_lon(1) - 0.5*dlon .or. lon > nav_lon(ni_p) + 0.5*dlon) profile_pass_qc = .false.
      if (pos_qc /= 1) profile_pass_qc = .false.
      if (prof_potm_qc /= 1) profile_pass_qc = .false.
   end function profile_pass_qc

END MODULE mod_obs_en4_temp