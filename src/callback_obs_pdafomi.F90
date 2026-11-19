!> PDAF-OMI interfaces for SSH, EN4 and OSTIA observations.
!!
!! Keep type order identical in full initialization, observation operators,
!! and local initialization. EN4 packs temperature before salinity.
!! Disabled types still reset their flags at each analysis call.

!> Initialize the full observation vector on all filter ranks.
subroutine init_dim_obs_pdafomi(step, dim_obs)
   use obs_ssh_mgrid_pdafomi, only: assim_ssh_mgrid, &
                                  init_dim_obs_ssh_mgrid
   use obs_en4_pdafomi, only: init_dim_obs_en4
   use obs_ostia_pdafomi, only: init_dim_obs_ostia
   implicit none
   integer, intent(in) :: step
   integer, intent(out) :: dim_obs
   integer :: dim_ssh, dim_en4, dim_ostia

   dim_ssh = 0
   if (assim_ssh_mgrid) call init_dim_obs_ssh_mgrid(step, dim_ssh)
   call init_dim_obs_en4(step, dim_en4)
   call init_dim_obs_ostia(step, dim_ostia)
   dim_obs = dim_ssh + dim_en4 + dim_ostia
end subroutine init_dim_obs_pdafomi

!> Apply full observation operators in the initialized type order.
subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)
   use mod_kind_pdaf
   use obs_ssh_mgrid_pdafomi, only: assim_ssh_mgrid, obs_op_ssh_mgrid
   use obs_en4_pdafomi, only: obs_op_en4
   use obs_ostia_pdafomi, only: obs_op_ostia
   implicit none
   integer, intent(in) :: step, dim_p, dim_obs
   real(pwp), intent(in) :: state_p(dim_p)
   real(pwp), intent(inout) :: ostate(dim_obs)

   if (assim_ssh_mgrid) &
      call obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)
   call obs_op_en4(dim_p, dim_obs, state_p, ostate)
   call obs_op_ostia(dim_p, dim_obs, state_p, ostate)
end subroutine obs_op_pdafomi

!> Accumulate local observation dimensions around an analysis column.
subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)
   use nemo_pdaf, only: nwet
   use obs_ssh_mgrid_pdafomi, only: assim_ssh_mgrid, &
                                  init_dim_obs_l_ssh_mgrid
   use obs_en4_pdafomi, only: init_dim_obs_l_en4
   use obs_ostia_pdafomi, only: init_dim_obs_l_ostia
   implicit none
   integer, intent(in) :: domain_p, step, dim_obs
   integer, intent(out) :: dim_obs_l

   dim_obs_l = 0
   if (nwet <= 0) return
   if (assim_ssh_mgrid) &
      call init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_en4(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_ostia(domain_p, step, dim_obs, dim_obs_l)
end subroutine init_dim_obs_l_pdafomi
