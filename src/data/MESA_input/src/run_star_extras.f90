! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def

      implicit none

      ! x_logical_ctrl, x_integer_ctrl, x_ctrl variables
      logical :: do_mix_gradT, do_stop_at_carbondep, do_stop_at_corecollapse
      integer :: rolling_count

      ! Global variables -- helps with outputs, etc.
      real (dp) :: m_core, mass_PZ, delta_r_PZ, alpha_PZ, r_core, r_PZ, rho_core_top
      ! real (dp) :: X_ini
      integer :: k_PZ_top, k_PZ_bottom
      logical :: doing_DBP = .false., step2 = .false.

      ! Rolling history variables; need history of: dt and r_core (_PZ variable is core + PZ)
      real(dp), allocatable :: r_core_history(:), r_PZ_history(:), dt_history(:)
      integer :: history_count = 0
      real(dp) :: avg_time = 0.0_dp
      real(dp) :: r_core_rolled, delta_r_PZ_rolled, m_core_rolled, mass_PZ_rolled, alpha_PZ_rolled

      !extra meshing controls
       real(dp) :: xtra_dist_below, xtra_dist_above_ov, xtra_dist_above_bv, xtra_coeff_mesh

       namelist /xtra_coeff_cb/ &
          xtra_dist_below, &
          xtra_dist_above_ov, &
          xtra_dist_above_bv, &
          xtra_coeff_mesh


      ! these routines are called by the standard run_star check_model
      contains

      ! include 'standard_run_star_extras.inc'

      subroutine extras_controls(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          s% extras_startup => extras_startup
          s% extras_start_step => extras_start_step
          s% extras_check_model => extras_check_model
          s% extras_finish_step => extras_finish_step
          s% extras_after_evolve => extras_after_evolve
          s% how_many_extra_history_columns => how_many_extra_history_columns
          s% data_for_extra_history_columns => data_for_extra_history_columns
          s% how_many_extra_profile_columns => how_many_extra_profile_columns
          s% data_for_extra_profile_columns => data_for_extra_profile_columns

          s% how_many_extra_history_header_items => how_many_extra_history_header_items
          s% data_for_extra_history_header_items => data_for_extra_history_header_items
          s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
          s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

          s% other_overshooting_scheme => extended_convective_penetration
          s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
 
          ! Custom logical controls
          do_stop_at_carbondep = s% x_logical_ctrl(1)
          do_stop_at_corecollapse = s% x_logical_ctrl(2)
          do_mix_gradT = s% x_logical_ctrl(10)
          ! Add extra meshing
          s% use_other_mesh_delta_coeff_factor = s% x_logical_ctrl(11)
          call read_inlist_xtra_coeff_core_boundary(ierr) ! Read inlist
          if (ierr /= 0) return
          s% other_mesh_delta_coeff_factor => mesh_delta_coeff_core_boundary

          ! Custom integer controls
          rolling_count = s % x_integer_ctrl(1)
          if (rolling_count < 1) rolling_count = 1

          ! allocate DBCP rolling memory
          write(*,*) 'allocating'    
          allocate(r_core_history(rolling_count))
          allocate(r_PZ_history(rolling_count))
          allocate(dt_history(rolling_count))

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: restart
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return        
      end subroutine extras_startup


      integer function extras_start_step(id)
          integer, intent(in) :: id
          integer :: ierr, i
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          doing_DBP = .false.
          extras_start_step = 0

          if (history_count < rolling_count) then
              history_count = history_count + 1
          end if
          do i = rolling_count, 2, -1
              r_core_history(i) = r_core_history(i-1)
              r_PZ_history(i) = r_PZ_history(i-1)
              dt_history(i) = dt_history(i-1)
          end do
          r_core_history(1) = 0
          r_PZ_history(1) = 0
          dt_history(1) = s % dt
          avg_time = 0
          do i = 1, history_count
                avg_time = avg_time + dt_history(i)
          end do
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          logical :: do_retry
          integer k
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_check_model = keep_going

          ! Flag PZ as anonymous_mixing
          if (doing_DBP) then
            do k=s%nz, 1, -1
                if (s% r(k) < r_PZ .and. s% r(k) > r_core_history(1)) then
                    s%mixing_type(k) = anonymous_mixing
                end if
            end do
          endif

          do_retry = .false.

          ! by default, indicate where (in the code) MESA terminated
          if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_history_columns = 12
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_history_column_name) :: names(n)
          real(dp) :: vals(n)
          real(dp) :: r_cb
          integer :: k
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          call star_eval_conv_bdy_r(s, 1, r_cb, ierr)

          ! TODO: update this
          if (.not. doing_DBP) then
              mass_PZ=0
              delta_r_PZ=0
              alpha_PZ=0

              if (s% mixing_type(s% nz) /= convective_mixing) then
                  r_core = 0
                  m_core = 0
                  rho_core_top = 0
              else
                  call star_eval_conv_bdy_k(s, 1, k, ierr)
                  r_core = r_cb
                  m_core = s%m(k)
                  rho_core_top = s%rho(k)
              endif
          endif

          names(1) = 'm_core'
          names(2) = 'mass_pen_zone'
          names(3) = 'delta_r_pen_zone'
          names(4) = 'alpha_pen_zone'
          names(5) = 'r_core'
          names(6) = 'rho_core_top_pen'
          names(7) = 'r_cb'
          names(8)  = 'roll_m_core'
          names(9)  = 'roll_mass_pen_zone'
          names(10) = 'roll_delta_r_pen_zone'
          names(11) = 'roll_alpha_pen_zone'
          names(12) = 'roll_r_core'

          delta_r_PZ_rolled = r_PZ - r_core_rolled
          m_core_rolled = 0.0_dp
          do k = s% nz, 1, -1
            if (s% r(k) >= r_core_rolled .and. m_core_rolled .eq. 0.0_dp) then
                m_core_rolled = s % m(k)
                alpha_PZ_rolled = delta_r_PZ_rolled / s % scale_height(k)
            else if (s % r(k) >= r_PZ) then
                mass_PZ_rolled = s % m(k) - m_core_rolled
                exit
            end if
          end do

!          write(*,*) "m_cores", m_core/Msun, m_core_rolled/Msun
!          write(*,*) "m_PZ", mass_PZ/Msun, mass_PZ_rolled/Msun
!          write(*,*) "alpha_PZ", alpha_PZ, alpha_PZ_rolled
!          write(*,*) "r_cores", r_core/Rsun, r_core_rolled/Rsun
!          write(*,*) "delta_r_PZ", delta_r_PZ/Rsun, delta_r_PZ_rolled/Rsun

          vals(1) = m_core/Msun
          vals(2) = mass_PZ/Msun
          vals(3) = delta_r_PZ/Rsun
          vals(4) = alpha_PZ
          vals(5) = r_core/Rsun
          vals(6) = rho_core_top
          vals(7) = r_cb/Rsun
          vals(8) = m_core_rolled / Msun
          vals(9) = mass_PZ_rolled / Msun
          vals(10) = delta_r_PZ_rolled / Rsun
          vals(11) = alpha_PZ_rolled
          vals(12) = r_core_rolled / Rsun

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
          integer, intent(in) :: id, n, nz
          character (len=maxlen_profile_column_name) :: names(n)
          real(dp) :: vals(nz,n)
          integer, intent(out) :: ierr
          ! integer :: vals_nr
          type (star_info), pointer :: s
          integer :: k
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_history_column_name) :: names(n)
          real(dp) :: vals(n)
          type(star_info), pointer :: s
          integer, intent(out) :: ierr
          ierr = 0
          call star_ptr(id,s,ierr)
          if(ierr/=0) return

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_profile_column_name) :: names(n)
          real(dp) :: vals(n)
          type(star_info), pointer :: s
          integer, intent(out) :: ierr
          ierr = 0
          call star_ptr(id,s,ierr)
          if(ierr/=0) return

          ! here is an example for adding an extra profile header item
          ! also set how_many_extra_profile_header_items
          ! names(1) = 'mixing_length_alpha'
          ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s
          character (len=200) :: fname
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_finish_step = keep_going

          ! stop at carbon depletion
          if (do_stop_at_carbondep) then
             if ((s%xa(s%net_iso(io16), s%nz) >= 0.5 ) .and. (s%xa(s%net_iso(ic12), s%nz) <= 1e-5)) then
                write(*,*) "Carbon depletion"
                extras_finish_step = terminate
                write(fname, fmt="(a10)") 'C_depl.mod'
                call star_write_model(s% id, fname, ierr)
             end if
          end if
          ! stop at onset of core-collapse
          if (do_stop_at_corecollapse) then
             ! change net on the fly post C depletion
             if ((s%xa(s%net_iso(io16), s%nz) >= 0.5 ) .and. (s%xa(s%net_iso(ic12), s%nz) <= 1e-5)) then
                write(fname, fmt="(a10)") 'C_depl.mod'
                call star_write_model(s% id, fname, ierr)
                s% job% change_net = .true.
                s% job% change_initial_net = .true.
                s% job% new_net_name = "mesa_128.net"
                write(*,*) "Change net to ", s% job%new_net_name
                ! flip switch so we don't enter here again
                do_stop_at_corecollapse = .false.
             end if
          end if

          if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
      end subroutine extras_after_evolve

!!! CUSTOM

   subroutine other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
          real(dp) :: fraction, Peclet_number, diffusivity    ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
          integer :: k, i, indx, indx2, indx3
          logical, parameter :: DEBUG = .FALSE.
          real(dp), allocatable :: conv_vel(:)

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          allocate(conv_vel(s% nz))

          if (s%D_mix(1) .ne. s%D_mix(1)) return  ! To ignore iterations where Dmix and gradT are NaNs

          if (s%num_conv_boundaries .lt. 1) then  ! Is zero at initialisation of the run
          if (DEBUG) then
              write(*,*) 'runstarex_gradT: skip since there are no convective boundaries'
          end if
          return
          endif

          if (do_mix_gradT .eqv. .false.) then
              if (DEBUG) then
                  write(*,*) 'runstarex_gradT: skip since x_logical_ctrl(10) == false.'
              end if
              return
          end if

          do k= s%nz, 1, -1
              if (s%D_mix(k) <= s% min_D_mix) exit
              if (s%conv_vel(k) > 0) then
                conv_vel(k) = s%conv_vel(k)
              else if (k < s% nz) then
                ! Use conv_vel from last cell
                conv_vel(k) = conv_vel(k+1)
              else 
                conv_vel(k) = 0
              endif
              s% adjust_mlt_gradT_fraction(k) = calculate_peclet_fraction(s, k, conv_vel(k), ierr)
              !write(*,*) 'fraction, conv_vel', s% adjust_mlt_gradT_fraction(k), conv_vel(k)
          end do

      end subroutine other_adjust_mlt_gradT_fraction_Peclet

      function calculate_peclet_fraction(s, k, vc, ierr) result(fraction)
            type(star_info), pointer :: s
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            real(dp), intent(in) :: vc
            real(dp) :: fraction, Peclet_number, diffusivity    
            diffusivity = 16.0_dp * boltz_sigma * pow3(s% T(k)) / ( 3.0_dp * s% opacity(k) * pow2(s% rho(k)) * s% cp(k) )
            Peclet_number = vc * s% scale_height(k) * s% mixing_length_alpha / diffusivity

            if (Peclet_number >= 100.0_dp) then
                fraction = 1.0_dp
            else if (Peclet_number .le. 0.01_dp) then
                fraction = -1.0_dp
            else
                fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
            end if
      end function calculate_peclet_fraction

      subroutine extended_convective_penetration(id, i, j, k_a, k_b, D, vc, ierr)
          integer, intent(in) :: id, i, j
          integer, intent(out) :: k_a, k_b
          real(dp), intent(out), dimension(:) :: D, vc
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          logical, parameter :: DEBUG = .FALSE.
          real(dp) :: f, f0
          real(dp) :: D0, Delta0
          real(dp) :: w
          real(dp) :: factor
          real(dp) :: r_cb, Hp_cb
          real(dp) :: r_ob, D_ob, vc_ob
          logical  :: outward
          integer  :: dk, k, k_ob, indx, indx2, indx3
          real(dp) :: r, dr, r_step

          ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
          ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
          ! using the j'th set of overshoot parameters. The overshoot
          ! follows the extended convective penetration scheme description by Mathias
          ! Michielsen, "Probing the shape of the mixing profile and of the thermal
          ! structure at the convective core boundary through asteroseismology",
          ! A&A, 628, 76 (2019)

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if ((i .ne. 1) .or. (s%mixing_type(s%nz) .ne. convective_mixing)) then
              write(*,'(A,i2,A,i2)') 'ERROR: dissipation_balanced_penetration can only be used for core convection, &
                    &so the first convective boundary. The routine got called for convective boundary number ',i, &
                    &', and the mixing type in the core was', s%mixing_type(s%nz)
              ierr = -1
              return
          end if


          if (.not. step2) then         
              call dissipation_balanced_penetration(s, id) !, m_core, mass_PZ, delta_r_PZ, alpha_PZ, r_core, rho_core_top)
              end if
          ! Extract parameters
          ! alpha_PZ is distance from core boundary outward, so add f0 to it to make PZ zone reach that region              
          f = alpha_PZ + s%overshoot_f0(j) 

          ! extend of step function (a_ov)
          f0 = s%overshoot_f0(j)
          D0 = s%overshoot_D0(j)
          Delta0 = s%overshoot_Delta0(j)

          if (f < 0.0_dp .OR. f0 <= 0.0_dp) then
              write(*,*) 'ERROR: for extended convective penetration, must set f0 > 0, and f >= 0'
              write(*,*) 'see description of overshooting in star/defaults/control.defaults'
              ierr = -1
              return
          end if

          ! Evaluate convective boundary (_cb) parameters
          call star_eval_conv_bdy_r(s, i, r_cb, ierr)
          if (ierr /= 0) return


          if (s%r(s%nz-10) >= r_cb) then
             !Bail -- very small core.
            if (DEBUG) then
                write(*,*) 'core is not resolved; bailing on overshoot'
            end if
              ! Use no extra mixing.
              D(:)  = s%D_mix(:)
              vc(:) = s%conv_vel(:)
           
            ! Store history before bailing
            r_PZ_history(1) = r_cb
            r_core_rolled = r_core_history(1) * dt_history(1) / avg_time
            r_PZ = r_PZ_history(1) * dt_history(1) / avg_time
            do indx = 2, history_count
                !write(*,*) 'dt_history(indx), avg_time', indx, r_PZ_history(indx), dt_history(indx), avg_time
                r_PZ = r_PZ + r_PZ_history(indx) * dt_history(indx) / avg_time
                r_core_rolled = r_core_rolled + r_core_history(indx) * dt_history(indx) / avg_time
            end do
            return
          endif

          call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
          if (ierr /= 0) return

          ! Evaluate overshoot boundary (_ob) parameters
          call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
          if (ierr /= 0) return

          ! Loop over cell faces, adding overshoot until D <= overshoot_D_min
          outward = s%top_conv_bdy(i)

          if (outward) then
              k_a = s % nz
              k_b = 1
              dk = -1
          else
              k_a = k_ob+1
              k_b = s%nz
              dk = 1
          endif

          
          if (f > 0.0_dp) then
              r_step = f*Hp_cb
              r_PZ_history(1) = r_step + r_ob

              ! Get average outer radius of PZ
              r_core_rolled = r_core_history(1) * dt_history(1) / avg_time
              r_PZ = r_PZ_history(1) * dt_history(1) / avg_time
              do indx = 2, history_count
                  !write(*,*) 'dt_history(indx), avg_time', indx, r_PZ_history(indx), dt_history(indx), avg_time
                  r_PZ = r_PZ + r_PZ_history(indx) * dt_history(indx) / avg_time
                  r_core_rolled = r_core_rolled + r_core_history(indx) * dt_history(indx) / avg_time
              end do
              factor = 1.0_dp
          else
              factor = 0.0_dp
          endif


          face_loop : do k = k_a, k_b, dk 
              ! Evaluate the extended convective penetration factor
              r = s%r(k)
              if (r > r_PZ .and. r > r_PZ_history(1)) exit face_loop
              
              ! Store the diffusion coefficient and velocity
              D(k) = (D0 + Delta0*D_ob)*factor
              vc(k) = (D0/D_ob + Delta0)*vc_ob*factor
              if (s% D_mix(k) < D(k)) then !DBCP mixing.
                !write(*,*) 'overwriting D_mix with D'
                s% D_mix(k) = D(k)
                s%conv_vel(k) = vc(k)
                s% mixing_type(k) = anonymous_mixing
              else
                D(k) = s % D_mix(k)
                vc(k) = s % conv_vel(k)
              endif
              
              ! Check for early overshoot completion
              if (D(k) < s%overshoot_D_min) then
                  k_b = k
                  exit face_loop
              endif

          end do face_loop

          if (step2) then
              ! Free to continue evolution
              step2 = .false.
          else
              ! Adjust structure with new grads.
              step2 = .true.
              call star_set_vars(id, s% dt, ierr)
          endif

          if (DEBUG) then
              write(*,*) 'step exponential overshoot:'
              write(*,*) '  r_PZ, r_cb =', r_PZ, r_cb
              write(*,*) '  Hp_cb      =', Hp_cb
          end if

      end subroutine extended_convective_penetration


      subroutine dissipation_balanced_penetration(s, id)
         use eos_def
         use star_lib
         use kap_def
         type (star_info), pointer :: s
         integer, intent(in) :: id
         real(dp), parameter :: f = 0.86d0
         real(dp), parameter :: xi = 0.6d0
         integer :: k, j, ierr, indx
         real(dp) :: Lint, V_CZ, Favg, RHS, dr, h, dLint
         real(dp) :: r_cb

         ! If we ever hit this, set flag to true.
         doing_DBP = .true.
         V_CZ = 0d0
         Lint = 0d0

         call star_eval_conv_bdy_k(s, 1, k, ierr)
         call star_eval_conv_bdy_r(s, 1, r_cb, ierr)
         k_PZ_bottom = k
         r_core_history(1) = r_cb
         m_core = s%m(k)
         r_core = r_cb
         
         rho_core_top = s%rho(k)
         h = s%scale_height(k)

         ! prescription based on Jermyn A. et al (2022)  https://arxiv.org/pdf/2203.09525.pdf
         ! Equation A1 is used here.
         ! RHS refers to right-hand side of equation A1, and Lint the integrated
         ! luminosity on either left or right side of that equation

         ! Integrate over cells that are fully in CZ
         ! r and L_conv are face values, assume they change linear within cell
         do j=s%nz,k+1,-1
            if (j .eq. s%nz) then
                dr = s%r(j)
                Lint = Lint + s%L_conv(j)*0.5 * dr
            else
                dr = s%r(j) - s%r(j+1)
                Lint = Lint + (s%L_conv(j+1) + s%L_conv(j))*0.5 * dr
            endif
         end do

        ! Take cell that is partially convective
        ! convective part of cell k
        ! L_conv goes to 0 at edge of conv zone
         dr = r_cb - s%r(k+1)
         Lint = Lint + s%L_conv(k+1)*0.5 * dr

         ! Calculate target RHS
         V_CZ = 4d0/3d0 * pi * r_cb*r_cb*r_cb
         Favg = Lint / V_CZ
         RHS = (1d0 - f) * Lint
         Lint = 0d0

         ! Integrate over RZ until we find the edge of the PZ
         ! remainder of cell k (non-convective part)
         ! Do integration explicitely, f*xi*4*pi*r^2 is moved outside of integral
         dr = s%r(k) - r_cb
         dLint = xi * f * 4d0 * pi * (pow3(s%r(k))-pow3(r_cb))/3 * Favg + (s%L(k)*0.5 * (s%grada_face(k) / s%gradr(k) - 1d0)) * dr
         Lint = dLint

         ! If remainder of cell k would already satisfy Lint > RHS
         if (Lint > RHS) then
            dr = dr*(Lint - RHS)/dLint
            mass_PZ =  s%rho(k) * 4d0/3d0 * pi * (pow3(r_cb+dr) - pow3(r_cb)) !s%m(k) - m_core !only used for history output
            delta_r_PZ = dr
            alpha_PZ = delta_r_PZ / h
            k_PZ_top = k
            return
         end if
         ! Else calculate dL_int for each cell, untill the total integrated L > RHS
         do j=k-1,1,-1
            dr = s%r(j) - s%r(j+1)
            dLint = xi * f * 4d0 * pi * (pow3(s%r(j))-pow3(s%r(j+1)))/3 * Favg &
            + ( (s%L(j+1)*(s%grada_face(j+1) / s%gradr(j+1) - 1d0) +(s%L(j)*(s%grada_face(j) / s%gradr(j) - 1d0)) )*0.5 * dr)

            if (Lint + dLint > RHS) then
               dr = dr*(RHS - Lint)/dLint
               mass_PZ = s%m(j) - m_core !only used for history output
               delta_r_PZ = s%r(j+1)+dr - r_cb
               alpha_PZ = delta_r_PZ / h
               k_PZ_top = j
               return
            end if
            Lint = Lint + dLint
         end do

      end subroutine dissipation_balanced_penetration

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Other mesh routine
      subroutine mesh_delta_coeff_core_boundary(id, ierr)
          integer, intent(in) :: id
          ! real(dp), intent(in), dimension(:) :: eps_h, eps_he, eps_z
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          logical, parameter :: dbg = .false.
          integer :: k, k_max_N2comp
          real(dp) :: Hp, r_lower, r_upper_ov, r_upper_BV

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (xtra_coeff_mesh == 1d0 .or. xtra_coeff_mesh < 0) return
          if (s% mixing_type(s% nz) /= convective_mixing) return  ! only enable for convective cores

          r_upper_ov=-1
          r_upper_BV=-1
          k_max_N2comp = MAXLOC(s% brunt_n2_composition_term(:), 1)
          ! Find boundary of convective core, and go inwards by the specified distance (in Hp)
          do k = s% nz, 1, -1
              if (s% mixing_type(k) == convective_mixing) then
                  continue
              else
                  Hp = s% scale_height(k)  !s% P(k)/(s% rho(k)*s% grav(k))
                  r_lower = max (s% r(s%nz), s% r(k) - xtra_dist_below*Hp)
                  exit
              endif
          end do

          do k = s% nz, 1, -1
              ! Start increasing the mesh once closer than the given distance (in Hp) to the core boundary
              if (s%r(k) > r_lower) then
                  if (xtra_coeff_mesh < s% mesh_delta_coeff_factor(k)) then
                      s% mesh_delta_coeff_factor(k) = xtra_coeff_mesh
                  end if
              else
                  cycle
              endif

              ! Go up to the given distance past the overshoot boundary
              if (r_upper_ov<0 .and. s% mixing_type(k) /= overshoot_mixing .and. s% mixing_type(k) /= convective_mixing) then
                  if (xtra_dist_above_ov > 0) then
                      Hp = s% scale_height(k) !s% P(k)/(s% rho(k)*s% grav(k))
                      r_upper_ov = min(s% r(1), s% r(k) + xtra_dist_above_ov*Hp)
                  else
                      r_upper_ov = 0
                  end if
              end if

              ! Go up to the given distance past the order in magnitude decrease in BV composition term, outwards of its maximum
              if (r_upper_BV<0 .and. k < k_max_N2comp .and. (s% brunt_n2_composition_term(k)*10 < maxval(s% brunt_n2_composition_term(:))) ) then
                  if (xtra_dist_above_bv > 0) then
                      Hp = s% scale_height(k) !s% P(k)/(s% rho(k)*s% grav(k))
                      r_upper_BV = min(s% r(1), s% r(k) + xtra_dist_above_bv*Hp)
                  else
                      r_upper_BV = 0
                  end if
              end if

              ! Stop increasing mesh when further than the specified distance from both the overshoot boundary and BV composition peak
              if (s% r(k) > r_upper_ov .and. s% r(k) > r_upper_BV .and. r_upper_ov >=0 .and. r_upper_BV >=0) exit
          end do

      end subroutine mesh_delta_coeff_core_boundary
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine read_inlist_xtra_coeff_core_boundary(ierr)
          use utils_lib
          integer, intent(out) :: ierr
          character (len=256) :: filename, message
          integer :: unit

          filename = 'inlist_xtra_coeff_mesh'
          write(*,*) 'read inlist_xtra_coeff_mesh'

          ! set defaults
          xtra_dist_below = 0.1d0
          xtra_dist_above_ov = 0.1d0
          xtra_dist_above_bv = 0.1d0
          xtra_coeff_mesh = 0.15d0

          open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
          if (ierr /= 0) then
              write(*, *) 'Failed to open control namelist file ', trim(filename)
          else
              read(unit, nml=xtra_coeff_cb, iostat=ierr)
              close(unit)
              if (ierr /= 0) then
                  write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
                  write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
                  write(*, *)
                  open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
                  read(unit, nml=xtra_coeff_cb)
                  close(unit)
              end if
          end if

      end subroutine read_inlist_xtra_coeff_core_boundary

      end module run_star_extras
