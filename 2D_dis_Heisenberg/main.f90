! 2D Heisenberg model with disorder
!=============================================!
program XXZ_2d_model
!=============================================!
  use configuration
  use mpi
  use mpi_modules
  implicit none
  integer :: i, j, ierr, q, value
  logical :: status

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,id,ierr)
  call mpi_comm_size(mpi_comm_world,np,ierr)

  call read_parameters
  call init_vars
  call init_ran
  call init_lattice

  do while (no_finish<ndis)

98       call init_disorder
       if (finished) goto 88
       call init_qmc
       call qmc_eqm(istep)
       call qmc_run

       if (map.eq. 1) then
         call init_map(status)
         if (.not.status) goto 98
         call init_qmc
         call qmc_eqm(istep_map)
         call qmc_run_map

         call flip_map
         call init_qmc
         call qmc_eqm(istep_map)
         call qmc_run_map

         call end_map
       endif
       deallocate(opstring)
  enddo
88  call mpi_finalize(ierr)
 end program XXZ_2d_model
