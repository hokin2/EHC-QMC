

subroutine init_qmc
  use configuration
  implicit none
  call init_conf
  call init_meas
  call init_weights
  call init_vrtx
  call encode_states
endsubroutine init_qmc

subroutine qmc_eqm(steps)
  use configuration
  implicit none
  integer :: i, steps
  do i=1,steps
     call mcstep
     call adjstl
     if (mod(i,max(steps/20,1))==0) call adjnl
  end do
  write(*,'(1a,1i4.4,1a,1i5.5,1a,1i2.2,1a)') 'node-',id,' dis-',no_finish,' m-',mode ,' eqm  done'
endsubroutine qmc_eqm

subroutine qmc_run
  use configuration
  implicit none
  integer :: i, j
  do i=1,nbins
     call init_bin_meas
     do j=1,mstep
        call mcstep
        call meas_p
        call meas_e
        if (map.eq.1) call meas_e2
        if(mod(j,max(mstep/5,1))==0) then
           write(*,51) 'node-',id,' dis-',no_finish,' m-',mode ,' meas done ', dble(j)*100/(mstep), '% bin-',i
51 format(1a,1i4.4,1a,1i5.5,1a,1i2.2,1a,1f6.1,1a,1i4.4)
        endif
     enddo
     call record_entropy(mstep)
     call write_bin_entropy
     call record_energy(mstep)
     call write_bin_energy
     if (map.eq.1) call record_energy2(mstep) ! pass to findh3
  enddo
  call average_bins(nbins)
  call write_average
endsubroutine qmc_run

subroutine qmc_run_map
  use configuration
  implicit none
  integer :: i, j
  do i=1,nbins_map
    call init_bin_meas
     do j=1,mstep_map
        call mcstep
        call meas_e
        if(mod(j,max(mstep_map/5,1))==0) then
           write(*,51) 'node-',id,' dis-',no_finish,' m-',mode ,' meas done ', dble(j)*100/(mstep_map), '% bin-',i
51 format(1a,1i4.4,1a,1i5.5,1a,1i2.2,1a,1f6.1,1a,1i4.4)
        endif
     enddo
     call record_energy(mstep_map)
     ! call write_bin_energy
  enddo
  call average_bins(nbins_map)
  call write_energy
endsubroutine qmc_run_map
