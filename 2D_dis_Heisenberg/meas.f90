! standard output format in MPI_INTEGER
! do while (file_opens)
!    INQUIRE(FILE="xxx.dat", opened=file_opens)
!    if (file_opens) then
!      call sleep(1)
!    else
!       open(10,file='xxx.dat',status='unknown',position='append')
!       close(10)
!    endif
! enddo


subroutine init_meas
  use configuration
  implicit none
    QCM(:,:) = 0.0d0
    QCM_2(:,:) = 0.0d0
    energy_avg(:,:) = 0.0d0
    etot(:) = 0.0d0
    mag(:,:) = 0.0d0

    p_af1_avg(:)=0
    p_af2_avg(:)=0
    p_h_avg(:)=0
    p_max_avg(:)=0
    s_2_avg(:)=0
endsubroutine init_meas

subroutine init_bin_meas
  use configuration
  implicit none
  p_af1=0; p_af2=0; p_h=0; p_max=0; s_2=0;
  energy(:) = 0.0d0
  energy2(:,:) = 0.0d0
  endsubroutine init_bin_meas

subroutine average_bins(bins)
    use configuration
    implicit none
    integer :: bins
    real(8) :: c, c1
    c = 1/dble(bins)
    if (bins==1) then
      c1 = 0.0d0
    else
      c1 = 1/dsqrt(dble(bins-1))
    endif
    energy_avg(:,1) = energy_avg(:,1)*c
    energy_avg(:,2) = sqrt(energy_avg(:,2)*c - energy_avg(:,1)*energy_avg(:,1))*c1
    etot(1) = etot(1)*c
    etot(2) = sqrt(etot(2)*c - etot(1)*etot(1))*c1
    mag(:,1) = mag(:,1)*c
    mag(:,2) = sqrt(mag(:,2)*c - mag(:,1)*mag(:,1))*c1
    p_af1_avg(1) = p_af1_avg(1)*c
    p_af1_avg(2) = sqrt(p_af1_avg(2)*c - p_af1_avg(1)*p_af1_avg(1))*c1
    p_af2_avg(1) = p_af2_avg(1)*c
    p_af2_avg(2) = sqrt(p_af2_avg(2)*c - p_af2_avg(1)*p_af2_avg(1))*c1
    p_h_avg(1) = p_h_avg(1)*c
    p_h_avg(2) = sqrt(p_h_avg(2)*c - p_h_avg(1)*p_h_avg(1))*c1
    p_max_avg(1) = p_max_avg(1)*c
    p_max_avg(2) = sqrt(p_max_avg(2)*c - p_max_avg(1)*p_max_avg(1))*c1
    s_2_avg(1) = s_2_avg(1)*c
    s_2_avg(2) = sqrt(s_2_avg(2)*c - s_2_avg(1)*s_2_avg(1))*c1
  endsubroutine average_bins

  subroutine write_average
    use configuration
    implicit none
    call write_mag
    call write_entropy
    call write_energy
  endsubroutine write_average

!=======================================!
subroutine meas_p
!=======================================!
  use configuration
  implicit none

  integer :: i,o,op,b,s1,s2,l_i,l_q,q, q_2
  real(8) :: p_q, p_q_sum
  real(8) ::  rndm,p
  integer :: status
  integer(16) ::  val, max_count
  do i=1, qmax
    call h_table(i)%reset
  enddo
  no_states(:) = 0
  record_spn(:, :) = 0
  l_i = 1

  do q = 1, qmax
     do l_q = 1, l_len(q) !i is the propagation time
       op = opstring(l_i)
       call hist_spn(q)
       if (mod(op,2)==1) then !left off diagonal
          b=op/2
          s1=bond(1,b)
          s2=bond(2,b)
          spn(q,s1)= -spn(q,s1)
          spn(q,s2)= -spn(q,s2)
       endif

       l_i = l_i + 1
     end do
     if (q<qmax) l_i = l_i+1
  enddo
  if (l_i-1 .ne. l) then
    write(*,* ) 'sth wrong 1',l_i-1,l
    stop
  endif

  ! calculation of s_2
  p_q_sum = 0.0d0
  max_count = 0
  do i = 1,no_states(1)
    call h_table(1)%get_value(record_spn(1, i), val, status)
    p_q = float(val)
    if (val>max_count) then
      max_count = val
      code_max = record_spn(1, i)
    endif
    do q = 2, qmax
        call h_table(q)%get_value(record_spn(1, i), val, status)
        if (status.eq.-1) then
          val = 0
        endif
        p_q = p_q*float(val)
    end do
    ! later replace by Kahan’s summation
    p_q_sum = p_q_sum + p_q
  end do
  p_q_sum = p_q_sum/float(l_len(1))/float(l_len(2))
  s_2 = s_2 + p_q_sum

  ! count the most probable config
  do q = 1, qmax
    call h_table(q)%get_value(code_af1, val, status)
    if (status>-1) p_af1 = p_af1 + val
    call h_table(q)%get_value(code_af2, val, status)
    if (status>-1) p_af2 = p_af2 + val
    call h_table(q)%get_value(code_h, val, status)
    if (status>-1) p_h = p_h + val
    call h_table(q)%get_value(code_max, val, status)
    if (status>-1) p_max = p_max + val
  end do

end subroutine meas_p

!=======================================!
subroutine meas_e
!=======================================!
  use configuration
  implicit none

  integer :: i,o,op,b,s1,s2,l_i,l_q,q
  integer ::  val, status, max_count
  integer:: bond2(nb), sum_b, i_nn
  real(8):: prev_t(nn+2),t_h

  no_states(:) = 0
  l_i = 1

  do q = 1, qmax
      do i=1,nb
         bond2(i) = spn(q,bond(1,i))*spn(q,bond(2,i))
      enddo
      sum_b = sum(bond2)
      prev_t(:) = 0 ! to record the earliest tau before change
      t_h=0 ! the time for Hamiltonian
       do l_q = 1, l_len(q) !i is the propagation time
         op = opstring(l_i)
         if (op>0) then
           t_h = t_h+1
         if (mod(op,2)==1) then !off diagonal
           b=op/2
           s1=bond(1,b)+2
           s2=bond(2,b)+2
           ! diagponal
           energy(s1)=energy(s1)+(t_h-prev_t(s1))/dble(nh_q(q))*dble(spn(q,bond(1,b)))
           energy(s2)=energy(s2)+(t_h-prev_t(s2))/dble(nh_q(q))*dble(spn(q,bond(2,b)))
           energy(1)=energy(1)+(t_h-prev_t(1))*dble(sum_b)/dble(nh_q(q))
           ! off-diagonal
           energy(2)=energy(2)+1.

           prev_t(s1) = t_h
           prev_t(s2) = t_h
           prev_t(1) = t_h

            spn(q,bond(1,b))= -spn(q,bond(1,b))
            spn(q,bond(2,b))= -spn(q,bond(2,b))
            do i_nn=1,6
               bond2(bond_nn(i_nn,b))=-bond2(bond_nn(i_nn,b))
               sum_b=sum_b+2*bond2(bond_nn(i_nn,b))
            enddo
         endif
       endif
         l_i = l_i + 1
       end do
       ! add at the end of opstring
       do i=3,nn+2
          energy(i)=energy(i)+(t_h-prev_t(i))/dble(nh_q(q))*dble(spn(q,i-2))
       enddo
       energy(1)=energy(1)+(t_h-prev_t(1))*sum_b/dble(nh_q(q))

       if (q<qmax) l_i = l_i+1
    enddo
    if (l_i-1 .ne. l) then
      write(*,* ) 'sth wrong 1',l_i-1,l
      stop
    endif

end subroutine meas_e

!=======================================!
subroutine meas_e2
!=======================================!
  use configuration
  implicit none

  integer :: i,j, o,op,b,s1,s2,l_i,l_q,q, q_2,ss(2),k
  integer ::  val, status, max_count
  integer:: bond2(nb), sum_b, i_nn, prev_o
  real(8):: prev_t2(nn+2,nn+2),t_h

  no_states(:) = 0
  l_i = 1

  do q = 1, qmax
      do i=1,nb
         bond2(i) = spn(q,bond(1,i))*spn(q,bond(2,i))
      enddo
      sum_b = sum(bond2)

      prev_t2(:,:) = 0 ! to record the earliest tau before change
      t_h = 0 ! the time for Hamiltonian
      prev_o = 0

       do l_q = 1, l_len(q) !i is the propagation time
         op = opstring(l_i)
         if (op>0) then
           t_h = t_h+1
           o=mod(op,2)
         if (mod(op,2)==1) then !off diagonal
           b=op/2
           s1=bond(1,b)+2
           s2=bond(2,b)+2
           ! diagponal
           ! convariance matrix
           ! counting trick
           ! off-diagonal - off-diagonal
           if (prev_o==1) then
              energy2(2,2)=energy2(2,2)+dble(nh_q(q)-1)
           endif

           ! counting trick+ diagonal measurements
           ! off-diagonal - Jz diagonal
           energy2(2,1) = energy2(2,1)+dble(sum_b)
           energy2(1,2) = energy2(1,2)+dble(sum_b)
           ! off-diagonal - hz diagonal
           do j=3,nn+2
              energy2(2,j) = energy2(2,j)+dble(spn(q,j-2))
              energy2(j,2) = energy2(j,2)+dble(spn(q,j-2))
           enddo

           ! diagonal measurements
           ! Jz diagonal - Jz diagonal
           energy2(1,1) = energy2(1,1)+(t_h-prev_t2(1,1))*dble(sum_b*sum_b)/dble(nh_q(q))
           ! Jz diaongal - hz diagonal
           do j=3,nn+2
              energy2(1,j) = energy2(1,j)+(t_h-prev_t2(1,j))*dble(sum_b*spn(q,j-2))/dble(nh_q(q))
              energy2(j,1) = energy2(j,1)+(t_h-prev_t2(j,1))*dble(sum_b*spn(q,j-2))/dble(nh_q(q))
           enddo

           ! hz diagonal - hz diagoanl
           ! row
           do j=3,nn+2
              energy2(s1,j)= energy2(s1,j)+(t_h-prev_t2(s1,j))*dble(spn(q,s1-2)*spn(q,j-2))/dble(nh_q(q))
              energy2(s2,j)= energy2(s2,j)+(t_h-prev_t2(s2,j))*dble(spn(q,s2-2)*spn(q,j-2))/dble(nh_q(q))
           enddo
           ! column
           do j=3,nn+2
                 energy2(j,s1)= energy2(j,s1)+(t_h-prev_t2(j,s1))*dble(spn(q,s1-2)*spn(q,j-2))/dble(nh_q(q))
                 energy2(j,s2)= energy2(j,s2)+(t_h-prev_t2(j,s2))*dble(spn(q,s2-2)*spn(q,j-2))/dble(nh_q(q))
           enddo
           ! overlap btw row and column
           ss(1)=s1
           ss(2)=s2
           do j=1,2
              do k=1,2
energy2(ss(j),ss(k))= energy2(ss(j),ss(k))-(t_h-prev_t2(ss(j),ss(k)))*dble(spn(q,ss(j)-2)*spn(q,ss(k)-2))/dble(nh_q(q))
              enddo
           enddo

           prev_t2(:,s1)=t_h
           prev_t2(:,s2)=t_h
           prev_t2(s1,:)=t_h
           prev_t2(s2,:)=t_h
           prev_t2(1,:)=t_h
           prev_t2(:,1)=t_h

            spn(q,bond(1,b))= -spn(q,bond(1,b))
            spn(q,bond(2,b))= -spn(q,bond(2,b))
            do i_nn=1,6
               bond2(bond_nn(i_nn,b))=-bond2(bond_nn(i_nn,b))
               sum_b=sum_b+2*bond2(bond_nn(i_nn,b))
            enddo

          endif
       endif
         l_i = l_i + 1
         prev_o=o
       end do

        ! add at the end of opstring
        ! diagonal measurements
        ! Jz diagonal - Jz diagonal
        energy2(1,1) = energy2(1,1)+(t_h-prev_t2(1,1))*dble(sum_b*sum_b)/dble(nh_q(q))
        ! Jz diaongal - hz diagonal
        do j=3,nn+2
           energy2(1,j) = energy2(1,j)+(t_h-prev_t2(1,j))*dble(sum_b*spn(q,j-2))/dble(nh_q(q))
           energy2(j,1) = energy2(j,1)+(t_h-prev_t2(j,1))*dble(sum_b*spn(q,j-2))/dble(nh_q(q))
        enddo
        ! hz diagonal - hz diagoanl
        do i=3,nn+2
           do j=3,nn+2
              energy2(i,j)= energy2(i,j)+(t_h-prev_t2(i,j))*dble(spn(q,i-2)*spn(q,j-2))/dble(nh_q(q))
           enddo
        enddo

        ! count the contribution across boundary
        op=0
        i=1
        do while (op==0)
           op=opstring(i)
           if (op/=0) then
              o=mod(op,2)
              if ((op==1).and.(prev_o==1)) then
                 energy2(2,2)=energy2(2,2)+dble((nh_q(q)-1))
              endif
           endif
           i=i+1
        enddo

       if (q<qmax) l_i = l_i+1

    enddo
    if (l_i-1 .ne. l) then
      write(*,* ) 'sth wrong 1',l_i-1,l
      stop
    endif

end subroutine meas_e2

subroutine hist_spn(q)
  use configuration
  use hash_table
  use convert_state
  implicit none
  integer :: q, status, i
  integer(16):: key, value
  key = encode(spn(q,:))
  if (key<0) then
    write(*,*) spn(q,:),key
  endif

  call h_table(q)%get_value(key, value, status)
  if (status>-1) then
    call h_table(q)%store_value(key, value+1, status)
  else
   if (no_states(q)<mem_state) then
      value = 1
      call h_table(q)%store_value(key, value, status)
     no_states(q) = no_states(q)+1
     record_spn(q, no_states(q)) = key
   endif
  endif
  ! write(*,*) key,value
  ! write(*,*) h_table(q)%keys
  ! write(*,*) h_table(q)%vals
endsubroutine hist_spn



!========================================================!
subroutine record_energy(steps)
!========================================================!
  use configuration; implicit none

  integer :: i,steps,j,i2,j2

  real(8) :: de, v, e_avg
  real(8) :: vmax(nn+2), unit(nn+2)
  logical :: file_opens
  ! <H>
  energy = energy/dble(steps)/dble(qmax)
  mag(:,1) = mag(:,1)+ energy(3:nn+2)
  mag(:,2) = mag(:,2)+ energy(3:nn+2)*energy(3:nn+2)
  energy(1) = energy(1)/4.0
  energy(2) = -energy(2)/beta
  ! Jz already fliped negative in Jzz value
  energy(2) = energy(2)*jjz
  do i=3,nn+2
     energy(i) = energy(i)*hsites(i-2)/2.0
  enddo
  e_avg = sum(energy)/dble(nn)
  etot(1) = etot(1) + e_avg
  etot(2) = etot(2) + e_avg*e_avg

  energy_avg(:,1) = energy_avg(:,1) + energy
  energy_avg(:,2) = energy_avg(:,2) + energy*energy

end subroutine record_energy
!========================================================!

subroutine write_bin_energy
  use configuration
  implicit none
  integer :: i
  logical :: file_opens
  file_opens = .True.
   do while (file_opens)
      INQUIRE(FILE="out_bin_energy.dat", opened=file_opens)
      if (file_opens.eqv..False.) then
         open(10,file='out_bin_energy.dat',status='unknown',position='append')
         do i=1, nn+2
           write(10,'(I5,2X,I5,2X,f15.7)') no_finish, i, energy(i)
         enddo
         close(10)
      endif
   enddo
endsubroutine write_bin_energy


subroutine record_energy2(steps)
  use configuration
  implicit none
    integer :: i,steps,j,i2,j2
    real(8) :: de, v
    real(8) :: vmax(nn+2), unit(nn+2)
    logical :: file_opens

     de=beta*dble(steps)
     ! vmax : all the jmax , hmax to maintain postivity
     vmax(1)=amax(1,0)*dble(nn) !Jz
     vmax(2)=0             !J
     vmax(3:nn+2)=amax(1:nn,1)
     ! unit : the paramters in hamiltonian
     unit(1)=jjz
     unit(2)=1.0d0
     unit(3:nn+2)=hsites(1:nn)
     energy2(:,:) = energy2(:,:)/dble(qmax)
     ! <HH>
     ! Jxy-Jxy
     energy2(2,2) = energy2(2,2)/de/beta

     ! Jxy-Jz/h
     energy2(1,2) = -energy2(1,2)/de/4.0
     energy2(2,1) = -energy2(2,1)/de/4.0
     do i=3,nn+2
        energy2(i,2) = -energy2(i,2)/de*hsites(i-2)/2.0
        energy2(2,i) = -energy2(2,i)/de*hsites(i-2)/2.0
     enddo

     !Jz/h-Jz/h
     energy2(1,1)=energy2(1,1)/dble(steps)/16.0
     do i=3,nn+2
        energy2(i,1)=energy2(i,1)/dble(steps)*hsites(i-2)/8.0
        energy2(1,i)=energy2(1,i)/dble(steps)*hsites(i-2)/8.0
        do j=3,nn+2
           energy2(i,j)=energy2(i,j)/dble(steps)*hsites(i-2)*hsites(j-2)/4.0
        enddo
     enddo

     ! <HH> -<H><H>
     do i=1,nn+2
        do j=1,nn+2
           energy2(i,j) = (energy2(i,j)-energy(i)*energy(j))/unit(i)/unit(j)
        enddo
     enddo

     energy3=0.0d0
     do i=1,nn+2
        do j=1,nn+2
           if (i<3) then
              i2=1
           else
              i2=i-1
           endif
           if (j<3) then
              j2=1
           else
              j2=j-1
           endif
           energy3(i2,j2) = energy3(i2,j2) + energy2(i,j)
        enddo
     enddo

     do i=1,nn+1
        do j=1,i
           if (i/=j) then
              v=(energy3(i,j)+energy3(j,i))/2.0d0
              energy3(i,j)=v
              energy3(j,i)=v
           endif
        enddo
     enddo

     do i=1,nn+1
        do j=1,nn+1
           QCM(j,i) = QCM(j,i) + energy3(j,i)
           QCM_2(j,i) = QCM_2(j,i) + energy3(j,i)*energy3(j,i)
        enddo
     enddo

end subroutine record_energy2


subroutine record_entropy(steps)
  use configuration
  implicit none
  real(8):: c1 ,c2
  logical:: file_opens
  integer :: steps
  c1 = float(steps)
  c2 = float(steps)*float(l)
  p_af1 = p_af1/c2
  p_af2 = p_af2/c2
  p_h = p_h/c2
  p_max = p_max/c2
  s_2 = s_2/c1
  p_af1_avg(1) = p_af1_avg(1) + p_af1
  p_af1_avg(2) = p_af1_avg(2) + p_af1*p_af1
  p_af2_avg(1) = p_af2_avg(1) + p_af2
  p_af2_avg(2) = p_af2_avg(2) + p_af2*p_af2
  p_h_avg(1) = p_h_avg(1) + p_h
  p_h_avg(2) = p_h_avg(2) + p_h*p_h
  p_max_avg(1) = p_max_avg(1) + p_max
  p_max_avg(2) = p_max_avg(2) + p_max*p_max
  s_2_avg(1) = s_2_avg(1) + s_2
  s_2_avg(2) = s_2_avg(2) + s_2*s_2

endsubroutine record_entropy

subroutine write_bin_entropy
  use configuration
  implicit none
  logical :: file_opens
  file_opens =.True.
  do while (file_opens)
     INQUIRE(FILE="out_bin_entropy.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
        open(10,file='out_bin_entropy.dat',status='unknown',position='append')
        write(10,'(1a10, 1i10, 1i40)') 'h_st:', no_finish, code_h
        write(10,'(1a10, 1i10, 1i40)') 'max_st:',no_finish, code_max
        write(10,'(1a10, 1i10, 1e20.7)') 'p_2:', no_finish, s_2
        write(10,'(1a10, 1i10, 1e20.7)') 'p_AF1:', no_finish, p_af1
        write(10,'(1a10, 1i10, 1e20.7)') 'p_AF2:', no_finish, p_af2
        write(10,'(1a10, 1i10, 1e20.7)') 'p_h:',no_finish, p_h
        write(10,'(1a10, 1i10, 1e20.7)') 'p_max:',no_finish, p_max
        close(10)
     endif
  enddo
endsubroutine write_bin_entropy

!========================================================!
subroutine write_mag
!========================================================!
  use configuration
  implicit none
  logical:: file_opens
  integer:: i
  file_opens = .True.
   do while (file_opens)
      INQUIRE(FILE="out_mag_h0.dat", opened=file_opens)
      if (file_opens) then
        call sleep(1)
      else
         open(10,file='out_mag_h0.dat',status='unknown',position='append')
         do i=1, nn
           write(10,'(I5,2X,I5,2X,2f15.7)') no_finish, i, mag(i,1), mag(i,2)
         enddo
         close(10)
      endif
   enddo
endsubroutine write_mag

!========================================================!
subroutine write_entropy
!========================================================!
  use configuration
  implicit none
  logical:: file_opens
  integer:: i
  file_opens =.True.
  do while (file_opens)
     INQUIRE(FILE="out_entropy.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
        open(10,file='out_entropy.dat',status='unknown',position='append')
        write(10,'(1a10, 1i10, 1i40)') 'h_st:', no_finish, code_h
        write(10,'(1a10, 1i10, 1i40)') 'max_st:',no_finish, code_max
        write(10,'(1a10, 1i10, 2e20.7)') 'p_2:', no_finish, s_2_avg(1), s_2_avg(2)
        write(10,'(1a10, 1i10, 2e20.7)') 'p_AF1:', no_finish, p_af1_avg(1), p_af1_avg(2)
        write(10,'(1a10, 1i10, 2e20.7)') 'p_AF2:', no_finish, p_af2_avg(1), p_af2_avg(2)
        write(10,'(1a10, 1i10, 2e20.7)') 'p_h:',no_finish, p_h_avg(1), p_h_avg(2)
        write(10,'(1a10, 1i10, 2e20.7)') 'p_max:',no_finish, p_max_avg(1), p_max_avg(2)
        close(10)
     endif
  enddo
endsubroutine write_entropy

subroutine write_energy
  use configuration
  implicit none
  logical:: file_opens
  integer:: i
  file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_energy.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
        open(10,file='out_energy.dat',status='unknown',position='append')
        do i=1, nn+2
          write(10,'(I5,2X,I5,I5,2X,2f15.7)') no_finish, i, mode, energy_avg(i,1), energy_avg(i,2)
        enddo
        close(10)
     endif
  enddo

endsubroutine write_energy
