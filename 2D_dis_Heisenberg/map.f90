! for the excited state mapping

subroutine init_map(map_success)
  use configuration
  implicit none
  logical :: map_success
  qmax = 1
  hsites_save(:) = hsites(:)
  call find_H3(map_success)
  jjz= 1.0d0
  hsites(:) =h3(:)
  mode = 2
  deallocate(opstring)
endsubroutine init_map

subroutine flip_map
  use configuration
  implicit none
  min_e(:) = etot(:) ! record the e_min
  hsites(:) = -h3(:)
  jjz=-1.0d0
  mode = 3
  deallocate(opstring)
endsubroutine flip_map

subroutine end_map
  use configuration
  implicit none
  max_e(1) = -etot(1) ! record the e_max
  max_e(2) = etot(2)
  call find_eps
  write(*,30) 'node-',id,' dis-',no_finish,' finished'
  30 format(1a,1i4.4,1a,1i5.5, 1a)
endsubroutine end_map

!========================================================!
function find_dtyp(h_in,n_in) result(tmp)
!========================================================!
  use configuration
  implicit none
  integer:: n_in, i
  real(8):: h_in(n_in)
  real(8):: tmp
  tmp=0.0d0
  do i=1,n_in-1
     tmp = tmp+log(abs(h_in(i)-h_in(i+1)))
  enddo
  !tmp = tmp + abs(h_in(n_in)-h_in(1))
  tmp = tmp +log(abs(h_in(n_in)-h_in(1)))
  tmp = tmp/dble(n_in)
  tmp = exp(tmp)
endfunction find_dtyp
!========================================================!

!========================================================!
subroutine find_H3(map_success)
!========================================================!
  ! find the new Hamiltonian using covariance matrix
  use configuration
  implicit none
  integer:: i,j,info,n_space,lwork,maxoverlap1(1),maxoverlap2(1)
  real(8):: v
  !real(8):: tmpmat(nb+1,nb+1)
  real(8):: tmpval(nn+1)
  real(8):: wr(nn+1),wi(nn+1),vl(nn+1,nn+1),vr(nn+1,nn+1),work(4*(nn+1))
  real(8):: d_min,h3_tmp(nn),find_dtyp
  !real(8),allocatable:: d_rec(:)
  real(8):: trivial1(nn+1), trivial2(nn+1), dotprod(nn+1)
  logical :: found, file_opens, map_success
  real(8):: QCM_err_sum

  trivial1(1) = 0.0d0
  trivial1(2:nn+1) = 1.0d0

  trivial2(1) = 1.0d0
  trivial2(2:nn+1) = hsites_save(:)
  !1. average out different bins of covariance matrix:

  QCM_err_sum = 0
  QCM(:,:) = QCM(:,:)/dble(nbins)
  do i=1,nn+1
     do j=1,nn+1
        ! s.e. = (<x^2> - <x>^2)/sqrt(nbin-1)
        if (nbins==1) then
          QCM_err(j,i) = 0.0d0
        else
          QCM_err(j,i) = sqrt(QCM_2(j,i)/dble(nbins) - QCM(j,i)*QCM(j,i))/sqrt(dble(nbins-1))
        endif
        QCM_err_sum = QCM_err_sum + QCM_err(j,i)
      enddo
   enddo

   file_opens = .True.
   do while (file_opens)
      INQUIRE(FILE="out_cov_mat.dat", opened=file_opens)
      if (file_opens) then
        call sleep(1)
      else
        open(17,file='out_cov_mat.dat',status='unknown',position='append')
       do i=1,nn+1
          do j=1,nn+1
            write(17,77) no_finish, QCM(j,i), QCM_err(j,i)
          end do
        end do
          close(17)
          77    format(1i5, 2f15.7)
       endif
     end do

     !   open(17, file='cov_err.dat',status='unknown', position='append')
     !   write(17,'(f15.7,2X,f15.7)') QCM_err_sum, QCM_err_sum/dble((nn+1)*(nn+1))
     !   close(17)

  !2. Diagonalization
  call diasym(QCM,tmpval,nn+1)

  do i=1,nn+1
     vr(:,i)=QCM(:,i)
  enddo

file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_eigvals_cov_mat.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
      open(10,file='out_eigvals_cov_mat.dat',status='unknown',position='append')
      do i=1,nn+1
        write(10,70) no_finish, i, tmpval(i) !i-th eigenvalue after sorting
      enddo
70    format(I5,2X,I5,2X,f15.7)
      close(10)
    endif
  end do

  file_opens = .True.
    do while (file_opens)
       INQUIRE(FILE="out_eigvecs_cov_mat.dat", opened=file_opens)
       if (file_opens) then
         call sleep(1)
       else
        open(15,file='out_eigvecs_cov_mat.dat',status='unknown',position='append')
        do i=1,nn+1
              write(15,*) no_finish, i
              write(15,75) vr(:,i) !i-th eigenvector after sorting
        enddo
      75    format(f15.7)
            close(15)
        endif
    enddo

  do i=1,nn+1
    dotprod(i) = abs(sum(vr(:,i)*trivial1))
  enddo
  maxoverlap1 = maxloc(dotprod)

  file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_overlap.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
      open(10,file='out_overlap.dat',status='unknown',position='append')
      do i=1,nn+1
          write(10,61) no_finish,1,i, dotprod(i) ! eigevec(i) overlap with trivial1 eigenvector
      enddo
      close(10)
    endif
  enddo

  do i=1,nn+1
    if (i .ne. maxoverlap1(1)) then
      dotprod(i) = abs(sum(vr(:,i)*trivial2))
    else
      dotprod(i) = 0
    endif
  enddo
  maxoverlap2 = maxloc(dotprod)

  file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_overlap.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
      open(10,file='out_overlap.dat',status='unknown',position='append')
      do i=1,nn+1
          write(10,61) no_finish,2,i, dotprod(i) ! eigevec(i) overlap with trivial1 eigenvector
      enddo
      close(10)
    endif
  enddo
61    format(I5,2X,I5,2X,I5,2X,f15.7)

  ! select eigenvector as new disorder
  found=.False.
  i = 0
  do while (.not. found)
    i = i+1
    if ((i .ne. maxoverlap1(1)) .and. (i .ne. maxoverlap2(1))) then
        found = .True.
    endif
    if (found) then
        sign_ising = vr(1,i)
        h3 = vr(2:nn+1,i)/vr(1,i)
        d_h3 = find_dtyp(h3,nn)
        eig_3 = tmpval(i)
    endif
  enddo
  d_h0 = find_dtyp(hsites,nn)

  file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_3rd_eigval_cov_mat.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
      open(10,file='out_3rd_eigval_cov_mat.dat',status='unknown',position='append')
          write(10,72) no_finish, tmpval(3), maxoverlap1, maxoverlap2, eig_3 !3rd eigenvalue
    72    format(I5,2X,f15.7,2X,I5,2X,I5,2X,f15.7)
          close(10)
     endif
  enddo

  file_opens = .True.
  do while (file_opens)
     INQUIRE(FILE="out_hp.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
        open(10,file='out_hp.dat',status='unknown',position='append')
        do i=1,nn
          write(10,'(1i10, 1f10.5)') no_finish, h3(i)
        enddo
        close(10)
     endif
  enddo

  write(*,'(1a,1i4.4,1a,1i5.5,1a,1f8.3,1a,1f8.3)') 'node-',id,' dis-',no_finish,' from d=',d_h0 ,' map to d=',d_h3
  map_success = .True.
  if (d_h3>100.0) then
     write(*,'(1a,1i4.4,1a,1i5.5,1a,1f8.3)') 'node-',id,' dis-',no_finish,' exceed 201, h_avg=',d_h3
     file_opens=.True.
     do while (file_opens)
        INQUIRE(FILE="out_eps.dat", opened=file_opens)
        if (file_opens) then
          call sleep(1)
        else
           open(10,file='out_eps.dat',status='unknown',position='append')
           write(10,'(I5,2X,f15.7,2X,f15.7,2X,f15.7,2X,f15.7)') no_finish,-1.0,eig_3,sign_ising,d_h3
           close(10)
        endif
     enddo
     map_success = .False.
  endif
end subroutine find_H3
!========================================================!

!========================================================!
 subroutine diasym(a,eig,n)
!========================================================!
!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!========================================================!


!========================================================!
subroutine find_eps
!========================================================!
use configuration
implicit none
integer:: i
logical:: file_opens
energy_save(:) = energy_save(:)/dble(nbins)
! h3_e=0.0d0
do i=1,nn+2  ! 1: J 2: Jz 3to end:mag fields
   if (i>2) then
      h3_e(1) = h3_e(1) + energy_save(i)*h3(i-2)/hsites_save(i-2)
   else
      h3_e(1) = h3_e(1) + energy_save(i)
   endif
enddo
h3_e(1) = h3_e(1)/dble(nn)
file_opens=.True.
do while (file_opens)
   INQUIRE(FILE="out_h3e.dat", opened=file_opens)
   if (file_opens) then
     call sleep(1)
   else
      open(10,file='out_h3e.dat',status='unknown',position='append')
      write(10,20) no_finish, min_e(1),max_e(1),h3_e(1)
20    format(I5,2X,f15.7,2X,f15.7,2X,f15.7)
      close(10)
   endif
enddo

eps= (h3_e(1)-min_e(1))/(max_e(1)-min_e(1))

file_opens=.True.
do while (file_opens)
   INQUIRE(FILE="out_eps.dat", opened=file_opens)
   if (file_opens) then
     call sleep(1)
   else
      open(10,file='out_eps.dat',status='unknown',position='append')
      !write(10,'(I5,2X,f15.7,2X,f15.7,2X,f15.7,2X,f15.7,2X,I5,2X,I5)') no_finish,eps,eig_3,sign_ising,d_h3,maxoverlap1,maxoverlap2
      write(10,'(I5,2X,f15.7,2X,f15.7,2X,f15.7,2X,f15.7)') no_finish,eps,eig_3,sign_ising,d_h3
      close(10)
   endif
enddo
end subroutine find_eps
!========================================================!
