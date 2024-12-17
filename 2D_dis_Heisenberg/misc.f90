
!======================================================!
real(8) function rndm()
!======================================================!
  use configuration
  implicit none

  mzran=iir-kkr
  if (mzran.lt.0) mzran=mzran+2147483579
  iir=jjr
  jjr=kkr
  kkr=mzran
  nnr=69069*nnr+1013904243
  mzran=mzran+nnr
  rndm=.5+.23283064e-9*mzran
  return
end function rndm
!========================================================!

!======================================!
subroutine init_conf
!======================================!
  use configuration
  implicit none
  integer :: i, q, l_i, l_len_init
  !initialize spin as aligned with magnetic fields
  do i=1,nn
    if (hsites(i)>0) then
       spn(:,i)=-1
    else
       spn(:,i)=1
    endif
  enddo
  ! every replica with l_len_init steps and 1 gap
  l_len_init = 10
  l = l_len_init*qmax + qmax -1
  l_len(:) = l_len_init
  l_i = 1
  do q = 1, qmax
    l_start(q) = l_i
    l_i = l_i + l_len(q)
    if (q<qmax) l_i = l_i +1
  enddo
  if (l_i-1 .ne. l) then
    write(*,*) 'sth wrong', l_i-1, l, l_start
    stop
  endif
  allocate(opstring(l))
  opstring(:)=0
  nh_q(:) = 0

  nl=5
  lopers=0.d0
  nloops=0.d0

end subroutine init_conf

subroutine init_lattice
  use configuration
  implicit none

  integer :: i,j,is,ix,iy
  integer :: xy(2,nn),xyi(0:nx-1,0:ny-1), nnxy(4,nn)

  if (id==0) then
     open(10, file='latt_info.dat',status='unknown')
     write(10,*) 'nn sites'
  endif
  is=0
  do j=0,ny-1
     do i=0,nx-1
        is=is+1
        xy(1,is)=i
        xy(2,is)=j
        xyi(i,j)=is
        nnxy(1,is) = mod(is+nx-1,nn)+1
        nnxy(2,is) = is+1
        if (mod(nnxy(2,is),nx)==1) then
           nnxy(2,is) = nnxy(2,is)-nx
        endif
        nnxy(3,is) = is-nx
        if (nnxy(3,is)<1) then
           nnxy(3,is)=nnxy(3,is)+nn
        endif
        nnxy(4,is) = is-1
        if (mod(nnxy(4,is),nx)==0) then
           nnxy(4,is)=nnxy(4,is)+nx
        endif
        if (id==0) then
           write(10,*) 'is,nnxy',is,nnxy(:,is)
        endif
     enddo
  enddo

  ! two dimensions
  ! for calculation of energy by counting bonds
  if (id==0) then
     write(10,*) 'nn bonds'
  endif
  do i=1,nn
     ix=mod(xy(1,i)+1,nx)
     iy=xy(2,i)
     ! horizontal
     bond(1,i)=i
     bond(2,i)=xyi(ix,iy)
     bond_nn(1,i) = i+nn
     bond_nn(2,i) = nnxy(4,i)
     bond_nn(3,i) = nnxy(3,i)+nn
     bond_nn(4,i) = nnxy(3,xyi(ix,iy))+nn
     bond_nn(5,i) = xyi(ix,iy)
     bond_nn(6,i) = xyi(ix,iy)+nn

     ix=xy(1,i)
     iy=mod(xy(2,i)+1,ny)
     ! vertical
     bond(1,nn+i)=i
     bond(2,nn+i)=xyi(ix,iy)
     bond_nn(1,nn+i) = nnxy(3,i)+nn
     bond_nn(2,nn+i) = i
     bond_nn(3,nn+i) = xyi(ix,iy)
     bond_nn(4,nn+i) = xyi(ix,iy)+nn
     bond_nn(5,nn+i) = nnxy(4,xyi(ix,iy))
     bond_nn(6,nn+i) = nnxy(4,i)
     if (id==0) then
        write(10,*) 'is,bond_nn',i,bond_nn(:,i)
        write(10,*) 'is,bond_nn', i+nn,bond_nn(:,i+nn)
     endif
  enddo
  if (id==0) then
     close(10)
  endif
  write(*,'(1a, 1i4.4,1a,i4.4,1a)') 'node-',id,' done init_lattice'
end subroutine init_lattice

subroutine init_disorder
  use configuration
  implicit none
  integer :: i
  real(8) ::  rndm
  logical :: file_opens
  file_opens =.True.
  do while (file_opens)
     INQUIRE(FILE="no_finish.dat", opened=file_opens)
     ! if (file_opens) then
     !    call sleep(1)
     ! else
     call sleep(id)
     if (.NOT. file_opens) then
        open(10, file='no_finish.dat',status='old')
        read(10,*) no_finish
        close(10)
        no_finish=no_finish+1
        if (no_finish<=ndis) then
          open(10, file='no_finish.dat',status='unknown',action="write")
          write(10,*) no_finish
          write(*,51) 'node-',id,' dis-',no_finish,' started'
51 format(1a,1i4.4,1a,1i5.5,1a)

          close(10)
          finished = .False.
        else
          finished = .True.
        endif
     endif
  enddo

  ! disorder magnetic fields
  hsites=0.0d0
  do i=1,nn
     hsites(i) = (rndm()-0.5)*2.0*hhz
     ! read(15,*) hsites(i)
  enddo
  file_opens =.True.
  do while (file_opens)
     INQUIRE(FILE="out_h0.dat", opened=file_opens)
     if (file_opens) then
       call sleep(1)
     else
        open(10,file='out_h0.dat',status='unknown',position='append')
        do i=1,nn
          write(10,'(1i10, 1f10.5)') no_finish, hsites(i)
        enddo
        close(10)
     endif
  enddo
  qmax = qmax_save

  mode = 1
end subroutine init_disorder

!===========================================================!
subroutine init_weights
!===========================================================!
  use configuration
  implicit none

  integer :: s1,s2, i ,op, b

  amax=0.d0

  do op = 0,1
     do b = 1, nb
        do s1=-1,1,2
           do s2=-1,1,2
              if (op==0) then
                 awgt(s1,s2, b,op)=jjz*0.25*s1*s2+0.125*(s1*hsites(bond(1,b))+s2*hsites(bond(2,b))) !heisenberg term
              else
                 awgt(s1,s2, b,op)=0.5
              endif
              if (amax(b,op) < awgt(s1,s2, b,op)) amax(b,op)=awgt(s1,s2, b,op)
           enddo
        enddo
     enddo
     amax(:,0)=amax(:,0)+0.1! amax maximum adding weight
     do b = 1, nb
        do s1=-1,1,2
           do s2=-1,1,2
              awgt(s1,s2, b,op)=amax(b,op)-awgt(s1,s2, b,op)
              if (awgt(s1,s2,b,op)<1e-7) then
                 dwgt(s1,s2,b,op)=1e7
              else
                 dwgt(s1,s2, b,op)=1.d0/awgt(s1,s2, b,op)
              endif
           enddo
        enddo
     enddo
  enddo
  if (id==0) then
  open(10,file='wgt.dat',status='unknown',position='append')
  do b=1,nb
     do op=0,1
        write(10,*) 'b=',b,' op=',op
        write(10,'(9 f15.3)') awgt(:,:,b,op)
        write(10,'(9 f15.3)') dwgt(:,:,b,op)
     enddo
  enddo
  close(10)
  endif
end subroutine init_weights

module convert_state
contains
function encode(conf) result(key)
  use configuration
  implicit none
  integer(16):: i, conf(nn),key
  key = 0
  do i=1, nn
    key = key + (conf(i)+1)/2*2**(i-1)
  enddo
endfunction encode

function decode(key) result(conf)
  use configuration
  implicit none
  integer ::  i
  integer(16):: conf(nn), key, key2
  conf(:) = 0
  key2 = key
  do i = 1, nn
    conf(i) = mod(key2,2)
    key2 =key2/2
  enddo
endfunction decode
! Program test
!   use convert_state
!   use configuration
!   implicit none
!   integer, allocatable:: conf(:)
!
!   nn = 16
!   allocate(conf(nn))
!   conf = decode(42405)
!   write(*,*) conf
!   write(*,*) encode(conf)
!   deallocate(conf)
! endprogram test

endmodule convert_state

subroutine encode_states
  use configuration
  use convert_state
  implicit none
  integer:: i, j, is
  integer(16):: st_af1(nn), st_af2(nn), st_h(nn)

  is = 0
  do j=1,ny
    do i=1,nx
     is=is+1
     st_af1(is) = ((-1)**(i+j)+1)/2
     st_af2(is) = ((-1)**(i+j+1)+1)/2
     if (hsites(is)>0) then
       st_h(is) = 0
     else
       st_h(is) = 1
     endif
    end do
  end do
  code_af1 = encode(st_af1)
  code_af2 = encode(st_af2)
  code_h = encode(st_h)
endsubroutine encode_states
