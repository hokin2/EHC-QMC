! https://github.com/jannisteunissen/ffhash
module hash_table
#define FFH_KEY_TYPE integer*16
#define FFH_VAL_TYPE integer*16
#include "ffhash_inc.f90"
end module hash_table

!===============================
module configuration
!===============================
use hash_table
save

! #1
! define lattice related
integer:: nx, ny
! for ensemble two region A
integer:: nxsub, nysub

! no of sites
integer:: nn
!number of bonds 1D: *1  2D: *2 (= 2*nn)
integer:: nb
! parameters for model
!heisenberg,disorder strength, inverse temp
real(8) :: jjz,hhz,beta
! magnetic field for every sites
real(8),allocatable :: hsites(:), hsites_save(:)
! new Hamiltonian H3, disorder part, J/Jz assumed to be 1
real(8),allocatable :: h3(:)
! no of disorder configs
integer(8) :: ndis
! array of sites connected by each bond,
integer, allocatable :: bond(:,:)
! nearest neighbor of bond, for counting energy
integer, allocatable  :: bond_nn(:,:)

! #2
! ensemble method related (PRB 90, 125105 (2014), PRL 112, 057203 (2014))
! 1- glued ensemble 2- independent ensemble
integer :: ensemble
! the number of replica, how many beta
! read it from input
integer :: qmax, qmax_save
! global variables (replace nlist_)
integer :: q_
!string length
integer :: l
! list length for each part, M in the paper
integer, allocatable:: l_len(:)
! the start of the list for different sections
integer, allocatable :: l_start(:)
!** steps spent in each ensemble (not using Humeniuk and Roscilde)
!integer :: nstep(2,nn)
! label for subregion
integer(8),allocatable :: label(:)

! #3
! QMC related
integer:: nbins, mstep, istep
integer:: nbins_map, mstep_map, istep_map
! basis state for qmax replica
integer(16), allocatable :: spn(:,:)
integer, allocatable :: opstring(:)
!random number generator
integer :: iir,jjr,kkr,nnr,mzran
! MPI ----
! no of finished disorders
integer :: no_finish
integer :: np, id
logical :: finished
! maping to excited state
! do excited state mapping or not (1=do)
integer :: map
! temporary number of copy (being 1 in mapping phase)
integer :: qmax_c
! epsilon for the renormalized energy
! d_h0 for the old disorder strength
! d_h3 for the new disorder strength
! eig_3 for the third lowest eigenvalue
real(8):: eps, d_h0, d_h3, eig_3
! max_e for the max energy of mapped Hamitonian
! min_e for the min energy of mapped Hamitonian
! h3_e for the relative energy
! sign_ising for the sign of Ising term
real(8):: max_e(2), min_e(2), h3_e(2), sign_ising
! status of simulation
integer:: mode

! #4
! Measurement related
!number of non identity op, nhpart along l direction, nhpart2 for two regions
!** energy calculation, dont know whether use it or not
! integer :: nh
integer,allocatable :: nh_q(:)
!** phase -1^bond number, related to stiffness(stagger)
!ph(0:nn)
!number of loop to apply in each step
integer :: nl
! count how many loop applied in equilibrium state
real(8) :: nloops
! count how many vertex involved in loopupdate in equilibrium state
real(8) :: lopers
! histogram
type(ffh_t),allocatable:: h_table(:)
! max number of record
integer :: mem_state
! new states generated
integer,allocatable :: no_states(:)
integer(16),allocatable :: record_spn(:,:)
! probability of matching
real(8) :: p_af1, p_af2, p_h, p_max
real(8) :: p_af1_avg(2), p_af2_avg(2), p_h_avg(2), p_max_avg(2)
! encode no_states
integer(16) :: code_af1, code_af2, code_h, code_max
! Renyi entropy
real(8) :: s_2, s_inf
real(8) :: s_2_avg(2)
! energy
real(8), allocatable :: energy(:), energy_save(:), energy_avg(:,:)
! energy*energy
real(8), allocatable :: energy2(:,:), energy3(:,:)
! Hamiltonian and Hamiltonian correlation
real(8), allocatable :: QCM(:,:), QCM_1(:,:), QCM_2(:,:), QCM_err(:,:)
! total energy
real(8) :: etot(2)
! magnetization
real(8), allocatable :: mag(:,:)


! #5
! vertex related
! number of vertex (6*nb)(for directed loop update=6*nb, for normal loop update=6)
integer:: nvx
!weight od adding/ deleting diagonal op
!real(8) :: awgt(-1:1,-1:1),dwgt(-1:1,-1:1)
real(8), allocatable:: awgt(:,:,:,:),dwgt(:,:,:,:)
!constant added to ham for positive definiteness
real(8), allocatable :: amax(:,:)
!vxoper takes in vertex no and returns op type,
integer, allocatable :: vxoper(:)
!vxcode takes in op type, spin state at btm site 1 and btm site 2 and returns vx no.
integer, allocatable :: vxcode(:,:,:,:)
!legs 2 vertex, uses the spin state at the 4 legs to return vx no.
integer, allocatable :: legvx(:,:,:,:,:,:)
! vertex 2 leg, takes in vextex leg number and vertex type to return spin state at that leg
integer, allocatable :: vxleg(:,:)
!equation table and number of eq
integer, allocatable :: etab(:,:,:),neqs(:)
!vertex weight, using awgt
real(8), allocatable:: vxw(:)
!for directed loop update
integer, allocatable:: vxbnd(:), vtyp(:)
! input (legin, leg-out, vertex no) return new vertex no
integer, allocatable :: vxnew(:,:,:)
! input (legin, leg-out, vertex no) return new vertex weight
real(8), allocatable :: vxprb(:,:,:)

end module configuration

module mpi_modules
  use configuration
  use mpi
  implicit none
!#include "mpif.h"
!===
contains
!===
  !==============================================================!
  subroutine read_parameters
  !==============================================================!
    implicit none
    integer :: ierr
    logical :: file_exists

    write(*,'(1a, 1i4.4,1a,i4.4,1a)') 'node-',id,' initialized among ',np,' nodes'

    if (id.eq.0) then
      open(10,file='read.in',status='old')
      read(10,*) nx, ny
      read(10,*) jjz, hhz, beta
      read(10,*) nbins, mstep, istep, ndis
      read(10,*) ensemble, qmax
      ! for ensemble 2
      read(10,*) nxsub, nysub
      ! for memory in hash
      read(10,*) mem_state
      ! for mapping to excited state or not (1=do)
      read(10,*) map
      ! qmc for mapping
      read(10,*) nbins_map, mstep_map, istep_map
      close(10)

    endif



    CALL MPI_BCAST(nx ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ny ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(jjz ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(hhz ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(beta ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(nbins ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(mstep ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(istep ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ndis ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(ensemble ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(qmax ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(nxsub ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nysub ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(map ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(mem_state ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(nbins_map ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(mstep_map ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(istep_map ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (id.eq.0) then
      INQUIRE(FILE="no_finish.dat", EXIST=file_exists)
      if (file_exists) then
         open(10, file='no_finish.dat',status='old')
         read(10,*) no_finish
         close(10)
      else
         open(10, file='no_finish.dat',status='new')
         write(10,*) 0
         close(10)
         no_finish=0
      endif
    endif
    CALL MPI_BCAST(no_finish ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    write(*,'(1a, 1i4.4,1a,i4.4,1a)') 'node-',id,' done read_parameters'
    write(*,'(1a, 1i4.4,1a,i4.4,1a,i4.4,1a)') 'node-',id,' done ',no_finish, '/',ndis,' disorders'
  end subroutine read_parameters

  !======================================================!
  subroutine init_ran
  !======================================================!
    use configuration;   implicit none
    integer :: is,js,ks,ls ,ierr
    real(8) ::    rndm
    if (id.eq.0) then
    open(10,file='rand.in',status='old')
    read(10,*)is
    read(10,*)js
    read(10,*)ks
    read(10,*)ls
    close(10)
    endif
    CALL MPI_BCAST(is ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(js ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ks ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ls ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    iir=1+abs(is)+id
    jjr=1+abs(js)
    kkr=1+abs(ks)
    nnr=ls
    if (id.eq.0) then
      open(10,file='rand.in',status='unknown')
      write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
      write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
      write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
      write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
      close(10)
    endif

  end subroutine init_ran
  !=======================================================!
  !==============================================================!
  subroutine init_vars
  !==============================================================!
    use configuration
    use hash_table
    implicit none
    integer :: i
    nn = nx*ny
    nb = 2*nn
    nvx = 6*nb

    ! #1
    ! define lattice related
    allocate(hsites(nn))
    allocate(hsites_save(nn))
    if (map.eq.1) allocate(h3(nn))
    ! #2
    ! ensemble method related (PRB 90, 125105 (2014), PRL 112, 057203 (2014))
    allocate(l_len(qmax),l_start(qmax))
    ! label for subregion
    allocate(label(nn))
    qmax_save = qmax

    ! #3
    ! QMC related
    ! basis state
    allocate(spn(qmax,nn))
    allocate(bond(2,nb))
    allocate(bond_nn(6,nb))

    ! #4
    ! Measurement related
    allocate(nh_q(qmax))
    allocate(h_table(qmax))
    allocate(no_states(qmax))
    no_states(:) = 0
    allocate(record_spn(qmax,mem_state))
    record_spn(:,:) = 0
    ! energy
    allocate(energy(nn+2), energy_save(nn+2), energy_avg(nn+2,2))
    ! energy*energy
    allocate(energy2(nn+2,nn+2),energy3(nn+1,nn+1))
    ! quantum covariance matrix
    allocate(QCM(nn+1,nn+1), QCM_2(nn+1,nn+1), QCM_err(nn+1,nn+1))
    ! magentization
    allocate(mag(nn,2))

    ! #5
    ! vertex related
    allocate(awgt(-1:1,-1:1, nb,0:1),dwgt(-1:1,-1:1, nb,0:1))
    !constant added to ham for positive definiteness
    allocate(amax(nb,0:1))
    !vxoper takes in vertex no and returns op type,
    allocate(vxoper(nvx))
    !vxcode takes in op type, spin state at btm site 1 and btm site 2 and returns vx no.
    allocate(vxcode(0:1,-1:1,-1:1,nb))
    !legs 2 vertex, uses the spin state at the 4 legs to return vx no.
    allocate(legvx(-1:1,-1:1,-1:1,-1:1,nb,0:0))
    ! vertex 2 leg, takes in vextex leg number and vertex type to return spin state at that leg
    allocate(vxleg(0:3,nvx))
    !equation table and number of eq
    allocate(etab(0:5,3,2*nvx),neqs(2*nvx))
    !vertex weight, using awgt
    allocate(vxw(nvx))
    !for directed loop update
    allocate(vxbnd(nvx), vtyp(nvx))
    ! input (legin, leg-out, vertex no) return new vertex no
    allocate(vxnew(0:3,0:3,nvx))
    ! input (legin, leg-out, vertex no) return new vertex weight
    allocate(vxprb(0:3,0:3,nvx))

  end subroutine init_vars
end module
