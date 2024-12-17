!===============================
module configuration
!===============================
save

integer, parameter :: nx = 16
integer, parameter :: ny = 1
integer, parameter :: nn = nx*ny
integer, parameter :: nb = nn!2*nn !number of bonds
integer, parameter :: nvx = 6 !6*nb ! number of vertex
integer, parameter :: ntau = 200 ! unequal green function

real(8) :: jjp,jjz,hhz,beta !heisenberg, inverse temp

integer :: l !string length
integer :: mloop !number of loops
integer :: nh !number of non identity op
integer :: spn(nn),bond(2,nb),ph(0:nn) !basis state, array of sites connected by each bond, phase -1^bond number, related to stiffness(stagger)
real(8) :: magsite(0:nn), magsite2(0:nn, 0:nn) !magnetization for each site
integer, allocatable :: opstring(:) !

!real(8) :: awgt(-1:1,-1:1),dwgt(-1:1,-1:1) !weight od adding/ deleting diagonal op
real(8) :: awgt(-1:1,-1:1),dwgt(-1:1,-1:1) !weight od adding/ deleting diagonal op
real(8) :: amax !constant added to ham for positive definiteness

integer :: vxoper(nvx),vxcode(0:1,-1:1,-1:1) !vxoper takes in vertex no and returns op type, vxcode takes in op type, spin state at btm site 1 and btm site 2 and returns vx no.
integer:: vxbnd(nvx), vtyp(nvx)
integer :: legvx(-1:1,-1:1,-1:1,-1:1) !legs 2 vertex, uses the spin state at the 4 legs to return vx no.
integer :: vxleg(0:3,nvx) ! vertex 2 leg, takes in vextex leg number and vertex type to return spin state at that leg
integer :: etab(0:5,3,2*nvx),neqs(2*nvx) !equation table and number of eq
real(8) :: vxw(nvx) !vertex weight, using awgt

integer :: vxnew(0:3,0:3,nvx)! input (legin, leg-out, vertex no) return new vertex no
real(8) :: vxprb(0:3,0:3,nvx)! input (legin, leg-out, vertex no) return new vertex weight

integer :: nl !number of loop
real(8) :: nloops!Number of maximum loops
real(8) :: lopers !size of loop

real(8) :: etot,avu,avk,umag,sxu,ssa,sxa,rho(2) !measurables

integer :: iir,jjr,kkr,nnr,mzran !random number generator
real(8) ::    rr

! Hamiltonian and Hamiltonian correlation
real(8) :: pattern1(2*nb)
real(8) :: pattern2(2*nb, 2*nb),pattern3(2*nb, 2*nb),pattern3_sym(2*nb, 2*nb)
real(8) :: QCM(2*nb, 2*nb)

integer:: nbins
real(8):: energy_rec(nb+2)
real(8):: p3_tmp(2*nb,2*nb)
integer:: modes

real(8):: accept(4),tries(4)
end  module configuration
!===============================







!=============================================!
program XXZ_2d_model
!=============================================!
  use configuration; implicit none

  integer :: i,j,mstep,istep
  
  real(8):: eigval(2*nn)
  real(8) :: eigvec(2*nn,2*nn)

  open(10,file='read.in',status='old')
  read(10,*)jjp,jjz,hhz,beta,nbins,istep,mstep
  !mstep is number of monte carlo measurements in each bin, istep is equilibration step
  write(*,'(1a,1i3,1a,1i3,1a,1i3)') 'L= ',nx,' * ',ny,'  nvx= ',nvx
  call initran
  


     !1. simulating original Hamiltonian H0 to get covariance matrix
     modes=1   
      
     call lattice
     call initconf(1) !initial state
     call pvect0
     call initvrtx

     lopers=0.d0
     nloops=0.d0
     mloop=100*l

     do i=1,istep
        if (mod(i,istep/10)==0) then
           !write(*,*) i,istep
        endif
        call mcstep
        call adjstl
        if (mod(i,istep/20)==0) then
           call adjnl
           !write(*,*)'l,nh,nl=',l,nh,nl
           call errchk
           !write(*,*)spn(:)
        endif
     enddo
     
     
     write(*,*)'l,nh,nl=',l,nh,nl
     write(*,*)'Equilibrium done'
     
     
     do i=1,nbins
        call zerodata
        do j=1,mstep
           call mcstep
           !write(*,*)'done mcstep',l,nh,nl
           if(mod(j,mstep/10)==0) then 
              call errchk
           endif
        enddo
        call writeres(mstep,1)
        
        write(*,*)'completed bin', i 
     enddo
     
     write(*,*)'Measurement done'
     write(*,*) 'Covariant matrix calculation done'
     
     deallocate(opstring)
     
     
	
     stop




     !2. Diagonalization of cavariant matrix - Real symmetric matrix
	 call diasym(QCM,eigval,2*nb)

	 open(15,file='eigvals_cov_mat.dat',status='unknown')
     do i=1,2*nb
         write(15,65) eigval(i) !i-th eigenvalue
     enddo
65   format(f15.7)
     close(15)
		

	 do i=1,2*nb
		 eigvec(:,i)=QCM(:,i)
	 enddo
   
     open(25,file='eigvecs_cov_mat.dat',status='unknown')
     do i=1,2*nb
         write(25,75) eigvec(:,i) !i-th eigenvector 
     enddo
75   format(f15.7)
     close(25)
 

end program XXZ_2d_model
!===================================!





!==============================================================!
subroutine lattice
!==============================================================!
  use configuration; implicit none

  integer :: i,j,is,ix,iy
  integer :: xy(2,nn),xyi(0:nx-1,0:ny-1)
  real(8) :: rndm

  is=0
  do j=0,ny-1
  do i=0,nx-1
     is=is+1
     xy(1,is)=i
     xy(2,is)=j
     xyi(i,j)=is
     ph(is)=(-1)**(i+j)
  enddo
  enddo

! 1D
  do i=1,nn
     ix=mod(xy(1,i)+1,nx)
     iy=xy(2,i)
     bond(1,i)=i
     bond(2,i)=xyi(ix,iy)
  enddo

end subroutine lattice



!==========================================================!
subroutine adjstl !adjust string length
!==========================================================!
  use configuration; implicit none

  integer :: i,p1,p2,dl,l1
  integer,allocatable :: tstring(:)
  real(8) :: r,rndm

  dl=nh/4
  if (nh+dl<l) return
  l1=l+dl

  allocate(tstring(l))
  tstring(:)=opstring(:)

  deallocate(opstring)
  allocate(opstring(l1))
  do i=1,l
     opstring(i)=tstring(i)
  enddo
  do i=l+1,l1
     opstring(i)=0
  enddo
  r=dble(nh)/dble(l1)
  do i=l1,1,-1
     if(opstring(i) /= 0) then
        p1=i
        exit
     endif
  enddo
  outer: do p2=l1,1,-1
     if (p2 <= p1) then
        exit outer
     elseif (rndm() <= r) then
        opstring(p2)=opstring(p1)
        opstring(p1)=0
        innner: do i=p1-1,1,-1
           if (opstring(i) /= 0) then
              p1=i
              cycle outer
           endif
        enddo innner
        exit outer
     endif
  enddo outer

  l=l1
  mloop=100*l

end subroutine adjstl

!=============================================================!
subroutine errchk
!============================================================!
  use configuration; implicit none

  integer :: sconf(nn),b,o,op,i!,s1,s2
    
  sconf(:)=spn(:)
  
  do i=1,l
     op=opstring(i)
     b=op/2; o=mod(op,2)      
     if (op.eq.0) cycle
     if (b.lt.1.or.b.gt.nb) then
        open(12,file='log.txt',status='unknown',position='append')
        write(12,*)'wrong bond: ',b
        close(12)
        stop
     endif
     if(o == 1) then
        sconf(bond(1,b))=-sconf(bond(1,b))
        sconf(bond(1,b))=-sconf(bond(1,b))
     endif
  enddo
  do i=1,nn
     if (spn(i).ne.sconf(i)) then
        open(12,file='log.txt',status='unknown',position='append')
        write(12,*)'wrong propagation of spin state'
        close(12)
        stop
     endif
  enddo

end subroutine errchk

!===========================================================!
subroutine pvect0
!===========================================================!
  use configuration; implicit none

  integer :: s1,s2
  real(8) :: z

  z=dble(2*2) !coordination number
  amax=0.d0

  do s1=-1,1,2
  do s2=-1,1,2
     awgt(s1,s2)=jjz*0.25*s1*s2 !heisenberg term
     awgt(s1,s2)=awgt(s1,s2)-(1/z)*hhz*0.5*(s1+s2) !zeeman term 0.5 is for spin half
     if (amax < awgt(s1,s2)) amax=awgt(s1,s2)
  enddo
  enddo

  amax=amax+1.d0! amax maximum adding weight
  do s1=-1,1,2
  do s2=-1,1,2
     awgt(s1,s2)=amax-awgt(s1,s2)
     dwgt(s1,s2)=1.d0/awgt(s1,s2)
  enddo
  enddo

end subroutine pvect0
!==========================================================!


!===========================!
subroutine writeres(nmsr,mode)
!===========================!
  use configuration; implicit none

  integer :: i,nmsr,j,i2,j2, mode,k
  integer :: k1=(nb/2)+1

  real(8) :: denom, KE2, PE2, KE_PE, PE_KE
  real(8) :: vmax(2*nn), unit(2*nn)
  
  !to symmetrise the covariance matrix
  real(8) :: sum_pattern1_1,sum_pattern1_2, avg_pattern1_1,avg_pattern1_2
  real(8) :: sym_pattern3(3*nb)
  integer :: sym_counter(3*nb) 


  avu=avu/(dble(nn)*beta*dble(nmsr))
  avu=-avu+amax
  avk=-avk/(dble(nn)*beta*dble(nmsr))
  
  denom=(0.5*beta)**2 
      

     KE2=0.0
     PE2=0.0
	 do i=(1),nb
		 KE2 = KE2 + pattern1(i)
         PE2 = PE2 + pattern1(nb+i)
     enddo
     
     PE2=PE2/(dble(nn)*beta*dble(nmsr))
	 PE2=-PE2+amax
	 KE2=-KE2/(dble(nn)*beta*dble(nmsr))
	 
	
	sum_pattern1_1 = 0.d0
	sum_pattern1_2 = 0.d0
	 do i=1,nb
         pattern1(i) = - 1.0*( pattern1(i)/(beta*dble(nmsr)) )
         pattern1(nb+i) = amax - ( pattern1(nb+i)/(beta*dble(nmsr))) 
         
         sum_pattern1_1 = sum_pattern1_1 + pattern1(i) 
         sum_pattern1_2 = sum_pattern1_2 + pattern1(nb+i)
     enddo
	 avg_pattern1_1 = sum_pattern1_1/nb
	 avg_pattern1_2 = sum_pattern1_2/nb


	 do i=1,nb
         pattern1(i)    = avg_pattern1_1
         pattern1(nb+i) = avg_pattern1_2
     enddo

	 open(28,file='pattern1.dat',status='unknown',position='append')
     do i=1,2*nb
         write(28,83) pattern1(i)
     enddo
83   format(f15.7)     
     close(28)

if (mode==1) then
     open(10,file='enr.dat',status='unknown',position='append')!position='append'
     write(10,20)avu,avk, PE2, KE2 !, (PE2+KE2), sum_pattern  !total energy, average potential, average kinetic
20   format(4f18.10)
     close(10)

          
     
     
open(27,file='pattern2.dat',status='unknown',position='append')
!73   format(3f15.7)     
!     close(27)     
     
     
     do i=1,nb 
        do j=1,nb
           write(27,73) (pattern1(i)*pattern1(j))  
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb           
           write(27,73) (pattern1(i)*pattern1(j+nb)) 
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb           
           write(27,73) (pattern1(i+nb)*pattern1(j)) 
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb           
           write(27,73) (pattern1(i+nb)*pattern1(j+nb)) 
        enddo
     enddo 
     
73   format(f15.10)     
     close(27)  
     
     
     
     
     !NEW
     do i=1,nb
        do j=1,nb
        
			!!upper diagonal block
			pattern3(i,j) = (pattern2(i,j)/(nb*denom*dble(nmsr)))*(nb/4.0) 
			
			!!left off-diag block
			pattern3(i,j+nb) = -(-amax*pattern1(i) - (pattern2(i,j+nb)/(nb*denom*dble(nmsr)))*(nb/4.0) )
						
			!!right off-diag block
			pattern3(i+nb,j) = -(-amax*pattern1(j) - (pattern2(i+nb,j)/(nb*denom*dble(nmsr)))*(nb/4.0) ) 
			
			!!lower digonal  block
			pattern3(i+nb,j+nb) = (pattern2(i+nb,j+nb)/(nb*denom*dble(nmsr)))*(nb/4.0) &
			                      +amax*pattern1(i+nb) + amax*pattern1(j+nb) - amax*amax						  
			
        enddo
     enddo 
      
      
      
     !symmetrising pattern3 
     sym_pattern3(:) = 0.d0 
     sym_counter(:)  = 0
    
     do i=1,nb   
        do j=1,(nb+1-i)
			!symmetrise upper diagonal block
			sym_pattern3(i) = sym_pattern3(i) + pattern3(j+i-1,j)
			sym_counter(i)  = sym_counter(i) + 1  
			
			!symmetrise upper off-diagonal block
			sym_pattern3(i+nb) = sym_pattern3(i+nb) + pattern3(j+i+nb-1,j)
			sym_counter(i+nb)  = sym_counter(i+nb) + 1  
			
			!symmetrise lower diagonal block
			sym_pattern3(i+2*nb) = sym_pattern3(i+2*nb) + pattern3(j+i+nb-1,j+nb)  
			sym_counter(i+2*nb)  = sym_counter(i+2*nb) + 1 					
        enddo
        sym_pattern3(i) = sym_pattern3(i)/sym_counter(i)        
        sym_pattern3(i+nb) = sym_pattern3(i+nb)/sym_counter(i+nb)        
        sym_pattern3(i+2*nb) = sym_pattern3(i+2*nb)/sym_counter(i+2*nb)
     enddo 
    
    
    do i=1,nb 
        do j=i,nb
        
            k = mod(j-i,nb)+1
            
            if (k>k1) then
			    k= k1 - (k-k1)
		    endif
            
			pattern3(i,j) = sym_pattern3(k)        
            pattern3(j,i) = pattern3(i,j)
            
          
            pattern3(i+nb,j) = sym_pattern3(k+nb)       
            pattern3(j+nb,i) = pattern3(i+nb,j)
            
            pattern3(i+nb,j+nb) = sym_pattern3(k+2*nb)        
            pattern3(j+nb,i+nb) = pattern3(i+nb,j+nb)
            
        enddo
     enddo
     
     
     
     do i=1,nb 
        do j=1,nb
            pattern3(i,j+nb) = pattern3(i+nb,j)
        enddo
     enddo
    
     
open(37,file='pattern3.dat',status='unknown',position='append')
         
     do i=1,nb 
        do j=1,nb
           write(37,103) pattern3(i,j)  
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb
           write(37,103) pattern3(i,j+nb)       
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb
           write(37,103) pattern3(i+nb,j)   
        enddo
     enddo
     
     do i=1,nb 
        do j=1,nb
           write(37,103) pattern3(i+nb,j+nb)        
        enddo
     enddo
     
103   format(f15.10)     
     close(37)
      
       
      
      
      
     do i=1,nb
        do j=1,nb
        
			!!upper diagonal block
			QCM(i,j) = pattern3(i,j) - (pattern1(i)*pattern1(j))  
			
			!!left off-diag block
			QCM(i,j+nb) = pattern3(i,j+nb) - (pattern1(i)*pattern1(j+nb)) 
			
			!!right off-diag block
			!QCM(i+nb,j) = pattern3(i+nb,j) - (pattern1(i+nb)*pattern1(j)) 
			QCM(i+nb,j) = QCM(i,j+nb)
			
			!!lower digonal  block
			QCM(i+nb,j+nb) = pattern3(i+nb,j+nb) - (pattern1(i+nb)*pattern1(j+nb))			   
	
        enddo
     enddo
        

     open(23,file='cov_mat.dat',status='unknown',position='append')
     do i=1,2*nb
        do j=1,2*nb
           write(23,63) QCM(i,j)
        enddo
     enddo
63   format(f15.7)     
     close(23)

  
  endif

end subroutine writeres
!=======================!





!==================!
subroutine zerodata
!==================!
  use configuration; implicit none

  etot=0.d0
  avu=0.d0
  avk=0.d0
  umag=0.d0
  sxu=0.d0
  ssa=0.d0
  sxa=0.d0
  rho(:)=0.d0 
  magsite(:) = 0.d0
  magsite2(:,:) = 0.d0
  pattern1(:)=0.0d0
  pattern2(:,:)=0.0d0
  pattern3(:,:)=0.0d0
  QCM(:,:)=0.0d0
  pattern3_sym(:,:)=0.0d0
end subroutine zerodata
!======================!
!======================================!
subroutine initconf(mode)
!======================================!
  use configuration; implicit none
  
  integer :: i,mode
  do i=1,nn
     spn(i)=(-1)**i
  enddo

  l=20
  allocate(opstring(l))
  opstring(:)=0
  nh=0
  nl=5
  tries=0.0
  accept=0.0
end subroutine initconf
!=====================================!

!======================================================!
subroutine initran
!======================================================!
  use configuration;   implicit none

  integer :: is,js,ks,ls
  real(8) ::    rndm

  open(10,file='rand.in',status='old')
  read(10,*)is
  read(10,*)js
  read(10,*)ks
  read(10,*)ls
  close(10)
  iir=1+abs(is)
  jjr=1+abs(js)
  kkr=1+abs(ks)
  nnr=ls
  open(10,file='rand.in',status='unknown')
  write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
  write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
  write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
  write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
  close(10)

end subroutine initran
!=======================================================!

!======================================================!
real(8) function rndm()
!======================================================!
  use configuration; implicit none

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


