!=====================================!
subroutine mcstep
!=====================================!
  use configuration; implicit none
  
  integer :: passed

10 call dupdate

   call updloop(passed)
   
   if (passed==0) then
     !write(*,*) 'not passed'
     goto 10
   endif
  !write(*,*) 'passed'

   call meas_CM
end subroutine mcstep



!=======================================!
subroutine meas_CM
!=======================================!
  use configuration; implicit none

  integer :: i,j,o,op,d,b,b1,b2,nbs,s1,s2,ss1,ss2,nh1,nms
  integer :: last(1:nb),bsum(1:nb), bloczz(1:nb), blocxx(1:nb)
  integer :: bcorxx(1:nb,1:nb), bcorxz(1:nb,1:nb)
  integer :: bcorzx(1:nb,1:nb), bcorzz(1:nb,1:nb)
  
  real(8) :: scc(0:nn/2),dcc(0:nn/2),ssum(1:nn),dsum(1:nn)
  real(8) :: sc1(0:nn/2),dc1(0:nn/2),bc1(0:nn/2)
  real(8) :: sx1(0:nn/2),dx1(0:nn/2),bx1(0:nn/2)
  real(8) :: dl1,dl2,nth1,nth2

  bcorxx(:,:)=0
  bcorxz(:,:)=0
  bcorzx(:,:)=0
  bcorzz(:,:)=0
  
  
  bloczz(:)=0
  blocxx(:)=0

  b1=0
  b2=0
  nh1=0
  nh=0


  do i=1,l
     op=opstring(i)

     if (op.ne.0) then
     
        nh = nh+1   ! prefac for 4-spin correlator
        
        s1=op/2
        o=mod(op,2)
         
        if (o.eq.0) then       ! counting non-empty string
        
	       bloczz(s1)= bloczz(s1) +1
        
           if (b1.ne.0) then ! last op was diagonal
              bcorzz(s1,b1)=bcorzz(s1,b1)+1
           elseif (b2.ne.0) then
              bcorzx(s1,b2)=bcorzx(s1,b2)+1
           else
              if (nh1.ne.0) then
                 write(*,*)'something wrong. either b1 or b2 has to be non-zero - 0', b1, b2
                 stop
              endif
           endif
              
           nh1=nh1+1
           b1=s1
           b2=0
           
           !write(*,*) s1, b1
        else
        
		   blocxx(s1)= blocxx(s1) +1
           
           if (b1.ne.0) then
              bcorxz(s1,b1)=bcorxz(s1,b1)+1
           elseif (b2.ne.0) then
              bcorxx(s1,b2)=bcorxx(s1,b2)+1
           else
              if (nh1.ne.0) then
                 write(*,*)'something wrong. either b1 or b2 has to be non-zero - 1', b1, b2
                 stop
              endif
           endif
           
           nh1=nh1+1    
           b1=0
           b2=s1
           
         endif
      endif
   enddo
   
   
	! count the contribution across boundary
	 i=1
	 if(opstring(i).ne.0) then
		op=opstring(i)
		s1=opstring(i)/2
		o=mod(op,2)
		if (o.eq.0) then
				
		  if (b1.ne.0) then ! last op was diagonal
			 bcorzz(s1,b1)=bcorzz(s1,b1)+1
		  elseif (b2.ne.0) then
			 bcorxz(s1,b2)=bcorxz(s1,b2)+1
		  else
			 if (nh1.ne.0) then
				write(*,*)'something wrong. either b1 or b2 has to be non-zero -2' , b1, b2
				stop
			 endif
		  endif
		  		  
	    else
	    	    
		  if (b1.ne.0) then
			 bcorzx(s1,b1)=bcorzx(s1,b1)+1
		  elseif (b2.ne.0) then
			 bcorxx(s1,b2)=bcorxx(s1,b2)+1
		  else
			 if (nh1.ne.0) then
				write(*,*)'something wrong. either b1 or b2 has to be non-zero -3' , b1, b2
				stop
			 endif
		  endif
		   
		endif
	 endif

   
   
   do i=1, nb
	   pattern1(i)       = pattern1(i)    + blocxx(i)
	   pattern1(i+nb)    = pattern1(i+nb) + bloczz(i)
   enddo
   
  
   do i=1, nb         
      do j=1, nb        
         pattern2(i,j)       = pattern2(i,j)       + (nh-1)*bcorxx(i,j) 
         pattern2(i+nb,j)    = pattern2(i+nb,j)    + (nh-1)*bcorxz(i,j) 
         pattern2(i,j+nb)    = pattern2(i,j+nb)    + (nh-1)*bcorzx(i,j) 
         pattern2(i+nb,j+nb) = pattern2(i+nb,j+nb) + (nh-1)*bcorzz(i,j)
      enddo
   enddo

end subroutine meas_CM
!=====================================!




!=======================================!
subroutine dupdate
!=======================================!
  use configuration; implicit none

  integer :: i,o,op,b,s1,s2,ss1,ss2,nh1,last,nu,nk!,j
  integer :: t1,t2,cur(2),su,sa,nms,ph1,ph2
  real(8) :: addp,delp,ssum,sstr,ssus,dl
  real(8) ::  rndm,p
  
  
  !write(*,*) 'diag update'

  su=0
  sa=0
  nms=0
  cur(:)=0
  su=sum(spn) !sum of total mag
  do i=1,nn
     sa=sa+ph(i)*spn(i) ! staggered mag
  enddo

  su=su/2
  sa=sa/2
  ssum=0.d0
  sstr=dble(sa)**2       !staggered mag squared             
  nu=0
  nk=0

  nh1=0
  last=0
  do i=1,l!i is the propagation time
     op=opstring(i)
     if (op == 0) then !0 means I  ! from empty string to diagonal string
        b=min(int(rndm()*dble(nb))+1,nb) ! select the bond
        ss1=spn(bond(1,b)); ss2=spn(bond(2,b))
        addp=beta*dble(nb)/dble(l-nh) 
        p=awgt(ss1,ss2)*addp !probability of changing to diagonal
        !write(*,*)b,bond(1,b),bond(2,b),ss1,ss2,p
        tries(1)=tries(1)+1
        if (p > 1.d0) then
           opstring(i)=2*b
           nh=nh+1
           accept(1)=accept(1)+1
        elseif (rndm() < p) then
           opstring(i)=2*b
           nh=nh+1
           accept(1)=accept(1)+1
        endif
     elseif (mod(op,2)==0) then     !means op is diagonal op       
        o=mod(op,2)
        b=op/2 !obtain the bond number
        ss1=spn(bond(1,b)); ss2=spn(bond(2,b)) !obtain states at the sites
        delp=dble(l-nh+1)/(beta*dble(nb))
        p=dwgt(ss1,ss2)*delp !probability of changing to I
        tries(2)=tries(2)+1
        if (p > 1.d0) then
           opstring(i)=0
           nh=nh-1
           accept(2)=accept(2)+1
        elseif (rndm() < p) then
           opstring(i)=0
           nh=nh-1
           accept(2)=accept(2)+1
        endif
        nu=nu+1 !diagonal op contribute to potential energy ** might need to add unit
        nh1=nh1+1
     elseif (mod(op,2)==1) then !left off diagonal
        b=op/2;  o=mod(op,2)
        nk=nk+1 !off diagonal contributes to kinetic energy
        dl=dble(nh1-last)
        ssum=ssum+dl*dble(sa)
        sstr=sstr+dl*dble(sa)**2
        s1=bond(1,b);  s2=bond(2,b)
        ss1=spn(s1);   ss2=spn(s2) !spin at site 1 and 2
        ph1=ph(s1);   ph2=ph(s2)
        spn(s1)= -spn(s1);  spn(s2)= -spn(s2)
        t1=spn(s1);    t2=spn(s2)
        sa=sa+ph1*(t1-ss1)+ph2*(t2-ss2)
        if (b<nn) then !check for x or y direction
           cur(1)=cur(1) + spn(s2)-spn(s1)
        else
           cur(2)=cur(2) + spn(s2)-spn(s1)
        endif
        last=nh1
     endif
  enddo
  dl=dble(nh1-last)
  ssum=ssum+dl*dble(sa)
  sstr=sstr+dl*dble(sa)**2 !

  if (nh1.ne.0) then
     ssus=(ssum**2)/(dble(nh1)*dble(nh1+1)*dble(nn)) ! beta/(nh+1)/nh
     ssus=ssus+sstr/(dble(nh+1)*dble(nh1+1)*dble(nn)) !beta/(nh+1)^2
  else
     ssus=sstr/dble(nh+1)/dble(nn)
  endif

  etot=etot+dble(nh)
  avu=avu+dble(nu)
  avk=avk+dble(nk)
  
  umag=umag+dble(su)**2/dble(nn)/dble(nn) ! S(0,0), maximum 1/4
  sxu=sxu+dble(beta)*(dble(su)**2)/dble(nn) ! uniform susceptibility

  ssa=ssa+sstr/dble(nh1+1)/dble(nn)/dble(nn) ! S(pi,pi), maximum 1/4
  sxa=sxa+dble(beta)*ssus/dble(nn) ! staggered susceptibility
  rho(:)=rho(:) + dble(cur(:))**2/dble(nn) ! spin stiffness
  
  
  !call meas

end subroutine dupdate

!========================================!
subroutine updloop(passed)
!========================================!
  use configuration; implicit none

  integer,allocatable :: vert(:),link(:), blocation(:)
  integer :: frst(nn),last(nn)

  integer :: i,j,o,p,p0,p1,vx,vp,ic,ic0,oc!,op
  integer :: i0,ii,i1,nv,nv1,n4,ml
  real(8) :: rndm,r
  integer :: passed

  integer :: s0,s1,ss0,ss1,b,vx0
  
  !write(*,*) 'loop update'

!------ make linked list ------!

  !allocate(vert(0:(l-1))) ! temporary storage for the loop linkage
  !allocate(link(0:(4*l-1))) !link is what leg of the vertex is linked to
  allocate(vert(0:(nh-1)))
  allocate(link(0:(4*nh-1)))
  allocate(blocation(0:(nh-1))) ! the bond location for non-identity opearation 
!  write(*,*) 'l',l,'nh',nh
  frst(:)=-1 !first level an op occurs at the site
  last(:)=-1 !laat level an operator acting on the site
  ii=0
  i0=0
  i1=1
  vert(:)=0
  link(:)=0
  blocation(:)=0
  do p=1,l !along the opertor string
     if (opstring(p)/=0) then ! exclude identity
        o=mod(opstring(p),2) !diag or off diag op
        b=opstring(p)/2 !bond of the op
        blocation(ii) = b
        s0=bond(1,b) ! site number 1 for b
        s1=bond(2,b) ! site number 2 for b
        ss0=spn(s0) !spin state of the site before the action of the operator
        ss1=spn(s1) !spin state 
        if (o==0) then !diagonal
           vert(ii)=legvx(ss0,ss1,ss0,ss1)
           !vert(ii)=legvx(ss0,ss1,ss0,ss1,b,0)
        else           !off-diagonal
           spn(s0)=-spn(s0); spn(s1)=-spn(s1)
           vert(ii)=legvx(ss0,ss1,spn(s0),spn(s1))
           !vert(ii)=legvx(ss0,ss1,spn(s0),spn(s1),b,0)
        endif
        p0=last(s0) !last p that had an op acting on s0
        p1=last(s1)
        if (p0/=-1) then ! site 1 check, tend to link to others
           link(p0)=i0
           link(i0)=p0
        else
           frst(s0)=i0
        endif
        if (p1/=-1) then ! site 2 check, tend to link to others
           link(p1)=i1
           link(i1)=p1
        else
           frst(s1)=i1
        endif
        last(s0)=i0+2 ! the top left leg
        last(s1)=i1+2 ! the top right leg
        ii=ii+1 !this counts the number of vertices vert(ii) givea ith vertex. we have this to track non I op
        i0=i0+4 !4 becuse 4 legs per vertex ! bottom left for the next non-identity vertex
        i1=i1+4 ! bottom right for the next non-identity vertex
     endif
  enddo
  ! for the top of the last vertex and the bottom of the first vertex to link
  do s0=1,nn ! for linkage across the boundary
     i0=frst(s0)
     if (i0/=-1) then ! if site s0 got some vertex along p direction
        p0=last(s0)   ! along p direction pick up the last vertex and the first vertex to link up
        link(p0)=i0
        link(i0)=p0
     endif
  enddo
!write(*,*) ii, nh
!------ loop update ----!

!     dn0=1: worm-head is a
!     dn0=-1: worm-head is a^+

  !write(*,*) 'start vert',vert
  n4=4*nh  ! no of non-identity operators should be equal to i0
  ml=100*l
  nv=0
  do i=1,nl ! nl time of randomly picking the leg from vertexes
     nv1=0
     p0=min(int(rndm()*n4),n4-1) !choose a random leg
     p1=p0
     vx0=p0/4 !vx0 is ii the vertex number that the leg belongs
     ic0=mod(p0,4) !in- leg for starting
     do j=1,ml !ml is maximum loop length
        tries(3)=tries(3)+1
        vp=p1/4
        vx=vert(vp)
       !if (vx>-1) then 
       !    write(*,*) vx
       !    write(*,*) vert
       ! endif
        ic=mod(p1,4) ! in-leg for p0
        r=rndm()
        !write(*,*) 'b',blocation(vp),vx,vp,vert(vp)
        !do o=0,1
           do oc=0,3             ! propose the new output leg
           !write(*,*) ic, oc,vx,o,r,vxprb(oc,ic,vx, blocation(vp),o)
              !write(*,*) ic,oc,vx,vxprb(oc,ic,vx)
              if (r.le.vxprb(oc,ic,vx)) then
                 goto 10
              endif
           enddo
        !enddo
        !o=1
        oc=3 !oc is confirm 3 if the above fails
!        vert(vp)=vxnew(oc,ic,vx)
10      vert(vp)=vxnew(oc,ic,vx)
        if (oc/=ic) then
           accept(3)=accept(3)+1
        endif
        !write(*,*) 'change',vp,vx,vert(vp)
!        do o=0,nh-1
!           if ((vert(o)<1) .or. (vert(o)>10)) then
!           write(*,*) 'step', vert
!           stop
!           endif
!        enddo
        !if (vert(vp)==0) write(*,*) 'got zero'
        !write(*,*) oc,ic,vx,o,vert(vp)
        p1=4*vp+oc ! leave at the output 
        nv=nv+1 !number of vertices visited
        if (p1==p0) goto 20 ! back to p0 leg
        p1=link(p1) 
        if (p1==p0) goto 20 ! back to p0 leg
     enddo
     passed=0
     return !return to diagonal update if number of loops more than ml
20   lopers=lopers+dble(nv)
  enddo
  !write(*,*) 'vert',vert
  j=0
  do i=1,l ! update the opstring from vert(temporary storage of vertexes)
     if (opstring(i) /= 0) then
        opstring(i)=2*(opstring(i)/2)+vxoper(vert(j))
        j=j+1
     endif
  enddo

  do i=1,nn ! update the basis
     if (frst(i) /= -1) then
        ic=mod(frst(i),2)
        vp=frst(i)/4
        spn(i)=vxleg(ic,vert(vp))
     endif
  enddo

  nloops=nloops+DBLE(nl)
  passed=1

  deallocate(vert)
  deallocate(link)
  deallocate(blocation)
end subroutine updloop
!====================================!
!====================================!
subroutine adjnl
!====================================!
  use configuration; implicit none

  integer :: nl1

  lopers=lopers/nloops
  write(*,*) 'lopers',lopers
  nl1=1+int(DBLE(2*l)/lopers)
  nl=(nl+nl1)/2
  lopers=0.d0
  nloops=0.d0
  
end subroutine adjnl
!===================================!
