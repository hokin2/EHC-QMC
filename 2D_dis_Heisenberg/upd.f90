! subroutines updating the config
!=====================================!
subroutine mcstep
!=====================================!
  use configuration
  implicit none
  integer :: pass

10 call dupdate
  call updloop(pass)
 if (pass==0) goto 10

end subroutine mcstep

!=======================================!
subroutine dupdate
!=======================================!
  use configuration
  implicit none

  integer :: i,o,op,b,s1,s2,l_i,l_q,q
  real(8) :: addp,delp
  real(8) ::  rndm,p
  l_i = 1
  do q = 1, qmax
     do l_q = 1, l_len(q) !i is the propagation time
       op = opstring(l_i)
       if (op == 0) then !0 means I  ! from empty string to diagonal string
          b=min(int(rndm()*dble(nb))+1,nb) ! select the bond
          s1=spn(q,bond(1,b))
          s2=spn(q,bond(2,b))
          addp=beta*dble(nb)/dble(l_len(q)-nh_q(q))
          p=awgt(s1,s2,b,0)*addp !probability of changing to diagonal
          !write(*,*)b,bond(1,b),bond(2,b),ss1,ss2,p
          if (p > 1.d0) then ! to save power
             opstring(l_i)=2*b
             nh_q(q)=nh_q(q)+1
          elseif (rndm() < p) then
             opstring(l_i)=2*b
             nh_q(q)=nh_q(q)+1
          endif
       elseif (mod(op,2)==0) then     !means op is diagonal op
          o=mod(op,2)
          b=op/2 !obtain the bond number
          s1=spn(q,bond(1,b))
          s2=spn(q,bond(2,b)) !obtain states at the sites
          delp=dble(l_len(q)-nh_q(q)+1)/(beta*dble(nb))
          p=dwgt(s1,s2,b,0)*delp !probability of changing to I
          if (p > 1.d0) then ! to save power
             opstring(l_i)=0
             nh_q(q)=nh_q(q)-1
          elseif (rndm() < p) then
             opstring(l_i)=0
             nh_q(q)=nh_q(q)-1
          endif
       elseif (mod(op,2)==1) then !left off diagonal
          b=op/2
          o=mod(op,2)
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

end subroutine dupdate

!==========================================================!
subroutine adjstl !adjust string length
!==========================================================!
  use configuration; implicit none

  integer :: p1,p2, q,l_i, dl
  integer,allocatable :: opstr_q(:), opstr_new(:),opstr_q_new(:)
  real(8) :: r,rndm
  integer :: l_len_new(qmax),l_start_new(qmax), l_new
  logical :: extend(qmax)

  ! calculate the new length of string
  l_i = 1
  do q = 1, qmax
    l_start_new(q) = l_i
    dl=nh_q(q)/4
    if (nh_q(q)+dl<l_len(q)) then
      extend(q) = .False.
      l_len_new(q) = l_len(q)
    else
      extend(q) = .True.
      l_len_new(q)=l_len(q)+dl
    endif
    l_i = l_i + l_len_new(q) + 1
  enddo
  l_new = sum(l_len_new) + qmax -1
  if (l_new .ne. l_i-2) then
    write(*,*) 'sth wrong2', l_i-2, l_new
    stop
  endif

  ! randomize the empty part
  ! opstring -> opstr_q -> opstr_q_new -> opstr_new
  allocate(opstr_new(l_new))
  opstr_new = 0
  do q = 1, qmax
    if (extend(q)) then
      allocate(opstr_q(l_len(q)))
      allocate(opstr_q_new(l_len_new(q)))
      opstr_q(:) = opstring(l_start(q):l_start(q)+l_len(q)-1)
      do l_i = 1, l_len(q)
         opstr_q_new(l_i) = opstr_q(l_i)
      enddo
      do l_i = l_len(q)+1, l_len_new(q)
         opstr_q_new(l_i) = 0
      enddo

      ! find the ratio
      r=dble(nh_q(q))/dble(l_len_new(q))
      ! find the first non-empty
      do l_i = l_len_new(q),1,-1
         if (opstr_q_new(l_i) /= 0) then
            p1=l_i
            exit
         endif
      enddo
      outer: do p2=l_len_new(q),1,-1
         if (p2 <= p1) then
            exit outer
         elseif (rndm() <= r) then
            opstr_q_new(p2)=opstr_q_new(p1)
            opstr_q_new(p1)=0
            innner: do l_i=p1-1,1,-1
               if (opstr_q_new(l_i) /= 0) then
                  p1=l_i
                  cycle outer
               endif
            enddo innner
            exit outer
         endif
      enddo outer
      opstr_new(l_start_new(q):l_start_new(q)+l_len_new(q)-1) = opstr_q_new(:)
      deallocate(opstr_q_new, opstr_q)
    else
      opstr_new(l_start_new(q):l_start_new(q)+l_len_new(q)-1) = opstring(l_start(q):l_start(q)+l_len(q)-1)
    endif
  enddo

  deallocate(opstring)
  allocate(opstring(l_new))
  opstring(:) = opstr_new(:)
  l_len = l_len_new(:)
  l_start = l_start_new(:)
  l = l_new
  deallocate(opstr_new)
end subroutine adjstl

!========================================!
subroutine updloop(pass)
!========================================!
  use configuration
  implicit none

  integer,allocatable :: vert(:),link(:)
  integer :: frst(nn),last(nn)

  integer :: i,j,o,p,p0,p1,vx,vp,ic,ic0,oc, q!,op
  integer :: i0,ii,i1,nv,nv1,n4,ml
  real(8) :: rndm,r
  integer :: pass

  integer :: s0,s1,ss0,ss1,b,vx0

  do q = 1, qmax
  !------ make linked list ------!

    !allocate(vert(0:(l-1))) ! temporary storage for the loop linkage
    !allocate(link(0:(4*l-1))) !link is what leg of the vertex is linked to
    allocate(vert(0:(nh_q(q)-1)))
    allocate(link(0:(4*nh_q(q)-1)))
    !  write(*,*) 'l',l,'nh',nh
    frst(:)=-1 !first level an op occurs at the site
    last(:)=-1 !laat level an operator acting on the site
    ii=0
    i0=0
    i1=1
    vert(:)=0
    link(:)=0
    do p=l_start(q),l_start(q)+l_len(q)-1 !along the opertor string
       if (opstring(p)/=0) then ! exclude identity
          o=mod(opstring(p),2) !diag or off diag op
          b=opstring(p)/2 !bond of the op
          s0=bond(1,b) ! site number 1 for b
          s1=bond(2,b) ! site number 2 for b
          ss0=spn(q,s0) !spin state of the site before the action of the operator
          ss1=spn(q,s1) !spin state
          if (o==0) then !diagonal
             vert(ii)=legvx(ss0,ss1,ss0,ss1,b,0)
          else           !off-diagonal
             spn(q,s0)=-spn(q,s0)
             spn(q,s1)=-spn(q,s1)
             vert(ii)=legvx(ss0,ss1,spn(q,s0),spn(q,s1),b,0)
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

  !------ loop update ----!
    n4=4*nh_q(q)  ! no of non-identity operators should be equal to i0
    ml=100*l_len(q)
    nv=0
    do i=1,nl ! nl time of randomly picking the leg from vertexes
       nv1=0
       p0=min(int(rndm()*n4),n4-1) !choose a random leg
       p1=p0
       vx0=p0/4 !vx0 is ii the vertex number that the leg belongs
       ic0=mod(p0,4) !in- leg for starting
       do j=1,ml !ml is maximum loop length
          vp=p1/4
          vx=vert(vp)
          ic=mod(p1,4) ! in-leg for p0
          r=rndm()
          do oc=0,3             ! propose the new output leg
            if (r.le.vxprb(oc,ic,vx)) then
               goto 10
            endif
          enddo
          oc=3 !oc is confirm 3 if the above fails
  10      vert(vp)=vxnew(oc,ic,vx)
          p1=4*vp+oc ! leave at the output
          nv=nv+1 !number of vertices visited
          if (p1==p0) goto 20 ! back to p0 leg
          p1=link(p1)
          if (p1==p0) goto 20 ! back to p0 leg
       enddo
       pass=0
       return !return to diagonal update if number of loops more than ml
  20   lopers=lopers+dble(nv)
    enddo
    !write(*,*) 'vert',vert
    j=0
    do i=l_start(q),l_start(q)+l_len(q)-1 ! update the opstring from vert(temporary storage of vertexes)
       if (opstring(i) /= 0) then
          opstring(i)=2*(opstring(i)/2)+vxoper(vert(j))
          j=j+1
       endif
    enddo

    do i=1,nn ! update the basis
       if (frst(i) /= -1) then
          ic=mod(frst(i),2)
          vp=frst(i)/4
          spn(q,i)=vxleg(ic,vert(vp))
       endif
    enddo
    nloops=nloops+DBLE(nl)
    pass=1
    deallocate(vert)
    deallocate(link)
end do
end subroutine updloop
!====================================!

!====================================!
subroutine adjnl
!====================================!
  use configuration
  implicit none
  integer :: nl1
  lopers=lopers/nloops
  nl1=1+int(DBLE(2*l)/lopers)
  nl=(nl+nl1)/2
  lopers=0.d0
  nloops=0.d0

end subroutine adjnl
!===================================!
