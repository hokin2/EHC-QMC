!==========================================!
subroutine init_vrtx
!==========================================!
  use configuration; implicit none

!  integer :: i,iv,s1,s2,t1,t2,lo,li,ie,ne

  call vtxwgts !for directed loop
!  call eqtable
  call prtable

!  call prtablea !what we know


!      open(unit=10,file='vtx.dat',status='unknown')
!      do iv =1,nvx
!      write(10,*)'vx=',iv
!      do li=1,4
!      write(10,*)'   li=',li
!      do lo=1,4
!         write(10,10)lo,vxnew(lo,li,iv),vxprb(lo,li,iv)
!      enddo
!      write(10,*)'============================================='
!      enddo
!      write(10,*)'*******************************************'
!      enddo
! 10   format(2i3,' ',f8.4)
!      close(10)

end subroutine init_vrtx

subroutine prtable
!=======================================!
!     solution b (bounce minimization) vertex probabilities.
!=============================================================!
  use configuration; implicit none

  integer :: j,k,ij,vxn,vx,ic,oc,st(0:3),vt,op,b
  real(8) ::  legwgt(0:3,0:3,nvx)

  integer :: jn(0:3),jo(0:3),jv(0:3)
  real(8) :: fac,jw(0:3),mwgt(3,3)

  integer :: jord(3),kord(3),icn,nj,nk,vxk
  jord =(/1,2,1/)
  kord =(/2,3,2/)

  legwgt=0.d0 ! legwgt(0:3,0:3,nvx)=0.d0 !
  vxnew=0 ! vxnew(0:3,0:3,nvx)=0 !
  vxprb=0.d0 ! vxprb(0:3,0:3,nvx)=0.d0 !

  !Create vxnew links
  do vx=1,nvx
     vt=vtyp(vx)
     op=vxoper(vx)
     b=vxbnd(vx)
     do ic=0,3
        do oc=0,3
           st(0)=vxleg(0,vx)
           st(1)=vxleg(1,vx)
           st(2)=vxleg(2,vx)
           st(3)=vxleg(3,vx)
           st(ic)=-st(ic)
           st(oc)=-st(oc)
           vxn=legvx(st(0),st(1),st(2),st(3),b,vt)
           if (vxn.ne.0) then
              vxnew(oc,ic,vx)=vxn
           endif
        enddo
     enddo
  enddo

  !Choose a vertex and in channel
  do vx=1,nvx
     do ic=0,3
        !Now we have list of vertices in Directed Loop equations
        !Initialize vertex order
        j=0
        do oc=0,3
           vxn=vxnew(oc,ic,vx)
           if (vxn.ne.0) then
              j=j+1
              jn(j)=j
              jo(j)=oc
              jv(j)=vxn
              jw(j)=vxw(vxn)
           endif
        enddo
        if (j.ne.3) stop 'DL FAIL 1'
       !Three moves to relative order
        do ij=1,3
           j=jord(ij)
           k=kord(ij)
           if (jw(j).lt.jw(k)) then
              jn(0)=jn(j)
              jn(j)=jn(k)
              jn(k)=jn(0)
              jw(0)=jw(j)
              jw(j)=jw(k)
              jw(k)=jw(0)
           endif
        enddo
        !Implement Eqs. 22 and 23 of Sylju\r{a}sen
        mwgt(1,1)=max(0.d0,jw(1)-jw(2)-jw(3))
        mwgt(1,2)=min(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)))
        mwgt(1,3)=min(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)))
        mwgt(2,1)=min(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)))
        mwgt(2,2)=0.d0
        mwgt(2,3)=max(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)))
        mwgt(3,1)=min(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)))
        mwgt(3,2)=max(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)))
        mwgt(3,3)=0.d0
        !Finally, assign weights to Directed Loop vertex elements
        do j=1,3
           nj=jn(j)
           icn=jo(nj)
           vxn=jv(nj)
           do k=1,3
              nk=jn(k)
              vxk=jv(nk)
              ! Cycle through to find correct out channel
              ! (the one that changes vxj to vxk)
              do oc=0,3
                 if (vxnew(oc,icn,vxn).eq.vxk) exit
                 if (oc.eq.3) stop 'DL FAIL 2'
              enddo
              legwgt(oc,icn,vxn)=mwgt(j,k)
           enddo
        enddo
     enddo
  enddo

  !Convert matrix weights to probabilities
  do vx=1,nvx
     fac=1.d0/vxw(vx)
     do ic=0,3
        do oc=0,3
           vxprb(oc,ic,vx)=fac*legwgt(oc,ic,vx)
        enddo
     enddo
  enddo

  !Create cumulative probabilities
  do vx=1,nvx
     do ic=0,3
        do oc=1,3
           vxprb(oc,ic,vx)=vxprb(oc,ic,vx)+vxprb(oc-1,ic,vx)
        enddo
     enddo
  enddo

  !Normalize cumulative probabilities
  do vx=1,nvx
     do ic=0,3
        do oc=0,3
           if (vxprb(3,ic,vx).lt.1.d-6) then
              vxprb(oc,ic,vx)=0.d0
           else
              vxprb(oc,ic,vx)=vxprb(oc,ic,vx)/vxprb(3,ic,vx)
              if (vxprb(oc,ic,vx).lt.1.d-6) vxprb(oc,ic,vx)=-1.d0
           endif
        enddo
     enddo
  enddo

end subroutine prtable

!==========================================================!
subroutine vtxwgts
!==========================================================!
  use configuration; implicit none

  integer :: s1,s2,t1,t2,op,vx
  integer :: b,vt

  legvx(:,:,:,:,:,:)=0

  vx=0
  do b=1,nb
     vt=0
     do op=0,1
     do s2=-1,1,2
     do s1=-1,1,2
        if((op == 0) .or. ((op == 1) .and. (s1 /= s2))) then
           vx=vx+1
           !write(*,*) vx,op,s1,s2
           vtyp(vx)=vt
           vxoper(vx)=op
           vxbnd(vx)=b
           vxcode(op,s1,s2,b)=vx
           if (op == 0) then
              t1=s1; t2=s2
           else
              t1=-s1; t2=-s2
           endif
           vxleg(0,vx)=s1
           vxleg(1,vx)=s2
           vxleg(2,vx)=t1
           vxleg(3,vx)=t2
           legvx(s1,s2,t1,t2,b,vt)=vx
           if (op.eq.0) then
              vxw(vx)=awgt(s1,s2,b,op)
           else
              vxw(vx)=0.5
           endif
        else
           vxcode(op,s1,s2,b)=0
        endif
     enddo
     enddo
     enddo
  enddo

end subroutine vtxwgts
!==========================================================!
