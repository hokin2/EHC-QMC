!==========================================!
subroutine initvrtx
!==========================================!
  use configuration; implicit none

  integer :: i,iv,s1,s2,t1,t2,lo,li,ie,ne

  call vtxwgts !for directed loop
  call eqtable
  call prtable
  
!  call prtablea !what we know

      open(unit=10,file='vxw.dat',status='unknown')
      do i=1,nvx
         write(10,*)i,vxw(i)
      enddo
      close(10)

      open(unit=10,file='etab.dat',status='unknown')
      do ie=1,12
         do ne=1,3
            write(10,*)'vx=',etab(4,ne,ie),'  li=',etab(5,ne,ie)
            write(10,*)etab(0,ne,ie),etab(1,ne,ie),etab(2,ne,ie),etab(3,ne,ie)
         enddo
         write(10,*)'=============================='
      enddo
      close(10)

end subroutine initvrtx
!=======================================!
subroutine prtable
!=======================================!
!     solution b (bounce minimization) vertex probabilities.
!=============================================================!
  use configuration; implicit none

  integer :: i,j,iv,j1,j2,j3,jj,ie,li,lo,vxold,vx,ic,oc
  real(8) ::  w1,w2,w3,aa,bb,cc

  vxnew(:,:,:)=0
  vxprb(:,:,:)=0.d0

  do ie=1,2*nvx
     if (neqs(ie).eq.2) then
        j1=1
        j2=2
        w1=vxw(etab(4,j1,ie))
        w2=vxw(etab(4,j2,ie))
        if (w2.lt.w1) then
           j1=2
           j2=1
           w1=vxw(etab(4,j1,ie))
           w2=vxw(etab(4,j2,ie))
        endif
        li=etab(5,j1,ie)
        etab(li,j1,ie)=0
        vxold=etab(4,j1,ie)
        do lo=0,3
           if (etab(lo,j1,ie).ne.0) then
              vxnew(lo,li,vxold)=etab(lo,j1,ie)
              vxprb(lo,li,vxold)=1.d0
           endif
        enddo
        li=etab(5,j2,ie)
        vxold=etab(4,j2,ie)
        do lo=0,3
           if (etab(lo,j2,ie).eq.etab(4,j1,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j2,ie)
              vxprb(lo,li,vxold)=w1/w2
           elseif (etab(lo,j2,ie).ne.0) then
              vxnew(lo,li,vxold)=etab(lo,j2,ie)
              vxprb(lo,li,vxold)=(w2-w1)/w2
           endif
        enddo
     elseif (neqs(ie).eq.3) then
        j1=1
        j2=2
        j3=3
        w1=vxw(etab(4,j1,ie))
        w2=vxw(etab(4,j2,ie))
        w3=vxw(etab(4,j3,ie))
        
        if (w2.lt.w1) then
           j1=2
           j2=1
           w1=vxw(etab(4,j1,ie))
           w2=vxw(etab(4,j2,ie))
        endif
        
        if (w3.lt.w2) then
           jj=j2
           j2=j3
           j3=jj
           w2=vxw(etab(4,j2,ie))
           w3=vxw(etab(4,j3,ie))
        endif
        if (w2.lt.w1) then
           jj=j1
           j1=j2
           j2=jj
           w1=vxw(etab(4,j1,ie))
           w2=vxw(etab(4,j2,ie))
        endif
        
        aa=0.5d0*(w1+w2-w3)
        if (aa.lt.0.d0) aa=0.d0
        bb=w1-aa
        cc=w2-aa
        
        li=etab(5,j1,ie)
        etab(li,j1,ie)=0
        vxold=etab(4,j1,ie)
        do lo=0,3
           if (etab(lo,j1,ie).eq.etab(4,j2,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j1,ie)
              vxprb(lo,li,vxold)=aa/w1
           elseif (etab(lo,j1,ie).eq.etab(4,j3,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j1,ie)
              vxprb(lo,li,vxold)=bb/w1
           endif
        enddo
        li=etab(5,j2,ie)
        etab(li,j2,ie)=0
        vxold=etab(4,j2,ie)
        do lo=0,3
           if (etab(lo,j2,ie).eq.etab(4,j1,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j2,ie)
              vxprb(lo,li,vxold)=aa/w2
           elseif (etab(lo,j2,ie).ne.0) then
              vxnew(lo,li,vxold)=etab(lo,j2,ie)
              vxprb(lo,li,vxold)=cc/w2
           endif
        enddo
        li=etab(5,j3,ie)
        vxold=etab(4,j3,ie)
        do lo=0,3
           if (etab(lo,j3,ie).eq.etab(4,j1,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j3,ie)
              vxprb(lo,li,vxold)=bb/w3
           elseif (etab(lo,j3,ie).eq.etab(4,j2,ie)) then
              vxnew(lo,li,vxold)=etab(lo,j3,ie)
              vxprb(lo,li,vxold)=cc/w3
           elseif (etab(lo,j3,ie).ne.0) then
              vxnew(lo,li,vxold)=etab(lo,j3,ie)
              vxprb(lo,li,vxold)=(w3-w2-w1+2.d0*aa)/w3
           endif
        enddo
     endif
  enddo

  do iv=1,nvx
  do li=0,3
  do lo=1,3
     vxprb(lo,li,iv)=vxprb(lo,li,iv)+vxprb(lo-1,li,iv)
  enddo
  enddo
  enddo
  
  do iv=1,nvx
     do li=0,3
        vxnew(li,li,iv)=iv
     enddo
  enddo
  
  do vx=1,nvx
  do ic=0,3
  do oc=0,3
     if (vxprb(oc,ic,vx).lt.0.00001) vxprb(oc,ic,vx)=-1.
  enddo
  enddo
  enddo
           
end subroutine prtable
!=============================================!
subroutine eqtable ! equation table 
!=============================================!
  use configuration; implicit none

  integer :: i,j,k,li,l1,l2,iv,jv,ie,ns,iiv,ll1,ll2,ic,oc,vx,vxn,dir
  integer :: ftab(0:3,0:3,nvx),st(0:3),stak(2,2*nvx)
  integer :: eset(0:3,nvx)

  ftab(:,:,:)=0
    ! book keeping, how the in-leg and out-leg affect vertex
  do vx=1,nvx
     do ic=0,3
     do oc=0,3
        st(0)=vxleg(0,vx)
        st(1)=vxleg(1,vx)
        st(2)=vxleg(2,vx)
        st(3)=vxleg(3,vx)
        st(ic)=-st(ic)
        st(oc)=-st(oc)
        vxn=legvx(st(0),st(1),st(2),st(3))
        if (vxn.ne.0) then
           ftab(oc,ic,vx)=vxn
        endif
     enddo
     enddo
  enddo

  open(10,file='ftab.dat',status='unknown')
  do vx=1,nvx
  do li=0,3
     write(10,*)vxleg(li,vx)
     write(10,*)ftab(0,li,vx),ftab(1,li,vx),ftab(2,li,vx),ftab(3,li,vx)
  enddo
  write(10,*)'=============================='
  enddo
  close(10)


  eset(:,:)=0
  etab(:,:,:)=0

  ie=0
  do iv=1,nvx
     do l1=0,3
        if (eset(l1,iv).eq.0) then
           ie=ie+1
           eset(l1,iv)=ie
           ns=1
           stak(1,ns)=iv
           stak(2,ns)=l1
40         if (ns.ne.0) then
              iiv=stak(1,ns)
              ll1=stak(2,ns)
              ns=ns-1
              do ll2=0,3
                 if (ll1.ne.ll2.and.ftab(ll2,ll1,iiv).ne.0) then
                    if (eset(ll2,ftab(ll2,ll1,iiv)).eq.0) then
                       eset(ll2,ftab(ll2,ll1,iiv))=ie
                       ns=ns+1
                       stak(1,ns)=ftab(ll2,ll1,iiv)
                       stak(2,ns)=ll2
                    endif
                 endif
              enddo
              goto 40
           endif
        endif
     enddo
  enddo

  do ie=1,2*nvx
     neqs(ie)=0
     do iv=1,nvx
        do l1=0,3
           if (eset(l1,iv).eq.ie) then
              neqs(ie)=neqs(ie)+1
              etab(0,neqs(ie),ie)=ftab(0,l1,iv)
              etab(1,neqs(ie),ie)=ftab(1,l1,iv)
              etab(2,neqs(ie),ie)=ftab(2,l1,iv)
              etab(3,neqs(ie),ie)=ftab(3,l1,iv)
              etab(4,neqs(ie),ie)=iv
              etab(5,neqs(ie),ie)=l1
           endif
        enddo
     enddo
  enddo

end subroutine eqtable

!==========================================================!
subroutine vtxwgts
!==========================================================!
  use configuration; implicit none

  integer :: s1,s2,t1,t2,op,vx,ic,oc,vxn,st(0:3)
  real(8) :: vxwgt(nvx)


  legvx(:,:,:,:)=0 !given the leg states can obtain vertex order type

  vx=0
  do op=0,1
     do s2=-1,1,2 
     do s1=-1,1,2
        if((op == 0) .or. ((op == 1) .and. (s1 /= s2))) then
           write(*,*)op,s1,s2
           vx=vx+1
           vxoper(vx)=op
           vxcode(op,s1,s2)=vx
           if (op == 0) then
              t1=s1; t2=s2 !t1 t2 are states after action of the op
           else
              t1=-s1; t2=-s2
           endif
           vxleg(0,vx)=s1 !vxleg gives the leg state of the vertex
           vxleg(1,vx)=s2
           vxleg(2,vx)=t1
           vxleg(3,vx)=t2
           legvx(s1,s2,t1,t2)=vx
           if (op.eq.0) then
              vxw(vx)=awgt(s1,s2)
           else
              vxw(vx)=0.5
           endif
        else
           vxcode(op,s1,s2)=0 !given op s1 s2 what is the vertex number
        endif
     enddo
     enddo
  enddo

end subroutine vtxwgts
!==========================================================!
!==========================================================!


