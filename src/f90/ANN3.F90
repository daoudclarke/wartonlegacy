       program ann3
!      third generaslized MPQS program
!      computation of kernel of matrix quad sieve 66 digit nos.        
!      for use by quadratic sieve program and MPQS       
       
       
       dimension indic(40000),iv(40000)
       dimension ic(40000),jpiv(40000),kct(1000)
       dimension nar(40000),littr(20)
       dimension iharr(1000),larr(1000,20),n(60)
       integer,allocatable,dimension (:,:)::iqi,iqi2
       open(unit=7,file='annpar',access='direct',form=&
       'formatted',recl=24,status='old')
       read(7,6,rec=1)irecnn,kkll,nzz
       narc9=kkll+91
6      format(i8,i8,i8)       
       open(unit=1,file='annf4',access='direct',form=&
       'formatted',recl=narc9,status='old')
!       nzz=kkll+100
       if (kkll.gt.2000)goto 700
       nunit=500
       
       goto 701
       
700    nunit=1000000/kkll
       
701    ntem=kkll/nunit 
      nzz=(ntem+2)*nunit 
      nardiv=ntem+1 
      if (nzz.gt.irecnn)goto 703
      nardiv=nardiv+1
      nzz=nzz+nunit

703   print *,'nardiv',nardiv
      
      
      allocate (iqi(nunit,kkll)) 
      allocate (iqi2(nunit,kkll)) 
      goto 702
!      nunit=100 
!      nzz=kkll+100 
!      ntem=kkll/nunit 
!      nardiv=ntem+1 
       
702   print *,'nzz',nzz
      write(7,6,rec=1)irecnn,kkll,nzz




       
       open(unit=2,file='annf5',access='direct',form=&
       'formatted',recl=nzz,status='old')
       read (1,1,rec=20)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
       print *,kkb,(nar(jf),jf=1,300)
       do i=1,kkll
       nar(i)=0
       end do
       kkb=0
       ihitn=0
       do i=1,20
       littr(i)=0
       end do

       print *,'irecnn',irecnn,'nzz',nzz,'nardiv',nardiv,'nunit',nunit
       print *,'ntem',ntem,'kkll',kkll
       




       if (irecnn.ge.nzz)goto 74




       do i=irecnn+1,nzz
       write (1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
       end do
       
       
       
       
74     mm=nzz
       nn=kkll
       do i=1,mm
       indic(i)=0
       end do
       do i=1,nn
       ic(i)=0
       end do
       ir=0
       
       do nnn=1,nardiv
       ibeg=(nnn-1)*nunit+1
       iend =ibeg+nunit-1
       nk1=1
       do k=ibeg,iend
       print *,'k=',k,'ir=',ir
       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       iqi(nk1,1)=kkb
       do jk=1,nn-1
       iqi(nk1,jk+1)=nar(jk)
       end do
       nk1=nk1+1
       end do
       do jf=1,nunit
       kct(jf)=0
       end do
       do k=1,nunit
       do j=1,nn
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
       
       
       
       
       
1      format(i1,i10,20i4,50000i1)       
       
       
34     jpiv(k)=j    
       
              
       

       do jj=1,nn
       if(jj.eq.j)goto 38
       if (iqi(k,jj).eq.0)goto 38
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,nunit
       itt = iqi(i,jj) + iqi(i,j) 
       iqi(i,jj) =mod(itt,2)
       
       

36     end do
38     end do  
       knd=(nnn-1)*nunit+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     ir =ir +1
       
       
       
       do is =1,mm
       iv(is) = 0
       end do
       knd=(nnn-1)*nunit+k
       iv(knd) =1
       do is =1,mm

       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is) =iqi(k,jv)
42     end do
       print *,'vector found','ir',ir,'ivs',(iv(nq),nq=1,mm)
       write(2,2,rec=ir)(iv(nq),nq=1,mm)
2      format(50000i1)     
       
       if (ir.lt.60)goto 50
       close(unit=1)
       close(unit=2)
       print *, 'ir=',ir
       stop
50     end do
       if (nnn.ne.1)goto 501
       
       
       
501    if (nnn.eq.nardiv+1)goto 100
       
       do nnn2=nnn+1,nardiv+1
       print *,'nnn2=',nnn2
       
       ibeg=(nnn2-1)*nunit+1
       iend=ibeg+nunit-1
       nk2=1
       do k=ibeg,iend
       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       
       
       iharr(nk2)=ihitn
       
       
       do jf=1,20
       larr(nk2,jf)=littr(jf)
       end do
       
       
       iqi(nk2,1)=kkb
       do jf=1,nn-1
       iqi(nk2,jf+1)=nar(jf)
       
       end do
       
       nk2=nk2+1 
       end do
       
       do k=1,nunit
       if (kct(k).eq.0)goto 601
       jnd=jpiv(k)
       do jf=1,kct(k)
       jnd2=iqi2(k,jf)
       
       do k2=1,nunit
       itt=iqi(k2,jnd2)+iqi(k2,jnd)
       iqi(k2,jnd2)=mod(itt,2)
       end do
       end do
       
       

601    end do
       
       
       
       nnd=(nnn2-1)*nunit
       do k=1,nunit
       kkb=iqi(k,1)
       ihitn=iharr(k)
       do jf=1,20
       littr(jf)=larr(k,jf)
       end do
       
       
       do jf=1,nn-1
       nar(jf)=iqi(k,jf+1)
       end do
       ivan=nnd+k
       write(1,1,rec=ivan)kkb,ihitn,(littr(jf),jf=1,20),&
       (nar(jk),jk=1,nn-1)
       end do
       end do
       end do
100    print *,'ir=',ir
       close(unit=1)
       close(unit=2)



       
       
       end

       

