       program bobnfs5
! Generalisation of bobxt, uses bobgpx3        
       
       
!      third generaslized exponent recovery program
!      computation of kernel of matrix from quadratic double sieving        
       
       dimension indic(40000),iv(40000)
       dimension ic(40000),jpiv(40000),kct(1000)
       dimension nar(40000),littr(20),mpiv(20000)
       dimension ihkia(1000),ihkib(1000),n(60),jfak(30),jfreq(30)
       dimension ibase(50),jbase(50),itarg(50),jtarg(50)
       dimension ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10)
       dimension m1(10),m2(10)
       integer,allocatable,dimension (:,:)::iqi,iqi2,jqi2
       iig=5
       ia=47
       ib=43
       call subbw6(ia,ib,ivv)
       ipd=2

       open(unit=7,file='nfspar',access='direct',form=&
       'formatted',recl=488,status='old')
       read (7,4001,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),&
       (ia4(jf),jf=1,10),(ia3(jf),jf=1,10),(ia2(jf),jf=1,10),&
       (ia1(jf),jf=1,10),(ia0(jf),jf=1,10),(m1(jf),jf=1,10),&
       (m2(jf),jf=1,10),irecnn,klim,kkll,izz4,kmx1,izz6
4001   format (30i4,80i4,6i8)       

1006   format (60i4)       
       print *,'m1',(m1(jf),jf=1,m1(2)+2)
       print *,'m2',(m2(jf),jf=1,m2(2)+2)
       print *,'n',(n(jf),jf=1,30)
       print *,'irecnn',irecnn,'kkll',kkll,'klim',klim,'kmx1',kmx1,'mm',izz6
!       stop
       narc9=kkll+13

       
       open(unit=1,file='nfsf4',access='direct',form=&
       'formatted',recl=narc9,status='old')

       nn=kkll

       
60     a=a

              
!       if (kkll.gt.200000)goto 700
!       nunit=200
       
       
700    nunit=600000/kkll
              
       ntem=(kkll)/nunit 
      nzz=(ntem+2)*nunit 
      nardiv=ntem+1 
      
      
      
      allocate (iqi(nunit,kkll)) 
      allocate (iqi2(nunit,kkll)) 
      allocate (jqi2(nunit,kkll))
      goto 702
 
       
702   print *,'nzz',nzz,'nardiv',nardiv
      

       open(unit=2,file='nkernel',access='direct',form=&
       'formatted',recl=nzz,status='old')
           
       
       

!       read(1,1,rec=3)kia,kib,isgn,(nar(ii),ii=1,kkll)
       

       do i=1,kkll 
       nar(i)=0
       end do
       
       kia=0
       kib=0
       isgn=0
       

       print *,'irecnn',irecnn,'nzz',nzz,'nardiv',nardiv,'nunit',nunit
       print *,'ntem',ntem,'kkll',kkll
       




       if (irecnn.ge.nzz)goto 74




       do i=irecnn+1,nzz
       write(1,1,rec=i)kia,kib,isgn,(nar(ii),ii=1,kkll)
       
!       write (1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
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
       read(1,1,rec=k)kia,kib,isgn,(nar(ii),ii=1,nn-1)
       
       iqi(nk1,1)=isgn
       do k1=1,nn-1
       iqi(nk1,k1+1)=nar(k1)
       end do
       
       

       print *,'k=',k,'ir=',ir
!       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
!       iqi(nk1,1)=mod(littr(3),ipd)
!       do jk=1,nn-1
!       iqi(nk1,jk+1)=mod(nar(jk),ipd)
!       end do
       nk1=nk1+1
!       print *,'nk1',nk1,'iend',iend,'littr',(littr(jd),jd=1,20),'k',k
       end do
       do jf=1,nunit
       kct(jf)=0
       jpiv(jf)=0
       end do
       do k=1,nunit
       do j=1,nn
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
       
       
       
       
       
1      format(i6,i6,i1,40000i1)       
       
       
34     jpiv(k)=j    
       ia=ipd
       ib=iqi(k,j)
       call subbw6(ia,ib,ivv)
       mpiv(k)=ivv       
       mul=ivv
       do i=k,nunit
       iqi(i,j)=(-1)*iqi(i,j)*mul
       iqi(i,j)=mod(iqi(i,j),ipd)
       if (iqi(i,j).ge.0)goto 33
       iqi(i,j)=iqi(i,j)+ipd
33     end do

       do jj=1,nn
       if(jj.eq.j)goto 38
       if (iqi(k,jj).eq.0)goto 38
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       jqi2(k,ktem)=iqi(k,jj)
       mul=iqi(k,jj)
       do i=k,nunit
       itt = iqi(i,jj) + iqi(i,j)*mul 
       iqi(i,jj) =mod(itt,ipd)
       if (iqi(i,jj).ge.0)goto 36
       iqi(i,jj)=iqi(i,jj)+ipd

36     end do
38     end do  
       knd=(nnn-1)*nunit+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     a=a

       
       
       
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
!       if (iv(1).eq.0)goto 610
       ir=ir+1
       print *,'vector found','ir',ir,'ivs',(iv(nq),nq=1,mm),'knd',knd
       write(2,2,rec=ir)(iv(nq),nq=1,mm)
2      format(50000i1)     
       
       if (ir.lt.60)goto 50
       goto 100
610    a=a
       goto 50

50     end do
       if (nnn.ne.1)goto 501
       
       
       
501    if (nnn.eq.nardiv+1)goto 100
       
       do nnn2=nnn+1,nardiv+1
       print *,'nnn2=',nnn2,'nnn',nnn,'ir',ir
       
       ibeg=(nnn2-1)*nunit+1
       iend=ibeg+nunit-1
       nk2=1
       do k=ibeg,iend
       read(1,1,rec=k)kia,kib,isgn,(nar(ii),ii=1,nn-1)
       
!       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       
       
!       iharr(nk2)=ihitn
       
       
!       do jf=1,20
!       larr(nk2,jf)=littr(jf)
!       end do
        ihkia(nk2)=kia
        ihkib(nk2)=kib
        
!       iqi(nk2,1)=mod(littr(3),ipd)
       
       iqi(nk2,1)=isgn
       do k1=1,nn-1
       iqi(nk2,k1+1)=nar(k1)
       end do
       
!       do jf=1,nn-1
!       iqi(nk2,jf+1)=mod(nar(jf),ipd)
       
!       end do
       
       nk2=nk2+1 
       end do
       
       do k=1,nunit
!       if (kct(k).eq.0)goto 601
       if (jpiv(k).eq.0)goto 601
       jnd=jpiv(k)
!       do jf=1,kct(k)
       mul=mpiv(k)
       do ii=1,nunit
       iqi(ii,jnd)=(-1)*iqi(ii,jnd)*mul
       iqi(ii,jnd)=mod(iqi(ii,jnd),ipd)
       if (iqi(ii,jnd).ge.0)goto 710
       iqi(ii,jnd)=iqi(ii,jnd)+ipd
710    end do       
       if (kct(k).eq.0)goto 601
       do jf=1,kct(k)
       
       jnd2=iqi2(k,jf)
       mul=jqi2(k,jf)
       do k2=1,nunit
       itt=iqi(k2,jnd2)+iqi(k2,jnd)*mul
       iqi(k2,jnd2)=mod(itt,ipd)
       if (iqi(k2,jnd2).ge.0)goto 720
       iqi(k2,jnd2)=iqi(k2,jnd2)+ipd
720    end do
       end do
       
       

601    end do
       
       
       
       nnd=(nnn2-1)*nunit
       do k=1,nunit
!       kkb=iqi(k,1)
!       ihitn=iharr(k)
!       larr(k,3)=iqi(k,1)
       kia=ihkia(k)
       kib=ihkib(k)


!       do jf=1,20
!       littr(jf)=larr(k,jf)
!       end do
        isgn=iqi(k,1)
        do jf=1,nn-1
        nar(jf)=iqi(k,jf+1)
        end do
        

       
!       do jf=1,nn-1
!       nar(jf)=iqi(k,jf+1)
!       end do
       ivan=nnd+k
       write(1,1,rec=ivan)kia,kib,isgn,(nar(ii),ii=1,nn-1)
       
!       write(1,1,rec=ivan)kkb,ihitn,(littr(jf),jf=1,20),&
!       (nar(jk),jk=1,nn-1)
       end do
       end do
       end do
100    print *,'ir=',ir
       close(unit=1)
       close(unit=2)
!       n(58)=mm
!       n(60)=ir
       
       write (7,4001,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),&
       (ia4(jf),jf=1,10),(ia3(jf),jf=1,10),(ia2(jf),jf=1,10),&
       (ia1(jf),jf=1,10),(ia0(jf),jf=1,10),(m1(jf),jf=1,10),&
       (m2(jf),jf=1,10),irecnn,klim,kkll,izz4,kmx1,mm
       
!       write (9,1006,rec=1)(n(jf),jf=1,60)
       
       close (unit=7)
       return
       
       end

       subroutine subbw6(ia,ib,ivv)
       ib=mod(ib,ia)
       if (ib.ge.0)goto 10
       ib=ia+ib
10     iu=1
       id=ia
       if (ib.eq.0)goto 888
       iv1=0
       iv3=ib
815    if (iv3.eq.0)goto 830
       iqq=int(id/iv3)
       it3=id-iqq*iv3
       it1=iu-iqq*iv1
       iu=iv1
       id=iv3
       iv1=it1
       iv3=it3
       goto 815
830    iv=(id-ia*iu)/ib
       if (iu.le.0)goto 870
       iv=iv*(1-ia)
       iv=mod(iv,ia)
       if (iv.ge.0)goto 870
       iv=iv+ia
870    a=a
       if (id.gt.1)goto 890
!       print *,'inverse',iv
       goto 890
888    iv=ib
890    ivv=iv
       
       end


