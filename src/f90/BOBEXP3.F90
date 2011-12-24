       program bobexp3
!      third generaslized exponent recovery program
!      computation of kernel of matrix         
       
       
       
       dimension indic(40000),iv(40000)
       dimension ic(40000),jpiv(40000),kct(1000)
       dimension nar(40000),littr(20),mpiv(20000)
       dimension iharr(1000),larr(1000,20),n(60),jfak(30),jfreq(30)
       dimension ibase(50),jbase(50),itarg(50),jtarg(50)
       integer,allocatable,dimension (:,:)::iqi,iqi2,jqi2
       ia=47
       ib=43
       call subbw6(ia,ib,ivv)
       ipd=2
       open(unit=7,file='exppar',access='direct',form=&
       'formatted',recl=344,status='old')
       open (unit=9,file='expparn',access='direct',form=&
       'formatted',recl=240,status='old')
1006   format (60i4)       
       read (9,1006,rec=1)(n(jf),jf=1,60)
       read(7,6,rec=1)irecnn,kkll,ijpow,nfak,(ibase(jf),jf=1,20),&
       (jbase(jf),jf=1,20),(itarg(jf),jf=1,20),(jtarg(jf),jf=1,20)

       narc9=kkll*4+91
6      format(i6,i6,i8,i4,20i4,20i4,20i4,20i4)       
       open(unit=1,file='sbexp1',access='direct',form=&
       'formatted',recl=narc9,status='old')
!       nzz=kkll+100
       if (kkll.gt.200000)goto 700
       nunit=200
       goto 701
       
700    nunit=600000/kkll
       
701    ntem=kkll/nunit 
      nzz=(ntem+2)*nunit 
      nardiv=ntem+1 
      
      
      
      allocate (iqi(nunit,kkll)) 
      allocate (iqi2(nunit,kkll)) 
      allocate (jqi2(nunit,kkll))
      goto 702
!      nunit=100 
!      nzz=kkll+100 
!      ntem=kkll/nunit 
!      nardiv=ntem+1 
       
702   print *,'nzz',nzz
!      write(7,6,rec=1)irecnn,kkll,nzz
      write (7,6,rec=1)irecnn,kkll,ijpow,nzz,(ibase(jf),jf=1,20),&
      (jbase(jf),jf=1,20),(itarg(jf),jf=1,20),(jtarg(jf),jf=1,20)
      nzz2=4*nzz

       
       open(unit=2,file='shexp1',access='direct',form=&
       'formatted',recl=nzz2,status='old')
       read (1,1,rec=20)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
       print *,kkb,(nar(jf),jf=1,300),'ihitn',ihitn
       print *,'littr',(littr(jf),jf=1,20)
      do ii=1,20 
!      read (2,2,rec=ii)(iv(jf),jf=1,nzz)
!       print *,'iv1',iv(1)
       end do
!       stop
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
!       ijsum=0
!      read (2,2,rec=2)(iv(jf),jf=1,mm)
!      print *,'iv1',iv(1)
!       do ii=1,mm
!      read (1,1,rec=ii)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
!      if (ii.ne.1)goto 741
!      print *,'nar4',nar(4)
!741   a=a      
!      ijsum=ijsum+nar(4)*iv(ii)
       
!       ijsum=mod(ijsum,ipd)
!       end do
!      print *,'ijsum',ijsum
!       stop
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
!       print *,'nk1',nk1,'iend',iend,'littr',(littr(jd),jd=1,20),'k',k
       end do
       do jf=1,nunit
       kct(jf)=0
       end do
       do k=1,nunit
       do j=1,nn
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
       
       
       
       
       
1      format(i1,i10,20i4,50000i4)       
       
       
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
       if (iv(1).eq.0)goto 610
       ir=ir+1
       print *,'vector found','ir',ir,'ivs',(iv(nq),nq=1,mm),'knd',knd
       write(2,2,rec=ir)(iv(nq),nq=1,mm)
2      format(50000i4)     
       
       if (ir.lt.60)goto 50
       goto 100
610    a=a
       goto 50

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
!       if (kct(k).eq.0)goto 601
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
       n(58)=mm
       n(60)=ir
       
       write (9,1006,rec=1)(n(jf),jf=1,60)
       close (unit=9)
       close (unit=7)
       
       
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
       print *,'inverse',iv
       goto 890
888    iv=ib
890    ivv=iv
       return
       end


