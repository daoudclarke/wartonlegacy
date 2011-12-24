       program bobker3
!      program to obtain kernels of matrices in z single -precision
!      of small matrices. works right to left         
       
       
       
       dimension indic(40000),iv(40000)
       dimension ic(40000),jpiv(40000),kct(1000),izzy(10)
       dimension nar(200,200),littr(20),mpiv(20000),nused(200)
       dimension iharr(1000),larr(1000,20),n(60),jfak(30),jfreq(30)
       dimension ibase(50),jbase(50),itarg(50),jtarg(50),iset(2000)
       dimension iqi(200,200),iqi2(200,200),jqi2(200,200),irno(60)
       
       
       iig=5
       ia=-43
       ib=11
       call subbw2(5,1,igcd2)
       print *,'igcd2',igcd2
       
       
       call subbw6(ia,ib,ivv,iu)
       print *,'ivv',ivv,'iu',iu
       

       ipd=2
       nar(1,1)=8
       nar(1,2)=5
       nar(1,3)=3
       nar(1,4)=1
       nar(1,5)=2
       nar(1,6)=3
       nar(1,7)=5
       nar(1,8)=9
       nar(1,9)=4
       nar(2,1)=9
       nar(2,2)=10
       nar(2,3)=7
       nar(2,4)=12
       nar(2,5)=4
       nar(2,6)=5
       nar(2,7)=7
       nar(2,8)=9
       nar(2,9)=6
       nar(3,1)=9
       nar(3,2)=10
       nar(3,3)=11
       nar(3,4)=12
       nar(3,5)=6
       nar(3,6)=7
       nar(3,7)=4
       nar(3,8)=14
       nar(3,9)=8
       nar(4,1)=13
       nar(4,2)=14
       nar(4,3)=15
       nar(4,4)=4
       nar(4,5)=8
       nar(4,6)=9
       nar(4,7)=6
       nar(4,8)=19
       nar(4,9)=5
       nar(5,1)=17
       nar(5,2)=18
       nar(5,3)=19
       nar(5,4)=3
       nar(5,5)=10
       nar(5,6)=6
       nar(5,7)=8
       nar(5,8)=24
       nar(5,9)=7
       nar(6,1)=7
       nar(6,2)=8
       nar(6,3)=10
       nar(6,4)=11
       nar(6,5)=13
       nar(6,6)=14
       nar(6,7)=16
       nar(6,8)=35
       nar(6,9)=15
       nar(7,1)=8
       nar(7,2)=9

       nar(7,3)=13
       nar(7,4)=7
       nar(7,5)=15
       nar(7,6)=16
       nar(7,7)=18
       nar(7,8)=35
       nar(7,9)=17
       nar(8,1)=9
       nar(8,2)=10
       nar(8,3)=7
       nar(8,4)=8
       nar(8,5)=17
       nar(8,6)=18
       nar(8,7)=15
       nar(8,8)=40
       nar(8,9)=19
       nar(9,1)=10
       nar(9,2)=11
       nar(9,3)=8
       nar(9,4)=9
       nar(9,5)=19
       nar(9,6)=20
       nar(9,7)=17
       nar(9,8)=45
       nar(9,9)=16
       nar(10,1)=11
       nar(10,2)=12
       nar(10,3)=14
       nar(10,4)=15
       nar(10,5)=21
       nar(10,6)=22
       nar(10,7)=24
       nar(10,8)=55
       nar(10,9)=23
       nar(11,1)=12
       nar(11,2)=13
       nar(11,3)=15
       nar(11,4)=11
       nar(11,5)=23
       nar(11,6)=24
       nar(11,7)=26
       nar(11,8)=55
       nar(11,9)=25
       nar(12,1)=13
       nar(12,2)=14
       nar(12,3)=11
       nar(12,4)=12
       nar(12,5)=25
       nar(12,6)=26
       nar(12,7)=23
       nar(12,8)=60
       nar(12,9)=27
       nar(13,1)=14
       nar(13,2)=15
       nar(13,3)=12
       nar(13,4)=13
       nar(13,5)=27
       nar(13,6)=28
       nar(13,7)=25
       nar(13,8)=65
       nar(13,9)=24
       nar(14,1)=15
       nar(14,2)=11
       nar(14,3)=13
       nar(14,4)=14
       nar(14,5)=29
       nar(14,6)=25
       nar(14,7)=27
       nar(14,8)=70
       nar(14,9)=26

       mm=130
       nn=120
       nzz=nn
       nunit=mm
       nzz2=4*mm
!       open (unit=2,file='nexexp2',access='direct',form=&
!       'formatted',recl=nzz2,status='old')
      open(unit=1,file='zfile2',access='direct',form=&
      'formatted',recl=213,status='old') 
       read (1,499,rec=40)kia,kib,isgn,(nused(i),i=1,70),&
       (nar(1,j),j=1,120),(izzy(j1),j1=1,10)
      print *,(nar(1,jf),jf=1,120)
      



        
      nardiv=1 
      kmx1=70
      kmx2=120
499   format (i6,i6,1i1,70i1,120i1,10i1)      
       
       
702   print *,'nzz',nzz,'nardiv',nardiv,'nzz2',nzz2
      

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
!       do k=1,2
       print *,'k=',k,'ir=',ir
       read (1,499,rec=k)kia,kib,isgn,(nused(i),i=1,kmx1),&
       (nar(nk1,j),j=1,kmx2),(izzy(j1),j1=1,10)
       iqi(nk1,1)=nar(nk1,1)
       do jk=1,nn-1
       iqi(nk1,jk+1)=nar(nk1,jk+1)
       end do
       nk1=nk1+1
!       print *,'nk1',nk1,'iend',iend,'littr',(littr(jd),jd=1,20),'k',k
       end do
       do jf=1,nunit
       kct(jf)=0
       end do
       do k=1,nunit
       issw=0
       do j=nn,1,-1
       if (iqi(k,j).eq.0)goto 1000
       if (ic(j).ne.0)goto 1000
       issw=1
       if (abs(iqi(k,j)).eq.1)goto 34
1000   end do
       if (issw.eq.0)goto 40
       imin=0
       ibmin=0
       irmin=100000000
       print *,'k',k,'iqi',(iqi(k,jf),jf=1,nn),'ic',(ic(jf),jf=1,nn)
       do j=nn,2,-1
!       print *,'funnyj',j
       if (iqi(k,j).eq.0)goto 101
       if (ic(j).ne.0)goto 101
       jnnn=j
       ibig=abs(iqi(k,j))
       
       do jb=j-1,1,-1
       print *,'j',j,'jb',jb,'k',k,'iqi',iqi(k,jb)
       if (iqi(k,jb).eq.0)goto 102
       if (ic(jb).ne.0)goto 102
       
       little=abs(iqi(k,jb))
       print *,'big',ibig,'little',little
       call subbw2(ibig,little,igcd2)
       print *,'j',j,'jb',jb,'k',k,'igcd2',igcd2
       
       if (igcd2.eq.1)goto 111
       if (igcd2.ge.irmin)goto 102
       imin=j
       ibmin=jb
       irmin=igcd2
       print *,'imin',imin,'ibmin',ibmin,'k',k,'irmin',irmin
102    end do
101    end do
       print *,'imin',imin
       j=jnnn

       if (imin.eq.0)goto 160
       print *,'k',k,'irmin',irmin
       ia=iqi(k,imin)/irmin
       ib=iqi(k,ibmin)/irmin
       
       goto 1111
160    do jf=2,nn
       if ((iqi(k,jf).ne.0).and.(ic(jf).eq.0))goto 34
       
161    end do
       j=1
       goto 34
111    imin=j
       ibmin=jb
       ia=iqi(k,imin)
       ib=iqi(k,ibmin)
       
1111   call subbw6(ia,ib,ivv,iu)
       print *,'ia',ia,'ib',ib,'ivv',ivv,'iu',iu
       do kkv=1,nunit
       iqi(kkv,imin)=iqi(kkv,imin)*iu+iqi(kkv,ibmin)*ivv
       end do
       little=1
       do kkv=1,nunit
       if (iqi(kkv,imin).eq.0)goto 120
       ibig=abs(iqi(kkv,imin))
       call subbw2(ibig,little,igcd2)
       litttle=igcd2
       if (igcd2.eq.1)goto 122
120    end do
       
       do kkv=1,nunit
       iqi(kkv,imin)=iqi(kkv,imin)/igcd2
       end do
122    a=a
       j=imin



       goto 34
110    a=a

       
       
       
       
       
1      format(i1,i10,20i4,50000i4)       
       
       
34     jpiv(k)=j    
       
       do jj=1,nn
       if(jj.eq.j)goto 38
       if (ic(jj).ne.0)goto 38
       if (iqi(k,jj).eq.0)goto 38
!       kct(k)=kct(k)+1
!       ktem=kct(k)
!       iqi2(k,ktem)=jj
!       jqi2(k,ktem)=iqi(k,jj)
       mul=iqi(k,jj)/iqi(k,j)
       mul=mul*(-1)
       do i=k,nunit
       itt = iqi(i,jj) + iqi(i,j)*mul 
       iqi(i,jj) =itt
!       if (iqi(i,jj).ge.0)goto 36
!       iqi(i,jj)=iqi(i,jj)+ipd

36     end do
       little=1
       do kkv=k,nunit
       if (iqi(kkv,jj).eq.0)goto 130
       ibig=abs(iqi(kkv,jj))
       call subbw2(ibig,little,igcd2)
       litttle=igcd2
       if (igcd2.eq.1)goto 132
130    end do
       
       do kkv=k,nunit
       iqi(kkv,jj)=iqi(kkv,jj)/igcd2
       end do
132    a=a






38     end do  
       knd=(nnn-1)*nunit+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     a=a
       issw=0
       
       
       
!       do is =1,mm
!       iv(is) = 0
!       end do
       knd=(nnn-1)*nunit+k
       ir=ir+1
!       iv(knd) =1
       irno(ir)=knd
       

!       do is =1,mm

!       jv = indic(is)
!       if (jv.eq.0)goto 42
!       iv(is) =iqi(k,jv)
!42     end do
!       if (iv(1).eq.0)goto 610
!       ir=ir+1
       print *,'vector found','ir',ir,'knd',knd
       
!       write (2,2,rec=ir)(iv(nq),nq=1,mm)
2      format(50000i4)     
       
       if (ir.lt.60)goto 50
       goto 100
610    a=a
       goto 50

50     end do
       end do
       
100    print *,'ir=',ir,'ic',(ic(jf),jf=1,nn),'indic',(indic(jf),jf=1,14)
       do i=1,irno(ir)
!       do i=1,irno(ir)
       print *,'iqi',(iqi(i,jf),jf=1,nn)
       end do
      close(unit=1) 
!       do k=1,1
       do k=1,ir
       do k1=1,irno(ir)
       iv(k1)=0
       end do
       ileng=irno(k)-k+1
       indk=0
       do kk=1,irno(k)-1
       if (indic(kk).eq.0)goto 140
       indk=indk+1
       iset(indk)=kk
       

140    end do 
       indk=indk+1
       iset(indk)=irno(k)
       iresn=iset(ileng)
       iresn2=iset(ileng-1)
!       print *,'iresn',iresn,'iresn2',iresn2,'ileng',ileng
!       print *,'iset',(iset(jf),jf=1,ileng)
!       do kk2=1,ileng-1
       ind=indic(iresn2)
!       print *,'ind',ind
       
       iv(iresn)=iqi(iresn,ind)
       iv(iresn2)=iqi(iresn2,ind)
       ibig=abs(iv(iresn))
       little=abs(iv(iresn2))
       call subbw2(ibig,little,igcd2)
       itemp=iv(iresn)/igcd2
       iv(iresn2)=iv(iresn2)/igcd2
       iv(iresn)=iv(iresn2)
       iv(iresn2)=(itemp*iv(iresn2))/iv(iresn2)*(-1)
       do kk2=1,ileng-2
       
       ind2=iset(ileng-kk2-1)
       im=indic(ind2)
       ires=0
!       print *,'im',im



       do jz=ileng-kk2,ileng
       jj=iset(jz)
!       print *,'jj',jj,'firiv',iv(jj),'ind2',ind2,'iqi',iqi(ind2,im)
       ires=ires+iv(jj)*iqi(jj,im)
       end do
       ibig=abs(iqi(ind2,im))
       little=abs(ires)
       print *,'ibig',ibig,'little',little
       if (little.eq.0)goto 142
       call subbw2(ibig,little,igcd2)
!       print *,'igcd2',igcd2
       ires=ires/igcd2
       iv(ind2)=iqi(ind2,im)/igcd2
       do jz=ileng-kk2,ileng
       jj=iset(jz)
!       print *,'jj',jj,'iv',iv(jj)
       iv(jj)=iv(jj)*iv(ind2)
       end do
       iv(ind2)=(iv(ind2)*ires)/iv(ind2)*(-1)
       goto 143
142    iv(ind2)=0
143    end do
       
       print *,'iv',(iv(jf),jf=1,irno(ir))
       end do
       isum=0
       do i=1,irno(ir)
       isum=isum+iqi(i,2)*iv(i)
!       print *,'isumind',isum,'i',i,'iqi',iqi(i,4),'iv',iv(i)
       if (isum.eq.0)goto 99
!       stop
99     end do
       print *,'isum',isum
       
       
       
       stop
! end of copying
       
       iresn=irno(1)
       
       ind=indic(iresn-1)
       iv(iresn)=iqi(iresn,ind)
       iv(iresn-1)=iqi(iresn-1,ind)
       ibig=abs(iv(iresn))
       little=abs(iv(iresn-1))
       call subbw2(ibig,little,igcd2)
       itemp=iv(iresn)/igcd2
       iv(iresn-1)=iv(iresn-1)/igcd2
       iv(iresn)=iv(iresn-1)
       iv(iresn-1)=(itemp*iv(iresn-1))/iv(iresn-1)*(-1)
       do j=1,iresn-2
       ind=iresn-j
       ind2=ind-1
       im=indic(ind2)
       ires=0
       do jj=ind,iresn
       ires=ires+iv(jj)*iqi(jj,im)
       end do
       ibig=abs(iqi(ind2,im))
       little=abs(ires)
       call subbw2(ibig,little,igcd2)
       ires=ires/igcd2
       iv(ind2)=iqi(ind2,im)/igcd2
       do jj=ind,iresn
       iv(jj)=iv(jj)*iv(ind2)
       end do
       iv(ind2)=(iv(ind2)*ires)/iv(ind2)*(-1)
       end do
       
       print *,'iv',(iv(jf),jf=1,irno(1))
       
       end

       subroutine subbw6(ia,ib,ivv,iu)
!       ib=mod(ib,ia)
!       if (ib.ge.0)goto 10
!       ib=ia+ib
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
       goto 870
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

       subroutine subbw2(ibig,little,igcd2)
918    itemp=int(ibig/little)
       irem1=ibig-itemp*little
       if (irem1.eq.0)goto 940
       ibig=little
       little=irem1
       goto 918
940    igcd2=little
       return
       end
