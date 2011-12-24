      program bopbexq3
!     finds kernel of integer matrix reducing basis as it proceeds
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000),irrr(2000)
      common mnum(50)
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      dimension jfak(30),jfreq(30),isumt(50),itt(50),iuu(100)
      dimension itempar(100),ibase(50),jbase(50),itarg(50),jtarg(50)
      open (unit=3,file='exppar',access='direct',form=&
      'formatted',recl=344,status='old')
      read (3,1003,rec=1)nn,mm,ijpow,nfak,(ibase(jf),jf=1,20),&
      (jbase(jf),jf=1,20),(itarg(jf),jf=1,20),(jtarg(jf),jf=1,20)
1003  format (i6,i6,i8,i4,20i4,20i4,20i4,20i4)      
      print *,'nn',nn,'mm',mm
      
      mmr=10*mm+91
      nnr=100*nn
      llr=2*nnr
      kkb=0
      ihitn=0
      do jf=1,20
      littr(jf)=0
      end do
      
      open(unit=1,file='zbexp1',access='direct',form=& 
      'formatted',recl=mmr,status='old')
      
      open(unit=2,file='zhexp1',access='direct',form=&
      'formatted',recl=nnr,status='old')
     
     


!      close(unit=2,status='delete') 
      open (unit=4,file='lexp1',access='direct',form=&
      'formatted',recl=llr,status='old')
      close (unit=1,status='delete') 
      close (unit=2,status='delete')
      close (unit=4,status='delete')
      stop



!      read (1,1001,rec=6)kkb,ihitn,(littr(jk),jk=1,20),(ibbr(jf),jf=1,mm)
!      print *,'ibbr 6',(ibbr(jf),jf=1,20)
      
!      goto 950
      do j=1,nn

      do i=1,nn 
      ihhl(i,1)=0
      ihhl(i,2)=0
      ihhl(i,3)=0
      end do
      ihhl(j,1)=0
      ihhl(j,2)=1
      ihhl(j,3)=1
      do ii=4,25
      ihhl(i,ii)=0
      end do
      
      
      
      
      
      write (2,1002,rec=j)((ihhl(ir,jf),jf=1,25),ir=1,nn)
1001  format(i1,i10,20i4,1400i10)
1002  format(1400(25i4))      
      end do
950   k=2
      kmax=1
      idd(1,1)=0
      idd(1,2)=1
      idd(1,3)=1
      ittt=0
      llam=0
      isumtt=0
      itt(1)=0
      itt(2)=0
      lamdar(1,1)=0
      lamdar(1,2)=0
      read (1,1001,rec=1)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      read (1,1001,rec=k)kkb,ihitn,(littr(jk),jk=1,20),(ibbr(jf),jf=1,mm)
      print *,'ok1','ibbl15',ibbl(5)
      do i=1,mm
      marr(1)=0
      if (ibbl(i).ge.0)goto 804
      marr(1)=1
804   a=a
      if (abs(ibbl(i)).ge.10000)goto 900
      marr(2)=1
      if (ibbl(i).ne.0)goto 8041
      marr(2)=0


8041  marr(3)=abs(ibbl(i))
      mbarr(1)=0
      goto 902
900   marr(2)=2      
      marr(3)=int(abs(ibbl(i))/10000)
      marr(4)=abs(ibbl(i))-marr(3)*10000
      mbarr(1)=0
      
902   if (ibbr(i).ge.0)goto 806
      mbarr(1)=1
806   a=a
      if (abs(ibbr(i)).ge.10000)goto 904

      mbarr(2)=1      
      if (ibbr(i).ne.0)goto 8042
      mbarr(2)=0
      
8042  mbarr(3)=abs(ibbr(i))
      goto 906
904   mbarr(2)=2
      mbarr(3)=int(abs(ibbr(i))/10000)
      mbarr(4)=abs(ibbr(i))-mbarr(3)*10000

906   call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,lamdar(1,2)+2
      kbarr(jf)=lamdar(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      lamdar(1,jf)=kcarr(jf)
      end do
      if (i.gt.15)goto 861
      print *,'i',i,'ibbli',ibbl(i),'ibbri',ibbr(i),'llam',llam
861   ittt=ittt+ibbl(i)*ibbl(i)
      isumtt=isumtt+ibbr(i)*ibbr(i)
      llam=llam+ibbl(i)*ibbr(i)
      end do
      ind=lamdar(1,2)+3
      do jf=ind,50
      lamdar(1,jf)=0
      end do
      do ii=2,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=2)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      print *,'ittt',ittt,'isumtt',isumtt,'lamdar1',(lamdar(1,jf),&
      jf=1,lamdar(1,2)+2),'llam',llam
      
      itt(1)=0
      itt(2)=1
      itt(3)=ittt
      if (ittt.eq.0)goto 86
      if (itt(2).eq.0)goto 86
      do jf=1,itt(2)+2
      idd(2,jf)=itt(jf)
      end do

      iff(1)=1
      goto 90
86    idd(2,1)=0
      idd(2,2)=1
      idd(2,3)=1
      itt(2)=0
      iff(1)=0

90    a=a      
      
      
      
      
      print *,'ok2','ibbr25',ibbr(5)
      
      
      print *,'ok211','lamdar1',(lamdar(1,jk),jk=1,50),'isumtt',isumtt
      print *,'iff1',iff(1)
      if (iff(1).eq.0)goto 700
      marr(1)=0
      marr(2)=1
      marr(3)=isumtt
      
      if (isumtt.ne.0)goto 901
      marr(2)=0
!      do jf=1,isumt(2)+2
!      marr(jf)=isumt(jf)
!      end do
901   do jf=1,idd(2,2)+2
      mbarr(jf)=idd(2,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      idd(3,jf)=mcarr(jf)
      end do
      do jf=1,lamdar(1,2)+2
      marr(jf)=lamdar(1,jf)
      mbarr(jf)=lamdar(1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,idd(3,2)+2
      karr(jf)=idd(3,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      idd(3,jf)=kcarr(jf)
      end do
      print *,'idd3',(idd(3,jf),jf=1,idd(3,2)+2)






100   if (idd(3,2).eq.0)goto 106
      iff(2)=1
      goto 108
700   idd(3,1)=0
      idd(3,2)=1
      idd(3,3)=isumtt

!700   do jf=1,isumt(2)+2
!      idd(3,jf)=isumt(jf)
!      end do
      
      goto 100

106   iff(2)=0
!      idd(3)=idd(2)
      do jf=1,idd(2,2)+2
      idd(3,jf)=idd(2,jf)
      end do
108   if (iff(1).eq.0)goto 119
      lpar=1
      print *,'ok2115'
      call subred(k,lpar,mm,nn)
      print *,'ok212'
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=ind,50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok22'



1004  format(1400(50i4))      
      
      
      
      
      if (iff(1).eq.0)goto 119
      if (iff(2).eq.0)goto 118
      goto 119
118   call subswap(k,mm,nn,kmax)
119   kmax=2
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      
      k=3
      print *,'whats going on here'
!      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)


123   if (k.le.kmax)goto 300
      print *,'ok23',' k',k,'ibbr',(ibbr(jk),jk=1,20)
      kmax=k
!      lamdar(k,1)=0
      read (1,1001,rec=1)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      print *,'read first ok'
      read (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      print *,'k',k,'ibbr5',ibbr(5),'ibbl5',ibbl(5)
      lamdar(1,1)=0
      lamdar(1,2)=0
      do i=1,mm
      marr(1)=0
      if (ibbl(i).ge.0)goto 860
      marr(1)=1
860   if (abs(ibbl(i)).ge.10000)goto 960
      marr(2)=1
      marr(3)=abs(ibbl(i))
      mbarr(1)=0
      goto 962
960   marr(2)=2
      marr(3)=int(abs(ibbl(i))/10000)
      marr(4)=abs(ibbl(i))-marr(3)*10000
      mbarr(1)=0
962   if (ibbr(i).ge.0)goto 862      
      mbarr(1)=1
862   if (abs(ibbr(i)).ge.10000)goto 964      
      mbarr(2)=1
      mbarr(3)=abs(ibbr(i))
      goto 966
964   mbarr(2)=2      
      mbarr(3)=int(abs(ibbr(i))/10000)
      mbarr(4)=abs(ibbr(i))-mbarr(3)*10000
966   call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,lamdar(1,2)+2
      kbarr(jf)=lamdar(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      lamdar(1,jf)=kcarr(jf)
      end do
      end do






!      do i=1,mm
!      lamdar(k,1)=lamdar(k,1)+ibbr(i)*ibbl(i)
!      end do
!      print *,'lamdark1',lamdar(k,1)
      do j=2,k 
       
      
      if ((iff(j).eq.0).and.(j.lt.k))goto 152
      iuu(1)=0
      iuu(2)=0
      print *,'j',j,'and here is where itbombs out'
      read (1,1001,rec=j)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      print *,'ok3','j',j
      do i=1,mm
      
!      iuu=iuu+ibbl(i)*ibbr(i)
      marr(1)=0
      if (ibbl(i).ge.0)goto 808
      marr(1)=1
808   if (abs(ibbl(i)).ge.10000)goto 920


      marr(2)=1
      if (ibbl(i).ne.0)goto 8084
      marr(2)=0
      
8084  marr(3)=abs(ibbl(i))
      mbarr(1)=0
      goto 922
920   marr(2)=2
      marr(3)=int(abs(ibbl(i))/10000)
      marr(4)=abs(ibbl(i))-marr(3)*10000
      mbarr(1)=0



922   if (ibbr(i).ge.0)goto 810
      mbarr(1)=1
810   if (abs(ibbr(i)).ge.10000)goto 930

      

      mbarr(2)=1
      if (ibbr(i).ne.0)goto 8085
      mbarr(2)=0
      
8085  mbarr(3)=abs(ibbr(i))
      goto 932
930   mbarr(2)=2
      mbarr(3)=int(abs(ibbr(i))/10000)
      mbarr(4)=abs(ibbr(i))-mbarr(3)*10000
      
932   call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,iuu(2)+2
      kbarr(jf)=iuu(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iuu(jf)=kcarr(jf)
      end do
      
      end do
      print *,'about to stop','j=',j ,'k',k,'startiuu',&
      
      (iuu(jf),jf=1,iuu(2)+2)
      
      if (j.eq.k)goto 159
      read (4,1004,rec=j)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      
      goto 160
159   do ii=1,nn
      do jf=1,50
      lamdar2(ii,jf)=lamdar(ii,jf)
      end do
      end do
      goto 160

!152   lamdar(k,j)=0
152   lamdar(j,1)=0      
      lamdar(j,2)=0
      
      
      
      
      goto 200
160   do i=1,j-1
      if (iff(i).eq.0)goto 166
      do jf=1,lamdar(i,2)+2
      marr(jf)=lamdar(i,jf)
      end do
      do jf=1,lamdar2(i,2)+2
      mbarr(jf)=lamdar2(i,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      itempar(jf)=mcarr(jf)
      end do
      do jf=1,iuu(2)+2
      marr(jf)=iuu(jf)
      end do
      do jf=1,idd(i+1,2)+2
      mbarr(jf)=idd(i+1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,itempar(2)+2
      kbarr(jf)=itempar(jf)
      end do
      call mpadd(1)
      print *,'karr',(karr(jf),jf=1,karr(2)+2),'kbarr',(kbarr(jf),&
      jf=1,kbarr(2)+2)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,idd(i,2)+2
      mbarr(jf)=idd(i,jf)
      end do
      call mendiv
      if (mcarr(2).ne.0)goto 811
      do jf=1,mdarr(2)+2
      iuu(jf)=mdarr(jf)
      end do
      print *,'inter iuu',(iuu(jk),jk=1,iuu(2)+2),'iddi',(idd(i,jk),&
      jk=1,idd(i,2)+2),'marr',(marr(jk),jk=1,marr(2)+2),'mbarr',&
      (mbarr(jk),jk=1,mbarr(2)+2)
      goto 166
811   print *,'iuu in error',' num',(marr(jk),jk=1,marr(2)+2),'denom',&
      (mbarr(jk),jk=1,mbarr(2)+2)
      stop
166   end do



!      iuu=iuu*idd(i+1)-lamdar(k,i)*lamdar(j,i)
!      print *,'iuu',iuu,'iddi',idd(i)
!      iuu=iuu/idd(i)
!166   end do
      print *,'iuu',(iuu(jf),jf=1,iuu(2)+2)
      if (j.lt.k)goto 190
      if (iuu(2).eq.0)goto 180
      do jf=1,iuu(2)+2
      idd(k+1,jf)=iuu(jf)
      end do
      
!      idd(k+1)=iuu
      iff(k)=1
      goto 200
!180   idd(k+1)=idd(k)
180   do jf=1,idd(k,2)+2      
      idd(k+1,jf)=idd(k,jf)
      end do
      iff(k)=0
      goto 200
!190   lamdar(k,j)=iuu
190   do jf=1,iuu(2)+2
      lamdar(j,jf)=iuu(jf)
      end do


200   end do
      

300   if (iff(k-1).eq.0)goto 320
      lpar=k-1
      print *,'lok300 lpar',lpar
      call subred(k,lpar,mm,nn)
      if (iff(k).eq.0)goto 310
      goto 320
310   call subswap(k,mm,nn,kmax)
      if (k-1.gt.2)goto 318
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=ind,50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      
      k=2
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      print *,'stop k2'
      goto 300
318   a=a
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=ind,50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      
     k=k-1
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      print *,'stop k-1',' k' ,k
      goto 300
320   ll=k-2
340   if (ll.lt.1)goto 350
      lpar=ll
      call subred(k,lpar,mm,nn)
      ll=ll-1
      goto 340

      

350    a=a      
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=ind,50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      print *,'i think problem is here',' k', k 
      
      k=k+1
!      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,nn)
      if (k.le.nn)goto 123
      read (2,1002,rec=1)((ihhl(ir,jf),jf=1,25),ir=1,nn)
      
      print *,'vector1',((ihhl(ir,jf),jf=1,25),ir=1,nn)
      close (unit=1)
      close(unit=2)
      close(unit=3)
      end

      subroutine subred(k,lpar,mm,nn)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000),irrr(2000)
      common mnum(50)
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      dimension itempar(50),iqq(50)
      do jf=1,lamdar(lpar,2)+2
      karr(jf)=lamdar(lpar,jf)
      kbarr(jf)=lamdar(lpar,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      itempar(jf)=kcarr(jf)
      end do
      do ii=2,idd(lpar+1,2)+2
      if (itempar(ii).lt.idd(lpar+1,ii))goto 490
      if (itempar(ii).gt.idd(lpar+1,ii))goto 100
      end do
      
!      if (abs(2*lamdar(k,lpar)).le.idd(lpar+1))goto 490
100   a=a      
      
!      iqq=int((2*lamdar(k,lpar)+idd(lpar+1))/(2*idd(lpar+1)))
      read (2,1002,rec=lpar)((ihhl(ir,jf),jf=1,25),ir=1,nn)
      read (2,1002,rec=k)((ihhr(ir,jf),jf=1,25),ir=1,nn)
      read (1,1001,rec=lpar)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      read (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
1001  format(i1,i10,20i4,1400i10)
1002  format(1400(25i4))
1004  format (1400(50i4))      
      print *,'okred1','idd2',(idd(2,jf),jf=1,idd(2,2)+2),'lpar',lpar
      kkb=0
      ihitn=0
!      do i=1,nn
!      ihhr(i)=ihhr(i)-iqq*ihhl(i)
!      end do
!      do i=1,mm
!      ibbr(i)=ibbr(i)-iqq*ibbl(i)
!      end do
!      write (2,1002,rec=k)(ihhr(jf),jf=1,nn)
!      write (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
!      iqq=int((2*lamdar(k,lpar)+idd(lpar+1))/(2*idd(lpar+1)))
      do jf=1,lamdar(lpar,2)+2
      karr(jf)=lamdar(lpar,jf)
      kbarr(jf)=lamdar(lpar,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,idd(lpar+1,2)+2
      kbarr(jf)=idd(lpar+1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,idd(lpar+1,2)+2
      karr(jf)=idd(lpar+1,jf)
      kbarr(jf)=idd(lpar+1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      mbarr(jf)=kcarr(jf)
      end do
      call mendiv
      if (mdarr(1).eq.0)goto 10
      do jf=1,mdarr(2)+2
      karr(jf)=mdarr(jf)
      end do
      kbarr(1)=1
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iqq(jf)=kcarr(jf)
      end do
      goto 11
10    do jf=1,mdarr(2)+2
      iqq(jf)=mdarr(jf)
      end do
11    print *,'k',k,'lpar',lpar,'iqq',(iqq(jf),jf=1,iqq(2)+2)
      do jf=1,iqq(2)+2
      marr(jf)=iqq(jf)
      end do
      do jf=1,idd(lpar+1,2)+2
      mbarr(jf)=idd(lpar+1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,lamdar(lpar,2)+2
      karr(jf)=lamdar(lpar,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      lamdar(lpar,jf)=kcarr(jf)
      end do
      if (lpar.eq.1)goto 109
      read (4,1004,rec=lpar)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      
      do ii=1,mm
      do jf=1,lamdar2(ii,2)+2
      marr(jf)=lamdar2(ii,jf)
      end do
      do jf=1,iqq(2)+2
      mbarr(jf)=iqq(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,lamdar(ii,2)+2
      karr(jf)=lamdar(ii,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      lamdar(ii,jf)=kcarr(jf)
      end do
      if (ii.eq.lpar-1)goto 109
      end do
      
      
      
109   if (iqq(2).gt.2)goto 179
      if (iqq(2).eq.0)goto 1221
      if (iqq(2).eq.1)goto 110
      iqm=iqq(3)*10000+iqq(4)
      goto 112
110   iqm=iqq(3)
112   if (iqq(1).eq.1)goto 120
      goto 122
120   iqm=(-1)*iqm
122   goto 124
1221  iqm=0
      goto 124
179   print *,'iqq too large length=',iqq(2)
      stop
124   a=a
!      do i=1,nn
!      ihhr(i)=ihhr(i)-iqm*ihhl(i)
!      end do
      ibmax=0
      do i=1,mm
      ibbr(i)=ibbr(i)-iqm*ibbl(i)
      if (abs(ibbr(i)).le.ibmax)goto 500
      ibmax=abs(ibbr(i))
500   end do
!      write (2,1002,rec=k)(ihhr(jf),jf=1,nn)
      write (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      print *,'okred2',' k',k,'iqm',iqm ,'ibmax',ibmax
      do i=1,nn
      do jf=1,iqq(2)+2
      marr(jf)=iqq(jf)
      end do
      do jf=1,ihhl(i,2)+2
      mbarr(jf)=ihhl(i,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,ihhr(i,2)+2
      karr(jf)=ihhr(i,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ihhr(i,jf)=kcarr(jf)
      end do
      end do
      do ii=1,nn
      ind=ihhr(ii,2)+3
      do jf=ind,25
      ihhr(ii,jf)=0
      end do
      end do
      write (2,1002,rec=k)((ihhr(ir,jf),jf=1,25),ir=1,nn)
      print *,'okred 3',' k=',k
      





      






490   a=a
      print *,'okred 4',' k',k,'lpar',lpar
      return
      end


      subroutine subswap(k,mm,nn,kmax)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000),irrr(2000)
      common mnum(50)
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      dimension itempar(100),lamda(50),itt(50)
      
      


1001  format (i1,i10,20i4,1400i10) 
1002  format (1400(25i4))
      kkb=0
      ihitn=0
      read (2,1002,rec=k)((ihhr(ir,jf),jf=1,25),ir=1,nn)
      read (2,1002,rec=k-1)((ihhl(ir,jf),jf=1,25),ir=1,nn)
      read (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      read (1,1001,rec=k-1)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      write (2,1002,rec=k)((ihhl(ir,jf),jf=1,25),ir=1,nn)
      write(2,1002,rec=k-1)((ihhr(ir,jf),jf=1,25),ir=1,nn)
      write (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      write (1,1001,rec=k-1)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      print *,'ok swap0 k=',k
      if (k.lt.3)goto 522
      read (4,1004,rec=k-1)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok swap 1',' k', k
      
1004  format(1400(50i4))      
      if (k.lt.3)goto 522
      do j=1,k-2
      do jf=1,lamdar(j,2)+2
      itempar(jf)=lamdar(j,jf)
      end do
      do jf=1,lamdar2(j,2)+2
      lamdar(j,jf)=lamdar2(j,jf)
      end do
      do jf=1,itempar(2)+2
      lamdar2(j,jf)=itempar(jf)
      end do
      end do
      
      do ii=1,k-2
      ind=lamdar2(ii,2)+3
      do jf=ind,50
      lamdar2(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k-1)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok swap 2 k',k




!      itemp=lamdar(k,j)
!      lamdar(k,j)=lamdar(k-1,j)
!      lamdar(k-1,j)=itemp
!      end do
! 522   lamda=lamdar(k,k-1)
522   do jf=1,lamdar(k-1,2)+2      
      lamda(jf)=lamdar(k-1,jf)
      end do
      
      if (lamda(2).eq.0)goto 560
      if (k+1.gt.kmax)goto 537
      do i=k+1,kmax
      do jf=1,lamda(2)+2
      marr(jf)=lamda(jf)
      end do
      read (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok swap 3 i=',i
      do jf=1,lamdar2(k-1,2)+2
      mbarr(jf)=lamdar2(k-1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,idd(k,2)+2
      mbarr(jf)=idd(k,jf)
      end do
      call mendiv
      if (mcarr(2).ne.0)goto 811
      do jf=1,mdarr(2)+2
      lamdar2(k-1,jf)=mdarr(jf)
      end do
      write (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok swap4 wr i=',i
      goto 812
811   print *,'error in subswap type 1'
      stop
812   end do
!      lamdar(i,k-1)=(lamda*lamdar(i,k-1))/idd(k)
!      end do
       do jf=1,idd(k+1,2)+2
       itt(jf)=idd(k+1,jf)
       end do
       do jf=1,lamda(2)+2
       marr(jf)=lamda(jf)
       mbarr(jf)=lamda(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,idd(k,2)+2
       mbarr(jf)=idd(k,jf)
       end do
       call mendiv
       if (mcarr(2).ne.0)goto 813
       do jf=1,mdarr(2)+2
       idd(k,jf)=mdarr(jf)
       end do
       goto 814
813    print *,'error in subswap type 2'
       stop
814    a=a
       do jf=1,idd(k,2)+2
       idd(k+1,jf)=idd(k,jf)
       end do

!      itt=idd(k+1)
!      idd(k)=(lamda*lamda)/idd(k)
!      idd(k+1)=idd(k)
537    if (k+2.gt.kmax)goto 545       
       do i=k+2,kmax
       do j=k+1,i-1

!      do j=k+1,kmax-1
!      do i=j+1,kmax
      read (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'ok swap 5 re i=',i
      do jf=1,lamdar2(j,2)+2
      marr(jf)=lamdar2(j,jf)
      end do
      do jf=1,idd(k,2)+2
      mbarr(jf)=idd(k,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,itt(2)+2
      mbarr(jf)=itt(jf)
      end do
      call mendiv
      if (mcarr(2).ne.0)goto 815
      do jf=1,mdarr(2)+2
      lamdar2(j,jf)=mdarr(jf)
      end do
      write (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      
      goto 816
815   print *,'swapsub error type 3'
      stop
!      lamdar(i,j)=(lamdar(i,j)*idd(k))/itt
816   end do
      end do
545   if (k+1.gt.kmax)goto 590      
      do j=k+1,kmax
      do jf=1,idd(j+1,2)+2
      marr(jf)=idd(j+1,jf)
      end do
      do jf=1,idd(k,2)+2
      mbarr(jf)=idd(k,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,itt(2)+2
      mbarr(jf)=itt(jf)
      end do
      call mendiv
      if (mcarr(2).ne.0)goto 817
      do jf=1,mdarr(2)+2
      idd(j+1,jf)=mdarr(jf)
      end do
      goto 818
817   print *,'swapsub error type 4'
      stop



!      idd(j+1)=(idd(j+1)*idd(k))/itt
818   end do
      goto 590
560   do jf=1,idd(k-1,2)+2
      idd(k,jf)=idd(k-1,jf)
      end do


!560   idd(k)=idd(k-1)
      itemp=iff(k)
      iff(k)=iff(k-1)
      iff(k-1)=itemp
!      lamdar(k,k-1)=0
      lamdar(k-1,1)=0
      lamdar(k-1,2)=0
      if (k+1.gt.kmax)goto 590
      do i=k+1,kmax
      read (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn)
      print *,'okswap 6 re i=',i
      do jf=1,lamdar2(k-1,2)+2
      lamdar2(k,jf)=lamdar2(k-1,jf)
      lamdar2(k-1,jf)=0
      end do
!      lamdar(i,k)=lamdar(i,k-1)
!      lamdar(i,k-1)=0
!      end do
      
      write (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,nn) 
      print *,'ok swap7 wr i=',i
      end do
590   return


      end



      
      
      subroutine menmul
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000),irrr(2000)
      common mnum(50)
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      




      if ((marr(2).eq.0).or.(mbarr(2).eq.0))goto 9
      ilen=marr(2)
      ilen2=mbarr(2)
      do jf=1,ilen
      karr(jf)=marr(jf+2)
      end do
      do jf=1,ilen2
      kbarr(jf)=mbarr(jf+2)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      mcarr(jf+2)=kcarr(jf)
      end do
      mcarr(2)=ilen3
      mcarr(1)=mod(marr(1)+mbarr(1),2)
      goto 10
9     mcarr(1)=0
      mcarr(2)=0
10    return
      end
      subroutine mendiv
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000),irrr(2000)
      common mnum(50)
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      
      if (mbarr(2).eq.0)goto 11
      if (marr(2).eq.0)goto 9
      ilen=marr(2)
      ilen2=mbarr(2)
      do jf=1,ilen
      karr(jf)=marr(jf+2)
      end do
      do jf=1,ilen2
      kbarr(jf)=mbarr(jf+2)
      end do
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 10
      do jf=1,irlen
      mcarr(jf+2)=irrr(jf)
      end do
      mcarr(2)=irlen
      mcarr(1)=mod(marr(1)+mbarr(1),2)
6     if (icont.eq.0)goto 12
      do jf=1,icont
      mdarr(jf+2)=ipqt(jf)
      end do
      mdarr(2)=icont
      mdarr(1)=mod(marr(1)+mbarr(1),2)
      goto 15
11    print *,'halted: attempted division by zero, routine mendiv'
      stop
10    mcarr(1)=0
      mcarr(2)=0
      goto 6
9     mcarr(1)=0
      mcarr(2)=0
12    mdarr(1)=0
      mdarr(2)=0
15    return
      end





      

      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000)
      common irrr(2000)
      common mnum(50)
      
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      
      do i=1,ilen+ilen2
      kcarr(i) =0
      
      end do
      do i =1,ilen
      do j=1,ilen2
      itemp =karr(i) *kbarr(j)
      
      itemp2=int(itemp/10000)
      irem1 =itemp-itemp2 *10000
      
      
      kcarr(i+j)=kcarr(i+j) +irem1
      kcarr(i+j-1)=kcarr(i+j-1) +itemp2
      
      if (kcarr(i+j).lt.10000)goto 20
      kcarr(i+j)=kcarr(i+j)-10000
      kcarr(i+j-1) =kcarr(i+j-1)+1
20    do k =1,i+j-2
      if (kcarr(i+j-k).lt.10000)goto 22
      kcarr(i+j-k)=kcarr(i+j-k) -10000
      kcarr(i+j-k-1)=kcarr(i+j-k-1)+1
      end do
22    end do
      end do
      ilen3 =ilen +ilen2
      if(kcarr(1).ne.0)goto 100
      ilen3 =ilen3-1
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do
      
      
100   return
      end

      subroutine mpdiv(ilen,ilen2,irlen,icont,iswq)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000)
      common irrr(2000)
      common mnum(50)
      
      
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      
      
      
      
      
      dimension kdum(2000),isub(2000)
      if (kbarr(1).eq.0)goto 906
      if (ilen2.eq.0)goto 906
      do i =1,ilen
      kdum(i)=karr(i)
      ipqt(i)=0
      end do
      icont =0
      icont2 =0
      icong =0
      iswm=0
      iswq=1
      isw2=0
      ind2 =ilen2
      if(ilen2.ne.1)goto 8
      kbarr(2)=0
      ilen2=2
      ind2=2
      ilen=ilen+1
      karr(ilen)=0
      kdum(ilen)=0
      ipqt(ilen)=0
      isw2=1
8     if (ilen2.gt.ilen)goto 905
      if (ilen2.lt.ilen)goto 10
      if (kbarr(1).gt.karr(1))goto 905
      do i=1,ilen2
      if(kbarr(i).gt.karr(i))goto 905
      if(kbarr(i).lt.karr(i))goto 10
      end do
10    kbap=kbarr(1)
      lop =ilen +1 -ilen2
      ll=1
201   do i =1,ilen2
      if (kbarr(i).eq.kdum(ll+i-1))goto 202
      if (kbarr(i).gt.kdum(ll +i-1))goto 203
      if(i.gt.2)goto 2021
      goto 20
202   end do
2021  icont =ind2
      iapd =1
      goto 31

203   iapd =int((kdum(ll)*10000 +kdum(ll+1))/(kbap+1))
      goto 30
20    kbig = kdum(ll) *10000 +kdum(ll+1)
      kbap =kbap*10000 +kbarr(2) +1
      icont =ind2
      
      
      iapd=int(kbig/kbap)
      goto 31
      
      
30    icont =ind2 +1
      


31    ipqt(icont) =ipqt(icont)+iapd
      icont2 =icont2 +1
      do i=icont,2,-1
      if (ipqt(i).lt.10000)goto 35
      ipqt(i)=ipqt(i)-10000
      ipqt(i-1)=ipqt(i-1) +1
      end do

35    isub(1)=0
      do i=1,ilen2
      isub(i+1)= kbarr(i)*iapd
      end do
      do i=1,ilen2 
      ii=ilen2 +2 -i
      if (isub(ii).lt.10000)goto 37
      itemp =int(isub(ii)/10000)
      irem1 =isub(ii)-itemp *10000
      isub(ii) =irem1
      isub(ii-1) =isub(ii-1) +itemp
37    end do
      if(isub(1).ne.0)goto 38
      if (kdum(ll).lt.isub(2))goto 38
      ilens =ilen2
      do i =1,ilens
      isub(i) =isub(i+1)
      end do
      goto 39

      
38    ilens =ilen2 +1
39    do i =1,ilens
      ii = ilens +1 -i
      ig=icont-ilens
      
      kdum(ii+ig) =kdum(ii+ig) -isub(ii)
      if (kdum(ii+ig) .ge.0)goto 40
      kdum(ii+ig) =kdum(ii+ig) +10000
      kdum(ii+ig-1) =kdum(ii+ig-1) -1
40    end do
      icong =icong +1
      
      
      
   
402   do i =1,ilen
      if(kdum(i).ne.0)goto 70
      end do
      goto 90
      
150   kbap = kbarr(1)
      goto 201

70    if (ilen2 +i-1.gt.ilen)goto 90
      ll =i
      
      ind2 =ll +ilen2-1   
      if (ind2.lt.ilen)goto 150
      do i=1,ilen2
      if (kbarr(i).eq.kdum(ll+i-1))goto 71
      if (kbarr(i).gt.kdum(ll+i-1))goto 90
      if (kbarr(i).lt.kdum(ll+i-1))goto 150
71    end do
      goto 150
90    do i=1,ilen
      if (kdum(i).ne.0)goto 92
      end do
      irlen =0
      goto 921
92    irlen=ilen +1 -i
      do j=1,irlen
      irrr(j)=kdum(i+j-1)
      end do
921   do i=1,ilen
      if (ipqt(i).ne.0)goto 901
      end do
      icont =0
      goto 902
901   icont =ilen +1 -i
      do j=1,icont
      ipqt(j)=ipqt(j+i-1)
      end do
902   if (isw2.eq.0)goto 910
      ilen=ilen-1
      ilen2 =ilen2-1
      if (irlen.eq.0)goto 910
      irlen=irlen-1
      goto 910
905   do i=1,ilen
      irrr(i)=karr(i)
      end do
      irlen=ilen
      icont=0
      if (isw2.eq.1)goto 902
      goto 910
906   print *,'halted:attempted division by zero'
      stop
      
910   return      
      end





      subroutine mpadd(isora)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000)
      common irrr(2000)
      common mnum(50)
      
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      
      if(karr(2).eq.0)goto 50
      if(kbarr(2).eq.0)goto 55
      if(kbarr(2).gt.karr(2))goto 2
      kmx=0
      goto 4
2     kmx=1
4     if((karr(1).eq.kbarr(1)).and.(isora.eq.0))goto 30
      if((karr(1).ne.kbarr(1)).and.(isora.eq.1))goto 30
      if(karr(2).ne.kbarr(2))goto 22
      do i=3,karr(2) +2
      if(karr(i).lt.kbarr(i))goto 5
      if(karr(i).gt.kbarr(i))goto 24
      end do
      goto 26
22    if(kmx.eq.1)goto 5

24    kcarr(1)=karr(1)
      do i=2,karr(2)+2
      
      
      kcarr(i)=karr(i)
      end do
      
      
      ii=kbarr(2)+2
      do i=karr(2)+2,3,-1
      kcarr(i)=kcarr(i)-kbarr(ii)
      
      
      if(kcarr(i).ge.0)goto 6
      kcarr(i)=kcarr(i)+10000
      kcarr(i-1)=kcarr(i-1)-1
6     ii =ii-1
      if(ii.lt.3)goto 17
      end do
7     do i =3,kcarr(2)+2
      if(kcarr(i).ne.0)goto 8
      end do
8     kcarr(2) =kcarr(2) +3-i
      do j=3,kcarr(2) +2
      kcarr(j)=kcarr(j+i-3)
      end do
      goto 100
5     do i=2,kbarr(2)+2
      kcarr(i)=kbarr(i)
      end do
      kcarr(1)=mod(kbarr(1)+isora,2)
      ii =karr(2)+2
      do i =kbarr(2)+2,3,-1
      kcarr(i)=kcarr(i)-karr(ii)
      
      
      if(kcarr(i).ge.0)goto 16
      kcarr(i)=kcarr(i)+10000
      kcarr(i-1)=kcarr(i-1)-1
16    ii=ii-1
      if(ii.lt.3)goto 17
      end do
17    jj=i-1 
      if(jj.lt.3)goto 7
      do j=jj,3,-1
      if (kcarr(j).ge.0)goto 7
      kcarr(j)=kcarr(j)+10000
      kcarr(j-1)=kcarr(j-1)-1
18    end do
26    kcarr(1)=0
      kcarr(2)=0
      goto 100
30    if(kmx.eq.1)goto 40
      do i=3,karr(2)+2
      kcarr(i+2)=karr(i)
      
      end do
      kcarr(1)=karr(1)
      kcarr(2)=karr(2)+2
      kcarr(3)=0
      kcarr(4)=0
      ii =kbarr(2)+2
      do i=kcarr(2)+2,3,-1
      kcarr(i)=kcarr(i)+kbarr(ii)
      if(kcarr(i).lt.10000)goto 32
      kcarr(i)=kcarr(i)-10000
      kcarr(i-1)=kcarr(i-1)+1
32    ii =ii-1
      if(ii.lt.3)goto 34
      end do
34    jj=i-1
      do j=jj,3,-1
      if(kcarr(j).lt.10000)goto 36
      kcarr(j)=kcarr(j)-10000
      kcarr(j-1)=kcarr(j-1)+1
      end do
36    do i=3,kcarr(2)+2
      if(kcarr(i).ne.0)goto 38
      end do
38    kcarr(2)=kcarr(2)+3-i
      
      if(i.eq.3)goto 100
      do j =3,kcarr(2)+2
      kcarr(j)=kcarr(j+i-3)
      end do
      goto 100
40    do i=1,kbarr(2)+2
      kcarr(i+2)=kbarr(i)
      end do
      kcarr(1)=karr(1)
      kcarr(2)=kbarr(2)+2

      kcarr(3)=0
      kcarr(4)=0
      ii=karr(2)+2
      do i=kcarr(2)+2,1,-1
      kcarr(i)=kcarr(i) +karr(ii)
      if(kcarr(i).lt.10000)goto 42
      kcarr(i)=kcarr(i)-10000
      kcarr(i-1)=kcarr(i-1)+1
42    ii=ii-1
      if(ii.lt.3)goto 34
      end do
      goto 34

50    do i=2,kbarr(2)+2
      kcarr(i)=kbarr(i)
      end do
      kcarr(1)=mod(kbarr(1) +isora,2)
      if (kcarr(2).ne.0)goto 100
      kcarr(1)=0
      goto 100
55    do i=1,karr(2)+2
      kcarr(i)=karr(i)
      end do
100   return
      end


      

      subroutine subgcd(ibig,little,igcd2)
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000)
      common irrr(2000)
      common mnum(50)
      
      
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)






918   itemp=int(ibig/little)
      irem1=ibig -itemp *little
      if (irem1.eq.0)goto 940
      ibig=little
      little =irem1
      goto 918
940   igcd2=little
      return
      end


      subroutine mpgcd
      common karr(2000),kbarr(2000),kcarr(2000),ipqt(2000)
      common irrr(2000)
      common mnum(50)
      
      
      common kara(2000),karb(2000),kard(2000),karp(2000),karv(2000)
      
      common iarq(2)
      
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      
      dimension karu(2000),karv1(2000),karv3(2000),karqq(2000)
      dimension kart3(2000),kart1(2000)
      
      
      
      karu(1)=0
      karu(2)=1
      karu(3)=1
      do i=1,kara(2)+2
      kard(i)=kara(i)
      end do
      
      if (karb(2).eq.0)goto 888
3     karv1(2)=0
      karv1(1)=0
      do i=1,karb(2) +2
      karv3(i)=karb(i)
      end do
815   if (karv3(2).eq.0)goto 830
      
      
      
      do i=3,kard(2)+2
      karr(i-2)=kard(i)
      end do
      ilen=kard(2)
      do i=3,karv3(2)+2
      kbarr(i-2)=karv3(i)
      end do
      ilen2=karv3(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (icont.eq.0)goto 4
      do i=1,icont
      karqq(i+2)=ipqt(i)
      end do
4     karqq(2)=icont
      karqq(1)=karv3(1)
      
      if (irlen.eq.0)goto 41
      do i=1,irlen
      kart3(i+2)=irrr(i)
      end do
      kart3(2)=irlen
      kart3(1)=0
      
      
      
      goto 51
41    kart3(2)=0
      kart3(1)=0
51    if (karv1(2).eq.0)goto 6
      if (karqq(2).eq.0)goto 6
      do i=3,karv1(2) +2
      kbarr(i-2)=karv1(i)
      end do
      do i=3,karqq(2)+2
      karr(i-2)=karqq(i)
      end do
      ilen=karqq(2)
      ilen2=karv1(2)
      call mpmul(ilen,ilen2,ilen3)
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(2)=ilen3
      kbarr(1)=mod(karqq(1) +karv1(1),2)
      do i=1,karu(2)+2
      karr(i)=karu(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2)+2
      kart1(i)=kcarr(i)
      end do
      
      
      goto 7
6     do i=1,karu(2)+2
      kart1(i)=karu(i)
      end do
      karu(1)=0
      karu(2)=0
      goto 8
7     do i=1,karv1(2)+2
      karu(i)=karv1(i)
      end do
8     do i=1,karv3(2)+2
      kard(i)=karv3(i)
      end do
      do i=1,kart1(2)+2
      karv1(i)=kart1(i)
      end do
      do i=1,kart3(2)+2
      karv3(i)=kart3(i)
      end do
      goto 815
      
830   if ((kara(2).eq.0).or.(karu(2).eq.0))goto 9   
      do i=3,kara(2)+2
      karr(i-2)=kara(i)
      end do
      ilen=kara(2)
      do i=3,karu(2)+2
      kbarr(i-2)=karu(i)
      end do
      ilen2=karu(2)
      call mpmul(ilen,ilen2,ilen3)
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(2)=ilen3
      do i=1,kard(2)+2
      karr(i)=kard(i)
      end do
      kbarr(1)=mod(kara(1)+karu(1),2)
      call mpadd(1)
      do i=3,kcarr(2)+2
      karr(i-2)=kcarr(i)
      end do
      ksgn=kcarr(1)
      ilen=kcarr(2)
      goto 10
9     do i=3,kard(2)+2
      karr(i-2)=kard(i)
      end do
      ilen =kard(2)
      ksgn =kard(1)

10    do i=3,karb(2) +2    
      kbarr(i-2)=karb(i)
      end do
      ilen2=karb(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      karv(2)=icont
      if (icont.eq.0)goto 11
      do i=1,icont
      karv(i+2)=ipqt(i)
      end do
      karv(1)=mod(ksgn +karb(1),2)
      


      if ((karu(2).eq.0).or.(karu(1).eq.1))goto 890
      karr(1)=0
      karr(2)=1
      karr(3)=1
      do i=1,kara(2) +2
      kbarr(i)=kara(i)
      end do
      call mpadd(1)
      ksgn=kcarr(1)
      ilen2=kcarr(2)
      do i=3,kcarr(2)+2
      kbarr(i-2)=kcarr(i)
      end do
      
      do i=3,karv(2)+2
      karr(i-2)=karv(i)
      end do
      ilen=karv(2)
      call mpmul(ilen,ilen2,ilen3)
      karv(1)=mod(ksgn +karv(1),2)
      karv(2)=ilen3
      do i=1,ilen3
      karr(i)=kcarr(i)
      karv(i+2)=kcarr(i)
      end do
      do i=3,kara(2)+2
      kbarr(i-2)=kara(i)
      end do
      ilen2=kara(2)
      ilen=ilen3
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (karv(1).eq.0)goto 12
      do i=1,irlen
      kbarr(i+2)=irrr(i)
      end do
      kbarr(1)=0
      kbarr(2)=irlen
      do i=1,kara(2)+2
      karr(i)=kara(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2)+2
      karv(i)=kcarr(i)
      end do
      goto 890
12    do i=1,irlen
      karv(i+2)=irrr(i)
      end do
      karv(2)=irlen
      goto 890
11    karv(1)=0
      karv(2)=0
      goto 890
888   do i=1,karb(2)+2
      karv(i)=karb(i)
      end do
890   karv(1)=0
      return
      end

      


      
      
      
