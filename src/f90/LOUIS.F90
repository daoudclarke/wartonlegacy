      program bopbexp3
!     finds kernel of integer matrix reducing basis as it proceeds
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      dimension jfak(30),jfreq(30),isumt(50),itt(50),iuu(100)
      open (unit=3,file='exppar',access='direct',form=&
      'formatted',recl=256,status='old')
      read (3,1003,rec=1)nn,mm,(jfak(jf),jf=1,30),(jfreq(jf),jf=1,30)
1003  format (i6,i6,i4,30i4,30i4)      
      mmr=9*mm+91
      nnr=100*nn
      kkb=0
      ihitn=0
      do jf=1,20
      littr(jf)=0
      end do
      
      open(unit=1,file='bexp1',access='direct',form=& 
      'formatted',recl=mmr,status='old')
      open(unit=2,file='hexp1',access='direct',form=&
      'formatted',recl=nnr,status='old')
!      close(unit=2,status='delete') 
      open (unit=4,file='lexp1',access='direct',form=&
      'formatted',recl=280000,status='new')
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
1001  format(i1,i10,20i4,1400i9)
1002  format(1400(25i4))      
      end do
      k=2
      kmax=1
      idd(1,1)=0
      idd(1,2)=1
      idd(1,3)=1
      itt(1)=0
      itt(2)=0
      
      
      read (1,1001,rec=1)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      print *,'ok1','ibbl15',ibbl(5)
      do i=1,mm
      marr(1)=0
      if (ibbl(i).ge.0)goto 804
      marr(1)=1
804   a=a
      if (abs(ibbl(i)).ge.10000)goto 900
      marr(2)=1
      marr(3)=abs(ibbl(i))
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
      mbarr(3)=abs(ibbr(i))
      goto 906
904   mbarr(2)=2
      mbarr(3)=int(abs(ibbr(i))/10000)
      mbarr(4)=abs(ibbr(i))-mbarr(3)*10000

906   call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,itt(2)+2
      kbarr(jf)=itt(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      itt(jf)=kcarr(jf)
      end do
      end do


!      itt=itt+ibbl(i)*ibbl(i)
!      end do
      if (itt(2).eq.0)goto 86
      do jf=1,itt(2)+2
      idd(2,jf)=itt(jf)
      end do
!      idd(2)=itt
      iff(1)=1
      goto 90
86    idd(2,1)=0
      idd(2,2)=1
      idd(2,3)=1
      iff(1)=0
! 90    lamdar(2,1)=0
90    lamdar(1,1)=0      
      lamdar(1,2)=0
      
      
      isumt(1)=0
      isumt(2)=0
      read (1,1001,rec=2)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      print *,'ok2','ibbr25',ibbr(5)
      
      do i=1,mm
!      lamdar(2,1)=lamdar(2,1)+ibbl(i)*ibbr(i)
!      isumt=isumt+ibbr(i)*ibbr(i)
      marr(1)=0
      if (ibbl(i).ge.0)goto 800
      marr(1)=1
800   if (abs(ibbl(i)).ge.10000)goto 908      
      marr(2)=1
      marr(3)=abs(ibbl(i))
      mbarr(1)=0
      goto 910
908   marr(2)=2
      marr(3)=int(abs(ibbl(i))/10000)
      marr(4)=abs(ibbl(i))-marr(3)*10000

910   mbarr(1)=0      
      if (ibbr(i).ge.0)goto 802
      mbarr(1)=1
802   if (abs(ibbr(i)).ge.10000)goto 912



      mbarr(2)=1
      mbarr(3)=abs(ibbr(i))
      goto 914
912   mbarr(2)=2
      mbarr(3)=int(abs(ibbr(i))/10000)
      mbarr(4)=abs(ibbr(i))-mbarr(3)*10000




914   call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,isumt(2)+2
      kbarr(jf)=isumt(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      isumt(jf)=kcarr(jf)
      end do
      do jf=1,lamdar(1,2)+2
      kbarr(jf)=lamdar(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      lamdar(1,jf)=kcarr(jf)
      end do







      end do
      if (iff(1).eq.0)goto 700
      do jf=1,isumt(2)+2
      marr(jf)=isumt(jf)
      end do
      do jf=1,idd(2,2)+2
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






!      idd(3)=isumt*idd(2)-lamdar(2,1)*lamdar(2,1)
100   if (idd(3,2).eq.0)goto 106
      iff(2)=1
      goto 108
700   do jf=1,isumt(2)+2
      idd(3,jf)=isumt(jf)
      end do
      
      goto 100

106   iff(2)=0
!      idd(3)=idd(2)
      do jf=1,idd(2,2)+2
      idd(3,jf)=idd(2,jf)
      end do
108   if (iff(1).eq.0)goto 119
      lpar=1
      call subred(k,lpar,mm,nn)
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=lamdar(ii,ind),50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      



1004  format(1400(50i4))      
      
      
      
      
      if (iff(1).eq.0)goto 119
      if (iff(2).eq.0)goto 118
      goto 119
118   call subswap(k,mm,nn,kmax)
119   kmax=2
      k=3
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)


123   if (k.le.kmax)goto 300
      kmax=k
      lamdar(k,1)=0
      read (1,1001,rec=1)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      read (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
      print *,'k',k,'ibbr5',ibbr(5),'ibbl5',ibbl(5)
      do i=1,mm
      lamdar(k,1)=lamdar(k,1)+ibbr(i)*ibbl(i)
      end do
      print *,'lamdark1',lamdar(k,1)
      do j=2,k 
      if ((iff(j).eq.0).and.(j.lt.k))goto 152
      iuu(1)=0
      iuu(2)=0

      read (1,1001,rec=j)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      print *,'ok3','j',j
      do i=1,mm
!      iuu=iuu+ibbl(i)*ibbr(i)
      marr(1)=0
      if (ibbl(i).ge.0)goto 808
      marr(1)=1
808   if (abs(ibbl(i)).ge.10000)goto 920


      marr(2)=1
      marr(3)=abs(ibbl(i))
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
      mbarr(3)=abs(ibbr(i))
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
      read (4,1004,rec=j)((lamdar2(ir,jf),jf=1,50),ir=1400)
      
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
      mbarr(jf)=lamdar(i,jf)
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
      goto 166
811   print *,'iuu in error'
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
      call subred(k,lpar,mm,nn)
      if (iff(k).eq.0)goto 310
      goto 320
310   call subswap(k,mm,nn,kmax)
      if (k-1.gt.2)goto 318
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=lamdar(ii,ind),50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      
      k=2
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      goto 300
318   a=a
      do ii=1,k-1
      ind=lamdar(ii,2)+3
      do jf=lamdar(ii,ind),50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      
     k=k-1
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
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
      do jf=lamdar(ii,ind),50
      lamdar(ii,jf)=0
      end do
      end do
      do ii=k,nn
      do jf=1,50
      lamdar(ii,jf)=0
      end do
      end do
      write (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      
      
      k=k+1
      read (4,1004,rec=k)((lamdar(ir,jf),jf=1,50),ir=1,1400)
      if (k.le.nn)goto 123
      read (2,1002,rec=1)(ihhl(jf),jf=1,nn)
      
      print *,'vector1',(ihhl(jf),jf=1,nn)
      close (unit=1)
      close(unit=2)
      close(unit=3)
      end

      subroutine subred(k,lpar,mm,nn)
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
      do ii=2,idd(ipar+1,2)+2
      if (itempar(ii).lt.idd(lpar+1,ii))goto 490
      if (itempar(ii).gt.idd(lpar+1,ii))goto 100
      end do
      
!      if (abs(2*lamdar(k,lpar)).le.idd(lpar+1))goto 490
100   a=a      
      
!      iqq=int((2*lamdar(k,lpar)+idd(lpar+1))/(2*idd(lpar+1)))
      read (2,1002,rec=lpar)(ihhl(jf),jf=1,nn)
      read (2,1002,rec=k)(ihhr(jf),jf=1,nn)
      read (1,1001,rec=lpar)kkb,ihitn,(littr(j2),j2=1,20),(ibbl(jf),jf=1,mm)
      read (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
1001  format(i1,i10,20i4,1400i9)
1002  format(1400(25i4))
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
11    a=a
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
      if (iqq(2).gt.2)goto 179
      if (iqq(2).eq.1)goto 110
      iqm=iqq(3)*10000+iqq(4)
      goto 112
110   iqm=iqq(3)
112   if (iqq(1).eq.1)goto 120
      goto 122
120   iqm=(-1)*iqm
122   goto 124
179   print *,'iqq too large length=',iqq(2)
      stop
124   a=a
!      do i=1,nn
!      ihhr(i)=ihhr(i)-iqm*ihhl(i)
!      end do
      do i=1,mm
      ibbr(i)=ibbr(i)-iqm*ibbl(i)
      end do
!      write (2,1002,rec=k)(ihhr(jf),jf=1,nn)
      write (1,1001,rec=k)kkb,ihitn,(littr(j2),j2=1,20),(ibbr(jf),jf=1,mm)
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
      ihhr(jf)=0
      end do
      write (2,1002,rec=k)((ihhr(ir,jf),jf=1,25),ir=1,nn)
      
      




!      lamdar(k,lpar)=lamdar(k,lpar)-iqq*idd(lpar+1)
      if (lpar.eq.1)goto 490
      read (4,1004,rec=lpar)((lamdar2(ir,jf),jf=1,50),ir=1,1400)
      do jf=1,iqq(2)+2
      marr(jf)=iqq(jf)
      end do
      
      do i=1,lpar-1
      do jf=1,lamdar2(i,2)+2
      mbarr(jf)=lamdar(i,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,lamdar(i,2)+2
      karr(jf)=lamdar(i,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      lamdar(i,jf)=kcarr(jf)
      end do
      





!lamdar(k,i)=lamdar(k,i)-iqq*lamdar(lpar,i)
      end do
490   return
      end


      subroutine subswap(k,mm,nn,kmax)
      common ihhl(1400,25),ihhr(1400,25),ibbl(1400),ibbr(1400)
      common lamdar(1400,50),iff(1400),idd(1400,50),littr(20)
      common lamdar2(1400,50)
      dimension itempar(100),lamda(50),itt(50)
1001  format (i1,i10,20i4,1400i9) 
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
      if (k.lt.3)goto 522
      read (4,1004,rec=k-1)((lamdar2(ir,jf),jf=1,50),ir=1,1400)
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
      if (k.lt.3)goto 522
      do ii=1,k-2
      ind=lamdar2(ii,2)+3
      do jf=ind,50
      lamdar2(ii,jf)=0
      end do
      end do
      write (4,1004,k-1)((lamdar2(ir,jf),jf=1,50),ir=1,1400)





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
      read (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,1400)
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
      write (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,1400)
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
       mbarr(jf)=idd(jf)
       end do
       call mendiv
       if (mcarr(2),ne.0)goto 813
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
      read (4,1004,rec=i)((lamdar2(ir,jf),jf=1,50),ir=1,1400)
      do jf=1,lamdar(j,2)+2
      marr(jf)=lamdar(j,jf)
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
      lamdar(j,jf)=mdarr(jf)
      end do
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
560   idd(k)=idd(k-1)
      itemp=iff(k)
      iff(k)=iff(k-1)
      iff(k-1)=itemp
!      lamdar(k,k-1)=0
      lamdar(k-1,1)=0
      lamdar(k-1,2)=0
      if (k+1.gt.kmax)goto 590
      do i=k+1,kmax
      lamdar(i,k)=lamdar(i,k-1)
      lamdar(i,k-1)=0
      end do

590   return


      end



