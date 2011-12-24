      
      program bobmodj2
!     program for computing modular invariants ex "bobtrial"
      
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
      dimension n(50)
      dimension numb1(3000),numb2(50)
      dimension ipre(3000),icurr(3000)
      dimension ipisum(3,420),iptop(3,420)
      dimension ipi(420),itop1(520),itop2(520),ifack(420),iesum(420)
      dimension iiq(2,420),idel(2,2,1000),itab(400),iffy(2,420)
      dimension jay(2,420),ieesum(2,400),itemsy(420),ifsum(420)
      dimension igsum(420),itemst(2,420),ittop(2,420),inorm(2,1000)
      dimension iden(420),ntop(2,420),ieesum2(2,1000)
      msgn=0
      idis=715
      iia=1
      iib=-1
      ipi(1)=0
      ipi(2)=14
      ipi(3)=3
      ipi(4)=1415
      ipi(5)=9265
      ipi(6)=3589
      ipi(7)=7932
      ipi(8)=3846
      ipi(9)=2643
      ipi(10)=3832
      ipi(11)=7950
      ipi(12)=2884
      ipi(13)=1971
      ipi(14)=6939
      ipi(15)=9375
      ipi(16)=1100
      
      do jf=2,201
      karr(jf)=0
      end do
      karr(1)=1
      ilen=201
      kbarr(1)=18 
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iptop(1,jf+2)=ipqt(jf)
      end do
      iptop(1,1)=0
      iptop(1,2)=icont
      kbarr(1)=57 
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iptop(2,jf+2)=ipqt(jf)
      end do
      iptop(2,1)=0
      iptop(2,2)=icont
      kbarr(1)=239
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iptop(3,jf+2)=ipqt(jf)
      end do
      iptop(3,1)=0
      iptop(3,2)=icont
      do k=1,3
      do jf=1,iptop(1,2)+2
      ipisum(k,jf)=iptop(k,jf)
      end do
      end do
      do k=1,3
      do jf=3,iptop(k,2)+2
      karr(jf-2)=iptop(k,jf)
      kbarr(jf-2)=iptop(k,jf)
      end do
      ilen=iptop(k,2)
      ilen2=iptop(k,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      do jf=1,ilen3
      itop2(jf+2)=kcarr(jf)
      end do
      itop2(1)=0
      itop2(2)=ilen3
      do i=1,200
      do jf=3,iptop(k,2)+2
      karr(jf-2)=iptop(k,jf)
      end do
      ilen=iptop(k,2)
      do jf=3,itop2(2)+2
      kbarr(jf-2)=itop2(jf)
      end do
      ilen2=itop2(2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 500
      ilen3=ilen3-200
      do jf=1,ilen3
      iptop(k,jf+2)=kcarr(jf)
      end do
      iptop(k,1)=0
      iptop(k,2)=ilen3
      do jf=3,iptop(k,2)+2
      karr(jf-2)=iptop(k,jf)
      end do
      do jf=iptop(k,2)+1,iptop(k,2)+200
      karr(jf)=0
      end do


      ilen=iptop(k,2)+200
      kbarr(1)=2*i+1
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      icont=icont-200
      do jf=1,icont
      kbarr(jf+2)=ipqt(jf)
      end do
      kbarr(1)=mod(i,2)
      kbarr(2)=icont
      do jf=1,ipisum(k,2)+2
      karr(jf)=ipisum(k,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ipisum(k,jf)=kcarr(jf)
      end do
      end do
500   end do
      ipi(1)=0
      ipi(2)=0
      do jf=3,ipisum(1,2)+2
      karr(jf-2)=ipisum(1,jf)
      end do
      ilen=ipisum(1,2)
      kbarr(1)=48
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilen3
      do jf=1,ipi(2)+2
      karr(jf)=ipi(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ipi(jf)=kcarr(jf)
      end do
      do jf=3,ipisum(2,2)+2
      karr(jf-2)=ipisum(2,jf)
      end do
      ilen=ipisum(2,2)
      kbarr(1)=32
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilen3
      do jf=1,ipi(2)+2
      karr(jf)=ipi(jf)
      end do
      call mpadd(0) 
      do jf=1,kcarr(2)+2
      ipi(jf)=kcarr(jf)
      end do
      do jf=3,ipisum(3,2)+2
      karr(jf-2)=ipisum(3,jf)
      end do
      ilen=ipisum(3,2)
      kbarr(1)=20
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilen3
      do jf=1,ipi(2)+2
      karr(jf)=ipi(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      ipi(jf)=kcarr(jf)
      end do
      print *,'pi=',ipi(3),'.',(ipi(jf),jf=4,ipi(2)+2)
      






















!     compute rational integer square root
11112 numb1(1)=0
      numb1(2)=401
      numb1(3)=idis
      do jf=4,403
      numb1(jf)=0
      end do
      
      do jf=1,numb1(2)+2
      ipre(jf)=numb1(jf)
      end do
720   do jf=3,numb1(2)+2
      karr(jf-2)=numb1(jf)
      end do
      ilen=numb1(2)
      do jf=3,ipre(2)+2
      kbarr(jf-2)=ipre(jf)
      end do
      ilen2=ipre(2)
      
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      do jf=1,icont
      kbarr(jf+2)=ipqt(jf)
      end do
      kbarr(2)=icont
      kbarr(1)=0
      do jf=1,ipre(2)+2
      karr(jf)=ipre(jf)
      end do
      call mpadd(0)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      ilen2=1
      kbarr(1)=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      icurr(jf+2)=ipqt(jf)
      end do
      icurr(2)=icont
      icurr(1)=0
      
      
      if (icurr(2).lt.ipre(2))goto 730
      do jf=3,ipre(2)+2
      if (icurr(jf).lt.ipre(jf))goto 730
      end do
      goto 732
730   do jf=1,icurr(2)+2
      ipre(jf)=icurr(jf)
      end do
      goto 720
732   print *,'int. root',(icurr(jf),jf=1,icurr(2)+2)
      
      do jf=3,icurr(2)+2
      karr(jf-2)=icurr(jf)
      kbarr(jf-2)=icurr(jf)
      end do
      ilen=icurr(2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      print *,'chsq',(kcarr(jf),jf=1,ilen3)
      if (ilen3.ne.numb1(2))goto 740
      do jf=1,ilen3
      if (kcarr(jf).ne.numb1(jf+2))goto 740
      end do
      goto 742
740   print *,'not perfect square'
      print *,'len',icurr(2)
      
      
      
      ilen=ipi(2)
      ilen2=icurr(2)
      do jf=3,ipi(2)+2
      karr(jf-2)=ipi(jf)
      end do
      do jf=3,icurr(2)+2
      kbarr(jf-2)=icurr(jf)
      end do
      call mpmul(ilen,ilen2,ilen3)
      
      ilen3=ilen3-200
      
      do jf=1,ilen3
      itop1(jf+2)=kcarr(jf)
      itemsy(jf+2)=kcarr(jf)
      end do
      itemsy(1)=0
      itemsy(2)=ilen3
      itop1(1)=0
      itop1(2)=ilen3
!     provision for iia.ne.1      
      do jf=3,itop1(2)+2
      karr(jf-2)=itop1(jf)
      end do
      ilen=itop1(2)
      kbarr(1)=abs(iia)
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      itop1(jf+2)=ipqt(jf)
      itemsy(jf+2)=ipqt(jf)
      end do
      itop1(2)=icont
      itemsy(2)=icont


      
      
      itop2(1)=0
      itop2(2)=201
      itop2(3)=1
      do jf=4,203
      itop2(jf)=0
      iesum(jf)=0
      end do
      
      ifack(1)=0
      ifack(2)=1
      ifack(3)=1
      iesum(1)=0
      iesum(2)=201
      iesum(3)=1
      do i=1,200
      do jf=1,itop1(2)+2
      marr(jf)=itop1(jf)
      end do
      do jf=1,itop2(2)+2
      mbarr(jf)=itop2(jf)
      end do
      
      call menmul
      
      mcarr(2)=mcarr(2)-200
      if (mcarr(2).lt.0)goto 1011
      goto 1012
1011  mcarr(2)=0      
1012  do jf=1,mcarr(2)+2
      itop2(jf)=mcarr(jf)
      end do
      
      
      
!      do jf=3,itop1(2)+2
!      karr(jf-2)=itop1(jf)
!      end do
!      do  jf=3,itop2(2)+2
!      kbarr(jf-2)=itop2(jf)
!      end do
!      ilen=itop1(2)
!      ilen2=itop2(2)
!      call mpmul(ilen,ilen2,ilen3)
!      ilen3=ilen3-200
!      inabc=1
!      if (ilen3.le.0)goto 1000
!      do jf=1,ilen3
!      itop2(jf+2)=kcarr(jf)
!      end do
!      itop2(1)=0
!      itop2(2)=ilen3
      ilen=1
      karr(1)=i
      ilen2=ifack(2)
      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      ifack(jf+2)=kcarr(jf)
      end do
      ifack(1)=0
      ifack(2)=ilen3
      do jf=3,itop2(2)+2
      karr(jf-2)=itop2(jf)
      end do
      do jf=itop2(2)+1,itop2(2)+200
      karr(jf)=0
      end do
      ilen=itop2(2)+200
      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      ilen2=ifack(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      inabc=100
!      if (icont.le.200)goto 1000
      if (icont.le.200)goto 621
      goto 622
621   kbarr(1)=0
      kbarr(2)=0
      goto 623
622   do jf=1,icont-200
      kbarr(jf+2)=ipqt(jf)
      end do
      kbarr(1)=0
623   kbarr(2)=icont-200
      do jf=1,iesum(2)+2
      karr(jf)=iesum(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iesum(jf)=kcarr(jf)
      end do
      end do
      print *,'tau',(iesum(jf),jf=1,iesum(2)+2)
      
      
      if (iia.lt.0)goto 100
!      iib=abs(iib)
!      iib=mod(iib,2)
      isgn=0
      goto 102
100   isgn=1      
102   if (iib.lt.0)goto 104      
      goto 106
104   isgn=isgn+1
106   msgn=mod(isgn,2)
      do jf=3,ipi(2)+2
      karr(jf-2)=ipi(jf)
      end do
      
      ilen=ipi(2)
      kbarr(1)=abs(iib)
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      kbarr(1)=abs(iia)
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      itemsy(jf+2)=ipqt(jf)
      end do
      itemsy(2)=icont
      itemsy(1)=msgn
      goto 107







      karr(1)=1500
      do jf=2,13
      karr(jf)=0
      end do
      ilen=13
      do jf=3,ipi(2)+2
      kbarr(jf-2)=ipi(jf)
      end do
      ilen2=14
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3-13
      itemsy(jf+2)=kcarr(jf)
      end do
      itemsy(1)=0
      itemsy(2)=ilen3-13
107   do jf=1,itemsy(2)+2
      marr(jf)=itemsy(jf)
      mbarr(jf)=itemsy(jf)
      end do
      call menmul
      mcarr(2)=mcarr(2)-200
      if (mcarr(2).lt.0)goto 1021
      goto 1022
1021  mcarr(2)=0
1022  do jf=1,mcarr(2)+2
      itop1(jf)=mcarr(jf)
      end do
! 107   do jf=3,itemsy(2)+2
!      kbarr(jf-2)=itemsy(jf)
!      karr(jf-2)=itemsy(jf)
!      end do
!      ilen=itemsy(2)
!      ilen2=itemsy(2)
!      call mpmul(ilen,ilen2,ilen3)
!      ilen3=ilen3-200
!      inabc=2
!      if (ilen3.le.0)goto 1000
!      do jf=1,ilen3
!      itop1(jf+2)=kcarr(jf)
!      end do
!      itop1(1)=0
!      itop1(2)=ilen3



      do jf=1,itemsy(2)+2
      itop2(jf)=itemsy(jf)
      ifsum(jf)=itemsy(jf)
      end do
      ifsum(1)=msgn
      ifack(1)=0
      ifack(2)=1
      ifack(3)=1
      do i=1,200
      do jf=3,itop1(2)+2
      karr(jf-2)=itop1(jf)
      end do
      ilen=itop1(2)
      do jf=3,itop2(2)+2
      kbarr(jf-2)=itop2(jf)
      end do
      ilen2=itop2(2)
      call mpmul(ilen,ilen2,ilen3)
      
      if (ilen3.le.200)goto 600
      ilen3=ilen3-200
      do jf=1,ilen3
      itop2(jf+2)=kcarr(jf)
      end do
      itop2(2)=ilen3
      itop2(1)=0
      karr(1)=2*i
      ilen=1
      kbarr(1)=karr(1)+1
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      ilen2=ifack(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      ifack(jf+2)=kcarr(jf)
      end do
      ifack(1)=0
      ifack(2)=ilen3
      print *,'ifack2',ifack(2),'i',i,'itop22',itop2(2)
      do jf=3,itop2(2)+2
      karr(jf-2)=itop2(jf)
      end do
      do jf=itop2(2)+1,itop2(2)+200
      karr(jf)=0
      end do
      
      ilen=itop2(2)+200
      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      ilen2=ifack(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      inabc=101
!      if (icont.le.200)goto 1000
      if (icont.le.200)goto 611
      goto 612
611   kbarr(1)=0
      kbarr(2)=0
      goto 613
612   do jf=1,icont-200
      kbarr(jf+2)=ipqt(jf)
      end do
      kbarr(1)=mod(i+msgn,2)
      kbarr(2)=icont-200
613   do jf=1,ifsum(2)+2
      karr(jf)=ifsum(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ifsum(jf)=kcarr(jf)
      end do
      end do
      
600   print *,'sine',(ifsum(jf),jf=1,ifsum(2)+2)
      
      do jf=1,itemsy(2)+2
      marr(jf)=itemsy(jf)
      mbarr(jf)=itemsy(jf)
      end do
      call menmul
      mcarr(2)=mcarr(2)-200
      if (mcarr(2).lt.0)goto 1031
      goto 1032
1031  mcarr(2)=0
1032  do jf=1,mcarr(2)+2
      itop1(jf)=mcarr(jf)
      end do


!      do jf=3,itemsy(2)+2
!      karr(jf-2)=itemsy(jf)
!      kbarr(jf-2)=itemsy(jf)
!      end do
!      ilen=itemsy(2)
!      ilen2=itemsy(2)
!      call mpmul(ilen,ilen2,ilen3)
!      ilen3=ilen3-200
!      inabc=3
!      if (ilen3.le.0)goto 1000
!      do jf=1,ilen3
!      itop1(jf+2)=kcarr(jf)
      
!      end do
      itop1(1)=0
      itop2(1)=0
      itop1(2)=mcarr(2)
      itop2(2)=201
      itop2(3)=1
      igsum(1)=0
      igsum(2)=201
      igsum(3)=1
      do jf=4,203
      igsum(jf)=0
      itop2(jf)=0
      end do
      ifack(1)=0
      ifack(2)=1
      ifack(3)=1
      do i=1,200
      do jf=3,itop1(2)+2
      karr(jf-2)=itop1(jf)
      end do
      ilen=itop1(2)
      do jf=3,itop2(2)+2
      kbarr(jf-2)=itop2(jf)
      end do
      ilen2=itop2(2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 602
      ilen3=ilen3-200
      do jf=1,ilen3
      itop2(jf+2)=kcarr(jf)
      end do
      itop2(1)=0
      itop2(2)=ilen3
      karr(1)=2*i
      ilen=1
      kbarr(1)=karr(1)-1
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      ilen2=ifack(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      ifack(jf+2)=kcarr(jf)
      end do
      ifack(1)=0
      ifack(2)=ilen3
      
      do jf=3,itop2(2)+2
      karr(jf-2)=itop2(jf)
      end do
      do jf=itop2(2)+1,itop2(2)+200
      karr(jf)=0
      end do
      ilen=itop2(2)+200

      do jf=3,ifack(2)+2
      kbarr(jf-2)=ifack(jf)
      end do
      ilen2=ifack(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      inabc=102
!      if (icont.le.200)goto 1000
      if (icont.le.200)goto 631
      goto 632
631   kbarr(1)=0      
      kbarr(2)=0
      goto 633
632   do jf=1,icont-200
      kbarr(jf+2)=ipqt(jf)
      end do
      kbarr(1)=mod(i,2)
      kbarr(2)=icont-200
633   do jf=1,igsum(2)+2
      karr(jf)=igsum(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      igsum(jf)=kcarr(jf)
      end do
      end do
602   print *,'cosine',(igsum(jf),jf=1,igsum(2)+2)
      

      
      
      karr(1)=1
      do jf=2,401
      karr(jf)=0
      end do
      ilen=401
      do jf=3,iesum(2)+2
      kbarr(jf-2)=iesum(jf)
      end do
      ilen2=iesum(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iiq(1,jf+2)=ipqt(jf)
      end do
      iiq(1,2)=icont
      iiq(1,1)=mod(msgn+1,2)
      do jf=3,iiq(1,2)+2
      karr(jf-2)=iiq(1,jf)
      end do
      ilen=iiq(1,2)
      do jf=3,igsum(2)+2
      kbarr(jf-2)=igsum(jf)
      end do
      ilen2=igsum(2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=4
      if (ilen3.le.0)goto 1000
      

      do jf=1,ilen3
      iiq(1,jf+2)=kcarr(jf)
      end do
      iiq(1,2)=ilen3
      iiq(1,1)=igsum(1)
      do jf=3,ifsum(2)+2
      kbarr(jf-2)=ifsum(jf)
      end do
      ilen2=ifsum(2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=5
      if (ilen3.le.0)goto 1000
      do jf=1,ilen3
      iiq(2,jf+2)=kcarr(jf)
      end do
      iiq(2,2)=ilen3
      iiq(2,1)=ifsum(1)
      print *,'qr',(iiq(1,jf),jf=1,iiq(1,2)+2)
      
      print *,'qim',(iiq(2,jf),jf=1,iiq(2,2)+2)
      print *,'qrlen',iiq(1,2),'qimlen',iiq(2,2)
      


      
      do i=1,200
      itab(i)=0
      end do
      print *,'okay1'
      do i=1,10
      ipow=(i*(3*i-1))/2
      jpow=(i*(3*i+1))/2

      
      if (mod(i,2).eq.0)goto 800
      
      itab(ipow)=1
      itab(jpow)=1
      goto 806
800   itab(ipow)=2 
      itab(jpow)=2
806   end do

         
      print *,'okay2'
      do k=1,2
      do jf=1,iiq(k,2)+2
      itemst(k,jf)=iiq(k,jf)
      end do
      end do
      do ibc=1,2
      ittop(1,1)=0
      ittop(1,2)=201
      ittop(1,3)=1
      do jf=4,203
      ittop(1,jf)=0
      end do
      ittop(2,1)=0
      ittop(2,2)=0

      
      ieesum(1,1)=0
      ieesum(1,2)=201
      ieesum(1,3)=1
      do jf=4,203
      
      ieesum(1,jf)=0
      end do
! insert 2 instructions here      
      ieesum(2,1)=0
      ieesum(2,2)=0
      do i=1,200
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      if ((iiq(1,2).eq.0).or.(ittop(1,2).eq.0))goto 110
      do jf=3,iiq(1,2)+2
      karr(jf-2)=iiq(1,jf)
      end do
      ilen=iiq(1,2)
      do jf=3,ittop(1,2)+2
      kbarr(jf-2)=ittop(1,jf)
      end do
      ilen2=ittop(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 110
      ilen3=ilen3-200
      

      do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(iiq(1,1)+ittop(1,1),2)
110   a=a
      print *,'i',i,'ilen3',ilen3
      if ((iiq(2,2).eq.0).or.(ittop(2,2).eq.0))goto 112
      do jf=3,iiq(2,2)+2
      karr(jf-2)=iiq(2,jf)
      end do
      ilen=iiq(2,2)
      do jf=3,ittop(2,2)+2
      kbarr(jf-2)=ittop(2,jf)
      end do
      ilen2=ittop(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 112
      ilen3=ilen3-200
      

      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(iiq(2,1)+ittop(2,1)+1,2)
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
112   a=a
      print *,'i',i,'ilen3',ilen3
      if ((iiq(1,2).eq.0).or.(ittop(2,2).eq.0))goto 114
      do jf=3,iiq(1,2)+2
      karr(jf-2)=iiq(1,jf)
      end do
      ilen=iiq(1,2)
      do jf=3,ittop(2,2)+2
      kbarr(jf-2)=ittop(2,jf)
      end do
      ilen2=ittop(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 114
      ilen3=ilen3-200
      print *,'i',i,'ilen3',ilen3
      do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(iiq(1,1)+ittop(2,1),2)
114   a=a
      print *,'i',i,'ilen3',ilen3
      if ((iiq(2,2).eq.0).or.(ittop(1,2).eq.0))goto 116
      do jf=3,iiq(2,2)+2
      karr(jf-2)=iiq(2,jf)
      end do
      ilen=iiq(2,2)
      do jf=3,ittop(1,2)+2
      kbarr(jf-2)=ittop(1,jf)
      end do
      ilen2=ittop(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 116
      ilen3=ilen3-200
      
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(iiq(2,1)+ittop(1,1),2)
      do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
116   print *,'sofar ok i',i,'ilen3',ilen3
      do k=1,2
      do jf=1,inorm(k,2)+2
      ittop(k,jf)=inorm(k,jf)
      end do
      end do
      print *,'still ok','sumlens',ieesum(1,2),ieesum(2,2),ittop(1,2)
      print *,'itab',itab(i),'i',i
      
      if (itab(i).eq.0)goto 120
      inj=mod(itab(i),2)
      do k=1,2
      print *,'ieesumlen',ieesum(k,2),'k',k
      do jf=1,ieesum(k,2)+2
      karr(jf)=ieesum(k,jf)
      end do
      print *,'inormlen',inorm(k,2),'k',k,'inj',inj
      do jf=1,inorm(k,2)+2
      kbarr(jf)=inorm(k,jf)
      
      end do
      call mpadd(inj)
      do jf=1,kcarr(2)+2
      ieesum(k,jf)=kcarr(jf)
      end do
      end do
120   end do
1201  print *,'ieesum real',(ieesum(1,jf),jf=1,ieesum(1,2)+2)
      
      
      
      do k=1,2
      do jf=1,ieesum(k,2)+2
      ieesum2(k,jf)=ieesum(k,jf)
      end do
      end do
      ie=16
      n1=24
      n1=n1-ie
125   if (ie.eq.1)goto 150
      ie=ie/2
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      if (ieesum2(1,2).eq.0)goto 128 
      do jf=3,ieesum2(1,2)+2
      karr(jf-2)=ieesum2(1,jf)
      kbarr(jf-2)=ieesum2(1,jf)
      end do
      ilen=ieesum2(1,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 128
      ilen3=ilen3-200
      
      do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=0
128   a=a
      print *,'1seclen3',ilen3
      if (ieesum2(2,2).eq.0)goto 132
      do jf=3,ieesum2(2,2)+2
      karr(jf-2)=ieesum2(2,jf)
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen=ieesum2(2,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 132
      ilen3=ilen3-200
      
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=1
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
132   a=a
      print *,'2seclen3',ilen3
      if ((ieesum2(1,2).eq.0).or.(ieesum2(2,2).eq.0))goto 134
      do jf=3,ieesum2(1,2)+2
      karr(jf-2)=ieesum2(1,jf)
      end do
      ilen=ieesum2(1,2)
      do jf=3,ieesum2(2,2)+2
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen2=ieesum2(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 134
      ilen3=ilen3-200
      
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      kbarr(1)=2
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(ieesum2(1,1)+ieesum2(2,1),2)
134   a=a
      print *,'3seclen3',ilen3
      do k=1,2
      
      do jf=1,inorm(k,2)+2
      ieesum2(k,jf)=inorm(k,jf)
      end do
      end do
      if (n1.lt.ie)goto 125
      n1=n1-ie
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      if ((ieesum(1,2).eq.0).or.(ieesum2(1,2).eq.0))goto 136
      do jf=3,ieesum(1,2)+2
      karr(jf-2)=ieesum(1,jf)
      end do
      ilen=ieesum(1,2)
      do jf=3,ieesum2(1,2)+2
      kbarr(jf-2)=ieesum2(1,jf)
      end do
      ilen2=ieesum2(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 136
      ilen3=ilen3-200
      
      do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(ieesum(1,1)+ieesum2(1,1),2)
136   a=a
      print *,'4seclen3',ilen3
      if ((ieesum(2,2).eq.0).or.(ieesum2(2,2).eq.0))goto 138
      do jf=3,ieesum(2,2)+2
      karr(jf-2)=ieesum(2,jf)
      end do
      ilen=ieesum(2,2)
      do jf=3,ieesum2(2,2)+2
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen2=ieesum2(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 138
      ilen3=ilen3-200
      
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(ieesum(2,1)+ieesum2(2,1)+1,2)
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
138   a=a
      print *,'5seclen3',ilen3
      if ((ieesum(1,2).eq.0).or.(ieesum2(2,2).eq.0))goto 140
      do jf=3,ieesum(1,2)+2
      karr(jf-2)=ieesum(1,jf)
      end do
      ilen=ieesum(1,2)
      do jf=3,ieesum2(2,2)+2
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen2=ieesum(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 140
      ilen3=ilen3-200
      
      do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(ieesum(1,1)+ieesum2(2,1),2)
140   a=a
      print *,'6seclen3',ilen3
      if ((ieesum(2,2).eq.0).or.(ieesum2(1,2).eq.0))goto 142
      do jf=3,ieesum(2,2)+2
      karr(jf-2)=ieesum(2,jf)
      end do
      ilen=ieesum(2,2)
      do jf=3,ieesum2(1,2)+2
      kbarr(jf-2)=ieesum2(1,jf)
      end do
      ilen2=ieesum2(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 142
      ilen3=ilen3-200
      
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(ieesum(2,1)+ieesum2(1,1),2)
      do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
142   a=a
      print *,'7seclen3',ilen3
      
      do k=1,2
      
      do jf=1,inorm(k,2)+2
      ieesum2(k,jf)=inorm(k,jf)
      end do
      end do
      goto 125

150   print *,'powering ok'
      print *,'ieesum powered',(ieesum2(1,jf),jf=1,ieesum2(1,2)+2)
      print *,'im ieesum pow',(ieesum2(2,jf),jf=1,ieesum2(2,2)+2)
      
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0 
      inorm(2,2)=0
      if ((iiq(1,2).eq.0).or.(ieesum2(1,2).eq.0))goto 152
      do jf=3,iiq(1,2)+2
      karr(jf-2)=iiq(1,jf)
      end do
      ilen=iiq(1,2)
      do jf=3,ieesum2(1,2)+2
      kbarr(jf-2)=ieesum2(1,jf)
      end do
      ilen2=ieesum2(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 152
      ilen3=ilen3-200
      
      do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(iiq(1,1)+ieesum2(1,1),2)
152   if ((iiq(2,2).eq.0).or.(ieesum2(2,2).eq.0))goto 154 
      do jf=3,iiq(2,2)+2
      karr(jf-2)=iiq(2,jf)
      end do
      ilen=iiq(2,2)
      do jf=3,ieesum2(2,2)+2
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen2=ieesum2(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 154
      ilen3=ilen3-200
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(iiq(2,1)+ieesum2(2,1)+1,2)
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
154   if ((iiq(1,2).eq.0).or.(ieesum2(2,2).eq.0))goto 156
      do jf=3,iiq(1,2)+2
      karr(jf-2)=iiq(1,jf)
      end do
      ilen=iiq(1,2)
      do jf=3,ieesum2(2,2)+2
      kbarr(jf-2)=ieesum2(2,jf)
      end do
      ilen2=ieesum2(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 156
      ilen3=ilen3-200
      do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(iiq(1,1)+ieesum2(2,1),2)
156   if ((iiq(2,2).eq.0).or.(ieesum2(1,2).eq.0))goto 158
      do jf=3,iiq(2,2)+2
      karr(jf-2)=iiq(2,jf)
      end do
      ilen=iiq(2,2)
      do jf=3,ieesum2(1,2)+2
      kbarr(jf-2)=ieesum2(1,jf)
      end do
      ilen2=ieesum2(1,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 158
      ilen3=ilen3-200
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(iiq(2,1)+ieesum2(1,1),2)
      do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
158   do k=1,2
      do jf=1,inorm(k,2)+2
      idel(ibc,k,jf)=inorm(k,jf)
      end do
      end do
      print *,'q times sum ok'
!     putting iq in itemst squared back into iq
      if (ibc.eq.2)goto 166
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      if (itemst(1,2).eq.0)goto 162
      do jf=3,itemst(1,2)+2
      karr(jf-2)=itemst(1,jf)
      kbarr(jf-2)=itemst(1,jf)
      end do
      ilen=itemst(1,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      print *,'what','ilen3',ilen3
      if (ilen3.le.200)goto 162
      ilen3=ilen3-200
      do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=0
162   if (itemst(2,2).eq.0)goto 164
      print *,'ok q1'
      do jf=3,itemst(2,2)+2
      karr(jf-2)=itemst(2,jf)
      kbarr(jf-2)=itemst(2,jf)
      end do
      ilen=itemst(2,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 164
      ilen3=ilen3-200
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=1
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
164   if ((itemst(1,2).eq.0).or.(itemst(2,2).eq.0))goto 166
      print *,'ok q2'
      do jf=3,itemst(1,2)+2
      karr(jf-2)=itemst(1,jf)
      end do
      ilen=itemst(1,2)
      do jf=3,itemst(2,2)+2
      kbarr(jf-2)=itemst(2,jf)
      end do
      ilen2=itemst(2,2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.le.200)goto 166
      ilen3=ilen3-200
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      kbarr(1)=2
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(itemst(1,1)+itemst(2,1),2)
166   do k=1,2
      print *,'crash',k
      do jf=1,inorm(k,2)+2
      iiq(k,jf)=inorm(k,jf)
      end do
      end do
      print *,'next cycle'
      end do
      print *,'okaynew'
      print *,'del1',(idel(1,1,jf),jf=1,idel(1,1,2)+2)
      print *,'im lens',idel(1,2,2),idel(2,2,2)
      
      print *,'del2',(idel(2,1,jf),jf=1,idel(2,1,2)+2)
      
      inorm(1,1)=0
      inorm(1,2)=0
      if ((idel(1,1,2).eq.0).or.(idel(1,2,2).eq.0))goto 168
      do jf=3,idel(1,1,2)+2
      karr(jf-2)=idel(1,1,jf)
      kbarr(jf-2)=idel(1,1,jf)
      end do
      ilen=idel(1,1,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=24
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 2401
      goto 2402
2401  inorm(1,1)=0
      inorm(1,2)=0
      goto 2403
2402  do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=0
2403  do jf=1,idel(1,2,2)+2
      marr(jf)=idel(1,2,jf)
      mbarr(jf)=idel(1,2,jf)
      end do
      call menmul
      mcarr(2)=mcarr(2)-200
      if (mcarr(2).lt.0)goto 2301
      goto 2302
2301  mcarr(2)=0       
2302  do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do



!      do jf=3,idel(1,2,2)+2
!      karr(jf-2)=idel(1,2,jf)
!      kbarr(jf-2)=idel(1,2,jf)
!      end do
!      ilen=idel(1,2,2)
!      ilen2=ilen
!      call mpmul(ilen,ilen2,ilen3)
!      ilen3=ilen3-200
!      inabc=23
!      if (ilen3.le.0)goto 1000
!      do jf=1,ilen3
!      karr(jf+2)=kcarr(jf)
!      end do
!      karr(2)=ilen3
!      karr(1)=0
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iden(jf)=kcarr(jf)
      end do
      goto 170
168   print *,'f den zero or num zero'
      stop
170   inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      print *,'oknew2'
      if ((idel(1,1,2).eq.0).or.(idel(2,1,2).eq.0))goto 172
      do jf=3,idel(1,1,2)+2
      karr(jf-2)=idel(1,1,jf)
      end do
      ilen=idel(1,1,2)
      do jf=3,idel(2,1,2)+2
      kbarr(jf-2)=idel(2,1,jf)
      end do
      ilen2=idel(2,1,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=22
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 2201
      goto 2202
2201  inorm(1,1)=0
      inorm(1,2)=0
      goto 172
2202  do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(idel(1,1,1)+idel(2,1,1),2)
172   if ((idel(1,2,2).eq.0).or.(idel(2,2,2).eq.0))goto 174
      do jf=1,idel(1,2,2)+2
      marr(jf)=idel(1,2,jf)
      end do
      do jf=1,idel(2,2,2)+2
      mbarr(jf)=idel(2,2,jf)
      end do
      call menmul
      mcarr(2)=mcarr(2)-200
      if (mcarr(2).lt.0)goto 2101
      goto 2102
2101  mcarr(2)=0
2102  do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do







!      do jf=3,idel(1,2,2)+2
!      karr(jf-2)=idel(1,2,jf)
!      end do
!      ilen=idel(1,2,2)
!      do jf=3,idel(2,2,2)+2
!      kbarr(jf-2)=idel(2,2,jf)
!      end do
!      ilen2=idel(2,2,2)
!      call mpmul(ilen,ilen2,ilen3)
!      ilen3=ilen3-200
!      inabc=21
!      if (ilen3.le.0)goto 1000
!      do jf=1,ilen3
!      karr(jf+2)=kcarr(jf)
!      end do
!      karr(2)=ilen3
!      karr(1)=mod(idel(1,2,1)+idel(2,2,1),2)
      do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
174   if ((idel(1,1,2).eq.0).or.(idel(2,2,2).eq.0))goto 176 
      do jf=3,idel(1,1,2)+2
      karr(jf-2)=idel(1,1,jf)
      end do
      ilen=idel(1,1,2)
      do jf=3,idel(2,2,2)+2
      kbarr(jf-2)=idel(2,2,jf)
      end do
      ilen2=idel(2,2,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=20
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 2001
      goto 2002
2001  inorm(2,1)=0      
      inorm(2,2)=0
      goto 176
2002  do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(idel(1,1,1)+idel(2,2,1),2)
176   if ((idel(1,2,2).eq.0).or.(idel(2,1,2).eq.0))goto 178
      do jf=3,idel(1,2,2)+2
      karr(jf-2)=idel(1,2,jf)
      end do
      ilen=idel(1,2,2)
      do jf=3,idel(2,1,2)+2
      kbarr(jf-2)=idel(2,1,jf)
      end do
      ilen2=idel(2,1,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=19
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 1901
      goto 1902
1901  karr(1)=0
      karr(2)=0
      goto 1903
1902  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(idel(1,2,1)+idel(2,1,1)+1,2)
1903  do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
178   do k=1,2
      do jf=1,inorm(k,2)+2
      ntop(k,jf)=inorm(k,jf)
      end do
      end do
      if (ntop(1,2).eq.0)goto 180
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      do jf=ntop(1,2)+1,ntop(1,2)+200
      karr(jf)=0
      end do
      ilen=ntop(1,2)+200
      do jf=3,iden(2)+2
      kbarr(jf-2)=iden(jf)
      end do
      ilen2=iden(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iffy(1,jf+2)=ipqt(jf)
      end do
      iffy(1,2)=icont
      iffy(1,1)=mod(ntop(1,1)+iden(1),2)
180   do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      do jf=ntop(2,2)+1,ntop(2,2)+200
      karr(jf)=0
      end do
      ilen=ntop(2,2)+200
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      iffy(2,jf+2)=ipqt(jf)
      end do
      iffy(2,2)=icont
      iffy(2,1)=mod(ntop(2,1)+iden(1),2)
      print *,'oknew3'
      print *,'iffyreal',(iffy(1,jf),jf=1,iffy(1,2)+2)
      print *,'iffyim',(iffy(2,jf),jf=1,iffy(2,2)+2)
      

      do jf=3,iffy(1,2)+2
      karr(jf-2)=iffy(1,jf)
      end do
      ilen=iffy(1,2)
      kbarr(1)=256
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=iffy(1,1)
      kbarr(1)=0
      kbarr(2)=201
      kbarr(3)=1
      do jf=4,203
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ntop(1,jf)=kcarr(jf)
      end do
      do jf=1,iffy(2,2)+2
      karr(jf-2)=iffy(2,jf)
      end do
      ilen=iffy(2,2)
      kbarr(1)=256
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      ntop(2,jf+2)=kcarr(jf)
      end do
      ntop(2,2)=ilen3
      ntop(2,1)=iffy(2,1)
      do k=1,2
      do jf=1,ntop(k,2)+2
      ittop(k,jf)=ntop(k,jf)
      end do
      end do
      do icc=1,2
      inorm(1,1)=0
      inorm(1,2)=0
      inorm(2,1)=0
      inorm(2,2)=0
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      ilen=ntop(1,2)
      do jf=3,ittop(1,2)+2
      kbarr(jf-2)=ittop(1,jf)
      end do
      ilen2=ittop(1,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=18
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 1801
      goto 1802
1801  inorm(1,1)=0       
      inorm(1,2)=0
      goto 1803
1802  do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(ntop(1,1)+ittop(1,1),2)
1803  do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      ilen=ntop(2,2)
      do jf=3,ittop(2,2)+2
      kbarr(jf-2)=ittop(2,jf)
      end do
      ilen2=ittop(2,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=30
      if (ilen3.le.0)goto 3001
      goto 3002
3001  karr(1)=0
      karr(2)=0
      goto 3003
!      if (ilen3.le.0)goto 1000
3002  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(ntop(2,1)+ittop(2,1)+1,2)
3003  do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      ilen=ntop(1,2)
      do jf=3,ittop(2,2)+2
      kbarr(jf-2)=ittop(2,jf)
      end do
      ilen2=ittop(2,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=31
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3101
      goto 3102
3101  inorm(2,1)=0
      inorm(2,2)=0
      goto 3103
3102  do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(ntop(1,1)+ittop(1,1),2)
3103  do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      ilen=ntop(2,2)
      do jf=3,ittop(1,2)+2
      kbarr(jf-2)=ittop(1,jf)
      end do
      ilen2=ittop(1,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=32
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3201
      goto 3202
3201  karr(1)=0 
      karr(2)=0
      goto 3203
3202  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(ntop(2,1)+ittop(1,1),2)
3203  do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
      do k=1,2
      do jf=1,inorm(k,2)+2
      ntop(k,jf)=inorm(k,jf)
      end do
      end do
      end do
      print *,'oknew4'
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      ilen=ntop(1,2)
      do jf=3,iffy(1,2)+2
      kbarr(jf-2)=iffy(1,jf)
      end do
      ilen2=iffy(1,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=33
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3301
      goto 3302
3301  inorm(1,1)=0
      inorm(1,2)=0
      goto 3303
3302  do jf=1,ilen3
      inorm(1,jf+2)=kcarr(jf)
      end do
      inorm(1,2)=ilen3
      inorm(1,1)=mod(ntop(1,1)+iffy(1,1),2)
3303  print *,'oknew42'
      do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      ilen=ntop(2,2)
      do jf=3,iffy(2,2)+2
      kbarr(jf-2)=iffy(2,jf)
      end do
      ilen2=iffy(2,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=34
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3401
      goto 3402
3401  karr(1)=0
      karr(2)=0
      goto 3403
3402  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=mod(ntop(2,1)+iffy(2,1),2)
3403  do jf=1,inorm(1,2)+2
      kbarr(jf)=inorm(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(1,jf)=kcarr(jf)
      end do
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      ilen=ntop(1,2)
      do jf=3,iffy(2,2)+2
      kbarr(jf-2)=iffy(2,jf)
      end do
      ilen2=iffy(2,2)
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=35
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3501
      goto 3502
3501  inorm(2,1)=0      
      inorm(2,2)=0
      goto 3503
3502  do jf=1,ilen3
      inorm(2,jf+2)=kcarr(jf)
      end do
      inorm(2,2)=ilen3
      inorm(2,1)=mod(ntop(1,1)+iffy(2,1)+1,2)
3503  print *,'oknew44'
      do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      ilen=ntop(2,2)
      if (ilen.eq.0)goto 3601
      do jf=3,iffy(1,2)+2
      kbarr(jf-2)=iffy(1,jf)
      end do
      ilen2=iffy(1,2)
      print *,'lens',ilen,ilen2
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=36
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3601
      goto 3602
3601  karr(1)=0      
      karr(2)=0
      goto 3603
3602  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      
      karr(1)=mod(ntop(2,1)+iffy(1,1),2)
3603  print *,'ilen3=',ilen3
      do jf=1,inorm(2,2)+2
      kbarr(jf)=inorm(2,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      inorm(2,jf)=kcarr(jf)
      end do
      do k=1,2
      do jf=1,inorm(k,2)+2
      ntop(k,jf)=inorm(k,jf)
      end do
      end do
      print *,'oknew5'
      do jf=3,iffy(1,2)+2
      karr(jf-2)=iffy(1,jf)
      kbarr(jf-2)=iffy(1,jf)
      end do
      ilen=iffy(1,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=37
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3701
      goto 3702
3701  iden(1)=0
      iden(2)=0
      goto 3703
3702  do jf=1,ilen3
      iden(jf+2)=kcarr(jf)
      end do
      iden(2)=ilen3
      iden(1)=0
3703  do jf=3,iffy(2,2)+2
      karr(jf-2)=iffy(2,jf)
      kbarr(jf-2)=iffy(2,jf)
      end do
      ilen=iffy(2,2)
      ilen2=ilen
      call mpmul(ilen,ilen2,ilen3)
      ilen3=ilen3-200
      inabc=38
!      if (ilen3.le.0)goto 1000
      if (ilen3.le.0)goto 3801
      goto 3802
3801  karr(1)=0    
      karr(2)=0
      goto 3803
3802  do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(2)=ilen3
      karr(1)=0
3803  do jf=1,iden(2)+2
      kbarr(jf)=iden(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iden(jf)=kcarr(jf)
      end do
      print *,'ntop1',(ntop(1,jf),jf=1,ntop(1,2)+2)
      do jf=3,ntop(1,2)+2
      karr(jf-2)=ntop(1,jf)
      end do
      do jf=ntop(1,2)+1,ntop(1,2)+200
      karr(jf)=0
      end do
      ilen=ntop(1,2)+200
      do jf=3,iden(2)+2
      kbarr(jf-2)=iden(jf)
      end do
      ilen2=iden(2)
      print *,'ilenend',ilen,ilen2,'ntops',ntop(1,3),ntop(1,4)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      print *,'endicont',icont
      do jf=1,icont 
      jay(1,jf+2)=ipqt(jf)
      end do
      jay(1,2)=icont
      jay(1,1)=mod(ntop(1,1)+iden(1),2)
      do jf=3,ntop(2,2)+2
      karr(jf-2)=ntop(2,jf)
      end do
      do jf=ntop(2,2)+1,ntop(2,2)+200
      karr(jf)=0
      end do
      ilen=ntop(2,2)+200
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      jay(2,jf+2)=ipqt(jf)
      end do
      jay(2,2)=icont
      jay(2,1)=mod(ntop(2,1)+iden(1),2)
      print *,'jreal=',(jay(1,jf),jf=1,jay(1,2)+2)
      
      print *,'jimag=',(jay(2,jf),jf=1,jay(2,2)+2)
      stop
742   print *,'square ok'


      
1000  print *,'error:length negative, inabc=',inabc
      print *,'itop22',itop2(2),'ifack2',ifack(2)






      
      
      end
      


   
      
      subroutine menmul
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000),irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
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
      

      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      
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
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      dimension kdum(3000),isub(3000)
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
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
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


      subroutine mpprime(icorp,inlen)
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      dimension ie(50),itwo(100),nnum(50)
      print *,'inlen',inlen
      print *,'number',(mnum(i),i=1,inlen)
      icong=0
      icont1=0
      icont2=0
      ie(1)=0
      
      ileng =1
      itwo(1)=2
      karr(1)=itwo(1)
      ilen=1
      kbarr(1)=itwo(1)
      ilen2=1
4     call mpmul(ilen,ilen2,ilen3)
      
      if(ilen3.lt.inlen)goto 22
      if(ilen3.gt.inlen)goto 20
      do i=1,ilen3
      if(mnum(i).lt.kcarr(i))goto 20
      if(mnum(i).gt.kcarr(i))goto 22
      end do
      goto 22
20    ie(2) =ilen3
      do i=1,ilen3
      ie(i+2)=kcarr(i)
      end do
      goto 30
22    do i=1,ilen3
      kbarr(i)=kcarr(i)
      end do
      ilen2=ilen3
      goto 4

30    do i=1,inlen
      nnum(i+2)=mnum(i)
      end do
      nnum(1)=0
      nnum(2)=inlen
      
      
      do i=1,ie(2)
      karr(i)=ie(i+2)
      end do
      ilen =ie(2)
      
      
      ilen2=1
      kbarr(1)=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do i=1,icont
      ie(i+2)=ipqt(i)
      end do
      ie(2)=icont
      
      

      ie(1)=0
      do i=1,nnum(2)+2
      karr(i)=nnum(i)
      end do
      do i=1,ie(2)+2
      kbarr(i)=ie(i)
      end do
      
      
      call mpadd(1)
      do i=1,kcarr(2)+2
      nnum(i)=kcarr(i)
      end do
31    if(ie(2).gt.1)goto 32
      
      
      if(ie(3).eq.1)goto 100
32    do i=1,ie(2)+2      
      karr(i)=ie(i+2)
      end do
      ilen=ie(2)
      ilen2=1
      kbarr(1)=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      ie(2)=icont
      do i=1,icont
      ie(i+2)=ipqt(i)
      end do
      ilen=ileng
      ilen2=ileng
34    do i=1,ilen
      karr(i)=itwo(i)
      end do
      icong=icong +1
      
      
      

36    do i=1,ilen2
      kbarr(i)=itwo(i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do i=1,ilen3
      karr(i)=kcarr(i)
      
      end do
      do i=1,inlen
      kbarr(i)=mnum(i)
      end do
      ilen2=inlen
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      icont1=icont1+1
      
      do i=1,irlen
      itwo(i)=irrr(i)
      end do
      ileng =irlen
      ilen =irlen
40    if (nnum(2).lt.ie(2))goto 31
      if(nnum(2).gt.ie(2))goto 42
      do i=3,nnum(2)+2
      if(nnum(i).lt.ie(i))goto 31
      if(nnum(i).gt.ie(i))goto 42
      end do
42    do i=1,nnum(2)+2
      karr(i)=nnum(i)
      end do
      do i=1,ie(2) +2
      kbarr(i)=ie(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2) +2
      nnum(i)=kcarr(i)
      end do
      
      
      
      do i=1,ilen
      karr(i)=itwo(i)
      end do
      ilen2=1
      kbarr(1)=2
      call mpmul(ilen,ilen2,ilen3)
      do i=1,ilen3
      karr(i)=kcarr(i)
      end do
      ilen=ilen3
      do i=1,inlen
      kbarr(i)=mnum(i)
      end do
      ilen2=inlen
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      icont2=icont2+1
      
      do i=1,irlen
      itwo(i)=irrr(i)
      end do
      ileng=irlen
      goto 31

100   print *,'number to be tested=',(mnum(i),i=1,inlen)
      if(ileng.ne.1)goto 110
      if (itwo(1).ne.2)goto 110
      print *,'number is prime'
      icorp=1
      goto 112
110   print *,'number is composite'
      
      
      
      
      icorp=0
112   return
      end

      subroutine subgcd(ibig,little,igcd2)
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)

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
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      dimension karu(3000),karv1(3000),karv3(3000),karqq(3000)
      dimension kart3(3000),kart1(3000)
      
      
      
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

      subroutine mpkron(k)
      common karr(3000),kbarr(3000),kcarr(3000),ipqt(3000)
      common irrr(3000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(3000),karb(3000),kard(3000),karp(3000),karv(3000)
      
      common iarq(2)
      common marr(2000),mbarr(2000),mcarr(2000),mdarr(2000)
      dimension karae(50),karde(50),karpe(50),karh(50),ktemp(50)
      
      
      
      
      do i=3,kard(2) +2
      karr(i-2)=kard(i)
      end do
      ksgn=kard(1)
      ilen =kard(2)
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.ne.0)goto 414
      do i=3,karp(2)+2
      karr(i-2)=karp(i)
      end do
      ilen =karp(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 512
414   do i=1,karp(2) +2
      karpe(i)=karp(i)
      end do
      iv =0
      do ii=1,100
      do i=3,karpe(2)+2
      karr(i-2)=karpe(i)
      end do
      ilen=karpe(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.ne.0)goto 432
      do i=1,icont
      karpe(i+2)=ipqt(i)
      end do
      karpe(2)=icont
      iv =iv+1
      end do
      print *,'large denominator?'
      stop
432   ive2 =mod(iv,2)
      if (ive2.eq.0)goto 444
      do i=3,kard(2)+2
      karr(i-2)=kard(i)
      kbarr(i-2)=kard(i)
      end do
      ilen =kard(2)
      ilen2 =kard(2)
      call mpmul(ilen,ilen2,ilen3)
      karr(1)=0
      karr(2)=ilen3
      do i=1,ilen3
      karr(i+2)=kcarr(i)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do i=3,kcarr(2) +2
      karr(i-2)=kcarr(i)
      end do
      ilen=kcarr(2)
      kbarr(1)=8
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      karae(1)=0
      karae(2)=icont
      do i=1,icont
      karae(i+2)=ipqt(i)
      end do
      do i=3,karae(2)+2
      karr(i-2)=karae(i)
      end do
      ilen=karae(2)
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 444
      k=-1
      goto 446
444   k=1
446   do i=1,kard(2)+2
      karde(i)=kard(i)
      end do
452   if (karde(2).eq.0)goto 510
      iv =0
      do ii=1,100
      do i=3,karde(2)+2
      karr(i-2)=karde(i)
      end do
      ilen =karde(2)
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.1)goto 470
      iv =iv+1
      do i=1,icont
      karde(i+2)=ipqt(i)
      end do
      karde(1)=0
      karde(2)=icont
      end do
      print *,'warning:numerator very large'
      stop
470   ive2 =mod(iv,2)
      if (ive2.eq.0)goto 486
      do i=3,karpe(2)+2
      karr(i-2)=karpe(i)
      kbarr(i-2)=karpe(i)
      end do
      ilen =karpe(2)
      ilen2=karpe(2)
      call mpmul(ilen,ilen2,ilen3)
      karr(1)=0
      karr(2)=ilen3
      do i=1,ilen3
      karr(i+2)=kcarr(i)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      if (kcarr(2).eq.0)goto 466
      do i=3,kcarr(2)+2
      karr(i-2)=kcarr(i)
      end do
      ilen =kcarr(2)
      kbarr(1)=8
      ilen2 =1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do i=1,icont
      karae(i+2)=ipqt(i)
      karr(i)=ipqt(i)
      end do
      karae(1)=0
      karae(2)=icont
      ilen =icont
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 466
      kmul =-1
      goto 468
466   kmul =1
468   k=kmul *k
486   do i=2,karde(2)+2
      karr(i)=karde(i)
      end do
      karr(1)=ksgn
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      if (kcarr(2).eq.0)goto 490
      do i=1,kcarr(2) +2
      ktemp(i)=kcarr(i)
      end do
      do i=1,karpe(2)+2
      karr(i)=karpe(i)
      end do
      call mpadd(1)
      if (kcarr(2).eq.0)goto 490
      ilen2=kcarr(2)
      do i=3,kcarr(2) +2
      kbarr(i-2)=kcarr(i)
      end do
      do i=3,ktemp(2)+2
      karr(i-2)=ktemp(i)
      end do
      ilen =ktemp(2)
      call mpmul(ilen,ilen2,ilen3)
      ilen =ilen3
      do i=1,ilen3
      karr(i)=kcarr(i)
      end do
      kbarr(1)=4
      
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      ilen=icont
      do i=1,icont
      karr(i)=ipqt(i)
      end do
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 490
      kmul =-1
      goto 492
490   kmul =1
492   k=kmul *k
      do i=2,karde(2)+2
      karh(i)=karde(i)
      end do
      karh(1)=0
      do i=3,karde(2) +2
      kbarr(i-2)=karde(i)
      end do
      ilen2=karde(2)
      do i=3,karpe(2)+2
      karr(i-2)=karpe(i)
      end do
      ilen =karpe(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      ksgn=0
      if (irlen.eq.0)goto 502
      do i=1,irlen
      karde(i+2)=irrr(i)
      end do
502   karde(2)=irlen
      kard(1)=ksgn
      do i=1,karh(2)+2 
      karpe(i)= karh(i)
      end do
      goto 452
510   if ((karpe(2).eq.1).and.(karpe(3).eq.1))goto 513
512   k=0
513   print *,'k',k
      return
      end





      
