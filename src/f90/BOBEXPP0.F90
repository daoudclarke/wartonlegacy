    

      program bobexpp0 
!     generalization 0f bobexp1
!     first program in exponent recovery suite  using large prime variant  

      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      dimension kr(40000,2),mm(20),mres(20),n(60),msqr(20)
      dimension mmsq(20),mstem(20),zlog(10000)
      dimension iaeq(50),ibeq(50),iceq(50),ninv(50),ixx(5)
      dimension ibbeq(50),kkm(40000,2),kcalc(40000,2),ity(50)
      dimension narr(2,50),nbarr(2,50),ncarr(2,50),numb1(60),ipre(60)
      dimension icurr(60),itab(40),jtab(40),jfreq(50),littr(20),npr(10)
      real krecarr(10000),lglm
!     open(unit=2,file='fquat',access='direct',form=&
!     'formatted',recl=410,status='old')
!     do ii=1,4
!     read (2,72222,rec=ii)irecnn,iconq,(narr(ii,jf),jf=1,50),&
!     (nbarr(ii,jk),jk=1,50)
!     print *,'irecnn=',irecnn,'iconq',iconq,'a=',(narr(ii,jf),jf=1,&
!     narr(ii,2)+2),'b=',(nbarr(ii,jk),jk=1,nbarr(ii,2)+2)
!     end do
!     close(unit=2)
!     stop
      do i=1,2000
      ibbp(i)=0
      nbbp(i)=0
      end do
      iccz=0
      iconr=0
      ispecp=0
      imatc=0
      lvim=500
      npr(1)=9973
      npr(2)=9967
      npr(3)=9949
      npr(4)=9941
      npr(5)=9931
      itab(1)=750
      itab(2)=800
      itab(3)=1000
      itab(4)=1100
      itab(5)=1250
      itab(6)=1350
      itab(7)=1350
      itab(8)=1350
      itab(9)=1500
      itab(10)=2000
      itab(11)=2200
      itab(12)=2500
      itab(13)=2800
      itab(14)=2600
      itab(15)=3300
      itab(16)=3500
      itab(17)=3800
      itab(18)=4100
      itab(19)=4400
      itab(20)=4700
      itab(21)=5000
      itab(22)=5400
      itab(23)=5800
      itab(24)=6300
      itab(25)=6800
      itab(26)=7300
      itab(27)=7750
      itab(28)=8000
      itab(29)=8000
      itab(30)=8100
      itab(31)=8300
      itab(32)=8400
      itab(33)=8600
      itab(34)=8800
      itab(35)=9200
      itab(36)=9400
      itab(37)=9600
      itab(38)=9800
      jtab(1)=112
      jtab(2)=112
      jtab(3)=111
      jtab(4)=111
      jtab(5)=110
      jtab(6)=110
      jtab(7)=109
      jtab(8)=109
      jtab(9)=108
      jtab(10)=108
      jtab(11)=107
      jtab(12)=107
      jtab(13)=106
      jtab(14)=106
      jtab(15)=105
      jtab(16)=105
      jtab(17)=103
      jtab(18)=103
      jtab(19)=102
      jtab(20)=101
      jtab(21)=100
      jtab(22)=99
      jtab(23)=99
      jtab(24)=98
      jtab(25)=97
      jtab(26)=97
      jtab(27)=96
      jtab(28)=95
      jtab(29)=94
      jtab(30)=93
      jtab(31)=92
      jtab(32)=91
      jtab(33)=91
      jtab(34)=91
      jtab(35)=91
      jtab(36)=91
      jtab(37)=91
      jtab(38)=91

      iprerec=0
      iconq=1
      ineq=1
      ncom(1)=0
      ncom(2)=2
      ncom(3)=4
      ncom(4)=878
      iprar(1)=0
      iprar(2)=3
      iprar(3)=12
      iprar(4)=3456
      iprar(5)=7891
!      call bigb5(noyes)
!      print *,'ans',(iansar(jf),jf=1,iansar(2)+2)
!      print *,'noyes',noyes
      
      do i=1,10000
      rel=i
      zlog(i)= log10(rel)
      end do
      do jf=1,10000
      krecarr(jf)=127.0
      end do
      iaeq(1)=0
      iaeq(2)=7
      iaeq(3)=4
      iaeq(4)=9382
      iaeq(5)=7158
      iaeq(6)=2993
      iaeq(7)=8271
      iaeq(8)=704
      iaeq(9)=1601
      ibeq(1)=0
      ibeq(2)=6
      ibeq(3)=8
      ibeq(4)=9626
      ibeq(5)=8029
      ibeq(6)=9769
      ibeq(7)=8900
      ibeq(8)=4215
      iceq(1)=1
      iceq(2)=10
      iceq(3)=24
      iceq(4)=6913
      iceq(5)=5791
      iceq(6)=4969
      iceq(7)=1353
      iceq(8)=7209
      iceq(9)=3181
      iceq(10)=7933
      iceq(11)=2287
      iceq(12)=9538
      do jf=3,ibeq(2)+2
      karr(jf-2)=ibeq(jf)
      end do
      ilen=ibeq(2)
      kbarr(1)=2
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      ibbeq(jf+2)=kcarr(jf)
      end do
      ibbeq(1)=0
      ibbeq(2)=ilen3
      
      
      irecpl=0
      ktim=1
      kkb=0
      n(1)=0
      n(2)=17
      n(3)=12
      n(4)=1932
      n(5)=6311
      n(6)=3702
      n(7)=1795
      n(8)=2261
      n(9)=8503
      n(10)=2733
      n(11)=9089
      n(12)=7271
      n(13)=3846
      n(14)=687
      n(15)=3927
      n(16)=2124
      n(17)=207
      n(18)=1596
      n(19)=5119
      do jf=1,60
      n(jf)=0
      end do
      print *,'length of number radix 10000'
      read  *,n(2)
      print *,'number'
      read *,(n(jf),jf=3,n(2)+2)
      do jf=1,n(2)+2
      ncom(jf)=n(jf)
      end do
      
      
      
      print *,'input target legth radix 10000'
      read *,jtarg(2)
      print *,'input target'
      read *,(jtarg(jf+2),jf=1,jtarg(2))
      jtarg(1)=0
      do jf=1,20
      jfakin(jf)=0
      end do
      
      print *,'length of base radix 10000'
      read *,jbase(2)
      print *,'input base here'
      read *,(jbase(jf),jf=3,jbase(2)+2)
      jbase(1)=0
      print *,'jbase1',(jbase(jf),jf=1,jbase(2)+2)
!      open(unit=7,file='expcode',access='direct',form=&
!      'formatted',recl=240,status='new')
!      write(7,56,rec=1)(n(jf),jf=1,60)
56    format(60i4)      
      do jf=1,n(2)+2
      numb1(jf)=n(jf)
      end do
      rel=n(3)
      vlg=log10(rel)
      llg=int(vlg)
!      indy=(n(2)-1)*4+llg-38
      indy=(n(2)-1)*4+llg-13
      if (indy.gt.38)goto 2
      if (indy.le.0)goto 58
      goto 59
2     print *,'number too large to handle within a reasonable timeframe'
      print *,'recommend using GRETA suite'
!      close(unit=7)
      stop
58    indy=1
                                                 
59    limprm=itab(indy)
      mimprm=jtab(indy)
      print *,'limprm',limprm,'mimprm',mimprm,'indy',indy
      
      limprm2=limprm+limprm+limprm/5
!     compute square root of n      
      do jf=1,numb1(2)+2
      ipre(jf)=numb1(jf)
      end do
40    do jf=3,numb1(2)+2   
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
      if (icurr(2).lt.ipre(2))goto 45
      do jf=3,ipre(2)+2
      if (icurr(jf).lt.ipre(jf))goto 45
      end do
      goto 50
45    do jf=1,icurr(2)+2
      ipre(jf)=icurr(jf)
      end do
      goto 40
50    print *,'int. root',(icurr(jf),jf=1,icurr(2)+2)
      
      do jf=3,n(2)+2
      karr(jf-2)=n(jf)
      end do
      ilen=n(2)
!      ilen2=2
!      kbarr(1)=100
!      kbarr(2)=0
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      do jf=1,icont
!      karr(jf)=ipqt(jf)
!      end do
!      ilen=icont
      kbarr(1)=500
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      mntes1(jf+2)=ipqt(jf)
      end do
      mntes1(1)=0
      mntes1(2)=icont
      kbarr(1)=2000
      kbarr(2)=0
      ilen2=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      mntes2(jf+2)=ipqt(jf)
      end do
      mntes2(1)=0
      mntes2(2)=icont
      print *,'mntes1',(mntes1(jf),jf=1,mntes1(2)+2)
      print *,'mntes2',(mntes2(jf),jf=1,mntes2(2)+2)
      
      
      
      
      
      
      
      
      
      mm(1)=0
      mm(2)=7
      mm(3)=170
      mm(4)=1734
      mm(5)=3658
      mm(6)=7625
      mm(7)=7726
      mm(8)=2588
      mm(9)=2735
      

      
      irecnn=0
      kkmx=0
      
      open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
      open(unit=4,file='texp1',access='sequential')
!     print *,'ibase',(ibase(jf),jf=1,ibase(2)+2)
      

      read(3,5,rec=1) (ipr(i),i=1,65000)
5     format(65000i6)      
      
      print *,'no of primes<800001=',ipr(1),ipr(2)
      close(3)
      irecnn=0
      kkb=0
      ihitn=0
      do ii=1,20
      littr(ii)=0
      end do
      do ii=1,50
      ity(ii)=0
      end do
      limprm2=limprm
      call sieve0(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2)
      
!      write (4,*)irecnn,kkb,ihitn,nfak,(littr(jf),jf=1,20)
!      write (4,*)irecnn,(jfak(jf),jfreq(jf),ity(jf),jf=1,nfak)
      
      limprm2=limprm
      print *,'limprm2',limprm2
      
!      rel=ipr(limprm2)
!      lglm=log10(rel) 
      kia=0
      kbb=0
      icdn=0
      iconq=0
      
      do jf=1,ibase(2)+2
      norma(jf)=ibase(jf)
      
      end do
      iconarr(1)=0
      iconarr(2)=1
      iconarr(3)=1
      do jf=4,20
      iconarr(jf)=0
      end do

77    do jf=1,iconarr(2)+2
      karr(jf)=iconarr(jf)
      
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iconarr(jf)=kcarr(jf)
      end do
!      print *,'norma',(norma(jf),jf=1,norma(2)+2)
      do jf=3,ibase(2)+2
      karr(jf-2)=ibase(jf)
      end do
      ilen=ibase(2)
      
      
      do jf=3,norma(2)+2      
      
      kbarr(jf-2)=norma(jf)
      end do
      ilen2=norma(2)
      call mpmul(ilen,ilen2,ilen3)
      if (ilen3.gt.n(2))goto 80
      if (ilen3.lt.n(2))goto 79
      do jf=1,ilen3
      if (kcarr(jf).gt.n(jf+2))goto 80
      if (kcarr(jf).lt.n(jf+2))goto 79
      end do
      print *,'input number not prime'
      stop
79    do jf=1,ilen3
      norma(jf+2)=kcarr(jf)
      end do
      norma(1)=0
      norma(2)=ilen3
      goto 77

80    do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
!      print *,'ilen',ilen,'ilen2',ilen2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 86
      do jf=1,irlen
      norma(jf+2)=irrr(jf)
!      karr(jf)=irrr(jf)
      end do
      norma(1)=0
      norma(2)=irlen
!      print *,'ok sofar irlen=',irlen
!      ilen=irlen
      goto 90
      if (norma(2).le.6)goto 90
      do kjf=1,2
      ilen2=1
      kbarr(1)=npr(kjf)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 90
      end do
      goto 77

90    a=a      
      call sieve(kia,kkb,irecnn,kkmax,icdn,iconq,limprm2,iconr,ispecp,imatc,&
      iccz,lvim)
      do jf=1,nfak
      if (jfakin(jf).eq.0)goto 77
      end do




      
      if (irecnn.gt.limprm+2)goto 730
!      if (irecnn.gt.limprm+iconq-1)goto 730
72222 format(i6,i4,50i4,50i4)      
!     end of big loop      
69001 a=a
      goto 77
86    print *,'input number is not prime'      
      stop
730   print *,'run completed',' no of records=',irecnn
      
!      close (unit=2)
     do i=nfak+1,30
     jfak(i)=0
     jfreq(i)=0
     end do
     do jf=ibase(2)+3,20
     ibase(jf)=0
     end do
     do jf=jbase(2)+3,20
     jbase(jf)=0
     end do
     do jf=itarg(2)+3,20
     itarg(jf)=0
     end do
     do jf=jtarg(2)+3,20
     jtarg(jf)=0
     end do
     
     
     open (unit=2,file='exppar',access='direct',form=& 
     'formatted',recl=344,status='old') 
1002 format (i6,i6,i8,i4,20i4,20i4,20i4,20i4)      
     write (2,1002,rec=1)irecnn,limprm,ijpow,nfak,(ibase(jf),jf=1,20),& 
     (jbase(jf),jf=1,20),(itarg(jf),jf=1,20),(jtarg(jf),jf=1,20)
     open (unit=6,file='expparn',access='direct',form=&
     'formatted',recl=240,status='old')
1006 format (60i4)
     write (6,1006,rec=1)(n(jf),jf=1,60)
     close (unit=6)
     close (unit=2) 
      
      
      close(unit=4)
      
      end
      subroutine sieve0(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50),jfreq(50)
      print *,'limprm2',limprm2
!     change parameter      
      lprx1=ipr(limprm2)                                      
      print *,'lprx1',lprx1,'jtarg',(jtarg(jf),jf=1,jtarg(2)+2)
      
      do jf=iconarr(2)+3,20
      iconarr(jf)=0
      end do
      larp(1)=int(lprx1/10000)
      larp(2)=lprx1 -larp(1) *10000
      ijpow=1
      do jf=1,jbase(2)+2
      ibase(jf)=jbase(jf)
      end do
      do jf=1,jtarg(2)+2
      itarg(jf)=jtarg(jf)
      end do
440   a=a      
      do jf=1,itarg(2)+2
      litt(jf)=itarg(jf)
      end do
      
      do i=1,50
      ity(i)=0
      iabp(i)=0
      iabpn(i)=0
      end do
      


!     finish of calculation of number to be sieved
      do i=1,litt(2)+2
      litz(i)=litt(i)
      litd(i)=litt(i)
      end do
      litd(1)=0
      
      iflim=10000
      iab=1
      icur=0
!      print *,'litd',(litd(i),i=1,litd(2) +2),'icdn',icdn,'irecnn',irecnn&
!      ,kkb
      do i=1,50

      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
      if (iab.eq.limprm2)goto 1500
      if (iab.eq.16)goto 400
      if (iab.eq.96)goto 410
      goto 421
400   do jf=2,litd(2)+2
      if (litd(jf).lt.mntes1(jf))goto 420
      if (litd(jf).gt.mntes1(jf))goto 1500
      end do
      goto 420
410   do jf=1,litd(2)+2
      if (litd(jf).lt.mntes2(jf))goto 420
      if (litd(jf).gt.mntes2(jf))goto 1500
      end do

!      if (gpr(iab).gt.90.0)goto 1010
420  a=a 
!     print *,'ok pre-abort',iab      
421   a=a
      iprx1=ipr(iab)
      
      iarp(1)=int(iprx1/10000)
      iarp(2)=iprx1-iarp(1) *10000
      ibsw=0
      
      
      
      if (iarp(1).gt.larp(1))goto 1500
      if (iarp(1).lt.larp(1))goto 1162
      if (iarp(2).gt.larp(2))goto 1500
1162  if (iarp(1).gt.ipfar(1))goto 1200
      if (iarp(1).lt.ipfar(1))goto 1166
      if (iarp(2).gt.ipfar(2))goto 1200


1166  do i=3,litd(2)+2
      karr(i-2)= litd(i)
      end do
      
      
      ilen=litd(2)
      if (iprx1.lt.10000)goto 1168
      ilen2 =2
      kbarr(1)=iarp(1)
      kbarr(2)=iarp(2)
      goto 1170
1168  kbarr(1)=iprx1
      ilen2=1
1170  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      
      if (icont.eq.0)goto 1300
      if(irlen.gt.0)goto 1010
      
      
      do i=1,icont
      litd(i+2)=ipqt(i)
      end do
      litd(2)=icont
      
      
      if(ibsw.eq.1)goto 1176
      ibsw=1
      icur=icur +1
      iabp(icur)=iprx1
      iabpn(icur)=1
      ity(icur)=1
      goto 1166
1176  iabpn(icur)=iabpn(icur)+1
      goto 1166
      


1200  if (litd(2).ne.iprex(2))goto 12001
      do i=3,litd(2)+2
      if (litd(i).ne.iprex(i))goto 12001
      end do
      goto 1204

 
12001 if (litd(2).gt.2)goto 2020
      if (litd(2).eq.1)goto 1166
      if (litd(3).lt.larp(1))goto 1166
      if (litd(3).gt.larp(1))goto 2020
      if (litd(4).le.larp(2))goto 1166
2020  do i=3,litd(2) +2
      mnum(i-2)=litd(i)
      end do
      inlen=litd(2)

      call mpprime(icorp,inlen)
      if (icorp.eq.1)goto 1500
      do i=1,litd(2)+2
      iprex(i)=litd(i)
      end do
1204  ipfar(1)=ipfar(1) +2
      goto 1166


      

1300  irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kkb,(iconarr(i),i=1,&
      iconarr(2)+2),'iconq',iconq
      if (icur.le.kkmx)goto 1302
      kkmx=icur

1302  do i=1,20
      normar(i)=0
      littr(i)=0
      end do
!      do i=1,litz(2) +2
!      littr(i)=litz(i)
!      end do
!      do i=1,norma(2)+2
!      normar(i)=norma(i)
!      end do
!      ihitn=icdn*10000+kia
      ihitn=0
!      print *,'nfak',nfak
      nfak=icur
      do jf=1,nfak
      jfak(jf)=iabp(jf)
      end do
      do jf=1,nfak
      jfreq(jf)=iabpn(jf)
      end do
      print *,'ijpow',ijpow
      print *,'nfak',nfak      
      print *,'jfaks',(jfak(jf),jf=1,nfak)
      print *,'jfreqs',(jfreq(jf),jf=1,nfak)
      
      
      
      
      
      write(4,*)irecnn,kkb,ihitn,icur,(iconarr(j1),j1=1,20)
      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      




!     note next instruction      
      if (irecnn.lt.190000)goto 1501
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      

1500  ijpow=ijpow+1
      ilen=jbase(2)
      do jf=3,jbase(2)+ 2
      karr(jf-2)=jbase(jf)
      end do
      ilen2=ibase(2)
      do jf=3,ibase(2)+2
      kbarr(jf-2)=ibase(jf)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ncom(2)+2
      kbarr(jf-2)=ncom(jf)
      end do
      ilen2=ncom(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,irlen
      ibase(jf+2)=irrr(jf)
      end do
      ibase(1)=0
      ibase(2)=irlen
      
      
      ilen=jtarg(2)
      do jf=3,jtarg(2)+2
      karr(jf-2)=jtarg(jf)
      end do
      ilen2=itarg(2)
      do jf=3,itarg(2)+2
      kbarr(jf-2)=itarg(jf)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      
      ilen=ilen3
      
      ilen2=ncom(2)
      do jf=3,ncom(2)+2
      kbarr(jf-2)=ncom(jf)
      end do
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,irlen
      itarg(jf+2)=irrr(jf)
      end do
      itarg(1)=0
      itarg(2)=irlen
!      print *,'ijpow',ijpow,'icont',icont,'itarg',(itarg(jf),jf=1,&
!      itarg(2)+2),'litd',(litd(jf),jf=1,litd(2)+2)
      if (icont.eq.0)goto 1500
      if (mod(ijpow,2).eq.0)goto 1500
      goto 440

      
      



1501  return      

      
      end
      
      subroutine sieve(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2,iconr,ispecp,&
      imatc,iccz,lvim)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50),iran(30),iranb(30)
      jfx1=norma(2)
      do jf=1,6
      iran(jf)=96
      iranb(jf)=16
      end do
      do jf=7,12
      iran(jf)=96
      iranb(jf)=16
      end do


!     change parameter      
      lprx1=ipr(limprm2)                                      
      do jf=iconarr(2)+3,20
      iconarr(jf)=0
      end do
      larp(1)=int(lprx1/10000)
      larp(2)=lprx1 -larp(1) *10000
            
      do jf=1,norma(2)+2
      litt(jf)=norma(jf)
      end do
      
      do i=1,50
      ity(i)=0
      iabp(i)=0
      iabpn(i)=0
      end do
      


!     finish of calculation of number to be sieved
      do i=1,litt(2)+2
      litz(i)=litt(i)
      litd(i)=litt(i)
      end do
      litd(1)=0
      
      iflim=10000
      iab=1
      icur=0
!      print *,'litd',(litd(i),i=1,litd(2) +2),'icdn',icdn,'irecnn',irecnn&
!      ,kkb
      do i=1,50

      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
      if (iab.eq.limprm2)goto 1490
      if (iab.eq.iranb(jfx1))goto 400
      if (iab.eq.iran(jfx1))goto 410
      goto 421
400   do jf=2,litd(2)+2
      if (litd(jf).lt.mntes1(jf))goto 420
      if (litd(jf).gt.mntes1(jf))goto 1500
      end do
      goto 420
410   do jf=1,litd(2)+2
      if (litd(jf).lt.mntes2(jf))goto 420
      if (litd(jf).gt.mntes2(jf))goto 1500
      end do

!      if (gpr(iab).gt.90.0)goto 1010
420  a=a 
!     print *,'ok pre-abort',iab      
421   a=a
      iprx1=ipr(iab)
      
      iarp(1)=int(iprx1/10000)
      iarp(2)=iprx1-iarp(1) *10000
      ibsw=0
      
      
      
      if (iarp(1).gt.larp(1))goto 1500
      if (iarp(1).lt.larp(1))goto 1162
      if (iarp(2).gt.larp(2))goto 1500
1162  if (iarp(1).gt.ipfar(1))goto 1200
      if (iarp(1).lt.ipfar(1))goto 1166
      if (iarp(2).gt.ipfar(2))goto 1200


1166  do i=3,litd(2)+2
      karr(i-2)= litd(i)
      end do
      
      
      ilen=litd(2)
      if (iprx1.lt.10000)goto 1168
      ilen2 =2
      kbarr(1)=iarp(1)
      kbarr(2)=iarp(2)
      goto 1170
1168  kbarr(1)=iprx1
      ilen2=1
1170  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      
      if (icont.eq.0)goto 1301
      if(irlen.gt.0)goto 1010
      
      
      do i=1,icont
      litd(i+2)=ipqt(i)
      end do
      litd(2)=icont
      
      
      if(ibsw.eq.1)goto 1176
      ibsw=1
      icur=icur +1
      iabp(icur)=iprx1
      iabpn(icur)=1
      ity(icur)=1
      goto 1166
1176  iabpn(icur)=iabpn(icur)+1
      goto 1166
      


1200  if (litd(2).ne.iprex(2))goto 12001
      do i=3,litd(2)+2
      if (litd(i).ne.iprex(i))goto 12001
      end do
      goto 1204

 
12001 if (litd(2).gt.2)goto 2020
      if (litd(2).eq.1)goto 1166
      if (litd(3).lt.larp(1))goto 1166
      if (litd(3).gt.larp(1))goto 2020
      if (litd(4).le.larp(2))goto 1166
2020  do i=3,litd(2) +2
      mnum(i-2)=litd(i)
      end do
      inlen=litd(2)

      call mpprime(icorp,inlen)
      if (icorp.eq.1)goto 1490
      do i=1,litd(2)+2
      iprex(i)=litd(i)
      end do
1204  ipfar(1)=ipfar(1) +2
      goto 1166
1490  if (litd(2).gt.2)goto 1500
      if (litd(3).gt.lvim)goto 1500
      icur=icur+1
      iabp(icur)=litd(3)*10000+litd(4)
      iabpn(icur)=iabpn(icur)+1
      ity(icur)=1
      if (iccz.eq.0)goto 1499
      do jf=1,iccz
      if (iabp(icur).ne.ibbp(jf))goto 1498
      iconr=iconr+1
      nbbp(iccz)=nbbp(iccz)+1
      if (nbbp(iccz).ne.2)goto 200
      iconr=iconr+1
      ispecp=ispecp+1
200   imatc=imatc+1
      goto 1300
1498  end do
1499  iccz=iccz+1
      ibbp(iccz)=iabp(icur)
      nbbp(iccz)=nbbp(iccz)+1
      print *,'ibbp',ibbp(iccz)
      goto 1300
1301  iconr=iconr+1      

1300  irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kkb,(iconarr(i),i=1,&
      iconarr(2)+2),'iconr',iconr,'ispecp',ispecp,'imatc',imatc
      if (icur.le.kkmx)goto 1302
      kkmx=icur

1302  do i=1,20
      normar(i)=0
      littr(i)=0
      end do
!      do i=1,litz(2) +2
!      littr(i)=litz(i)
!      end do
!      do i=1,norma(2)+2
!      normar(i)=norma(i)
!      end do
!      ihitn=icdn*10000+kia
      ihitn=0
!      print *,'nfak',nfak
      do ii=1,nfak      
      if (jfakin(ii).eq.1)goto 435
      do jj=1,icur
      if (iabp(jj).eq.jfak(ii))goto 430
      goto 434
430   jfakin(ii)=1      
      goto 436
434   end do
435   end do
436   a=a
      write(4,*)irecnn,kkb,ihitn,icur,(iconarr(j1),j1=1,20)
      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      
!     note next instruction      
      if (irecnn.lt.190000)goto 1500
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      

      

1500  return      
      end
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      character*(*) in,out
      out =in
      i =len(out)
      do while(out(i:i) ==' ')
      out(i:i)='*'
      i = i-1
      end do 
      return 
      end

      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      do i=1,200
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      dimension kdum(50),isub(50)
      
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtaerg(50)
      dimension ie(50),itwo(100),nnum(50)
      
      
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

   
100   if(ileng.ne.1)goto 110
      if (itwo(1).ne.2)goto 110
      
      icorp=1
      goto 112
110   a=a
      
      
      
      
      icorp=0
112   return
      end

      subroutine subgcd(ibig,little,igcd2)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),ntes1(50),ntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      dimension karu(50),karv1(50),karv3(50),karqq(50)
      dimension kart3(50),kart1(50)
      
      
      
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
      kart3(1)=0
      kart3(2)=irlen
      do i=1,irlen
      kart3(i+2)=irrr(i)
      end do
      
      
      
      goto 51
41    kart3(1)=0
      kart3(2)=0
51    if (karv1(2).eq.0)goto 6
      
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50),ibbp(2000),nbbp(2000)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
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
       subroutine bwq5(ip,ians)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50),ibbp(2000),nbbp(2000)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(40000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
       common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
       dimension itempz(2)
       print *,'ncom',(ncom(jf),jf=1,ncom(2)+2)
       
       
       do jf=3,ncom(2)+2
       karr(jf-2)=ncom(jf)
       end do
       ilen=ncom(2)
       if (ip.lt.10000)goto 1100
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 1102
1100   kbarr(1)=ip
       ilen2=1
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       if (irlen.eq.1)goto 1106
       iaa=irrr(1)*10000+irrr(2)
       goto 1108
1104   print *,'factor found early',ip
       stop
1106   iaa=irrr(1)


       
       
1108   ix = 0
       
       
       
       
       
10     iprecod =ip -1
       i=0
26     itemp = int(iprecod/2)
       irem1 =iprecod -itemp*2
       if(itemp.eq.0)goto 46
       if(irem1.gt.0)goto 40
       iprecod = itemp
       i =i+1
       if(i.lt.200)goto 26
40     iq =iprecod
       ie = i
       goto 48
46     iq =1
       ie = i
48     i =1
       print *,'qs',iq
       n = 1
52     n = n*607
54     itemp=int(n/1000)
56     irem1 =n-1000 *itemp
       n = irem1
       id = n
       if (ip.lt.10000)goto 5901
       karp(1)=0
       karp(2)=2
       karp(3)=int(ip/10000)
       karp(4)=ip-karp(3)*10000
       goto 5902
5901   karp(1)=0       
       karp(2)=1
       karp(3)=ip
5902   kard(1)=0
       kard(2)=1
       kard(3)=id

       call mpkron(k)
       
       print *,'k',k,'ip',ip
       if(k.eq.-1)goto 68
       
       
       i =i +1
       if(i.lt.1000)goto 52
68     ipn =iq
       print *,'k',k,'n',n,'ip',ip
       
       
5601   iaas = n
       call sub516(ibprod,iaas,ipn,ip)
       iz = ibprod
       print *,'zz',iz
       
       iy = iz
       ir = ie
       ipn = (iq-1)/2
       iaas = iaa
       call sub516(ibprod,iaas,ipn,ip)
       ix = ibprod
       print *,ix
       if (ix.lt.10000)goto 140
       karr(1)=int(ix/10000)
       karr(2)=ix-karr(1)*10000
       ilen=2
       goto 142
140    karr(1)=ix       
       ilen=1
142    if (iaa.lt.10000)goto 150       
       kbarr(1)=int(iaa/10000)
       kbarr(2)=iaa-kbarr(1)*10000
       ilen2=2
       goto 152
150    kbarr(1)=iaa
       ilen2=1
152    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 160
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 162
160    kbarr(1)=ip
       ilen2=1
162    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 169
       do jf=1,irlen
       karr(jf)=irrr(jf)
       itempz(jf)=irrr(jf)
       end do
       itempzl=irlen
       ilen=irlen
       goto 170
169    print *,'error type 2'
       stop
170    if (ix.lt.10000)goto 172
       kbarr(1)=int(ix/10000)
       kbarr(2)=ix-kbarr(1)*10000
       ilen2=2
       goto 174
172    kbarr(1)=ix
       ilen2=1
174    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 181
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 182
181    kbarr(1)=ip
       ilen2=1
182    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 169
       if (irlen.eq.1)goto 190
       ib=irrr(1)*10000+irrr(2)
       goto 192
190    ib=irrr(1)
192    if (itempzl.eq.1)goto 196
       ix=itempz(1)*10000+itempz(2)
       goto 100
196    ix=itempz(1)














       
       
       
       
100    irem1 =mod(ib,ip)
       
       if(irem1.eq.1)goto 200
       i=1
       m=1
110    ipow =2**m
       goto 700
112    irem2 =ibprod
       print *,'irem2',irem2,'ip=',ip
       if(irem2.eq.1)goto 130
       m =m+1
       goto 110
130    if(m.eq.ir)goto 180
       ipow =2**(ir-m-1)
       goto 800
134    it = ibprod
       print *,'it',it
       if (it.lt.10000)goto 2000
       karr(1)=int(it/10000)
       karr(2)=it-karr(1)*10000
       kbarr(1)=karr(1)
       kbarr(2)=karr(2)
       ilen=2
       ilen2=2
       goto 2002
2000   karr(1)=it       
       kbarr(1)=it
       ilen=1
       ilen2=1
2002   call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 2004
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000 
       ilen2=2
       goto 2006
2004   kbarr(1)=ip
       ilen2=1
2006   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2010
       iy=irrr(1)*10000+irrr(2)
       goto 2012
2008   print *,'error type 3'
       stop
2010   iy=irrr(1)
2012   ir=mod(m,ip)
       
       if (ix.lt.10000)goto 2014
       karr(1)=int(ix/10000)
       karr(2)=ix-karr(1)*10000
       ilen=2
       goto 2016
2014   karr(1)=ix
       ilen=1
2016   if (it.lt.10000)goto 2018
       kbarr(1)=int(it/10000)
       kbarr(2)=it-kbarr(1)*10000
       ilen2=2
       goto 2020
2018   kbarr(1)=it
       ilen2=1

2020   call mpmul(ilen,ilen2,ilen3)
       
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       
       if (ip.lt.10000)goto 2022
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 2024
2022   kbarr(1)=ip
       ilen2=1
2024   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2026
       ix=irrr(1)*10000+irrr(2)
       goto 2028
2026   ix=irrr(1)
2028   if (ib.lt.10000)goto 2030
       print *,'ix',ix
       karr(1)=int(ib/10000)
       karr(2)=ib-karr(1)*10000
       ilen=2
       goto 2032
2030   karr(1)=ib
       ilen=1
2032   if (iy.lt.10000)goto 2034
       kbarr(1)=int(iy/10000)
       kbarr(2)=iy-kbarr(1)*10000
       ilen2=2
       goto 2036
2034   kbarr(1)=iy
       ilen2=1
2036   call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 2038
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 2040
2038   kbarr(1)=ip
       ilen2=1
   

2040   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2042
       ib=irrr(1)*10000+irrr(2)
       
       goto 100
2042   ib=irrr(1)
       print *,'ib1',ib
       goto 100



       
       
       
180    print *,'no square root exists'
       
       stop
       goto 900
200    print *,'square root=',ix
       ians=ix
       
       
       goto 900
700    ipn =ipow
       if (ipn.gt.10000000)goto 950
       
       iaas =ib
       call sub516(ibprod,iaas,ipn,ip)
       goto 112
800    ipn =ipow
       iaas=iy
       call sub516(ibprod,iaas,ipn,ip)
       goto 134
900    goto 1000       
       print *,'loop too short'
       stop
950    ians=999999
1000   return       
       end
       subroutine sub400(id,ip,k)
       ide =int(id/2)
       ipe= int(ip/2)
       ide2 =id -ide *2
       ipe2 =ip -ipe*2
       if(ide2.eq.0)goto 414
       goto 416
414    if(ipe2.eq.0)goto 512
416    iv = 0
       ipe = ip
       ii = 0
419    print *,'pefirst',ipe
       ipe2 =int(ipe/2)
       ipe3 =ipe -ipe2 *2
       if(ipe3.eq.1)goto 432
       ipe =ipe2
       iv = iv+1
       ii =ii +1
       if(ii.lt.50)goto 419
432    ive = int(iv/2)
       ive2 =iv -ive*2
       if(ive2.eq.0)goto 450
       iae =(id **2 -1)/8
       iae1 =int(iae/2)
       iae2 =iae -iae1 *2
       iae3 =iae2 +2
       k=(-1)**iae3
       goto 451
450    k =1
451    ide = id
452    if(ide.eq.0)goto 510 
       ive = 0
       ii = 0
456    ide1 =int(ide/2)
       ide2 =ide -ide1 *2
       if(ide2.eq.1)goto 470
       iv = iv+1
       ide =ide1
       ii =ii+1
       if(ii.lt.50)goto 456
470    ive =int(iv/2)
       ive2=iv -ive *2
       if(ive2.eq.0)goto 486
       iae=(ipe **2 -1)/8
       iae2 =int(iae/2)
       iae3 =iae -iae2 *2
       iae3 = iae3 +2
       
       
       
       
       k = (-1)**iae3*k
486    iae2 =((ide -1)*(ipe-1))/4       
       iae3 =int(iae2/2)
       iae4 =iae2-iae3 *2
       iae4 =iae4 +2
       k =(-1) **iae4 *k
       ir=abs(ide)
       
       itemp=int(ipe/ir)
       ide =ipe -itemp*ir
       ipe = ir
       
       goto 452
510    if(ipe.eq.1)goto 513
512    k =0
513    print *,k,ipe
       return
       end
       
       subroutine sub516(ibprod,iaas,ipn,ip)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50),ibbp(2000),nbbp(2000)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(40000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
       common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
       iy=iaas
       n=ipn
       ipow=2
14     if (ipow.gt.n)goto 20
       ipow=ipow*2
       goto 14
20     ie=int(ipow/2)
       n1=n
       n1=n1-ie
26     if (ie.eq.1)goto 100
       ie=int(ie/2)
       if (iy.lt.10000)goto 200
       karr(1)=int(iy/10000)
       karr(2)=iy-karr(1)*10000
       ilen=2
       kbarr(1)=karr(1)
       kbarr(2)=karr(2)
       ilen2=2
       goto 202
200    karr(1)=iy       
       ilen=1
       kbarr(1)=iy
       ilen2=1

202    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 220
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 222
220    kbarr(1)=ip
       ilen2=1
222    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 230
       if (irlen.eq.1)goto 240
       iy=irrr(1)*10000+irrr(2)
       goto 242
230    print *,'error number 2'
       iy=0
240    iy=irrr(1)
242    if (n1.lt.ie)goto 26
       n1=n1-ie
       if (iy.lt.10000)goto 250
       karr(1)=int(iy/10000)
       karr(2)=iy-karr(1)*10000
       ilen=2
       goto 252
250    karr(1)=iy
       ilen=1
252    if (iaas.lt.10000)goto 260       
       kbarr(1)=int(iaas/10000)
       kbarr(2)=iaas-kbarr(1)*10000
       ilen2=2
       goto 262
260    kbarr(1)=iaas       
       ilen2=1
262    call mpmul(ilen,ilen2,ilen3)        
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 270
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 272
270    kbarr(1)=ip
       ilen2=1
272    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 230
       if (irlen.eq.1)goto 280
       iy=irrr(1)*10000+irrr(2)
       goto 26
280    iy=irrr(1)
       goto 26

       
       
100    ibprod=iy
       
620    return
       end
       










