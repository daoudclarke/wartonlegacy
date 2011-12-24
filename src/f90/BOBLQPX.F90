    

      program boblqpx 
!     generalization 1 :program uses linear sieve 

!     MPQS program applied to exponent recovery problem or discrete    
!     logarithm problem
!     smaller sieving interval
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      dimension kr(40000,2),mm(30),mres(30),n(60),msqr(20)
      dimension mmsq(50),mstem(20),zlog(10000),mmsq2(50)
      dimension iaeq(50),ibeq(50),iceq(50),ninv(50),ixx(10)
      dimension ibbeq(50),kkm(40000,2),kcalc(40000,2)
      dimension narr(2,50),nbarr(2,50),ncarr(2,50),numb1(60),ipre(60)
      dimension icurr(60),itab(40),jtab(40)
      dimension mr(40000),mkm(40000),itemp1(50),itemp2(50)
      dimension littr(20),iabp(50),iabpn(50),ity(50),icdarr(50)
      dimension matem(50),mntem(50),mainv(5),maans(5),mmrem(5),imm(50)
      real krecarr(10000),lglm,mrecarr(10000)
      ibbbsw=0
      irecnn=0
      icdn=0
! mm for 32 digit no.
      mm(1)=0
      mm(2)=5
      mm(3)=4969
      mm(4)=399
      mm(5)=2763
      mm(6)=8851
      mm(7)=8310
      mm(8)=4156
      mm(8)=4157
      do jf=1,mm(2)+2
      imm(jf)=mm(jf)
      end do
      iva=mm(2)+2
      ivan=mm(iva)
      if (mod(ivan,2).eq.1)goto 910
      do jf=1,mm(2)+2
      karr(jf)=mm(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      mm(jf)=kcarr(jf)
      end do
910   a=a      
      do jf=3,mm(2)+2
      karr(jf-2)=mm(jf)
      
      end do
      do jf=1,imm(2)+2
      kbarr(jf-2)=imm(jf)
      end do
      
      ilen=mm(2)
      ilen2=imm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      mmsq(jf+2)=kcarr(jf)
      end do
      mmsq(1)=0
      mmsq(2)=ilen3
      print *,'mmsq',(mmsq(jf),jf=1,mmsq(2)+2)
      


!      open(unit=4,file='gpx1',access='sequential')
      icon=0
      goto 11
10    read(4,*)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      read(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      icon=icon+1
      print *,'irecnn',irecnn,'littr',(littr(jf),jf=1,20)
      do i=1,icur
!      if (ity(i).eq.2)goto 12
      print *,'iabp',iabp(i),'iabpn',iabpn(i),'ity',ity(i)
12    end do
      if (icon.ne.2)goto 10
      stop
11    a=a
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
      do jf=1,20
      jfakin(jf)=0
      kjfakin(jf)=0
      end do
!     for 20 digit no itab=250
      itab(1)=150
!     for 28 digit no itab1=500 
!      itab(1) was 1000
!     for 32 digit no itab1=600      
!      itab(1)=500
      itab(1)=1000
      itab(2)=1800
      itab(3)=2000
      itab(4)=2100
      itab(5)=2250
      itab(6)=1350
      itab(7)=1500
      itab(8)=1600
      itab(9)=1750
! try 5000 for 48 digit no.
      itab(9)=7000
      itab(10)=1850
      itab(11)=2000
      itab(12)=2300
      itab(13)=2600
      itab(14)=2900
      itab(15)=3200
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
!     for 20 digit prime jtab1=116
      jtab(1)=116
!     for 28 digit prime jtab1=112
      jtab(1)=114
      jtab(2)=112
      jtab(3)=111
      jtab(4)=111
      jtab(5)=110
      jtab(6)=110
      jtab(7)=109
      jtab(8)=109
      jtab(9)=108
! try jtab9=100 for 48 digit no.
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
      call bigb5(noyes)
      print *,'ans',(iansar(jf),jf=1,iansar(2)+2)
      print *,'noyes',noyes
      
      do i=1,10000
      rel=i
      zlog(i)= log10(rel)
      end do
      do jf=1,10000
      krecarr(jf)=127.0
      mrecarr(jf)=127.0
      end do
      
      
      
      irecpl=0
      ktim=1
      kkb=0
      n(1)=0
      print *,'length of number radix 10000'
      read  *,n(2)
      print *,'number'
      read *,(n(jf),jf=3,n(2)+2)
      
      print *,'length of base radix 10000'
      read *,jbase(2)
      print *,'input base here'
      read *,(jbase(jf),jf=3,jbase(2)+2)
      jbase(1)=0
      print *,'input target length radix 10000'
      read *,jtarg(2)
      print *,'input target'
      read *,(jtarg(jf+2),jf=1,jtarg(2))
      jtarg(1)=0
!      open(unit=7,file='expparn',access='direct',form=&
!      'formatted',recl=240,status='old')
!      write(7,56,rec=1)(n(jf),jf=1,60)
56    format(60i4)      
      do jf=1,n(2)+2
      karr(jf)=n(jf)
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      numb1(jf)=kcarr(jf)
      end do
      rel=n(3)
      vlg=log10(rel)
      llg=int(vlg)
      indy=(n(2)-1)*4+llg-38
      if (n(2).lt.4)goto 62
      if (indy.gt.38)goto 2
      if (indy.le.0)goto 58
      goto 59
62    limprm=75
      mimprm=116
      indy=0
      goto 63
2     print *,'number too large to handle within a reasonable timeframe'
      print *,'recommend using GRETA suite'
      close(unit=7)
      stop
58    indy=1

59    limprm=itab(indy)
      mimprm=jtab(indy)
63    a=a      
      print *,'limprm',limprm,'mimprm',mimprm,'indy',indy
      




      do jf=1,n(2)+2
      ncom(jf)=n(jf)
      end do
      do jf=1,mmsq(2)+2
      karr(jf)=mmsq(jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      mmsq2(jf)=kcarr(jf)
      end do
      
      
      
      
      
      print *,'mm',(mm(jf),jf=1,mm(2)+2)
      print *,'mmsq',(mmsq(jf),jf=1,mmsq(2)+2)
      print *,'mmsq2',(mmsq2(jf),jf=1,mmsq2(2)+2)
      print *,'n',(n(jf),jf=1,n(2)+2)
      
      limprm2=limprm+limprm+limprm/4
!     for extra small numbers bypass square root phase      
      if (n(2).lt.4)goto 60
!     compute square root of 2 times n      
60    a=a      
      
      open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
!      open(unit=4,file='gpx1',access='sequential')

      

      read(3,5,rec=1) (ipr(i),i=1,65000)
5     format(65000i6)      
      
      print *,'no of primes<800001=',ipr(1),ipr(2)
      close(3)
      print *,'limprm2',limprm2
      rel=ipr(limprm)
      lglm=log10(rel) 
      
      






      
      
      
      iblim=0
      iulim=2000
      lprx1=100000
      icon=0
      isivl=2000
      
      
      
      
400   do i=1,limprm
      mr(i)=999999
      kr(i,1)=999999
      kr(i,2)=999999
      end do
      
      
      goto 4004
4002  print *,'factor found early=2'
      stop


4004  iva=n(2)+2
      goto 7000
      iccn=0
      karp(1)=0
      do i=3,limprm2
      do jf=1,n(2)+2
      kard(jf)=n(jf)
      end do
      print *,'kard',(kard(jf),jf=1,kard(2)+2)
      if (ipr(i).lt.10000)goto 500
      karp(1)=0
      karp(2)=2
      karp(3)=int(ipr(i)/10000)
      karp(4)=ipr(i)-karp(3)*10000
      goto 502
500   karp(2)=1
      karp(1)=0
      karp(3)=ipr(i)

502   call mpkron(k)
      if (k.ne.1)goto 590
      print *,'i',i,'ipr',ipr(i),'k',k
      ip=ipr(i)
5901  call bwq5(ip,ians)
      kr(i,1)=ians
      kr(i,2)=ipr(i)-ians
      iccn=iccn+1
      print *,'kr',kr(i,1),'ipr',ipr(i),'iccn',iccn
      if (iccn.eq.limprm)goto 512
590   end do      
      print *,'loop too short'
      stop
512   print *,'square root phase finished'
      print *,'kr',(kr(jf,1),jf=2,11) 
      print *,'kr2',(kr(jf,2),jf=2,11)
      print *,'ipr',(ipr(jf),jf=2,11)
7000 a=a      
      
      
      
      do jf=3,n(2)+2
      karr(jf-2)=n(jf)
      end do
      ilen=n(2)
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
      kr(2,1)=1
! new inst. year 2001
      goto 440

      do i=2,limprm
      
!      do kf=1,2
!      if (kr(i,kf).eq.999999)goto 5150
      do jf=3,mm(2)+2
      karr(jf-2)=mm(jf)
      end do
      ilen=mm(2)
      if (ipr(i).lt.10000)goto 5122
      kbarr(1)=int(ipr(i)/10000)
      kbarr(2)=ipr(i)-kbarr(1)*10000
      ilen2=2
      goto 5124
5122  kbarr(1)=ipr(i)
      ilen2=1
5124  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 5130
      if (irlen.eq.1)goto 5140
      itemp=irrr(1)*10000+irrr(2)
      goto 5142
5130  itemp=0
      goto 5142
5140  itemp=irrr(1)
      
5142  mr(i)=ipr(i)-itemp
      
5153  a=a

!      if (i.eq.2)goto 5150
!      if (kr(i,kf).eq.999999)goto 5150
!      kr(i,kf)=kr(i,kf)-itemp
!      if (kr(i,kf).ge.0)goto 5150
!      kr(i,kf)=kr(i,kf)+ipr(i)
!5150  end do
      end do
!      print *,'krs',(kr(jf,1),jf=3,20)
461   ivan=mm(2)+2
      ivan2=mm(ivan)
!      if (mod(ivan2,2).eq.1)goto 439
!      kr(2,1)=1
!      goto 440
!439   kr(2,1)=2
440   print *,'limprm2 at mid',limprm2
      iconq=0
900    inlen=mm(2)
!       goto 901
       
       do jf=3,mm(2)+2
       
       mnum(jf-2)=mm(jf)
       end do
       call mpprime(icorp,inlen)
       
       if (icorp.eq.1)goto 901
       
       do jf=1,mm(2)+2
       karr(jf)=mm(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2
       call mpadd(0)
       
       do jf=1,kcarr(2)+2
       mm(jf)=kcarr(jf)
!       imm(jf)=kcarr(jf)
       end do
       goto 900

901    print *,'prime for new poly',(mm(jf),jf=1,mm(2)+2)
       iconq=iconq+1
!       goto 902
       do jf=3,n(2)+2
       karr(jf-2)=n(jf)
       end do
       ilen=n(2)
       do jf=3,mm(2)+2
       kbarr(jf-2)=mm(jf)
       end do
       ilen2=mm(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       karr(jf+2)=ipqt(jf)
       end do
       karr(1)=0
       karr(2)=icont
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       imm(jf)=kcarr(jf)
       end do


      do jf=3,mm(2)+2
      karr(jf-2)=mm(jf)
      
      end do
      do jf=1,imm(2)+2
      kbarr(jf-2)=imm(jf)
      end do
      
      ilen=mm(2)
      ilen2=imm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      mmsq(jf+2)=kcarr(jf)
      end do
      mmsq(1)=0
      mmsq(2)=ilen3



      do jf=1,mmsq(2)+2
      karr(jf)=mmsq(jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      mmsq2(jf)=kcarr(jf)
      end do
      print *,'mm',(mm(jf),jf=1,mm(2)+2)
      print *,'imm',(imm(jf),jf=1,imm(2)+2)
      print *,'n',(n(jf),jf=1,n(2)+2)
      print *,'mmsq',(mmsq(jf),jf=1,mmsq(2)+2)
      print *,'mmsq2',(mmsq2(jf),jf=1,mmsq2(2)+2)
      if (iconq.lt.3)goto 902
!      stop






902    a=a






      do i=2,limprm
!       do i=2,11
!      do kf=1,2
!      if (kr(i,kf).eq.999999)goto 7150
       if (ipr(i).lt.10000)goto 7050
       iprar(1)=0
       iprar(2)=2
       iprar(3)=int(ipr(i)/10000)
       iprar(4)=ipr(i)-iprar(3)*10000
       goto 7051
7050   iprar(1)=0
       iprar(2)=1
       iprar(3)=ipr(i)
7051   a=a
       do jf=3,mm(2)+2
       karr(jf-2)=mm(jf)
       end do
       ilen=mm(2)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       if (irlen.eq.0)goto 8150
       do jf=1,irlen
       matem(jf+2)=irrr(jf)
       end do
       matem(1)=0
       matem(2)=irlen
       do jf=1,iprar(2)+2
       kara(jf)=iprar(jf)
       end do
       do jf=1,matem(2)+2
       karb(jf)=matem(jf)
       end do
       print *,'iprar',ipr(i)
       print *,'matem',(matem(jf),jf=1,matem(2)+2)
       
       call mpgcd
       do jf=1,karv(2)+2
       mainv(jf)=karv(jf)
       end do
       print *,'i',i,'mainv',(mainv(jf),jf=1,mainv(2)+2)
       do jf=3,mmsq2(2)+2
       karr(jf-2)=mmsq2(jf)
       end do
       ilen=mmsq2(2)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       if (irlen.eq.0)goto 7200
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       
       end do
       kbarr(1)=0
       kbarr(2)=irlen
       do jf=1,iprar(2)+2
       karr(jf)=iprar(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       mmrem(jf)=kcarr(jf)
       end do

       goto 7201
7200   kbarr(1)=0
       kbarr(2)=0
       mmrem(1)=0
       mmrem(2)=0
7201   a=a
       print *,'mmrem',(mmrem(jf),jf=1,mmrem(2)+2)
       do jf=1,mmrem(2)+2
       mntem(jf)=mmrem(jf)
       end do
       do jf=3,mainv(2)+2
       karr(jf-2)=mainv(jf)
       end do
       ilen=mainv(2)
       do jf=3,mntem(2)+2
       kbarr(jf-2)=mntem(jf)
       end do
       ilen2=mntem(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       if (irlen.eq.0)goto 7300
       do jf=1,irlen
       maans(jf+2)=irrr(jf)
       end do
       maans(1)=0
       maans(2)=irlen
       goto 7301
7300   maans(1)=0
       maans(2)=0
7301   a=a
              
       print *,'maans',(maans(jf),jf=1,maans(2)+2)



       
8210   isum=0
       do jf=3,maans(2)+2
       isum=isum+maans(jf)*10000**(maans(2)+2-jf)
       end do
       kr(i,1)=isum
       print *,'maansnew',(maans(jf),jf=1,maans(2)+2)
       print *,'kr i',kr(i,1)
!      if (i.eq.2)goto 8150
!      if (kr(i,kf).eq.999999)goto 8150
!      kr(i,kf)=kr(i,kf)-itemp
!      if (kr(i,kf).ge.0)goto 8150
!      kr(i,kf)=kr(i,kf)+ipr(i)
!8150  end do
8150   end do
       print *,'ipr',(ipr(jf),jf=2,11)
       print *,'kr',(kr(jf,1),jf=2,11)
       print *,'mmsq',(mmsq(jf),jf=1,mmsq(2)+2)
       print *,'mmsq2',(mmsq2(jf),jf=1,mmsq2(2)+2)
       
 

!      irecnn=0
!      kia=0
!      kkb=0
!      ihitn=0
      
      
!      icdn=0

      if (ibbbsw.eq.1)goto 451
      call sieve0(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2)
      ibbbsw=1
      goto 451
      
      do jz=2,limprm
      if (ipr(jz).lt.10000)goto 405
      iprar(1)=0
      iprar(2)=2
      iprar(3)=ipr(jz)/10000
      iprar(4)=ipr(jz)-iprar(3)*10000
      goto 401
405   iprar(1)=0
      iprar(2)=1
      iprar(3)=ipr(jz)
401   a=a
      do jf=3,n(2)+2
      karr(jf-2)=n(jf)
      end do
      ilen=n(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,irlen
      kbarr(jf+2)=irrr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=irlen
!  minus sieving option below      
      do jf=1,kbarr(2)+2
      itemp1(jf)=kbarr(jf)
      end do
      goto 430
      do jf=1,iprar(2)+2
      karr(jf)=iprar(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      itemp1(jf)=kcarr(jf)
      end do
430   a=a      
      isum=0
      do jf=3,itemp1(2)+2
      isum=isum+itemp1(jf)*10000**(itemp1(2)+2-jf)
      end do
      kr(jz,1)=isum
      end do
      print *,'kr',(kr(jf,1),jf=2,11) 
      print *,'kr2s',(kr(jf,2),jf=2,11)
      print *,'ipr',(ipr(jf),jf=2,11)
      stop
451   a=a      
      print *,'kr',(kr(jf,1),jf=2,11) 
!      print *,'kr2s',(kr(jf,2),jf=2,11)
      print *,'ipr',(ipr(jf),jf=2,11)
      print *,'mr',(mr(jf),jf=1,20)
      print *,'mm',(mm(jf),jf=1,mm(2)+2)
      









!     change parameter      
      do i=2,limprm
      if (kr(i,1).eq.999999)goto 592
      rel=ipr(i)
      gpr(i)=log10(rel)
      goto 594
592   gpr(i)=99.0      
594   end do      
      rel=ipr(limprm)
      lglm=log10(rel)
      
      
      print *,'i',i
      print *,'krs',(kr(jf,1),jf=1,50)
      print *,'krs2',(kr(jf,2),jf=1,50)
      
      print *,'end y stage'
!      open(unit=2,file='gpx2',access='direct',form=&
!      'formatted',recl=410,status='old')
      print *,'okk limprm',limprm
      
      
55555 a=a
      
      

!     start of big loop      
     goto 452 
      

      
      


!     setting up sieving for ax+b      
!      do jz=2,limprm
!      mr(jz)=ipr(jz)
!      end do
452   a=a      


!      icdn=0
     icdarr(1)=0 
     icdarr(2)=0


      
      
      
599   do jz=2,limprm
      
      if (gpr(jz).gt.90.0)goto 6501
      
      if (icdarr(2).gt.0)goto 5991
      markl=kr(jz,1)
      if (markl.gt.10000)goto 610
      goto 606
      
5991  markl=kkm(jz,1)
      if (markl.gt.10000)goto 610
      

      
      
606   krecarr(markl)=krecarr(markl)-gpr(jz)
608   markl=markl+ipr(jz)
      
      if (markl.gt.10000)goto 610
      krecarr(markl)=krecarr(markl)-gpr(jz)
      goto 608


610   kkm(jz,1)=markl-10000
      goto 6501
6101  if (ipr(jz).eq.2)goto 6501
      if (kr(jz,2).eq.999999)goto 6501
      if (icdarr(2).gt.0)goto 6141
      markr=kr(jz,2)
      if (markr.gt.10000)goto 6502
      goto 6142
6141  markr=kkm(jz,2)
612   if (markr.gt.10000)goto 6502
6142  krecarr(markr)=krecarr(markr)-gpr(jz)
      markr=markr+ipr(jz)
      goto 612
6502  kkm(jz,2)=markr-10000
6501  end do
      goto 668




      do jz=2,limprm
!      if ((gpr(jz).gt.90.0).and.(jz.gt.200))goto 210
      if (icdarr(2).gt.0)goto 211
      mbarkl=mr(jz)
      if (mbarkl.gt.10000)goto 216
      goto 213
211   mbarkl=mkm(jz)
      if (mbarkl.gt.10000)goto 216
213   mrecarr(mbarkl)=mrecarr(mbarkl)-gpr(jz)
215   mbarkl=mbarkl+ipr(jz)
      if (mbarkl.gt.10000)goto 216
      mrecarr(mbarkl)=mrecarr(mbarkl)-gpr(jz)
      goto 215
216   mkm(jz)=mbarkl-10000
210   end do
!      print *,'end of second sieve'







!     start of last phase
      
      
668   do kz=1,10000
      
      vimprm=mimprm -4.0
      vimprm3=mimprm
!      if ((krecarr(kz).gt.vimprm).or.(mrecarr(kz).gt.vimprm3))goto 690
      if (krecarr(kz).gt.vimprm)goto 690
!      print *,'krecarr',krecarr(kz)
!     compute accurately      
      do jf=1,icdarr(2)+2
      karr(jf)=icdarr(jf)
      end do
      if (kz.eq.10000)goto 410
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=kz
      goto 411
410   kbarr(1)=0
      kbarr(2)=2
      kbarr(3)=1
      kbarr(4)=0
411   call mpadd(0)
      do jf=1,kcarr(2)+2
      
      karr(jf)=kcarr(jf)
      end do
      do jf=1,imm(2)+2
      kbarr(jf)=imm(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      norma0(jf)=kcarr(jf)
      end do
      do jf=3,norma0(2)+2
      karr(jf-2)=norma0(jf)
      
      end do
      do jf=3,mm(2)+2
      kbarr(jf-2)=mm(jf)
      end do
      
      ilen=norma0(2)
      ilen2=mm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      norma(jf)=kcarr(jf)
      end do
      goto 454





! minus side of sieve here      
      norma0(1)=1
      karr(1)=1
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      norma(jf)=kcarr(jf)
      end do

454   a=a








      rel=zlog(norma(3))+4*norma(2)-4
      rel2=rel-127.0+krecarr(kz)
      rel=zlog(norma0(3))+4*norma0(2)-4
      rel3=rel-127.0 +mrecarr(kz)
!      print *,'rel3',rel3,'rel2',rel2,'lglm',lglm ,'mrec',mrecarr(kz)
      if (rel2.gt.lglm)goto 690
!      if ((rel2.gt.lglm).or.(rel3.gt.lglm))goto 690
      print *,'norma0',(norma0(jf),jf=1,norma0(2)+2),'limp2',limprm2
      kia=kz
      
      
      call sieve(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2,limprm)
      
      
      


690   end do
!      icdn=icdn+1
     do jf=1,icdarr(2)+2
     karr(jf)=icdarr(jf)
     end do
     kbarr(1)=0
     kbarr(2)=2
     kbarr(3)=1
     kbarr(4)=0
      call mpadd(0)
      do jf=1,kcarr(2)+2
      icdarr(jf)=kcarr(jf)
      end do


      do jf=1,10000
      krecarr(jf)=127.0
      mrecarr(jf)=127.0
      end do
      do jhj=1,nfak
!      if (jfakin(jhj).eq.0)goto 599
      end do
      do jhj=1,kfak
!      if (kjfakin(jhj).eq.0)goto 599
      end do
      if (irecnn.eq.400)goto 720
      if ((icdarr(2).eq.3).and.(icdarr(3).eq.2))goto 904
      goto 599
904   do jf=1,mm(2)+2
       karr(jf)=mm(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2
       call mpadd(0)
       
       do jf=1,kcarr(2)+2
       mm(jf)=kcarr(jf)
!       imm(jf)=kcarr(jf)
       end do
       





       goto 900










720   print *,'maximum no of primes=',kkmx
      print *,'irecnn=',irecnn
      
730   print *,'run completed',' no of records=',irecnn,'iconq',iconq

1002  format (i6,i6,i8,i4,20i4,20i4,20i4,20i4)      

      close (unit=1)
      close (unit=7)
      close (unit=2)
      
      
      
      
      close(unit=4)
      
      end
      subroutine sieve0(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      
      
      
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50),jfreq(20),itbase(50)
      ibbsw=0
!     change parameter      
      lprx1=ipr(limprm2)                                      
      
      
      
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
441   a=a





      
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
      jab=0
1010  iab=iab+1
      
      if (iab.eq.limprm2)goto 1500
      if (gpr(iab).gt.90.0)goto 1010
      jab=jab+1
      if (jab.eq.16)goto 400
      if (jab.eq.96)goto 410
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
420   a=a
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
      
      
!      if (icont.eq.0)goto 1300
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
      if ((icont.eq.1).and.(ipqt(1).eq.1))goto 1300
      goto 1166
1176  iabpn(icur)=iabpn(icur)+1
      if ((icont.eq.1).and.(ipqt(1).eq.1))goto 1300
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
1300  if (ibbsw.eq.1)goto 1301
      nfak=icur
      do jf=1,nfak
      jfak(jf)=iabp(jf)
      jfreq(jf)=iabpn(jf)
      end do
      print *,'ijpow',ijpow
      print *,'nfak',nfak
      print *,'jfaks',(jfak(jf),jf=1,nfak)
      print *,'jfreqs',(jfreq(jf),jf=1,nfak)
      goto 1600
      do jf=1,nfak
      jfreq(jf)=jfreq(jf)*ijpow2
      end do
      print *,'ijpow2',ijpow2
1301  kfak=icur
      do jf=1,kfak
      kfakin(jf)=iabp(jf)
      kfreq(jf)=iabpn(jf)
      end do
      do jf=1,nfak
      jfreq(jf)=jfreq(jf)*ijpow2
      iabp(jf)=jfak(jf)
      iabpn(jf)=jfreq(jf)
      end do
      
      print *,'ijpow2',ijpow2,'ijpow',ijpow
      print *,'ibase',(ibase(jf),jf=1,ibase(2)+2)
      


      irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kkb,(litt(i),i=1,litt(2) +2)&
      ,'iconq',iconq
      if (icur.le.kkmx)goto 1302
      kkmx=icur

1302  do i=1,20
      normar(i)=0
      littr(i)=0
      end do
!      do i=1,litz(2) +2
!      littr(i)=litz(i)
!      end do
      icur=nfak
      ihitn=0
!      ihitn=icdn*10000+kia
!      write(4,*)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
!      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      
!     note next instruction      
      if (irecnn.lt.190000)goto 1501
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      goto 1501
      


1500 if (ibbsw.eq.1)goto 1610
     ijpow=ijpow+1
     ilen=jbase(2)
     do jf=3,jbase(2)+2
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
     if (icont.eq.0)goto 1500
!     if (mod(ijpow,2).eq.0)goto 1500
     goto 440
1600 ibbsw=1
     ijpow2=1
     do jf=1,ibase(2)+2
     litt(jf)=ibase(jf)
     itbase(jf)=ibase(jf)
     end do
     goto 441
1610 do jf=3,ibase(2)+2 
     karr(jf-2)=ibase(jf)
     end do
     ilen=ibase(2)
     do jf=3,itbase(2)+2
     kbarr(jf-2)=itbase(jf)
     end do
     ilen2=itbase(2)
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
     ijpow2=ijpow2+1
     if (icont.eq.0)goto 1610
     if (mod(ijpow2,2).eq.0)goto 1610
     do jf=1,ibase(2)+2
     litt(jf)=ibase(jf)
     end do
     goto 441



1501 return






     end
      
      subroutine sieve(kia,kkb,irecnn,kkmx,icdn,iconq,limprm2,limprm)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      
      
      
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50),jabp(50),jabpn(50),jty(50)
      
!     change parameter      
      lprx1=ipr(limprm)                                      
      
      
      
      larp(1)=int(lprx1/10000)
      larp(2)=lprx1 -larp(1) *10000
      
      
      do i=1,50
      ity(i)=0
      iabp(i)=0
      iabpn(i)=0
      end do
      do jf=1,norma(2)+2
      litt(jf)=norma(jf)
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
      print *,'litd',(litd(i),i=1,litd(2) +2),'icdn',icdn,'irecnn',irecnn&
      ,'iconq',iconq
      do i=1,50

      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
      if (iab.eq.limprm)goto 1500
      if (gpr(iab).gt.90.0)goto 1010
       
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
1300  do i=1,norma0(2)+2
      litt(i)=norma0(i)
      litd(i)=norma0(i)
      end do
      print *,'norma0',(norma0(jf),jf=1,norma0(2)+2)
      
      lprx1=ipr(limprm)                                      
      
      
      
      larp(1)=int(lprx1/10000)
      larp(2)=lprx1 -larp(1) *10000
      
      indcur=icur
      litd(1)=0
      iab=1
      do i=1,50
      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1306  iab=iab+1
      iprx1=ipr(iab)
      if (iab.eq.limprm)goto 1500
      if ((gpr(iab).gt.90.0).and.(iab.gt.200))goto 1306
      iarp(1)=int(iprx1/10000)
      iarp(2)=iprx1-iarp(1)*10000
      ibsw=0
      if (iarp(1).gt.larp(1))goto 1500
      if (iarp(1).lt.larp(1))goto 1362
      if (iarp(2).gt.larp(2))goto 1500
1362  if (iarp(1).gt.ipfar(1))goto 1400
      if (iarp(1).lt.ipfar(1))goto 1366
      if (iarp(2).gt.ipfar(2))goto 1400
1366  do i=3,litd(2)+2
      karr(i-2)=litd(i)
      end do
      ilen=litd(2)
      if (iprx1.lt.10000)goto 1368
      kbarr(1)=iarp(1)
      kbarr(2)=iarp(2)
      ilen2=2
      goto 1370
1368  kbarr(1)=iprx1
      ilen2=1
1370  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (icont.eq.0)goto 1386
      if (irlen.gt.0)goto 1306
      do i=1,icont
      litd(i+2)=ipqt(i)
      end do
      litd(2)=icont
1382  if (ibsw.eq.1)goto 1384
      ibsw=1
      icur=icur+1
      if (icur.le.kkmx)goto 1383
      kkmx=icur
1383  if (icur.gt.50)goto 13861
      ity(icur)=2
      iabp(icur)=iprx1
      iabpn(icur)=1
      goto 1366
1384  iabpn(icur)=iabpn(icur)+1
      goto 1366
1386  a=a
      do jf=indcur+1,icur
      iabpn(jf)=iabpn(jf)+iabpn(jf)
      end do


      do ii=1,nfak
      if (jfakin(ii).eq.1)goto 435
      do jj=1,indcur
      if (iabp(jj).eq.jfak(ii))goto 430
      goto 434
430   jfakin(ii)=1
      goto 436
434   end do
435   end do
436   a=a
      do ii=1,kfak
      if (kjfakin(ii).eq.1)goto 445
      do jj=1,indcur
      if (iabp(jj).eq.kfakin(ii))goto 450
      goto 444
450   kjfakin(ii)=1
      
      goto 446
444   end do
445   end do
446   a=a 
      ind1=1
      ind2=1
      ind3=1
31    if (iabp(ind1).lt.kfakin(ind2))goto 30
      if (iabp(ind1).gt.kfakin(ind2))goto 35
      jabp(ind3)=iabp(ind1)
      jabpn(ind3)=iabpn(ind1)+kfreq(ind2)*(3-mod(irecnn+1,3))
      jty(ind3)=1
      ind1=ind1+1
      ind2=ind2+1
      ind3=ind3+1
      if (ind1.gt.indcur)goto 36
      if (ind2.gt.kfak)goto 37
      goto 31
30    jabp(ind3)=iabp(ind1)
      jabpn(ind3)=iabpn(ind1)
      jty(ind3)=1
      ind1=ind1+1
      ind3=ind3+1
      if (ind1.gt.indcur)goto 36
      goto 31
35    jabp(ind3)=kfakin(ind2)
      jabpn(ind3)=kfreq(ind2)*(3-mod(irecnn+1,3))
      jty(ind3)=1
      ind2=ind2+1
      ind3=ind3+1
      if (ind2.gt.kfak)goto 37
      goto 31
36    if (ind2.gt.kfak)goto 38
      do jf=ind2,kfak
      jabp(ind3)=kfakin(jf)
      jabpn(ind3)=kfreq(jf)*(3-mod(irecnn+1,3))
      jty(ind3)=1
      ind3=ind3+1
      end do
      goto 38
37    do jf=ind1,indcur
      jabp(ind3)=iabp(jf)
      jabpn(ind3)=iabpn(jf)
      jty(ind3)=1
      ind3=ind3+1
      end do
38    do jf=indcur+1,icur
      jabp(ind3)=iabp(jf)
      jabpn(ind3)=iabpn(jf)
      jty(ind3)=2
      ind3=ind3+1
      end do
      jcur=ind3-1








      irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kkb,(litt(i),i=1,litt(2) +2)&
      ,'iconq',iconq
      
      if (icur.le.kkmx)goto 1302
      kkmx=icur

1302  do i=1,20
      normar(i)=0
      littr(i)=0
      end do
      do i=1,litz(2) +2
      littr(i)=litz(i)
      end do
      do i=1,norma(2)+2
      normar(i)=norma(i)
      end do
      ihitn=icdn*10000+kia
!      write(4,*)irecnn,kkb,ihitn,jcur,(littr(j1),j1=1,20)
!      write(4,*)irecnn,(jabp(i),jabpn(i),jty(i),i=1,jcur)
      
!     note next instruction      
      if (irecnn.lt.190000)goto 1500
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      goto 1500
      

1400 if (litd(2).ne.iprex(2))goto 1402 
     do i=3,litd(2)+2
     if (litd(i).ne.iprex(i))goto 1402
     end do
     goto 1404
1402 if (litd(2).gt.2)goto 3020
     if (litd(2).eq.1)goto 1366
     if (litd(3).lt.larp(1))goto 1366
     if (litd(3).gt.larp(1))goto 3020
     if (litd(4).le.larp(2))goto 1366
3020 do i=3,litd(2)+2
     mnum(i-2)=litd(i)
     end do
     inlen=litd(2)
     call mpprime(icorp,inlen)
     if (icorp.eq.1)goto 1500
     do i=1,litd(2)+2
     iprex(i)=litd(i)
     end do
1404 ipfar(1)=ipfar(1)+2
     goto 1366
1500 return
     end


     

      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      
      
      
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
      
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(40000),norma0(50)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(40000),norma0(50)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
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
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(40000),norma0(50)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,ibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
       
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
       subroutine bigb5(noyes)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(40000),norma0(50)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
      common nfak,jfak(20),mntes1(50),mntes2(50),jfakin(20),iconarr(50)
      common ijpow,iibase(50),jbase(50),itarg(50),jtarg(50)
      common kfak,kfakin(20),kfreq(20),kjfakin(20)
       dimension ibase(50),ipow(50),ibst(50),nnum(50),ie(50)
       
       ibigsw=0
       
       do jf=3,ncom(2)+2
       karr(jf-2)=ncom(jf)
       end do
       ilen=ncom(2)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       print *,'iprar',(iprar(jf),jf=1,iprar(2)+2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1
       do jf=1,irlen
       ibase(jf+2)=irrr(jf)
       ibst(jf+2)=irrr(jf)
       ipow(jf+2)=irrr(jf)
       end do
       ibst(1)=0
       ibst(2)=irlen
       ibase(1)=0
       ibase(2)=irlen
       ipow(1)=0
       ipow(2)=irlen
       print *,'ipow',(ipow(jf),jf=1,ipow(2)+2)
       
       do jf=1,iprar(2)+2
       karr(jf)=iprar(jf)
       end do
       iconz=iprar(2)+2
       if (mod(iprar(iconz),4).ne.3)goto 2
       kbarr(1)=0 
       kbarr(2)=1
       kbarr(3)=1
       
       call mpadd(0)
15     do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       kbarr(1)=4
       ilen2=1
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       mnum(jf+2)=ipqt(jf)
       end do
       mnum(2)=icont
       mnum(1)=0
       
       goto 3
1      print *,'factor found early=',(iprar(jf),jf=1,iprar(2)+2)
       stop
2      kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=5
       call mpadd(1) 
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       kbarr(1)=8
       ilen2=1
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       print *,'icont',icont
       do jf=1,icont
       mnum(jf+2)=ipqt(jf)
       end do
       mnum(1)=0
       mnum(2)=icont
       
       do jf=3,ibase(2)+2
       karr(jf-2)=ibase(jf)
       end do
       ilen=ibase(2)
       kbarr(1)=4
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       ibase(jf+2)=irrr(jf)
       
       end do
       ibase(1)=0
       ibase(2)=irlen
       do jf=1,ibase(2)+2
       
       ipow(jf)=ibase(jf)
       end do



       
       
       
3      ie(1)=0
       ie(2)=1
       ie(3)=2
       karr(1)=2
       ilen=1
       inlen=mnum(2)
4      do jf=3,ie(2)+2
       kbarr(jf-2)=ie(jf)
       end do
       ilen2=ie(2)
       call mpmul(ilen,ilen2,ilen3)
       if (ilen3.lt.inlen)goto 22
       if (ilen3.gt.inlen)goto 20
       do i=1,ilen3
       if (mnum(i+2).lt.kcarr(i))goto 20
       if (mnum(i+2).gt.kcarr(i))goto 22
       end do
       goto 22
20     ie(2)=ilen3
       do i=1,ilen3
       ie(i+2)=kcarr(i)
       end do
       goto 30
22     do i=1,ilen3
       ie(i+2)=kcarr(i)
       end do
       ie(2)=ilen3
       goto 4
30     do i=1,inlen +2
       nnum(i)=mnum(i)
       end do
       do i=3,ie(2)+2
       karr(i-2)=ie(i)
       end do
       ilen=ie(2)
       kbarr(1)=2
       ilen2=1
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do i=1,icont
       ie(i+2)=ipqt(i)
       end do
       ie(1)=0
       ie(2)=icont
       
       do i=1,mnum(2)+2
       karr(i)=mnum(i)
       end do
       do i=1,ie(2)+2
       kbarr(i)=ie(i)
       end do
       call mpadd(1)
       do i=1,kcarr(2)+2
       nnum(i)=kcarr(i)
       end do
       print *,'nnum1',(nnum(jf),jf=1,nnum(2)+2)
31     if (ie(2).gt.1)goto 32
       if (ie(3).eq.1)goto 100
32     do i=3,ie(2)+2
       karr(i-2)=ie(i)
       end do
       ilen=ie(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       ie(2)=icont
       do i=1,icont
       ie(i+2)=ipqt(i)
       end do
       do jf=3,ipow(2)+2
       karr(jf-2)=ipow(jf)
       kbarr(jf-2)=ipow(jf)
       end do
       ilen=ipow(2)
       
       ilen2=ipow(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq) 
       do jf=1,irlen
       ipow(jf+2)=irrr(jf)
       end do
       ipow(2)=irlen
       ipow(1)=0
       if (nnum(2).lt.ie(2))goto 31
       if (nnum(2).gt.ie(2))goto 42
       do i=3,nnum(2)+2
       if (nnum(i).lt.ie(i))goto 31
       if (nnum(i).gt.ie(i))goto 42
       end do
42     do i=1,nnum(2)+2
       karr(i)=nnum(i)
       end do
       do i=1,ie(2)+2
       kbarr(i)=ie(i)
       end do
       call mpadd(1)
       do i=1,kcarr(2)+2
       nnum(i)=kcarr(i)
       end do
       do i=3,ipow(2)+2
       karr(i-2)=ipow(i)
       end do
       ilen=ipow(2)
       do i=3,ibase(2)+2
       kbarr(i-2)=ibase(i)
       end do
       ilen2=ibase(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       ipow(jf+2)=irrr(jf)
       end do
       ipow(1)=0
       ipow(2)=irlen
       goto 31
100    iconz=iprar(2)+2
       if (mod(iprar(iconz),4).eq.3)goto 110
       if (ibigsw.eq.1)goto 104
       print *,'ipow13',(ipow(jf),jf=1,ipow(2)+2)
       do jf=3,ipow(2)+2
       karr(jf-2)=ipow(jf)
       end do
       ilen=ipow(2)
       kbarr(1)=2
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ibst(2)+2
       kbarr(jf-2)=ibst(jf)
       end do
       ilen2=ibst(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       iansar(jf+2)=irrr(jf)
       end do
       iansar(1)=0
       iansar(2)=irlen
       goto 1041
1039   do jf=1,iprar(2)+2
       karr(jf)=iprar(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=3
       call mpadd(0)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       kbarr(1)=8
       ilen2=1
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       mnum(jf+2)=ipqt(jf)
       end do
       mnum(1)=0
       mnum(2)=icont
       do jf=1,ibst(2)+2
       ibase(jf)=ibst(jf)
       ipow(jf)=ibst(jf)
       end do
       
       print *,'ipowb',(ipow(jf),jf=1,ipow(2)+2)
       print *,'mnumb',(mnum(jf),jf=1,mnum(2)+2)
       goto 3
1041   do jf=3,iansar(2)+2
       karr(jf-2)=iansar(jf)
       kbarr(jf-2)=iansar(jf)
       end do
       ilen=iansar(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf+2)=kcarr(jf)
       end do
       karr(1)=0
       karr(2)=ilen3
       do jf=1,ibst(2)+2
       kbarr(jf)=ibst(jf)
       end do
       call mpadd(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       isgn=kcarr(1)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1141
       ibigsw=1
       goto 1039

       
104    print *,'bans',(ipow(jf),jf=1,ipow(2)+2)
       


110    do jf=1,ipow(2)+2
       iansar(jf)=ipow(jf)
       end do
       print *,'bigans',(iansar(jf),jf=1,iansar(2)+2)
       do jf=3,iansar(2)+2
       karr(jf-2)=iansar(jf)
       kbarr(jf-2)=iansar(jf)
       end do
       ilen=iansar(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf+2)=kcarr(jf)
       end do
       karr(1)=0
       karr(2)=ilen3
       do jf=1,ibst(2)+2
       kbarr(jf)=ibst(jf)
       end do
       call mpadd(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       isgn=kcarr(1)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.ne.0)goto 114 
1141   noyes=0
       
       goto 112
114    noyes=1
       print *,'irlen',irlen
       print *,'irrr',(irrr(jf),jf=1,irlen)











112    return
       end











