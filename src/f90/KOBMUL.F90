     

      program kobmul
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
      dimension idrv5(5),idrv4(5),idrv3(5),idrv2(5),idrv1(5)
      dimension match(2000),iv(3500),n(50),iabp(50),iabpn(50),ity(50)
      dimension numb1(20000),numb2(50),iquot1(20000),igcd(50)
      dimension ires(20000),ipre(20000),icurr(20000),iroot(20000)
      dimension mm1(10),mm2(10),mm3(10),mm4(10),iprod(90),littr(20)
      dimension ibm(10,90),ism(10,90),ibigm(90),mnrem(10,90),isumt(90)
      dimension iai(10,90),ixi(10,90),insum(90),nres(90),iresp(100)
      irecnn=0
      ipoin=0
      n(1)=0
      n(2)=13
      n(3)=6
      n(4)=8300
      n(5)=3888
      n(6)=1025
      n(7)=3257
      n(8)=1578
      n(9)=2089
      n(10)=3818
      n(11)=4804
      n(12)=3522
      n(13)=9365
      n(14)=2252
      n(15)=4763
      mm1(1)=0
      mm1(2)=7
      mm1(3)=2
      mm1(4)=6134
      mm1(5)=3430
      mm1(6)=7769
      mm1(7)=2474
      mm1(8)=7169
      mm1(9)=2876
      do jf=1,9
      mm2(jf)=mm1(jf)
      mm3(jf)=mm1(jf)
      mm4(jf)=mm1(jf)
      end do
      mm2(7)=2369
      mm3(7)=2454
      mm4(7)=2489

      modn=3
      mm=200
      open(unit=1,file='kits',access='sequential')
      open(unit=2,file='llrnel2',access='direct',form=&
      'formatted',recl=1000,status='old')
      read(2,2,rec=2)(iv(jf),jf=1,mm)
      
      do i=1,mm
      read(1, *)irecnn,kkb,ihitn,icur,(littr(jf),jf=1,20),&
      (numb1(jf),jf=1,20)
      read(1, *)irecnn,(iabp(jk),iabpn(jk),ity(jk),jk=1,icur)
      if (iv(i).eq.1)goto 20
      end do
      print *,'file error 1'
      stop
20    i1=i+1  
2     format(3164i5)      
      iroot(1)=0
      iroot(2)=1
      iroot(3)=1
      icc=1
      do i=i1,mm
      read(1, *)irecnn,kia,kib,icur,(littr(20),jf=1,20),&
      (numb2(jf),jf=1,20)
      
      read(1, *)irecnn,(iabp(jk),iabpn(jk),ity(jk),jk=1,icur)
      if(iv(i).eq.0)goto 99
      icc=icc+1
      do jk=1,numb1(2)+2
      kara(jk)=numb1(jk)
      end do
      do jk=1,numb2(2)+2
      karb(jk)=numb2(jk)
      end do
      call mpgcd
      do jk=1,kard(2)+2
      igcd(jk)=kard(jk)
      end do
      do jk=3,numb1(2)+2
      karr(jk-2)=numb1(jk)
      end do 
      ilen=numb1(2)
      ilen2=igcd(2)
      do jk=3,ilen2+2
      kbarr(jk-2)=igcd(jk)
      end do
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jk=1,icont
      iquot1(jk+2)=ipqt(jk)
      end do
      iquot1(1)=0
      iquot1(2)=icont
      do jk=3,numb2(2)+2
      karr(jk-2)=numb2(jk)
      end do
      
      ilen=numb2(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jk=1,icont
      
      kbarr(jk)=ipqt(jk)
      end do
      ilen2=icont
      do jk=3,iquot1(2)+2
      karr(jk-2)=iquot1(jk)
      end do
      ilen=iquot1(2)
      call mpmul(ilen,ilen2,ilen3)
      do jk=1,ilen3
      numb1(jk+2)=kcarr(jk)
      end do
      numb1(1)=0
      numb1(2)=ilen3
      do jk=3,iroot(2)+2
      karr(jk-2)=iroot(jk)
      end do
      ilen=iroot(2)
      do jk=3,igcd(2)+2
      kbarr(jk-2)=igcd(jk)
      end do
      ilen2=igcd(2)
      call mpmul(ilen,ilen2,ilen3)
      do jk=1,ilen3
      iroot(jk+2)=kcarr(jk)
      end do
      iroot(1)=0
      iroot(2)=ilen3
      
      print *,'numb1',(numb1(jf),jf=1,numb1(2)+2)
      print *,'iroot',(iroot(jf),jf=1,iroot(2)+2)
      print *,'i=',i,'icc=',icc
99    end do
999   print *,'numb1',(numb1(jf),jf=1,numb1(2)+2)
      
!     compute rational integer square root
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
      print *,'lenz',ilen,ilen2
      print *,'kars',(karr(jk),jk=1,ilen)
      print *,'kbars',(kbarr(jk),jk=1,ilen2)
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
      print *,'curr',(icurr(jk),jk=1,icurr(2)+2)
      
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
      stop
742   print *,'square ok'
      
      
      do jf=3,icurr(2)+2
      karr(jf-2)=icurr(jf)
      end do
      ilen=icurr(2)
      do jf=3,iroot(2)+2
      kbarr(jf-2)=iroot(jf)
      end do
      ilen2=iroot(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      
      
      print *,'norm=',(karr(jf),jf=1,ilen)
      kbarr(1)=157 
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      print *,'irrr1',irrr(1)
      
      
      
      print *,'icc=',icc
      modn=5
      ixi(1,1)=0
      ixi(1,2)=17
      ixi(1,3)=7
      ixi(1,4)=956
      ixi(1,5)=2654
      ixi(1,6)=6479
      ixi(1,7)=7906
      ixi(1,8)=5101
      ixi(1,9)=8044
      ixi(1,10)=9643
      ixi(1,11)=3072
      ixi(1,12)=1916
      ixi(1,13)=1680
      ixi(1,14)=1902
      ixi(1,15)=3208
      ixi(1,16)=5143
      ixi(1,17)=1618
      ixi(1,18)=5225
      ixi(1,19)=9586
      ixi(2,1)=0
      ixi(2,2)=9
      ixi(2,3)=990
      ixi(2,4)=8780
      ixi(2,5)=9448
      ixi(2,6)=7337
      ixi(2,7)=7405
      ixi(2,8)=2659
      ixi(2,9)=9707
      ixi(2,10)=9101
      ixi(2,11)=405
      ixi(3,1)=0
      ixi(3,2)=9
      ixi(3,3)=1918
      ixi(3,4)=9602
      ixi(3,5)=6417
      ixi(3,6)=7472
      ixi(3,7)=6378
      ixi(3,8)=476
      ixi(3,9)=234
      ixi(3,10)=1085
      ixi(3,11)=5674
      ixi(4,1)=0
      ixi(4,2)=10
      ixi(4,3)=27
      ixi(4,4)=8409
      ixi(4,5)=8370
      ixi(4,6)=9156
      ixi(4,7)=5087
      ixi(4,8)=3349
      ixi(4,9)=2919
      ixi(4,10)=1859
      ixi(4,11)=2123
      ixi(4,12)=5320









      ism(1,1)=0
      ism(1,2)=1
      ism(1,3)=199
      ism(2,1)=0
      ism(2,2)=1
      ism(2,3)=157
      ism(3,1)=0
      ism(3,2)=1
      ism(3,3)=167
      ism(4,1)=0
      ism(4,2)=1
      ism(4,3)=229
      ism(5,1)=0
      ism(5,2)=1
      ism(5,3)=257
      print *,'ixi4',(ixi(4,jf),jf=1,ixi(4,2)+2)
      
      n(1)=0
      n(2)=2
      n(3)=3005
      n(4)=3021
      goto 778
      ixi(1,1)=0
      ixi(1,2)=1
      ixi(1,3)=2
      ixi(2,1)=0
      ixi(2,2)=1
      ixi(2,3)=2
      ixi(3,1)=0
      ixi(3,2)=1
      ixi(3,3)=2
      ism(1,1)=0
      ism(1,2)=1
      ism(1,3)=3
      ism(2,1)=0
      ism(2,2)=1
      ism(2,3)=5
      ism(3,1)=0
      ism(3,2)=1
      ism(3,3)=7
      n(1)=0
      n(2)=1
      n(3)=143
      modn=3
      goto 777

      
778   do i=1,modn
      do jf=3,ism(i,2)+2
      karr(jf-2)=ism(i,jf)
      kbarr(jf-2)=ism(i,jf)
      end do
      ilen=ism(i,2)
      ilen2=ism(i,2)
      do j=1,5
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      kbarr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      ilen2=ilen3
      end do
      do jf=1,ilen3
      ism(i,jf+2)=kcarr(jf)
      end do
      ism(i,1)=0
      ism(i,2)=ilen3
      end do
      print *,'ism157',(ism(1,jf),jf=1,ism(1,2)+2)
      
      do i=2,modn
      do jf=1,ixi(1,2)+2
      ixi(i,jf)=ixi(1,jf)
      end do
      end do








          

      
      
      
777   do i=1,modn
      ibm(i,1)=0
      ibm(i,2)=1
      ibm(i,3)=1
      end do
      do i=1,modn
      do j=1,modn
      if (i.eq.j)goto 200
      do jf= 3,ibm(i,2)+2
      karr(jf-2)=ibm(i,jf)
      end do
      ilen=ibm(i,2)
      do jf=3,ism(j,2)+2
      kbarr(jf-2)=ism(j,jf)
      end do
      ilen2=ism(j,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ism(i,2)+2
      kbarr(jf-2)=ism(i,jf)
      end do
      ilen2=ism(i,2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 202
      do jf=1,irlen
      ibm(i,jf+2)=irrr(jf)
      end do

      
202   ibm(i,1)=0
      ibm(i,2)=irlen
200   end do
      end do
      do i=1,modn
      print *,'ibms',(ibm(i,jf),jf=1,ibm(i,2)+2)
      end do
      

!     compute ai mod mi
      do i=1,modn
      do jf=1,ism(i,2)+2
      kara(jf)=ism(i,jf)
      end do
      do jf=1,ibm(i,2)+2
      karb(jf)=ibm(i,jf)
      end do
      call mpgcd
      print *,'gcd',(kard(jf),jf=1,kard(2)+2)
      do jf=1,karv(2)+2
      iai(i,jf)=karv(jf)
      end do
      end do
      do i=1,modn
      print *,'iais',(iai(i,jf),jf=1,iai(i,2)+2)
      end do
      
!     compute prodm mod n 
      ibigm(1)=0
      ibigm(2)=1
      ibigm(3)=1
      do i=1,modn
      do jf=3,ism(i,2)+2
      karr(jf-2)=ism(i,jf)
      end do
      ilen=ism(i,2)
      do jf=3,ibigm(2)+2
      kbarr(jf-2)=ibigm(jf)
      end do
      ilen2=ibigm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 222
      do jf=1,irlen
      ibigm(jf+2)=irrr(jf)
      end do
      ibigm(1)=0
      ibigm(2)=irlen
      end do
      print *,'ibigm',(ibigm(jf),jf=1,ibigm(2)+2)
      
!     compute ibmi remainder n
      do i=1,modn
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      do jf=1,ism(i,2)+2
      karb(jf)=ism(i,jf)
      end do
      call mpgcd
      if (kard(2).ne.1)goto 220
      if (kard(3).ne.1)goto 220
      do jf=3,karv(2)+2
      kbarr(jf-2)=karv(jf)
      end do
      ilen2=karv(2)
      do jf=3,ibigm(2)+2
      karr(jf-2)=ibigm(jf)
      end do
      ilen=ibigm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 222
      do jf=1,irlen
      mnrem(i,jf+2)=irrr(jf)
      end do
      mnrem(i,1)=0
      mnrem(i,2)=irlen
      
      end do
      do i=1,modn
      print *,'mnrems',(mnrem(i,jf),jf=1,mnrem(i,2)+2)
      end do
      do jz=1,60
      isumt(jz)=0
      end do

!     compute rough numbers      
      isumt(1)=0
      isumt(2)=2
      isumt(3)=5000
      isumt(4)=0
      do i=1,modn
      do jf=3,iai(i,2)+2
      karr(jf-2)=iai(i,jf)
      end do
      ilen=iai(i,2)
      do jf=3,ixi(i,2)+2
      kbarr(jf-2)=ixi(i,jf)
      end do
      ilen2=ixi(i,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
!     note method of increase of number length      
      ilen=ilen3+2
      karr(ilen3+1)=0
      karr(ilen3+2)=0
      print *,'karrs',(karr(jf),jf=1,ilen)

      do jf=3,ism(i,2)+2
      kbarr(jf-2)=ism(i,jf)
      end do
      ilen2=ism(i,2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      karr(jf+2)=ipqt(jf)
      end do
      
      karr(1)=0
      karr(2)=icont
      do jf=1,isumt(2)+2
      kbarr(jf)=isumt(jf)
      end do
      
      call mpadd(0)
      do jf=1,kcarr(2)+2
      isumt(jf)=kcarr(jf)
      end do
      print *,'i',i,'isumt',(isumt(jf),jf=1,isumt(2)+2)
      
      end do
      isumt(2)=isumt(2)-2
      print *,'isumt',(isumt(jf),jf=1,isumt(2)+2)
      
!     final part
      insum(1)=0
      insum(2)=1
      insum(3)=0
      
      do i=1,modn
      do jf=3,iai(i,2)+2
      karr(jf-2)=iai(i,jf)
      end do 
      ilen=iai(i,2)
      do jf=3,mnrem(i,2)+2
      kbarr(jf-2)=mnrem(i,jf)
      end do
      ilen2=mnrem(i,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ixi(i,2)+2
      kbarr(jf-2)=ixi(i,jf)
      end do
      ilen2=ixi(i,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
      do jf=1,insum(2)+2
      kbarr(jf)=insum(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      insum(jf)=kcarr(jf)
      end do
      end do
      print *,'insum',(insum(jf),jf=1,insum(2)+2)
      
      do jf=3,isumt(2)+2
      karr(jf-2)=isumt(jf)
      end do
      ilen=isumt(2)
      do jf=3,ibigm(2)+2
      kbarr(jf-2)=ibigm(jf)
      end do
      ilen2=ibigm(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(2)=ilen3
      kbarr(1)=0
      do jf=1,insum(2)+2
      karr(jf)=insum(jf)
      end do
      call mpadd(1)
      if (kcarr(1).eq.0)goto 230
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      
230   do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 222
      do jf=1,irlen
      nres(jf+2)=irrr(jf)
      end do
210   nres(1)=0
      nres(2)=irlen
      goto 224
220   print *,'prime divisor=',(kard(jf),jf=1,kard(2)+2)
      stop
222   print *,'prime divisor=',(karr(jf),jf=1,karr(2)+2)
      stop
224   print *,'grand rem=',(nres(jf),jf=1,nres(2)+2)
      
      iresp(1)=0
      iresp(2)=4
      iresp(3)=60
      iresp(4)=2259
      iresp(5)=8966
      iresp(6)=7495
      do jf=1,iresp(2)+2 
      karr(jf)=iresp(jf)
      end do
      do jf=1,nres(2)+2
      kbarr(jf)=nres(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kara(jf)=kcarr(jf)
      end do
      do jf=1,n(2)+2
      karb(jf)=n(jf)
      end do
      call mpgcd
      print *,'gcd=',(kard(jf),jf=1,kard(2)+2)
      do jf=1,ixi(1,2)+2
      karr(jf)=ixi(1,jf)
      end do
      do jf=1,iresp(2)+2
      kbarr(jf)=iresp(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kara(jf)=kcarr(jf)
      end do
      do jf=1,n(2)+2
      karb(jf)=n(jf)
      end do
      call mpgcd
      print *,'gcd=',(kard(jf),jf=1,kard(2)+2)
      do jf=3,ixi(1,2)+2
      karr(jf-2)=ixi(1,jf)
      end do
      ilen=ixi(1,2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      print *,'rems',(irrr(jf),jf=1,irlen)
      do i=1,modn
      do jf=3,ibm(i,2)+2
      karr(jf-2)=ibm(i,jf)
      end do
      ilen=ibm(i,2)
      do jf=3,iai(i,2)+2
      kbarr(jf-2)=iai(i,jf)
      end do
      ilen2=iai(i,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,ism(i,2)+2
      kbarr(jf-2)=ism(i,jf)
      end do
      ilen2=ism(i,2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      print *,'irrr',(irrr(jf),jf=1,irlen)
      end do
      do i=1,modn
      print *,'i',i,(ixi(i,jf),jf=1,ixi(i,2)+2)
      end do








700   close(unit=2)
      close(unit=1)
      end
      


   
      
      
      subroutine addstar(in,out)

      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
      
      do i=1,30000
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
      dimension kdum(30000),isub(30000)
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      

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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      dimension karu(30000),karv1(30000),karv3(30000),karqq(30000)
      dimension kart3(30000),kart1(30000)
      
      
      
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
      common karr(30000),kbarr(30000),kcarr(30000),ipqt(30000)
      common irrr(30000)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(30000),karb(30000),kard(30000),karp(30000),karv(30000)
      
      common iarq(2)
      
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





      
