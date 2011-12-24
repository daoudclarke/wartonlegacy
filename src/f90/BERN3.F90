      
      program bern3 
!     first in generalised GNFS suite    
!     smaller sieving interval
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
      dimension ka3(50),ka2(50),ka1(50),kc(50),ks1(50),ks2(50)
      dimension ks3(50),ks4(50),ks12(50),ks13(50),ks14(50),krr(50)
      dimension kq(50),kpf(50),kqf(50),krf(50),ka(50),kb(50)
      dimension kp(50),inv2(50),inv3(50),inv4(50),inv8(50)
      dimension inv16(50),inv27(50),inv256(50),kcub(6) 
      dimension itempr(100),kicr(50),ians2(2,50)
      dimension khcr(50),icuba(6,2,50),iden(50),invcom(50),jnvcom(2,50)
      dimension iroota(2,50),kcapa(2,50),kdcapa(2,2,50),kcapb(2,2,50)
      dimension kcofc(2,2,50),kdis(2,2,50),ksol(2,2,2,50),idans(2,2,2,50)
      dimension isol(2,2,2,50),jinj(2)
      jinj(1)=0
      jinj(2)=0
      idegind=4
      goto 777
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=107
      ncom(2,1)=0
      ncom(2,2)=1
      ncom(2,3)=1

      ip(1)=0
      ip(2)=1
      ip(3)=367
      ipn(1)=0
      ipn(2)=1
      ipn(3)=3
      khsq(1)=0
      khsq(2)=1
      khsq(3)=197
!      call cornsq
!      stop
!      call cub5
!      stop
      iaas(1,1)=0
      iaas(1,2)=1
      iaas(1,3)=137
      iaas(2,1)=0
      iaas(2,2)=1
      iaas(2,3)=206
      call sub516
      print *,'icprod1',(icprod(1,jf),jf=1,icprod(1,2)+2)
      print *,'icprod2',(icprod(2,jf),jf=1,icprod(2,2)+2)
      stop
      
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=63
      ncom(2,1)=0
      ncom(2,2)=0
      ip(1)=0
      ip(2)=1
      ip(3)=131
!      call bwq5
!      stop
      idegind=4
      goto 777
      ip(1)=0
      ip(2)=1
      ip(3)=107
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=2
      ncom(2,1)=0
      ncom(2,2)=1
      ncom(2,3)=11
      call cub5
      stop
      
      goto 2001
      iaas(1,1)=0
      iaas(1,2)=1
      iaas(1,3)=4
      iaas(2,1)=0
      iaas(2,2)=0
      ipn(1)=0
      ipn(2)=1
      ipn(3)=4
      ip(1)=0
      ip(2)=1
      ip(3)=1009
      
      call sub516
      print *,'bprod1 ',(ibprod(1,jf),jf=1,ibprod(1,2)+2)
      print *,'bprod2',(ibprod(2,jf),jf=1,ibprod(2,2)+2)
      
      stop
2001  ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=7
      ncom(2,1)=0
      ncom(2,2)=1
      ncom(2,3)=5
      ip(1)=0
      ip(2)=2
      ip(3)=1000
      ip(4)=19
      goto 2000
      iacn(1,1)=0
      iacn(1,2)=1
      iacn(1,3)=563
      iacn(2,1)=0
      iacn(2,2)=1
      iacn(2,3)=494
      ibcn(1,1)=0
      ibcn(1,2)=1
      ibcn(1,3)=605
      ibcn(2,1)=0
      ibcn(2,2)=1
      ibcn(2,3)=600
      call sub1100
      stop
2000  call bwq5
      stop
      
      
      goto 777
      
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=17
      ipz=1009
      print *,'length radix 10000?'
      read *,ncom(1,2)
      print *,'number?'
      read *,(ncom(1,jf),jf=3,ncom(1,2)+2)
      print *,'modulus?'
      read *,ipz
      call cub5
!      print *,'ians',ians,'icubrn',icubrn,'sols',icubr1,icubr2,icubr3
      stop
777   ipz=131

      ip(1)=0
      ip(2)=1
      ip(3)=ipz
      iprar(1)=0
      iprar(2)=1
      iprar(3)=ipz
      
      ka3(1)=0
      ka3(2)=1
      ka3(3)=3
      ka2(1)=0
      ka2(2)=1
      ka2(3)=3
      ka1(1)=0
      ka1(2)=1
      ka1(3)=4
      kc(1)=0
      kc(2)=1
      kc(3)=20
      print *,'enter length radix 10000'
      read *,len
      
      
      print *,'enter coefficient 1  in standard form'
      read *,(ka3(jk),jk=1,len+2)
      print *,'enter coefficient 2 in standard form'
      read *,(ka2(jk),jk=1,len+2)
      print *,'enter coefficient 3 in standard form'
      read *,(ka1(jk),jk=1,len+2)
      print *,'enter constant coefficient in standard form'
      read *,(kc(jk),jk=1,len+2)
      print *,'enter prime modulus'
      read *,(ip(jk),jk=1,len+2)
      do jf=1,ip(2)+2
      iprar(jf)=ip(jf)
      end do


1     do jf=1,ip(2)+2
      kara(jf)=ip(jf)
      end do


2     karb(1)=0
      karb(2)=1
      karb(3)=2
      call mpgcd
      do jf=1,karv(2)+2
      inv2(jf)=karv(jf)
      end do
      karb(3)=3
      call mpgcd
      do jf=1,karv(2)+2
      inv3(jf)=karv(jf)
      end do
      karb(3)=4
      call mpgcd
      do jf=1,karv(2)+2
      inv4(jf)=karv(jf)
      end do
      karb(3)=8
      call mpgcd
      do jf=1,karv(2)+2
      inv8(jf)=karv(jf)
      end do
      karb(3)=16
      call mpgcd
      do jf=1,karv(2)+2
      inv16(jf)=karv(jf)
      end do
      karb(3)=27
      call mpgcd
      do jf=1,karv(2)+2
      inv27(jf)=karv(jf)
      end do
      karb(3)=256
      call mpgcd
      do jf=1,karv(2)+2
      inv256(jf)=karv(jf)
      end do
      print *,'inverses 4 8 27 256',inv4(3),inv8(3),inv27(3),inv256(3)
      print *,'inv16',inv16(3),'inv3',inv3(3),'inv2',inv2(3)
      
      if (ka3(2).eq.0)goto 3
      do jf=3,ka3(2)+2
      karr(jf-2)=ka3(jf)
      end do
      ilen=ka3(2)
      goto 4
3     ks1(1)=0
      ks1(2)=0
      goto 7
4     do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 3
      print *,'irrr',(irrr(jf),jf=1,irlen)
      
      do jf=1,irlen
      karr(jf+2)=irrr(jf)
      end do
      karr(2)=irlen
      karr(1)=mod(ka3(1)+1,2)
      if (karr(1).eq.0)goto 6
      do jf=1,iprar(2)+2
      kbarr(jf)=iprar(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ks1(jf)=kcarr(jf)
      end do
      goto 7
6     do jf=1,karr(2)+2
      ks1(jf)=karr(jf)
      end do
7     print *,'ks1',(ks1(jf),jf=1,ks1(2)+2)
      
      if (ka2(2).eq.0)goto 13
      do jf=3,ka2(2)+2
      karr(jf-2)=ka2(jf)
      end do
      ilen=ka2(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      goto 14
13    ks2(1)=0
      ks2(2)=0
      goto 17
14    call mpdiv(ilen,ilen2,irlen,icont,iswq) 
      if (irlen.eq.0)goto 13
      do jf=1,irlen
      ks2(jf+2)=irrr(jf)
      end do
      ks2(2)=irlen
      ks2(1)=ka2(1)
17    if (ka1(2).eq.0)goto 23      
      do jf=3,ka1(2)+2
      karr(jf-2)=ka1(jf)
      end do
      ilen=ka1(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      goto 24
23    ks3(1)=0
      ks3(2)=0
      goto 27
24    call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 23
      do jf=1,irlen
      karr(jf+2)=irrr(jf)
      end do
      karr(2)=irlen
      karr(1)=mod(ka1(1)+1,2)
      if (karr(1).eq.0)goto 26
      do jf=1,iprar(2)+2
      kbarr(jf)=iprar(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ks3(jf)=kcarr(jf)
      end do
      goto 27
26    do jf=1,karr(2)+2
      ks3(jf)=karr(jf)
      end do
27    if (kc(2).eq.0)goto 33
      do jf=3,kc(2)+2
      karr(jf-2)=kc(jf)
      end do
      ilen=kc(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      goto 34
33    ks4(1)=0
      ks4(2)=0
      goto 37
34    call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 33
      do jf=1,irlen
      ks4(jf+2)=irrr(jf)
      end do
      ks4(2)=irlen
      ks4(1)=kc(1)
37    if (ks1(2).eq.0)goto 43
      do jf=3,ks1(2)+2
      karr(jf-2)=ks1(jf)
      kbarr(jf-2)=ks1(jf)
      end do
      ilen=ks1(2)
      ilen2=ilen
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
      if (irlen.eq.0)goto 43
      do jf=1,irlen
      ks12(jf+2)=irrr(jf)
      end do
      ks12(1)=0
      ks12(2)=irlen
      goto 47
43    ks12(1)=0
      ks12(2)=0
47    if (ks12(2).eq.0)goto 53
      do jf=1,iprar(2)+2
      karr(jf)=iprar(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=3
      call mpadd(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,ks12(2)+2
      kbarr(jf-2)=ks12(jf)
      end do
      ilen2=ks12(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,inv8(2)+2
      kbarr(jf-2)=inv8(jf)
      end do
      ilen2=inv8(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
54    do jf=1,ks2(2)+2
      kbarr(jf)=ks2(jf)
      end do
      call mpadd(0)
      if (kcarr(2).eq.0)goto 55
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 55
      do jf=1,irlen
      kp(jf+2)=irrr(jf)
      end do
      kp(2)=irlen
      kp(1)=0
      goto 57
53    karr(1)=0
      karr(2)=0
      goto 54
55    kp(1)=0
      kp(2)=0
!     single precision from here onward
57    do jf=1,ks1(2)+2
      marr(jf)=ks1(jf)
      end do
      
      do jf=1,ks12(2)+2
      mbarr(jf)=ks12(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      ks13(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ks13(2)+2
      mbarr(jf)=ks13(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      
      do jf=1,inv8(2)+2
      mbarr(jf)=inv8(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kq(jf)=mcarr(jf)
      end do
      do jf=1,inv2(2)+2
      marr(jf)=inv2(jf)
      end do
      do jf=1,ks1(2)+2
      mbarr(jf)=ks1(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ks2(2)+2
      mbarr(jf)=ks2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
       
      do jf=1,kq(2)+2
      kbarr(jf)=kq(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,ks3(2)+2
      kbarr(jf)=ks3(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kq(jf)=mcarr(jf)
      end do

      do jf=1,ks1(2)+2
      marr(jf)=ks1(jf)
      end do
      do jf=1,ks13(2)+2
      mbarr(jf)=ks13(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      ks14(jf)=mcarr(jf)
      end do

      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=3
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,inv256(2)+2
      mbarr(jf)=inv256(jf)
      end do
      call menmul
      do jf=1,mcarr(jf)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ks14(2)+2
      mbarr(jf)=ks14(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krr(jf)=mcarr(jf)
      end do
      do jf=1,inv16(2)+2
      marr(jf)=inv16(jf)
      end do
      do jf=1,ks12(2)+2
      mbarr(jf)=ks12(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ks2(2)+2
      mbarr(jf)=ks2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,krr(2)+2
      kbarr(jf)=krr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krr(jf)=mcarr(jf)
      end do
      do jf=1,inv4(2)+2
      marr(jf)=inv4(jf)
      end do
      do jf=1,ks1(2)+2
      mbarr(jf)=ks1(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ks3(2)+2
      mbarr(jf)=ks3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,krr(2)+2
      karr(jf)=krr(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,ks4(2)+2
      kbarr(jf)=ks4(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krr(jf)=mcarr(jf)
      end do
      print *,'krr',(krr(jf),jf=1,krr(2)+2)
      print *,'kq',(kq(jf),jf=1,kq(2)+2)
      print *,'kp',(kp(jf),jf=1,kp(2)+2)
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do jf=1,kp(2)+2
      kbarr(jf)=kp(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kpf(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=4
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,krr(2)+2
      mbarr(jf)=krr(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kqf(jf)=mcarr(jf)
      end do
      do jf=1,kp(2)+2
      marr(jf)=kp(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=4
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,krr(2)+2
      mbarr(jf)=krr(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krf(jf)=mcarr(jf)
      end do
      do jf=1,kq(2)+2
      marr(jf)=kq(jf)
      mbarr(jf)=kq(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,krf(2)+2
      karr(jf)=krf(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krf(jf)=mcarr(jf)
      end do
      if (idegind.eq.4)goto 300
! to use this routine for solving cubic equations set 2nd coeff to kpf      
! 3rd coeff to kqf , const to krf AND minus 2nd to kp
      kp(1)=0
      kp(2)=1
      kp(3)=129
      kpf(1)=0
      kpf(2)=1
      kpf(3)=2
      kqf(1)=0
      kqf(2)=1
      kqf(3)=3
      krf(1)=0
      krf(2)=1
      krf(3)=77
300   do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,inv3(2)+2
      mbarr(jf)=inv3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      print *,'firka',(marr(jf),jf=1,marr(2)+2)
      
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      print *,'secka',(marr(jf),jf=1,marr(2)+2)
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,kqf(2)+2
      kbarr(jf)=kqf(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      ka(jf)=mcarr(jf)
      end do
      print *,'ka',(ka(jf),jf=1,ka(2)+2),'kpf',(kpf(jf),jf=1,kpf(2)+2)
      print *,'krf',(krf(jf),jf=1,krf(2)+2)
      print *,'kq',(kq(jf),jf=1,kq(2)+2)
      print *,'kqf',(kqf(jf),jf=1,kqf(2)+2)
      
      do jf=1,inv27(2)+2
      marr(jf)=inv27(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kb(jf)=mcarr(jf)
      end do
      print *,'firkb',(kb(jf),jf=1,kb(2)+2)
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,inv3(2)+2
      mbarr(jf)=inv3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kpf(2)+2
      mbarr(jf)=kpf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kqf(2)+2
      mbarr(jf)=kqf(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do

      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,kb(2)+2
      kbarr(jf)=kb(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,krf(2)+2
      kbarr(jf)=krf(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kb(jf)=mcarr(jf)
      end do
      print *,'kb',(kb(jf),jf=1,kb(2)+2)
      
     

      do jf=1,kb(2)+2
      marr(jf)=kb(jf)
      mbarr(jf)=kb(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      khsq(jf)=mcarr(jf)
      end do
      do jf=1,ka(2)+2
      marr(jf)=ka(jf)
      mbarr(jf)=ka(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ka(2)+2
      mbarr(jf)=ka(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,inv27(2)+2
      mbarr(jf)=inv27(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      itempr(jf)=mcarr(jf)
      end do
      do jf=1,khsq(2)+2
      marr(jf)=khsq(jf)
      end do
      do jf=1,inv4(2)+2
      mbarr(jf)=inv4(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,itempr(2)+2
      kbarr(jf)=itempr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      khsq(jf)=mcarr(jf)
      end do
      print *,'khsq',(khsq(jf),jf=1,khsq(2)+2)
!      if (khsq(2).eq.0)goto 200
      if (khsq(1).eq.0)goto 60
      do jf=1,khsq(2)+2
      karr(jf)=khsq(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      khsq(jf)=kcarr(jf)
      end do
      goto 60
60    a=a      
      do jf=1,kb(2)+2
      marr(jf)=kb(jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      khcr(jf)=mcarr(jf)
      end do
      khcr(1)=mod(khcr(1)+1,2)
      if (khcr(1).eq.0)goto 100
      do jf=1,khcr(2)+2
      karr(jf)=khcr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      khcr(jf)=kcarr(jf)
      end do
100   do jf=1,khcr(2)+2
      ncom(1,jf)=khcr(jf)
      end do
      ncom(2,1)=0
      ncom(2,2)=1
      ncom(2,3)=1
      print *,'khsq',(khsq(jf),jf=1,khsq(2)+2)
      print *,'khcr',(khcr(jf),jf=1,khcr(2)+2)
      
      call cub5
      
      ncub=ncubr
      if (ncubr.eq.0)goto 202
      do jbig=1,ncub
      do ibig=1,2
      do jf=1,icubar(jbig,ibig,2)+2
      icuba(jbig,ibig,jf)=icubar(jbig,ibig,jf)
      end do
      end do
      end do
202   a=a
      do jf=1,khcr(2)+2
      ncom(1,jf)=khcr(jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      ncom(2,jf)=kcarr(jf)
      end do
      call cub5
      
      if (ncubr.eq.0)goto 110
      goto 111
110   print *,'ncubr',ncubr,'ncub',ncub
      if (ncub+ncubr.eq.0)goto 999
      ncub=ncub+ncubr
      kkbig=1
      print *,'ncub',ncub
      goto 160
111   do jbig=1,ncubr
      do ibig=1,2
      do jf=1,icubar(jbig,ibig,2)+2
      icuba(ncub+jbig,ibig,jf)=icubar(jbig,ibig,jf)
      end do
      end do
      end do
      ncub=ncub+ncubr
      print *,'ncub',ncub
      kkbig=1
160   a=a
      
!     use of kcapa as subroutine      
      do ibig=1,2
      do jf=1,icuba(kkbig,ibig,2)+2
      kcapa(ibig,jf)=icuba(kkbig,ibig,jf)
      end do
      end do
      inxyz=1
      print *,'kcapa1',(kcapa(1,jf),jf=1,kcapa(1,2)+2)
      print *,'kcapa2',(kcapa(2,jf),jf=1,kcapa(2,2)+2)
      print *,'ka',(ka(jf),jf=1,ka(2)+2)
      print *,'kp',(kp(jf),jf=1,kp(2)+2)
      print *,'inv3',(inv3(jf),jf=1,inv3(2)+2),'kkbig',kkbig
      
      




120   do jf=1,kcapa(1,2)+2
      marr(jf)=kcapa(1,jf)
      mbarr(jf)=marr(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      iden(jf)=mcarr(jf)
      end do
      do jf=1,kcapa(2,2)+2
      marr(jf)=kcapa(2,jf)
      mbarr(jf)=marr(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,khsq(2)+2
      mbarr(jf)=khsq(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,iden(2)+2
      karr(jf)=iden(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      iden(jf)=mcarr(jf)
      end do
      if (iden(1).eq.0)goto 101
      do jf=1,iden(2)+2
      karr(jf)=iden(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iden(jf)=kcarr(jf)
      end do
101   a=a
      do jf=1,ip(2)+2
      karp(jf)=ip(jf)
      end do
      
      
      do jf=1,iden(2)+2
      karb(jf)=iden(jf)
      end do
      call mpgcd
      do jf=1,karv(2)+2
      invcom(jf)=karv(jf)
      end do
      do jf=1,invcom(2)+2
      marr(jf)=invcom(jf)
      end do
      do jf=1,kcapa(1,2)+2
      mbarr(jf)=kcapa(1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      jnvcom(1,jf)=mcarr(jf)
      end do
      do jf=1,invcom(2)+2
      marr(jf)=invcom(jf)
      end do
      do jf=1,kcapa(2,2)+2
      mbarr(jf)=kcapa(2,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      mcarr(1)=mod(mcarr(1)+1,2)
      if (mcarr(1).eq.0)goto 1021
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      jnvcom(2,jf)=kcarr(jf)
      end do
      goto 104
1021  do jf=1,mcarr(2)+2
      jnvcom(2,jf)=mcarr(jf)
      end do
104   a=a
      print *,'iden',(iden(jf),jf=1,iden(2)+2)
      print *,'invcom',(invcom(jf),jf=1,invcom(2)+2)
      print *,'jnvcom1',(jnvcom(1,jf),jf=1,jnvcom(1,2)+2)
      print *,'jnvcom2',(jnvcom(2,jf),jf=1,jnvcom(2,2)+2)
      




      
      
      
      if (inxyz.eq.1)goto 121
      if (inxyz.eq.2)goto 122
      
121   do ibig=1,2
      do jf=1,jnvcom(ibig,2)+2
      marr(jf)=jnvcom(ibig,jf)
      end do
      do jf=1,inv3(2)+2
      mbarr(jf)=inv3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ka(2)+2
      mbarr(jf)=ka(jf)
      end do
      call menmul
      print *,'prod',(mcarr(jk),jk=1,mcarr(2)+2),'ibig',ibig
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do

      do jf=1,icuba(kkbig,ibig,2)+2
      karr(jf)=icuba(kkbig,ibig,jf)
      end do
      call mpadd(1)
      print *,'sum',(kcarr(jk),jk=1,kcarr(2)+2),'ibig',ibig
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      print *,'quo',(mcarr(jk),jk=1,mcarr(2)+2),'ibig',ibig
      do jf=1,mcarr(2)+2
      iroota(ibig,jf)=mcarr(jf)
      end do
      print *,'ibig',ibig,'irt',(iroota(ibig,jk),jk=1,iroota(ibig,2)+2)
      end do
      
      do jf=1,kp(2)+2
      marr(jf)=kp(jf)
      end do
      do jf=1,inv3(2)+2
      mbarr(jf)=inv3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,iroota(1,2)+2
      kbarr(jf)=iroota(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      
      
      do jf=1,mcarr(2)+2
      iroota(1,jf)=mcarr(jf)
      end do
      if (iroota(2,1).eq.0)goto 103
      do jf=1,iroota(2,2)+2
      karr(jf)=iroota(2,jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iroota(2,jf)=kcarr(jf)
      end do
103   a=a
!     answers to cubic equation are in iroota       
      print *,'iroota',(iroota(1,jf),jf=1,iroota(1,2)+2)
      print *,'irootacom',(iroota(2,jf),jf=1,iroota(2,2)+2)
!      print *,'kkbig',kkbig,'kicr',(kicr(jf),jf=1,kicr(2)+2)
      print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
      print *,'icuba',(icuba(kkbig,1,jf),jf=1,icuba(kkbig,1,2)+2)
      print *,'icubac',(icuba(kkbig,2,jf),jf=1,icuba(kkbig,2,2)+2)
      print *,'khsq',(khsq(jf),jf=1,khsq(2)+2)
      print *,'khcr',(khcr(jf),jf=1,khcr(2)+2)
      
      
2012  a=a
      if (idegind.eq.4)goto 310
      if (kkbig.ge.ncub)goto 312
      if (iroota(2,2).eq.0)goto 312
      kkbig=kkbig+1
      goto 160
312   goto 999
310   if (iroota(2,2).eq.0)goto 314
      kkbig=kkbig+1
      if (kkbig.ge.ncub)goto 312
      goto 160


314   do jf=1,iroota(1,2)+2
      karr(jf)=iroota(1,jf)
      end do
      do jf=1,kp(2)+2
      kbarr(jf)=kp(jf)
      end do
      call mpadd(1)
      if (kcarr(1).eq.0)goto 102
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
102   do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2 
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      ncom(1,jf)=mcarr(jf)
      kard(jf)=mcarr(jf)
      end do
      do jf=1,iroota(2,2)+2
      ncom(2,jf)=iroota(2,jf)
      end do
      print *,'ncom1',(ncom(1,jf),jf=1,ncom(1,2)+2)
      print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
      print *,'inv3',(inv3(jf),jf=1,inv3(2)+2)
      print *,'ka',(ka(jf),jf=1,ka(2)+2)
      print *,'kp',(kp(jf),jf=1,kp(2)+2)
      print *,'icuba11',(icuba(kkbig,1,jf),jf=1,icuba(1,1,2)+2)
      print *,'icuba12',(icuba(kkbig,2,jf),jf=1,icuba(1,2,2)+2)
      print *,'invcom',(invcom(jf),jf=1,invcom(2)+2)
      print *,'jnvcom1',(jnvcom(1,jf),jf=1,jnvcom(1,2)+2)
      print *,'jnvcom2',(jnvcom(2,jf),jf=1,jnvcom(2,2)+2)
      
      do jf=1,ip(2)+2
      karp(jf)=ip(jf)
      end do
      call mpkron(k)
      if (k.eq.-1)goto 998

      
      call cornsq
      
      if (nsq.eq.0)goto 998
      print *,'ncub',ncub,'current end kkbig',kkbig
      print *,'ncom1',(ncom(1,jf),jf=1,ncom(1,2)+2)
      print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
      print *,'krr',(krr(jf),jf=1,krr(2)+2)
      print *,'isqr',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2)
      print *,'isqrc',(isqurar(1,2,jf),jf=1,isqurar(1,2,2)+2)
      
      do ibig=1,2
      do jf=1,isqurar(1,ibig,2)+2
      kdcapa(1,ibig,jf)=isqurar(1,ibig,jf)
      end do
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do ibig=1,2
      if (kdcapa(1,ibig,2).eq.0)goto 113
      do jf=1,kdcapa(1,ibig,2)+2
      kbarr(jf)=kdcapa(1,ibig,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kdcapa(2,ibig,jf)=kcarr(jf)
      end do
      goto 114
113   kdcapa(2,ibig,1)=0
      kdcapa(2,ibig,2)=0
114   end do
      
      do ibig=1,2
      do jf=1,kdcapa(1,ibig,2)+2
      kcapa(ibig,jf)=kdcapa(1,ibig,jf)
      end do
      end do
      iibig=1
      inxyz=2
      goto 120
122   print *,'jnvcom1',(jnvcom(1,jf),jf=1,jnvcom(1,2)+2)
      print *,'jnvcom2',(jnvcom(2,jf),jf=1,jnvcom(2,2)+2)
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      do ibig=1,2
      do jf=1,jnvcom(ibig,2)+2 
      marr(jf)=jnvcom(ibig,jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      jnvcom(ibig,jf)=mcarr(jf)
      end do
      end do
      
      
      nnsq=0
141   do ibig=1,2
      if (iibig.eq.2)goto 142
      do jf=1,jnvcom(ibig,2)+2
      marr(jf)=jnvcom(ibig,jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      
      call mendiv
      print *,'ok11'
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,kq(2)+2
      mbarr(jf)=kq(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      print *,'ok12'
      do jf=1,mcarr(2)+2
      kcapb(iibig,ibig,jf)=mcarr(jf)
      end do
      print *,'kcapb',(kcapb(iibig,ibig,jk),jk=1,kcapb(iibig,ibig,2)+2)
142   do jf=1,iroota(ibig,2)+2
      marr(jf)=iroota(ibig,jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,kcapb(iibig,ibig,2)+2
      kbarr(jf)=kcapb(iibig,ibig,jf)
      end do
      call mpadd(1)
      if (kcarr(1).eq.0)goto 124
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
124   do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kcofc(iibig,ibig,jf)=mcarr(jf)
      end do
      print *,ibig,'kcofc',(kcofc(iibig,ibig,jk),jk=1,kcofc(iibig,ibig,2)+2)
      end do
      ipn(1)=0
      ipn(2)=1
      ipn(3)=2
      do ibig=1,2
      do jf=1,kdcapa(iibig,ibig,2)+2
      iaas(ibig,jf)=kdcapa(iibig,ibig,jf)
      end do
      end do
      call sub516
!     keep square in icprod
      do ibig=1,2
      do jf=1,kcofc(iibig,ibig,2)+2
      marr(jf)=kcofc(iibig,ibig,jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=4
      call menmul
      
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,icprod(ibig,2)+2
      karr(jf)=icprod(ibig,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kcarr(jf)=mcarr(jf)
      end do

      if (kcarr(1).eq.0)goto 126
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
126   do jf=1,kcarr(2)+2 
      kdis(iibig,ibig,jf)=kcarr(jf)
      end do
      end do
      print *,'kdis1',(kdis(iibig,1,jk),jk=1,kdis(iibig,1,2)+2)
      print *,'kdisc',(kdis(iibig,2,jk),jk=1,kdis(iibig,2,2)+2)
      print *,'kdcapa1',(kdcapa(1,1,jk),jk=1,kdcapa(1,1,2)+2)
      print *,'kdcapa1c',(kdcapa(1,2,jk),jk=1,kdcapa(1,2,2)+2)
      print *,'kdcapa2',(kdcapa(2,1,jk),jk=1,kdcapa(2,1,2)+2)
      print *,'kdcapa2c',(kdcapa(2,2,jk),jk=1,kdcapa(2,2,2)+2)
      print *,'jnvcom',(jnvcom(1,jk),jk=1,jnvcom(1,2)+2)
      print *,'jnvcomc',(jnvcom(2,jk),jk=1,jnvcom(2,2)+2),'iibig',iibig
      
      
      
      
      jjbig=1
      do ibig=1,2
      do jf=1,kdis(iibig,ibig,2)+2
      ncom(ibig,jf)=kdis(iibig,ibig,jf)
      end do
      print *,'ncom',(ncom(ibig,jf),jf=1,ncom(ibig,2)+2),'kkbig',kkbig
      end do
      

      
!      if (iibig.eq.1)goto 144
      do jf=1,ncom(1,2)+2
      kard(jf)=ncom(1,jf)
      end do
      call mpkron(k)
      
      if (k.eq.-1)goto 144
      
      call cornsq
      
      nnsq=nnsq+nsq
      if (nsq.eq.0)goto 144
      jinj(iibig)=1
140   do ibig=1,2
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do jf=1,kdcapa(iibig,ibig,2)+2
      kbarr(jf)=kdcapa(iibig,ibig,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,isqurar(1,ibig,2)+2
      kbarr(jf)=isqurar(1,ibig,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      ksol(jjbig,iibig,ibig,jf)=mcarr(jf)
!      print *,'jjbig',jjbig,'iibig',iibig,'ibig',ibig
      end do
      end do
      jjbig=jjbig+1
      if (jjbig.eq.3)goto 135
!     put complement of isquare into isquare(isqurar)
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do ibig=1,2
      if (isqurar(1,ibig,2).eq.0)goto 132
      do jf=1,isqurar(1,ibig,2)+2
      kbarr(jf)=isqurar(1,ibig,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      isqurar(1,ibig,jf)=kcarr(jf)
      end do
132   end do
      goto 140
135  if (iibig.eq.2)goto 201
     
     do mmm=1,2
     do ibig=1,2 
!     print *,'mmm',mmm,'ibig',ibig,'ksol',(ksol(mmm,1,ibig,jk),jk=&
!     1,ksol(mmm,1,ibig,2)+2)
     end do 
     end do 
      
      if (inxyz.eq.3)goto 201
144   do ibig=1,2
      if (kcapb(1,ibig,2).eq.0)goto 128
      do jf=1,kcapb(1,ibig,2)+2
      kbarr(jf)=kcapb(1,ibig,jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kcapb(2,ibig,jf)=kcarr(jf)
      end do
      goto 129
128   kcapb(2,ibig,1)=0
      kcapb(2,ibig,2)=0
129   end do
      if (iibig.eq.2)goto 201
      iibig=2
      jjbig=1
      
      
      goto 141



      











201   if (nnsq.eq.0)goto 998
      print *,'nnsq',nnsq
      
      
      do i=1,2
      do j=1,2
      if (jinj(j).eq.0)goto 2010
!      do j=1,nnsq/2
      do k=1,2
      print *,'ksols',(ksol(i,j,k,jf),jf=1,&
      ksol(i,j,k,2)+2)
      end do
2010  end do
      end do
      
      if (nnsq.eq.0)goto 998
      do i=1,2
      do j=1,2
      if (jinj(j).eq.0)goto 2011
      do jf=1,inv4(2)+2
      marr(jf)=inv4(jf)
      end do
      do jf=1,ka3(2)+2
      mbarr(jf)=ka3(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,ksol(i,j,1,2)+2
      karr(jf)=ksol(i,j,1,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kcarr(jf)=mcarr(jf)
      end do
      if (mcarr(1).eq.0)goto 150
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
150   do jf=1,kcarr(2)+2
      isol(i,j,1,jf)=kcarr(jf)
      end do
      do jf=1,ksol(i,j,2,2)+2
      isol(i,j,2,jf)=ksol(i,j,2,jf)
      end do
      print *,'isols',(isol(i,j,1,jk),jk=1,isol(i,j,1,2)+2)
      print *,'isolcs',(isol(i,j,2,jk),jk=1,isol(i,j,2,2)+2)
2011  end do
      end do
!      print *,'kdis1',(kdis(iibig,1,jk),jk=1,kdis(iibig,1,2)+2)
      
!      print *,'kdcapa1',(kdcapa(1,1,jk),jk=1,kdcapa(1,1,2)+2)
      
!      print *,'kdcapa2',(kdcapa(2,1,jk),jk=1,kdcapa(2,1,2)+2)
      
!      print *,'jnvcom',(jnvcom(1,jk),jk=1,jnvcom(1,2)+2)
!      print *,'jnvcomc',(jnvcom(2,jk),jk=1,jnvcom(2,2)+2),'iibig',iibig
!      print *,'kcofc',(kcofc(2,1,jf),jf=1,kcofc(2,1,2)+2)
!      print *,'kcapb',(kcapb(2,1,jf),jf=1,kcapb(2,1,2)+2)
!      print *,'iroota',(iroota(1,jf),jf=1,iroota(1,2)+2)
!      print *,'kp',(kp(jf),jf=1,kp(2)+2)
!      print *,'krr',(krr(jf),jf=1,krr(2)+2)
!      print *,'kq',(kq(jf),jf=1,kq(2)+2),'kkbig',kkbig










      
      
      
      
      
      
      
      
      
      goto 1000
998   kkbig=kkbig+1
      if (kkbig.gt.ncub)goto 999
      goto 160








203   stop
999   print *,'no solutions mod',(ip(jf),jf=1,ip(2)+2)
1000  end
      
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
      character*(*) in,out
      out =in
      i =len(out)
      do while(out(i:i) ==' ')
      out(i:i)='*'
      i = i-1
      end do 
      return 
      end

      subroutine menmul
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
      
      if (mbarr(2).eq.0)goto 11
      if (marr(2).eq.0) goto 9
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
      mdarr(1)=mod (marr(1)+mbarr(1),2)
      goto 15

11    print *,'halted : attempted division by zero,routine mendiv'
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100),khsq(50)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(100),mbarr(100),mcarr(200),mdarr(100)
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
       subroutine bwq5
       
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100),khsq(50)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(100),mbarr(100),mcarr(200),mdarr(100)
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       iconz=1
       nn=1
       do jbig=1,2
       if (ncom(jbig,2).eq.0)goto 1104
       do jf=3,ncom(jbig,2)+2
       karr(jf-2)=ncom(jbig,jf)
       end do
       ilen=ncom(jbig,2)
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       do jf=1,irlen
       iaa(jbig,jf+2)=irrr(jf)
       end do
       iaa(jbig,1)=0
       iaa(jbig,2)=irlen
       goto 1103
1104   iaa(jbig,1)=0       
       iaa(jbig,2)=0
1103   end do
       

       
1108   do i=1,2
       do j=1,3
       ix(i,j)=0
       end do
       end do
       do jf=3,ip(2)+2
       karr(jf-2)=ip(jf)
       kbarr(jf-2)=ip(jf)
       end do
       ilen=ip(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf+2)=kcarr(jf)
       end do
       karr(1)=0
       karr(2)=ilen3
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       itot(jf)=kcarr(jf)
       iprecod(jf)=kcarr(jf)
       end do
       print *,'firprecod',(iprecod(jf),jf=1,iprecod(2)+2)
       
       do ibig=1,200
       do jf=3,iprecod(2)+2
       karr(jf-2)=iprecod(jf)
       end do
       ilen=iprecod(2)
       
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (icont.eq.0)goto 84
       if (irlen.ne.0)goto 84
       do jf=1,icont
       iprecod(jf+2)=ipqt(jf)
       end do
       iprecod(1)=0
       iprecod(2)=icont
       end do
       print *,'loop too short' 
       stop
84     ivan=iprecod(2)+2
       irem=mod(iprecod(ivan),2)
       print *,'precod',(iprecod(jf),jf=1,iprecod(2)+2),'irem',irem
       
       
       do jf=1,iprecod(2)+2
       karr(jf)=iprecod(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2-irem
       call mpadd(0)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ipn(jf+2)=ipqt(jf)
       end do
       ipn(1)=0
       ipn(2)=icont
       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 1600
98     do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ix(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2),'ix2',(ix(2,jf),jf=1,&
       ix(2,2)+2)
       print *,'ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ib(ibig,jf)=ibprod(ibig,jf)
       ibperm(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ixperm(ibig,jf)=ix(ibig,jf)
       end do
       end do
132    nn=nn*607        
       nn=mod(nn,5003)
       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
       
       nn=nn*607
       nn=mod(nn,5003)
       iaas(2,1)=0
       iaas(2,2)=1
       iaas(2,3)=nn
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       if (ibsw.eq.1)goto 1620
112    call sub516
       print *,'ok4'
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iz(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 172
       goto 174
172    if (ib(2,2).eq.0)goto 228
174    icon=1
173    do ibig=1,2
       do jf=1,ib(ibig,2)+2
       iacn(ibig,jf)=ib(ibig,jf)
       end do
       do jf=1,iz(ibig,2)+2
       ibcn(ibig,jf)=iz(ibig,jf)
       end do
       end do
       print *,'ok5',' icon',icon,'ib1',(ib(1,jf),jf=1,ib(1,2)+2)
       print *,'ok5','ib2',(ib(2,jf),jf=1,ib(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ib(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 192
       goto 194
192    if (ib(2,2).eq.0)goto 200
194    icon=icon+1
       if (icon.eq.50)goto 280
       goto 173
200    if (mod(icon,2).eq.1)goto 280
       do ibig=1,2
       do jf=1,iz(ibig,2)+2
       iaas(ibig,jf)=iz(ibig,jf)
       end do
       end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=icon/2
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ibcn(ibig,jf)=ix(ibig,jf)
       end do
       end do
       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
       ),jf=1,ibcn(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ix(ibig,jf)=icprod(ibig,jf)
       end do
       end do
228    print *,'square root=',(ix(1,jf),jf=1,ix(1,2)+2),'ipn',ipn(3)
       print *,'square root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       isqurar(1,ibig,jf)=ix(ibig,jf)
       end do
       end do
       nsq=2
       
       goto 300
280    print *,'no square root exists',' iconz',iconz
       do ibig=1,2
       do jf=1,ixperm(ibig,2)+2
       ix(ibig,jf)=ixperm(ibig,jf)
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2)
       
       
       do jf=1,ibperm(ibig,2)+2
       ib(ibig,jf)=ibperm(ibig,jf)
       end do
       end do
       iconz=iconz+1
       if (iconz.eq.30)goto 232
       goto 132
1600   ipn(1)=0
       ipn(2)=1
       ipn(3)=2
       ibsw=1
       goto 98
1620   ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       goto 112
232    nsq=0

300    return
       end


       
       subroutine sub400(id,ipz,k)
       ip=ipz
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
       
       
       
       
       
       
       
       subroutine cub5
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100),khsq(50)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(100),mbarr(100),mcarr(200),mdarr(100)
       dimension itempz(2),iaa(2,100),inv2(100),minusb(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       dimension nn(5)
       nscon=0
       ibsw=0
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       iconz=1
       nn(1)=0
       nn(2)=1
       nn(3)=1
       do jbig=1,2
       if (ncom(jbig,2).eq.0)goto 1104
       do jf=3,ncom(jbig,2)+2
       karr(jf-2)=ncom(jbig,jf)
       end do
       ilen=ncom(jbig,2)
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       do jf=1,irlen
       iaa(jbig,jf+2)=irrr(jf)
       end do
       iaa(jbig,1)=0
       iaa(jbig,2)=irlen
       goto 1103
1104   iaa(jbig,1)=0       
       iaa(jbig,2)=0
1103   end do
       

       
1108   do i=1,2
       do j=1,3
       ix(i,j)=0
       end do
       end do
       do jf=3,ip(2)+2
       karr(jf-2)=ip(jf)
       kbarr(jf-2)=ip(jf)
       end do
       ilen=ip(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf+2)=kcarr(jf)
       end do
       karr(1)=0
       karr(2)=ilen3
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       itot(jf)=kcarr(jf)
       iprecod(jf)=kcarr(jf)
       end do
       print *,'firprecod',(iprecod(jf),jf=1,iprecod(2)+2)
       
       do ibig=1,200
       do jf=3,iprecod(2)+2
       karr(jf-2)=iprecod(jf)
       end do
       ilen=iprecod(2)
       
       ilen2=1
       kbarr(1)=3
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (icont.eq.0)goto 84
       if (irlen.ne.0)goto 84
       do jf=1,icont
       iprecod(jf+2)=ipqt(jf)
       end do
       iprecod(1)=0
       iprecod(2)=icont
       end do
       print *,'loop too short' 
       stop
84     a=a
       irem=irrr(1)
       iremp=irem
       print *,'precod',(iprecod(jf),jf=1,iprecod(2)+2),'irem',irem
       
       
       do jf=1,iprecod(2)+2
       karr(jf)=iprecod(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=3-irem
       call mpadd(0)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       ilen2=1
       kbarr(1)=3
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ipn(jf+2)=ipqt(jf)
       end do
       ipn(1)=0
       ipn(2)=icont
!       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
!       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 1600
98     do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ix(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
!       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2),'ix2',(ix(2,jf),jf=1,&
!       ix(2,2)+2)
!       print *,'ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ib(ibig,jf)=ibprod(ibig,jf)
       ibperm(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ixperm(ibig,jf)=ix(ibig,jf)
       end do
       end do
132    izzc=1       
       
310    do jf=1,nn(2)+2       
       marr(jf)=nn(jf)
       end do
       mbarr(1)=0
       mbarr(2)=2
       mbarr(3)=4035
       mbarr(4)=3607
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=3
       mbarr(3)=1
       mbarr(4)=0
       mbarr(5)=0
       call mendiv
       do jf=1,mcarr(2)+2
       nn(jf)=mcarr(jf)
       iaas(izzc,jf)=mcarr(jf)
       end do
       if (izzc.eq.2)goto 312
       izzc=izzc+1
       goto 310
       

       
312    do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       if (ibsw.eq.1)goto 1620
112    call sub516
       print *,'ok4'
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iz(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 172
       goto 174
172    if (ib(2,2).eq.0)goto 228
174    icon=1
173    do ibig=1,2
       do jf=1,ib(ibig,2)+2
       iacn(ibig,jf)=ib(ibig,jf)
       end do
       do jf=1,iz(ibig,2)+2
       ibcn(ibig,jf)=iz(ibig,jf)
       end do
       end do
!       print *,'ok5',' icon',icon,'ib1',(ib(1,jf),jf=1,ib(1,2)+2)
!       print *,'ok5','ib2',(ib(2,jf),jf=1,ib(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ib(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 192
       goto 194
192    if (ib(2,2).eq.0)goto 200
194    icon=icon+1
       if (icon.eq.50)goto 280
       goto 173
200    if (mod(icon,3).ne.0)goto 280
       do ibig=1,2
       do jf=1,iz(ibig,2)+2
       iaas(ibig,jf)=iz(ibig,jf)
       end do
       end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=icon/3
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ibcn(ibig,jf)=ix(ibig,jf)
       end do
       end do
!       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
!       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
!       ),jf=1,ibcn(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ix(ibig,jf)=icprod(ibig,jf)
       end do
       end do
228    print *,'cube root=',(ix(1,jf),jf=1,ix(1,2)+2),'ipn',ipn(3)
       print *,'cube root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       icubar(1,ibig,jf)=ix(ibig,jf)
       end do
       end do
       if (iremp.eq.2)goto 1511
       
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       ncom(ibig,jf)=ix(ibig,jf)
       end do
       end do
       call bwq5
       if (nsq.ne.0)goto 1502
       print *,'nsq',nsq
       goto 1510
1502   do ibig=1,2
       do jf=1,isqurar(1,ibig,2)+2
       iaas(ibig,jf)=isqurar(1,ibig,jf)
       end do
       end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       call sub516
       print *,'icprod1',(icprod(1,jf),jf=1,icprod(1,2)+2)
       print *,'icprod2',(icprod(2,jf),jf=1,icprod(2,2)+2)
       
       print *,'iaa1',(iaa(1,jf),jf=1,iaa(1,2)+2)
       print *,'iaa2',(iaa(2,jf),jf=1,iaa(2,2)+2)
       print *,'khsq',(khsq(jf),jf=1,khsq(2)+2)
       print *,'icubar',(icubar(1,1,jf),jf=1,icubar(1,1,2)+2)
       print *,'icubarc',(icubar(1,2,jf),jf=1,icubar(1,2,2)+2)
       
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       if (icprod(ibig,jf).ne.iaa(ibig,jf))goto 1504
       end do
       end do
       goto 1505
1504   do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       do ibig=1,2
       if (isqurar(1,ibig,2).eq.0)goto 1506
       do jf=1,isqurar(1,ibig,2)+2
       kbarr(jf)=isqurar(1,ibig,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       iaas(ibig,jf)=kcarr(jf)
       icubar(1,ibig,jf)=kcarr(jf)
       end do
       goto 1507
       
1506   iaas(ibig,1)=0       
       iaas(ibig,2)=0
1507   end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       call sub516
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       if (icprod(ibig,jf).ne.iaa(ibig,jf))goto 1510
       end do
       end do
       goto 1511
1510   nscon=nscon+1       
       if (nscon.ge.20)goto 15110
       do jf=1,nn(2)+2
       karr(jf)=nn(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       nn(jf)=kcarr(jf)
       end do
       ibsw=0
       goto 1108

       
       
15110   print *,'major error'
        goto 232











       
       
1505   print *,'iremp=',iremp,'cube root=',(isqurar(1,1,jf),jf=1,&
       isqurar(1,1,2)+2)
       print *,'cube root comp=',(isqurar(1,2,jf),jf=1,&
       isqurar(1,2,2)+2)
1500   do ibig=1,2    
       do jf=1,isqurar(1,ibig,2)+2
       icubar(1,ibig,jf)=isqurar(1,ibig,jf)
       end do
       end do
1511   do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=3
       call mpadd(1)
       do jf=1,kcarr(2)+2
       ncom(1,jf)=kcarr(jf)
       end do
       ncom(2,1)=0
       ncom(2,2)=0
       call bwq5
       if (nsq.ne.0)goto 1512 
       ncubr=1
       print *,'only 1 square root'
       goto 300
1512   print *,'-3 root',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2) 
       print *,'-3 comp',(isqurar(1,2,jf),jf=1,isqurar(1,2,2)+2)
       do jf=3,ip(2)+2
       karr(jf-2)=ip(jf)
       end do
       ilen=ip(2)
       kbarr(1)=2
       ilen2=1
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
       inv2(jf)=kcarr(jf)
       end do
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       do ibig=1,2
       if (icubar(1,ibig,2).eq.0)goto 1514
       do jf=1,icubar(1,ibig,2)+2
       kbarr(jf)=icubar(1,ibig,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       minusb(ibig,jf)=kcarr(jf)
       end do
       goto 1515
1514   minusb(ibig,1)=0
       minusb(ibig,2)=0
1515   end do
       do ibig=1,2
       do jf=1,icubar(1,ibig,2)+2
       iacn(ibig,jf)=icubar(1,ibig,jf)
       end do
       do jf=1,isqurar(1,ibig,2)+2
       ibcn(ibig,jf)=isqurar(1,ibig,jf)
       end do
       end do
       call sub1100
       print *,'minusb',(minusb(1,jf),jf=1,minusb(1,2)+2)
       print *,'minusbcomp',(minusb(2,jf),jf=1,minusb(2,2)+2)
       print *,'inv2',(inv2(jf),jf=1,inv2(2)+2)
       
!  keep product of root(-3) and icubar in icprod       
       iconz2=2
1521   do ibig=1,2
       do jf=1,minusb(ibig,2)+2
       karr(jf)=minusb(ibig,jf)
       end do
       do jf=1,icprod(ibig,2)+2
       kbarr(jf)=icprod(ibig,jf)
       end do
       call mpadd(0)
       if (kcarr(2).eq.0)goto 1516
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,inv2(2)+2
       kbarr(jf-2)=inv2(jf)
       end do
       ilen2=inv2(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       print *,'ilen',ilen,'karrs',(karr(jf),jf=1,ilen)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       print *,'ilen2',ilen2,'iconz2',iconz2,'ibig',ibig
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       print *,'irlen',irlen,'kbarrs',(kbarr(jf),jf=1,ilen2)
       if (irlen.eq.0)goto 1516
       do jf=1,irlen
       icubar(iconz2,ibig,jf+2)=irrr(jf)
       end do
       icubar(iconz2,ibig,1)=0
       icubar(iconz2,ibig,2)=irlen
       goto 1517
1516   icubar(iconz2,ibig,1)=0
       icubar(iconz2,ibig,2)=0
1517   end do
       iconz2=iconz2+1
       if (iconz2.eq.4)goto 1520
! obtain complements of icprod and restore in icprod
       print *,'icubar',(icubar(2,1,jf),jf=1,icubar(2,1,2)+2)
       print *,'ipmiddle',(ip(jk),jk=1,ip(2)+2)
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       do ibig=1,2
       print *,'ip2mid',(ip(jk),jk=1,ip(2)+2)
       if (icprod(ibig,2).eq.0)goto 1518
       do jf=1,icprod(ibig,2)+2
       kbarr(jf)=icprod(ibig,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       icprod(ibig,jf)=kcarr(jf)
       end do
1518   end do 
       print *,'ipbottom',(ip(jk),jk=1,ip(2)+2)
       goto 1521
1520  print *,'2icubar',(icubar(2,1,jf),jf=1,icubar(2,1,2)+2) 
      print *,'2icubarcomp',(icubar(2,2,jf),jf=1,icubar(2,2,2)+2)
      print *,'3icubar',(icubar(3,1,jf),jf=1,icubar(3,1,2)+2)
      print *,'3icubarc',(icubar(3,2,jf),jf=1,icubar(3,2,2)+2)
      print *,'1icubar',(icubar(1,1,jf),jf=1,icubar(1,1,2)+2)
      print *,'1icubarc',(icubar(1,2,jf),jf=1,icubar(1,2,2)+2)
      print *,'iaa1',(iaa(1,jf),jf=1,iaa(1,2)+2)
      print *,'iaa2',(iaa(2,jf),jf=1,iaa(2,2)+2)
      ncubr=3


       
       
       goto 300
280    print *,'no cube root exists',' iconz',iconz
       do ibig=1,2
       do jf=1,ixperm(ibig,2)+2
       ix(ibig,jf)=ixperm(ibig,jf)
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2)
       
       
       do jf=1,ibperm(ibig,2)+2
       ib(ibig,jf)=ibperm(ibig,jf)
       end do
       end do
       iconz=iconz+1
       if (iconz.eq.300)goto 232
       goto 132
1600   ipn(1)=0
       ipn(2)=1
       ipn(3)=2
       ibsw=1
       goto 98
1620   ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       goto 112
232    ncubr=0
300    return
       end

       
       
          
       
       
       
       subroutine sub516
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100),khsq(50)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(100),mbarr(100),mcarr(200),mdarr(100)
       dimension ipn2(100),ie(100),igg(2,100)


       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 300
       goto 301
300    do ibig=1,2
       do jf=1,iaas(ibig,2)+2
       icprod(ibig,jf)=iaas(ibig,jf)
       ibprod(ibig,jf)=iaas(ibig,jf)
       end do
       end do
       goto 620

301    ie(1)=0
       ilen=1
       karr(1)=2
       ilen2=1
       kbarr(1)=2
4      call mpmul(ilen,ilen2,ilen3)
       if (ilen3.lt.ipn(2))goto 22
       if (ilen3.gt.ipn(2))goto 20
       do i=1,ilen3
       if (ipn(i+2).lt.kcarr(i))goto 20
       if (ipn(i+2).gt.kcarr(i))goto 22
       end do
       goto 22
20     ie(2)=ilen3       
       do i=1,ilen3
       ie(i+2)=kcarr(i)
       end do
       goto 30
22     do i=1,ilen3       
       kbarr(i)=kcarr(i)
       end do
       ilen2=ilen3
       goto 4
30     do i=1,ipn(2)+2
       ipn2(i)=ipn(i)
       end do
       do i=1,ie(2)
       karr(i)=ie(i+2)
       end do
       ilen=ie(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do i=1,icont
       ie(i+2)=ipqt(i)
       end do
       ie(2)=icont
       print *,'ie',(ie(jf),jf=1,ie(2)+2)
       
       do i=1,ipn2(2)+2
       karr(i)=ipn2(i)
       end do
       do i=1,ie(2)+2
       kbarr(i)=ie(i)
       end do
       call mpadd(1)
       do i=1,kcarr(2)+2
       ipn2(i)=kcarr(i)
       end do
       do ibig=1,2
       do jf=1,iaas(ibig,2)+2
       igg(ibig,jf)=iaas(ibig,jf)
       ibprod(ibig,jf)=iaas(ibig,jf)
       end do
       end do
31     if ((ie(2).eq.1).and.(ie(3).eq.1))goto 620
       do jf=3,ie(2)+2
       karr(jf-2)=ie(jf)
       end do
       ilen=ie(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ie(jf+2)=ipqt(jf)
       end do
       ie(2)=icont
!       print *,'ok2','ie',(ie(jf),jf=1,ie(2)+2)
       
!       if (ipn2(2).lt.ie(2))goto 31
!       if (ipn2(2).gt.ie(2))goto 42
!       do jf=3,ipn2(2)+2
!       if (ipn2(jf).lt.ie(jf))goto 31
!       if (ipn2(jf).gt.ie(jf))goto 42
!       end do
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       ibcn(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ibprod(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if (ipn2(2).lt.ie(2))goto 31
       if (ipn2(2).gt.ie(2))goto 42
       do jf=3,ipn2(2)+2
       if (ipn2(jf).lt.ie(jf))goto 31
       if (ipn2(jf).gt.ie(jf))goto 42
       end do
42     a=a       
       
       do jf=1,ipn2(2)+2
       karr(jf)=ipn2(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       ipn2(jf)=kcarr(jf)
       end do
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       do ibig=1,2
       do jf=1,igg(ibig,2)+2
       ibcn(ibig,jf)=igg(ibig,jf)
       end do
       end do
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ibprod(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       goto 31
620    print *,'ok3'
       return
       end


       subroutine sub1100
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50),
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       
       common kara(50),karb(50),kard(50),karp(50),karv(50)

       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100),khsq(50)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(100),mbarr(100),mcarr(200),mdarr(100)
       do jf=1,iacn(1,2)+2
       marr(jf)=iacn(1,jf)
       end do
       do jf=1,ibcn(1,2)+2
       mbarr(jf)=ibcn(1,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       icprod(1,jf)=mcarr(jf)
       end do
       do jf=1,iacn(2,2)+2
       marr(jf)=iacn(2,jf)
       end do
       do jf=1,ibcn(2,2)+2
       mbarr(jf)=ibcn(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ip(2)+2
       mbarr(jf)=ip(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,khsq(2)+2
       mbarr(jf)=khsq(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,icprod(1,2)+2
       kbarr(jf)=icprod(1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ip(2)+2
       mbarr(jf)=ip(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       icprod(1,jf)=mcarr(jf)
       end do
       
3      if ((iacn(1,2).eq.0).or.(ibcn(2,2).eq.0))goto 104
       do jf=3,iacn(1,2)+2
       karr(jf-2)=iacn(1,jf)
       end do
       ilen=iacn(1,2)
       do jf=3,ibcn(2,2)+2
       kbarr(jf-2)=ibcn(2,jf)
       end do
       ilen2=ibcn(2,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       icprod(2,jf+2)=kcarr(jf)
       end do
       icprod(2,1)=0
       icprod(2,2)=ilen3
       goto 141
104    icprod(2,1)=0
       icprod(2,2)=0
141    if ((iacn(2,2).eq.0).or.(ibcn(1,2).eq.0))goto 142       
       do jf=3,iacn(2,2)+2
       karr(jf-2)=iacn(2,jf)
       end do
       ilen=iacn(2,2)
       do jf=3,ibcn(1,2)+2
       kbarr(jf-2)=ibcn(1,jf)
       end do
       ilen2=ibcn(1,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=ilen3
       do jf=1,icprod(2,2)+2
       karr(jf)=icprod(2,jf)
       end do
       goto 151
142    do jf=1,icprod(2,2)+2
       kcarr(jf)=icprod(2,jf)
       end do
       goto 152
151    call mpadd(0)
152    if (kcarr(2).eq.0)goto 101
       isgn=kcarr(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 101
       if (isgn.eq.0)goto 102
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=isgn
       kbarr(2)=irlen
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       icprod(2,jf)=kcarr(jf)
       end do
       goto 200
101    icprod(2,1)=0
       icprod(2,2)=0
       goto 200
102    do jf=1,irlen
       icprod(2,jf+2)=irrr(jf)
       end do
       icprod(2,1)=isgn
       icprod(2,2)=irlen
200    a=a
! 200    print *,'icprod1',(icprod(1,jf),jf=1,icprod(1,2)+2)
!       print *,'icprod2',(icprod(2,jf),jf=1,icprod(2,2)+2)
!       print *,'iacn1',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',&
!       (iacn(2,jf),jf=1,iacn(2,2)+2)
!       print *,'ibcn1',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',&
!       (ibcn(2,jf),jf=1,ibcn(2,2)+2)
       
       return
       end
       subroutine cornsq
! field is p       
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100),khsq(50)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(100),mbarr(100),mcarr(200),mdarr(100)
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       iconz=1
       nn=1
       do jbig=1,2
       if (ncom(jbig,2).eq.0)goto 1104
       do jf=3,ncom(jbig,2)+2
       karr(jf-2)=ncom(jbig,jf)
       end do
       ilen=ncom(jbig,2)
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       do jf=1,irlen
       iaa(jbig,jf+2)=irrr(jf)
       end do
       iaa(jbig,1)=0
       iaa(jbig,2)=irlen
       goto 1103
1104   iaa(jbig,1)=0       
       iaa(jbig,2)=0
1103   end do
       

       
1108   do i=1,2
       do j=1,3
       ix(i,j)=0
       end do
       end do
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
!       kbarr(jf-2)=ip(jf)
       end do
!       ilen=ip(2)
!       ilen2=ilen
!       call mpmul(ilen,ilen2,ilen3)
!       do jf=1,ilen3
!       karr(jf+2)=kcarr(jf)
!       end do
!       karr(1)=0
!       karr(2)=ilen3
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       itot(jf)=kcarr(jf)
       iprecod(jf)=kcarr(jf)
       end do
       print *,'firprecod',(iprecod(jf),jf=1,iprecod(2)+2)
       
       do ibig=1,200
       do jf=3,iprecod(2)+2
       karr(jf-2)=iprecod(jf)
       end do
       ilen=iprecod(2)
       
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (icont.eq.0)goto 84
       if (irlen.ne.0)goto 84
       do jf=1,icont
       iprecod(jf+2)=ipqt(jf)
       end do
       iprecod(1)=0
       iprecod(2)=icont
       end do
       print *,'loop too short' 
       stop
84     ivan=iprecod(2)+2
       irem=mod(iprecod(ivan),2)
       print *,'precod',(iprecod(jf),jf=1,iprecod(2)+2),'irem',irem
       
       
       do jf=1,iprecod(2)+2
       karr(jf)=iprecod(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2-irem
       call mpadd(0)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ipn(jf+2)=ipqt(jf)
       end do
       ipn(1)=0
       ipn(2)=icont
       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 1600
98     do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ix(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2),'ix2',(ix(2,jf),jf=1,&
       ix(2,2)+2)
       print *,'ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ib(ibig,jf)=ibprod(ibig,jf)
       ibperm(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ixperm(ibig,jf)=ix(ibig,jf)
       end do
       end do
132    nn=nn*607        
       nn=mod(nn,5003)
       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
       iaas(2,1)=0
       iaas(2,2)=0
       kard(1)=0
       kard(2)=1
       kard(3)=nn
       do jf=1,ip(2)+2
       karp(jf)=ip(jf)
       end do
       call mpkron(k)
       if (k.ne.-1)goto 132

!       nn=nn*607
!       nn=mod(nn,5003)
!       iaas(2,1)=0
!       iaas(2,2)=1
!       iaas(2,3)=nn
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       if (ibsw.eq.1)goto 1620
112    call sub516
       print *,'ok4'
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iz(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 172
       goto 174
172    if (ib(2,2).eq.0)goto 228
174    icon=1
173    do ibig=1,2
       do jf=1,ib(ibig,2)+2
       iacn(ibig,jf)=ib(ibig,jf)
       end do
       do jf=1,iz(ibig,2)+2
       ibcn(ibig,jf)=iz(ibig,jf)
       end do
       end do
       print *,'ok5',' icon',icon,'ib1',(ib(1,jf),jf=1,ib(1,2)+2)
       print *,'ok5','ib2',(ib(2,jf),jf=1,ib(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ib(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 192
       goto 194
192    if (ib(2,2).eq.0)goto 200
194    icon=icon+1
       if (icon.eq.50)goto 280
       goto 173
200    if (mod(icon,2).eq.1)goto 280
       do ibig=1,2
       do jf=1,iz(ibig,2)+2
       iaas(ibig,jf)=iz(ibig,jf)
       end do
       end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=icon/2
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ibcn(ibig,jf)=ix(ibig,jf)
       end do
       end do
       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
       ),jf=1,ibcn(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ix(ibig,jf)=icprod(ibig,jf)
       end do
       end do
228    print *,'square root=',(ix(1,jf),jf=1,ix(1,2)+2),'ipn',ipn(3)
       print *,'square root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       isqurar(1,ibig,jf)=ix(ibig,jf)
       end do
       end do
       nsq=2
       
       goto 300
280    print *,'no square root exists',' iconz',iconz
       do ibig=1,2
       do jf=1,ixperm(ibig,2)+2
       ix(ibig,jf)=ixperm(ibig,jf)
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2)
       
       
       do jf=1,ibperm(ibig,2)+2
       ib(ibig,jf)=ibperm(ibig,jf)
       end do
       end do
       iconz=iconz+1
       if (iconz.eq.30)goto 232
       goto 132
1600   ipn(1)=0
       ipn(2)=1
       ipn(3)=2
       ibsw=1
       goto 98
1620   ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       goto 112
232    nsq=0

300    return
       end


       
       













































































































































