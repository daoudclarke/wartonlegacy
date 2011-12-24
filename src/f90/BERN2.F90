    

      program bern2 
!     multi-precisioning of bern1 and extending to complex roots      
!     first in generalised GNFS suite    
!     smaller sieving interval
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      dimension ka3(50),ka2(50),ka1(50),kc(50),ks1(50),ks2(50)
      dimension ks3(50),ks4(50),ks12(50)
      dimension kp(50),inv2(50),inv3(50),inv4(50),inv8(50)
      dimension inv16(50),inv27(50),inv256(50),kcub(6),kzans(30),ksol(4) 
      goto 777
      
      ncom(1)=0
      ncom(2)=1
      ncom(3)=17
      ip=1009
      print *,'length radix 10000?'
      read *,ncom(2)
      print *,'number?'
      read *,(ncom(jf),jf=3,ncom(2)+2)
      print *,'modulus?'
      read *,ip
      call cub5(ip,ians,icubrn,icubr1,icubr2,icubr3)
      print *,'ians',ians,'icubrn',icubrn,'sols',icubr1,icubr2,icubr3
      stop
777   ip=131



      iprar(1)=0
      iprar(2)=1
      iprar(3)=ip
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
      if (ip.lt.10000)goto 1 
      kara(1)=0
      kara(2)=2
      kara(3)=ip/10000
      kara(4)=ip-kara(3)*10000
      goto 2
1     kara(1)=0
      kara(2)=1
      kara(3)=ip
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
      print *,'inv16',inv16(3)
      
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
57    ks13=ks1(3)*ks12(3)
      ks13=mod(ks13,ip)
      print *,'ks1(3)',ks1(3),'ks2',ks2(3),'ks3',ks3(3),'ks4',ks4(3)
      kq=mod((ip-1)*ks13,ip)
      kq=mod(kq*inv8(3),ip)
      itemp=mod(inv2(3)*ks1(3),ip)
      itemp=mod(itemp*ks2(3),ip)
      kq=mod(kq+itemp-ks3(3),ip)
!      kq=(ip-1)*ks13*inv8(3)+inv2(3)*ks1(3)*ks2(3)
!      kq=kq-ks3(3)
!      kq=mod(kq,ip)
      ks14=mod(ks1(3)*ks13,ip)


      krr=mod((ip-3)*inv256(3),ip)
      krr=mod(krr*ks14,ip)
      itemp=mod(inv16(3)*ks12(3),ip)
      itemp=mod(itemp*ks2(3),ip)
      krr=mod(krr+itemp,ip)
      itemp=mod(inv4(3)*ks1(3),ip)
      itemp=mod(itemp*ks3(3),ip)
      
      krr=mod(krr-itemp+ks4(3),ip)
      print *,'krr',krr
      
      kpf=ip-kp(3)
      kqf=(ip-4)*krr
      kqf=mod(kqf,ip)
      krf=mod(4*kp(3)*krr,ip)
      krf=mod(krf-kq*kq,ip)
      ka=(ip-1)*inv3(3)
      ka=mod(ka,ip)
      ka=mod(ka*kpf,ip)
      ka=mod(ka*kpf+kqf,ip)
      print *,'ka',ka,'kpf',kpf,'krf',krf,'kq',kq
      
      kb=2*inv27(3)*kpf
      kb=mod(kb,ip)
      kb=mod(kb*kpf,ip)
      kb=mod(kb*kpf,ip)
      print *,'fir kb',kb
      itemp=mod((ip-1)*inv3(3),ip)
      itemp=mod(itemp*kpf,ip)
      itemp=mod(itemp*kqf,ip)
      kb=mod(itemp+kb+krf,ip)
      print *,'inc kb',kb,'itemp',itemp
      khsq=mod(kb*kb,ip)
      
      itemp=mod(ka*ka,ip)
      itemp=mod(itemp*ka,ip)
      itemp=mod(itemp*inv27(3),ip)
      khsq=mod(khsq*inv4(3)+itemp,ip)
      print *,'ks12',ks12(3),'ks13',ks13,'ks14',ks14,'kq',kq
      print *,'krr',krr,'kpf',kpf,'kqf',kqf,'krf',krf,'ka',ka
      print *,'kb',kb,'kp',kp(3),'inv3',inv3(3)
      karb(1)=0
      karb(2)=1
      karb(3)=75
      call mpgcd
      print *,'inv75',(karv(jf),jf=1,karv(2)+2)

      
      ncub=0
      icubrn=0
      icon=1
      print *,'khsq',khsq,'ka',ka
      
      
      if (khsq.eq.0)goto 200
      if (khsq.gt.0)goto 60
      khsq=khsq+ip
      goto 60
200   kicr=mod((ip-1)*kb*inv2(3),ip)
      if (kicr.eq.0)goto 999
      if (kicr.gt.0)goto 201
      kicr=kicr+ip
      goto 201

      
60    ncom(1)=0
      ncom(2)=1
      ncom(3)=khsq
      
      do jf=1,iprar(2)+2
      karp(jf)=iprar(jf)
      end do
      do jf=1,ncom(2)+2
      kard(jf)=ncom(jf)
      end do
      call mpkron(k)
      if (k.ne.1)goto 67
      call bwq5(ip,ians)
      print *,'ians',ians,'ka',ka
      ians2=ip-ians
      khcr=mod((ip-kb)*inv2(3)+ians,ip)
      if (khcr.eq.0)goto 202
      if (khcr.ge.0)goto 621
      khcr=khcr+ip
621   ncom(1)=0      
      ncom(2)=1
      ncom(3)=khcr
    
      call cub5(ip,ians,icubrn,icubr1,icubr2,icubr3)
      print *,'khcr',khcr
      print *,'roots',icubrn,icubr1,icubr2,icubr3,'khcr',khcr,'ka',ka
      
      ncub=0
      if (icubrn.eq.0)goto 202
      kcub(icon)=icubr1
      icon=icon+1
      if (icon.gt.icubrn)goto 62
      kcub(icon)=icubr2
      icon=icon+1
      if (icon.gt.icubrn)goto 62
      kcub(icon)=icubr3
      icon=icon+1


62    ncub=icon-1
      print *,'ians2',ians2
202   kicr=mod((ip-1)*kb*inv2(3)+ians2,ip)
      if (kicr.eq.0)goto 68
      if (kicr.gt.0)goto 66
      kicr=kicr+ip
66    print *,'kicr',kicr,'ka',ka
201   ncom(1)=0
      ncom(2)=1
      ncom(3)=kicr
      call cub5(ip,ians,icubrn,icubr1,icubr2,icubr3)
      print *,'roots2',icubrn,icubr1,icubr2,icubr3,'kicr',kicr,'ka',ka
67    ncub=ncub+icubrn
      if (ncub.eq.0)goto 999
      kcub(icon)=icubr1
      icon=icon+1
      if (icon.gt.ncub)goto 68
      kcub(icon)=icubr2
      icon=icon+1
      if (icon.gt.ncub)goto 68
      kcub(icon)=icubr3
68    print *,'no of roots=',ncub
      if (ncub.eq.0)goto 999
      print *,'cube roots',(kcub(jf),jf=1,ncub),'ka',ka
      
      jcon=1
      if (ncub.eq.1)goto 203
      do jf=jcon+1,ncub
      do kf=1,jcon
      if (kcub(jf).ne.kcub(kf))goto 204
      goto 205
204   end do
      jcon=jcon+1
      kcub(jcon)=kcub(jf)
205   end do
      ncub=jcon
203   kcon=0
      do nxx=1,ncub
      karb(1)=0
      karb(2)=1
      karb(3)=kcub(nxx)
      
      
      call mpgcd
      print *,'karv',karv(3)
      itemp=karv(3)
      itemp=mod(itemp*inv3(3),ip)
      itemp=mod(itemp*ka,ip)
      itemp=mod(kcub(nxx)-itemp,ip)
      print *,'itemp',itemp,'ka=',ka
      kuu=itemp
      print *,'kuu',kuu
      kzu=kuu+kp(3)*inv3(3)
      kzu=mod(kzu,ip)
      karb(3)=23
      call mpgcd
      print *,'inv23',karv(3)
      
      print *,'kzu',kzu,'kp',kp(3),'kq',kq,'krr',krr,'ka',ka,'kb',kb
      
      ncom(1)=0
      ncom(2)=1
      ncom(3)=mod(kzu-kp(3),ip)
      if (ncom(3).eq.0)goto 999
      if (ncom(3).gt.0)goto 70
      ncom(3)=ncom(3)+ip
70    kard(3)=ncom(3)
      goto 212
206   ians=0
      goto 221
207   ians=0
      goto 222
212   call mpkron(k)
      if (k.ne.1)goto 90
      call bwq5(ip,ians)
      print *,'ians',ians
      
      kcapa=ians
      kcapa2=ip-ians
      karb(3)=kcapa
      call mpgcd
      itemp=karv(3)
      
      itemp=mod(inv2(3)*itemp,ip)
      kcapb=mod(kq*itemp,ip)
      kcofc=mod(kzu*inv2(3)-kcapb,ip)
      kdis=mod(kcapa*kcapa-4*kcofc,ip)
      if (kdis.eq.0)goto 206
      if (kdis.gt.0)goto 72
      kdis=kdis+ip
72    ncom(1)=0
      ncom(2)=1
      ncom(3)=kdis
      kard(3)=ncom(3)
      call mpkron(k)
      if (k.ne.1)goto 74

      call bwq5(ip,ians)
      print *,'iansz',ians
      itemp=(ip-kcapa)+ians
      kcon=kcon+1
      kzans(kcon)=mod(itemp*inv2(3),ip)
      if (kzans(kcon).ge.0)goto 219
      kzans(kcon)=kzans(kcon)+ip
219   print *,'kzan1',kzans(kcon),'kcon=',kcon
      itemp=(ip-kcapa)+ip-ians
221   kcon=kcon+1
      kzans(kcon)=mod(itemp*inv2(3),ip)
      if (kzans(kcon).ge.0)goto 220
      kzans(kcon)=kzans(kcon)+ip
220   print *,'kzan2',kzans(kcon),'kcon=',kcon
      

      


74    kcapb2=ip-kcapb
      kcofc2=mod(kzu*inv2(3)-kcapb2,ip)
      kdis2=mod(kcapa2*kcapa2-4*kcofc2,ip)
      if (kdis2.eq.0)goto 207
      if (kdis2.gt.0)goto 76
      kdis2=kdis2+ip
76    ncom(1)=0
      ncom(2)=1
      ncom(3)=kdis2
      kard(3)=ncom(3)
      call mpkron(k)
      if (k.ne.1)goto 90
      call bwq5(ip,ians)
      print *,'iansz3',ians
      kcon=kcon+1
      itemp=(ip-kcapa2)+ians
      kzans(kcon)=mod(itemp*inv2(3),ip)
      if (kzans(kcon).ge.0)goto 229
      kzans(kcon)=kzans(kcon)+ip
229   print *,'kzan3',kzans(kcon),'kcapa2',kcapa2
222   kcon=kcon+1
      itemp=(ip-kcapa2)+ip-ians
      kzans(kcon)=mod(itemp*inv2(3),ip)
      if (kzans(kcon).ge.0)goto 230
      kzans(kcon)=kzans(kcon)+ip
230   print *,'kzans4',kzans(kcon),'kcon=',kcon

      print *,'kzu',kzu,'kp',kp(3),'kq',kq,'kr',krr

90    end do
      print *,'kcon',kcon
      
      if (kcon.eq.0)goto 999
      do jf=1,kcon
      print *,'i',jf,'sol',kzans(jf)
      end do
      jcon=1
      if (kcon.eq.1)goto 253
      do jf=jcon+1,kcon
      do kf=1,jcon
      if (kzans(jf).ne.kzans(kf))goto 254
      goto 255
254   end do      
      jcon=jcon+1
      kzans(jcon)=kzans(jf)
255   end do
253   print *,'jcon',jcon
      do jf=1,jcon
      print *,'kzans',kzans(jf)
      end do
      do jf=1,jcon
      ksol(jf)=mod(kzans(jf)-ka3(3)*inv4(3),ip)
      if (ksol(jf).ge.0)goto 260
      ksol(jf)=ksol(jf)+ip
260   print *,'jf',jf,'ksol',ksol(jf)
      end do

      
      
      
      
      
      
      
      
      
      
      
      
      goto 1000
999   print *,'no solutions mod=',ip

1000  end
      
      subroutine sieve(kia,kkb,irecnn,kkmx,icdn,iconq)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50)
      
!     change parameter      
      lprx1=ipr(13000)                                      
      
      
      
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
      do i=1,kcarr(2)+2
      litz(i)=litt(i)
      litd(i)=litt(i)
      end do
      litd(1)=0
      
      iflim=10000
      iab=1
      icur=0
      print *,'litd',(litd(i),i=1,litd(2) +2),'icdn',icdn,'irecnn',irecnn&
      ,kkb
      do i=1,50

      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
      if (iab.eq.13000)goto 1500
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


      

1300  irecnn=irecnn+1
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
      write(4,*)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      
!     note next instruction      
      if (irecnn.lt.19000)goto 1500
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      

1500  return
      

      
      end
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
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
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
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
       
       subroutine cub5(ip,ians,icubrn,icubr1,icubr2,icubr3)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       dimension itempz(5),ibar(5),izar(5),ixar(5),itar(5),iyar(5)
       dimension ibprodar(5),ianar(5),nstr(50),icubar(5),inv2(5),ksol(3)
       print *,'ncom',(ncom(jf),jf=1,ncom(2)+2)
       do jf=1,ncom(2)+2
       nstr(jf)=ncom(jf)
       end do
       
       icong=0
       iconz=0
       do jf=3,ncom(2)+2
       karr(jf-2)=ncom(jf)
       end do
       ilen=ncom(2)
       if (ip.lt.10000)goto 1100
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       iprar(1)=0
       iprar(2)=2
       iprar(3)=kbarr(1)
       iprar(4)=kbarr(2)
       
       goto 1102
1100   kbarr(1)=ip
       ilen2=1
       iprar(1)=0
       iprar(2)=1
       iprar(3)=ip
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       if (irlen.eq.1)goto 1106
       iaa=irrr(1)*10000+irrr(2)
       goto 1108
1104   print *,'factor found early',ip
       stop
1106   iaa=irrr(1)


       
       
1108   ix = 0
       itemp=ip/3
       iqrem=ip-3*itemp
       if (iqrem.eq.2)goto 1500
       iprecod=ip-iqrem
       
10     iprecod =ip -1
       i=0
26     itemp = int(iprecod/3)
       irem1 =iprecod -itemp*3
       if(itemp.eq.0)goto 46
       if(irem1.gt.0)goto 40
       iprecod = itemp
       i =i+1
       if(i.lt.200)goto 26
40     iq =iprecod
       ie = i
       goto 48
46     iq =iprecod
       ie = i
48     i =1
       print *,'qs',iq
       itemp=iq/3
       irem=iq-itemp*3
       ipn=(iq+3-irem)/3
       iprem=irem
       iaas=iaa
       call sub516(ibprod,iaas,ipn,ip)
       ix=ibprod
       ipn=iq
       iaas=iaa
       call sub516(ibprod,iaas,ipn,ip)
       ib=ibprod
       print *,'ib',ib,'ix',ix
       if (ib.lt.10000)goto 5401
       ibar(1)=0
       ibar(2)=2
       ibar(3)=ib/10000
       ibar(4)=ib-ibar(3)*10000
       goto 5402
5401   ibar(1)=0
       ibar(2)=1
       ibar(3)=ib
5402   ibperm=ibprod
       ixperm=ix
       if (ix.lt.10000)goto 5201
       ixar(1)=0
       ixar(2)=2
       ixar(3)=ix/10000
       ixar(4)=ix-ixar(3)*10000
       goto 5202
5201   ixar(1)=0
       ixar(2)=1
       ixar(3)=ix

       
5202   n = 1
52     n = n*607
54     itemp=int(n/1000)
56     irem1 =n-1000 *itemp
       n = irem1
       
       print *,'secib',ib,'ibperm',ibperm
       iaas=n
       ipn=iq
       call sub516(ibprod,iaas,ipn,ip)
       iz=ibprod
       print *,'iz',iz,'ie',ie,'ibperm',ibperm
       
       itemp=ib/ip
       irem=ib-ip*itemp
       if (irem.eq.1)goto 200
       icon=1
       ibbz=ib
       
       if (iz.lt.10000)goto 10801
       izar(1)=0
       izar(2)=2
       izar(3)=iz/10000
       izar(4)=iz-izar(3)*10000
       goto 10802
10801  izar(1)=0
       izar(2)=1
       izar(3)=iz
10802  a=a
       do jf=3,ibar(2)+2
       karr(jf-2)=ibar(jf)
       end do
       ilen=ibar(2)
       do jf=3,izar(2)+2
       kbarr(jf-2)=izar(jf)
       end do
       ilen2=izar(2)
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
       ibar(jf+2)=irrr(jf)
       end do
       ibar(1)=0
       ibar(2)=irlen

       isum=0
       do jf=irlen,1,-1
       isum=isum+irrr(jf)*10000**(irlen-jf)
       end do
       ib=isum




1080   if (icon.gt.14)goto 10803
!       print *,'ibfir',ib,'iz',iz,'ibperm',ibperm,'icon',icon
10803  if (ib.eq.1)goto 115
       icon=icon+1
       if (icon.eq.3**ie*500**(iqrem-1))goto 180
       
        
       

1074   do jf=3,izar(2)+2
       karr(jf-2)=izar(jf)
       end do
       ilen=izar(2)
       do jf=3,ibar(2)+2
       kbarr(jf-2)=ibar(jf)
       end do
       ilen2=ibar(2)
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
!       print *,'bresult',(irrr(jf),jf=1,irlen)
       do jf=1,irlen
       ibar(jf+2)=irrr(jf)
       end do
       ibar(1)=0
       ibar(2)=irlen
       isum=0
       do jf=irlen,1,-1
       isum=isum+irrr(jf)*10000**(irlen-jf)
       end do
       ib=isum
!       print *,'ibsum',ib
       goto 1080
115    itemp=icon/3   
       irem=icon-itemp*3
       
       
       if (irem.gt.0)goto 180
!       print *,'icon',icon,'iconz',iconz
       
       iaas=iz
       ipn=itemp
!       print *,'itemp',itemp,'iz',iz,'ibperm',ibperm
       
       call sub516(ibprod,iaas,ipn,ip)

       do jf=3,ixar(2)+2
       karr(jf-2)=ixar(jf)
       end do
       ilen=ixar(2)
       if (ibprod.lt.10000)goto 1151
       ibprodar(1)=0
       ibprodar(2)=2
       ibprodar(3)=ibprod/10000
       ibprodar(4)=ibprod-ibprodar(3)*10000
       goto 1152
1151   ibprodar(1)=0
       ibprodar(2)=1
       ibprodar(3)=ibprod
1152   do jf=3,ibprodar(2)+2
       kbarr(jf-2)=ibprodar(jf)
       end do
       ilen2=ibprodar(2)
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
       ixar(jf+2)=irrr(jf)
       end do
       ixar(1)=0
       ixar(2)=irlen
       isum=0
       do jf=irlen,1,-1
       isum=isum+irrr(jf)*10000**(irlen-jf)
       end do
       ix=isum
200    ipremy=3-iprem
       print *,'cube root=',ix,'qrem=',iqrem,'ind=',ipremy
       
       ians=ix
       if (ipremy.eq.1)goto 2301
       if (ix.lt.10000)goto 203
       ncom(1)=0
       ncom(2)=2
       ncom(3)=ix/10000
       ncom(4)=ix-ncom(3)*10000
       goto 204
203    ncom(1)=0       
       ncom(2)=1
       ncom(3)=ix
204    call bwq5(ip,ians)       
       print *,'ians',ians
2041   if (ians.lt.10000)goto 210
       ianar(1)=0
       ianar(2)=2
       ianar(3)=ians/10000
       ianar(4)=ians-ianar(3)*10000
       goto 212
210    ianar(1)=0
       ianar(2)=1
       ianar(3)=ians
212    do jf=3,ianar(2)+2
       karr(jf-2)=ianar(jf)
       kbarr(jf-2)=ianar(jf)
       end do
       ilen=ianar(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf+2)=kcarr(jf)
       end do
       karr(1)=0
       karr(2)=ilen3
       do jf=1,nstr(2)+2
       kbarr(jf)=nstr(jf)
       end do
       call mpadd(1)
       if (kcarr(2).eq.0)goto 230
       if (kcarr(1).eq.0)goto 220
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,iprar(2)+2
       kbarr(jf)=iprar(jf)
       end do
       call mpadd(0)
       if (kcarr(2).eq.0)goto 230
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
221    do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 230
       
       ians=ip-ians
       print *,'next time ians',ians
       goto 2041
       stop
220    do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       goto 221
230    print *,'required root=',ians
2301   icubr1=ians
       ix2=ip-3
       if (ix2.lt.10000)goto 205
       ncom(1)=0
       ncom(2)=2
       ncom(3)=ix2/10000
       ncom(4)=ix2-ncom(3)*10000
       goto 206
205    ncom(1)=0
       ncom(2)=1
       ncom(3)=ix2
206    do jf=1,iprar(2)+2
       karp(jf)=iprar(jf)
       end do
       do jf=1,ncom(2)+2 
       kard(jf)=ncom(jf)
       end do
       call mpkron(k)
       print *,'k',k
       if (k.eq.-1)goto 232
       call bwq5(ip,ians)
       print *,'sqrt-3',ians
       if (icubr1.lt.10000)goto 235
       icubar(1)=0
       icubar(2)=2
       icubar(3)=icubr1/10000
       icubar(4)=icubr1-icubar(3)*10000
       goto 236
235    icubar(1)=0
       icubar(2)=1
       icubar(3)=icubr1

       
236    do kk=1,2
       if (ians.lt.10000)goto 240
       ianar(1)=0
       ianar(2)=2
       ianar(3)=ians/10000
       ianar(4)=ians-ianar(3)*10000
       goto 241
240    ianar(1)=0
       ianar(2)=1
       ianar(3)=ians
241    do jf=3,icubar(2)+2
       karr(jf-2)=icubar(jf)
       end do
       ilen=icubar(2)
       do jf=3,ianar(2)+2
       kbarr(jf-2)=ianar(jf)
       end do
       ilen2=ianar(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=kk-1
       kbarr(2)=ilen3
       karr(1)=1
       do jf=2,icubar(2)+2
       karr(jf)=icubar(jf)
       end do
       
       
       call mpadd(0)
       do jf=1,kcarr(2)+2
       itempz(jf)=kcarr(jf)
       end do
       print *,'itempz',(itempz(jf),jf=1,itempz(2)+2)
       if (kk.eq.2)goto 242
       do jf=3,iprar(2)+2
       karr(jf-2)=iprar(jf)
       end do
       ilen=iprar(2)
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
       print *,'inv2',(inv2(jf),jf=1,inv2(2)+2)
242    do jf=3,itempz(2)+2
       karr(jf-2)=itempz(jf)
       end do
       ilen=itempz(2)
       do jf=3,inv2(2)+2
       kbarr(jf-2)=inv2(jf)
       end do
       ilen2=inv2(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       print *,'firkarr',(karr(jf),jf=1,ilen)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       karr(jf+2)=irrr(jf)
       end do
       karr(1)=itempz(1)
       karr(2)=irlen
       print *,'karrs',(karr(jf),jf=1,karr(2)+2)
       if (karr(1).eq.0)goto 243
       do jf=1,iprar(2)+2
       kbarr(jf)=iprar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
243    isum=0
       iling=karr(2)
       do jf=iling,1,-1
       isum=isum+karr(jf+2)*10000**(iling-jf)
       end do
       ksol(kk+1)=isum
       end do
       ksol(1)=icubr1
       icubrn=3
       icubr2=ksol(2)
       icubr3=ksol(3)
       goto 1000
232    icubrn=1
       goto 1000





202   stop










       

       
       
       
180    ipremy=3-iprem
       print *,'icong=',icong,'ind=',ipremy,'qrem',iqrem,'e',ie,'icon',icon
       print *,'icon',icon,'ib',ib,'iz',iz,'iconz',iconz
       
       if (iconz.eq.30)goto 196
183    ix=ixperm
       iconz=iconz+1
       ib=ibperm
       if (ib.lt.10000)goto 18301
       ibar(1)=0
       ibar(2)=2
       ibar(3)=ib/10000
       ibar(4)=ib-ibar(3)*10000
       goto 18302
18301  ibar(1)=0
       ibar(2)=1
       ibar(3)=ib
18302  if (iconz.eq.0)goto 5202
       goto 52
188    if (iqrem.eq.1)goto 196
       icong=icong+1
       if (icong.eq.3)goto 196
       iconz=0
       goto 183
196    icubrn=0
       
       goto 1000
1500   ipn=(ip+1)/6       
       iaas=iaa
       call sub516(ibprod,iaas,ipn,ip)
       ibprod2=ip-ibprod
       print *,'cube root=',ibprod,'or',ibprod2,'ind=1'
       ians=ibprod
       goto 2041
       stop
1000   return       
       end
       
       
       subroutine sub516(ibprod,iaas,ipn,ip)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       
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
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
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











