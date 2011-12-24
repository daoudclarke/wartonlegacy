      

      program bobnfs2
!     sieving for large numbers,sinle & double GNFS sieves     
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      dimension iprar(4),ifak(30),jfak(30),ifakn(30),jfakn(30)
      dimension iabp(50),iabpn(50),ity(50),littr(20),normar(20)
      dimension n(50)
      real krecarr,lgpr(10001),lglm 
      ktim=1
      
      
      do i=1,400000
      krecarr(i)=127.0
      kndx(i)=0
      end do
      irecnn=0
      kkmx=0
      
      open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
!      open(unit=4,file='nfsf2',access='sequential')
      
      read(3,5,rec=1) (ipr(i),i=1,65000)
5     format(65000i6)      
      print *,'ipr10000',ipr(10000),ipr(6500),'ipr2',ipr(2),ipr(3)
      print *,'no of primes<800001=',ipr(1)
      close(3)
! temporary use of files prefixed by x
      open (unit=2,file='xnfspar',access='direct',form=&
      'formatted',recl=488,status='old')


!      open (unit=2,file='nfspar',access='direct',form=&
!      'formatted',recl=488,status='old')
      read (2,4001,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),&
      (ia4(jf),jf=1,10),(ia3(jf),jf=1,10),(ia2(jf),jf=1,10),&
      (ia1(jf),jf=1,10),(ia0(jf),jf=1,10),(m1(jf),jf=1,10),&
      (m2(jf),jf=1,10),izz1,izz2,klim,izz4,izz5,izz6
4001  format (30i4,80i4,6i8)       
!      klim=3300
      rel=ipr(klim)
!      rel=ipr(4500)
      vlgls=log10(rel)
      lglm=2*vlgls
      print *,'lglm',lglm,'klim',klim,'izz1',izz1,'izz2',izz2,'izz4',izz4,&
      'izz5',izz5,'izz6',izz6
      print *,'m1',(m1(jf),jf=1,m1(2)+2),'m2',(m2(jf),jf=1,m2(2)+2)
      if (m2(3).ne.0)goto 731
      m2(3)=1
731   a=a      
      karr(1)=1
      do i=1,10000
      rel=i
      zlog(i)=log10(rel)
      end do
      karr(2)=3
      karr(3)=100
      karr(4)=0
      karr(5)=47
      kbarr(1)=1
      kbarr(2)=1
      kbarr(3)=100
      kbarr(4) =0
      kbarr(5)=54
      isora=1
      call mpadd(isora)
      print *,'add',(kcarr(i),i=1,kcarr(2)+2)
      
      
      
      kara(1)=0
      kara(2)=1
      kara(3)=14
      kara(4)=9839
      karb(1)=0
      karb(2)=1
      karb(3)=14
      call mpgcd
      print *,(karv(i),i=1,karv(2)+2)
      print *,'gcd=',(kard(i),i=1,kard(2) +2)
      
      kard(1)=0 
      kard(2)=2
      kard(3)=5
      kard(4)=3
      karp(1)=0
      karp(2)=2
      karp(3)=6
      karp(4)=7
      call mpkron(k)
      








      
      iblim=0
      iulim=2000
      lprx1=100000
      icon=0
      isivl=200000
      
      ia=1
      ib=1
400   ia=1
      ib=1

!      open (unit=2,file='pargkr',access='sequential',form=&
!      'unformatted')
!     goto 12344
!      read (2)kr
      kklim=7*(klim+300)
! temporary use of file prefixed by x
      open (unit=1,file='xnfsf1',access='direct',form=&
      'formatted',recl=kklim,status='old')


!      open (unit=1,file='nfsf1',access='direct',form=&
!      'formatted',recl=kklim,status='old')
      do iz2=1,9
      print *,'iz2',iz2
      iconff=0
      read (1,12347,rec=iz2)(kr(iz,iz2),iz=1,klim)
      do jf=1,klim
      if (kr(jf,iz2).eq.9999999)goto 730
      iconff=iconff+1
730   end do
      print *,'iz2',iz2,'no. of ideals',iconff


      end do
!      stop
!      print *,((kr(iz,iz2),iz=1,4500),iz2=1,2)
!      klim=2000
      print *,kr(2999,1),kr(2999,2),kr(3000,1),kr(3000,2)
      print *,'ipr2999',ipr(2999),'ipr3000',ipr(3000),ipr(3001),ipr(3002)
      print *,'kr3001',kr(3001,1),kr(3001,2),'kr3002',kr(3002,1),kr(3002,2)
      print *,'kr3003',kr(3003,1),kr(3003,2),'kr3004',kr(3004,1),kr(3004,2)
      print *,'kr3966',(kr(3966,jf),jf=1,10)
      print *,'kr4499',(kr(4499,jf),jf=1,10)
!      close (unit=2)
      
      icc=0
      do iz=1,2000 
!      do iz=1,4500
      if (kr(iz,1).eq.9999999)goto 12378
      icc=icc+1
!      print *,'iz',iz,'icc',icc
12378 end do      
      print *,'icc',icc,'klim',klim
      
      icc=0
      do iz=1,klim
      if ((kr(iz,1).ne.9999999).and.(kr(iz,2).eq.9999999))goto 12346
      goto 12377
12346 icc=icc+1
!      print *,'iz',iz,'kr12',kr(iz,1),kr(iz,2),'ipriz',ipr(iz)
12377 end do      
!      print *,'icc=',icc
!      print *,'ipr3000',ipr(3000),kr(3000,1),kr(3000,2),kr(3000,3)
      
      goto 12349
      
      
12344 do iz=1,klim
      do iz2=1,10
      kr(iz,iz2)=9999999
      end do
      end do
      
      
      
      
      kib=1
      ick=0
!     change parameter
      do iz=2,klim
!      print *,'iz',iz
      ikk=1
      do jf=3,ia0(2)+2
      karr(jf-2)=ia0(jf)
      end do
      ilen=ia0(2)
      if (ipr(iz).lt.10000)goto 8000
      kbarr(1)=int(ipr(iz)/10000)
      kbarr(2)=ipr(iz)-kbarr(1)*10000
      ilen2=2
      goto 8002
      
8000  kbarr(1)=ipr(iz)
      ilen2=1
8002  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.ne.0)goto 402
      kr(iz,ikk)=0
      ikk=ikk+1
402   do kia=1,ipr(iz)-1
5005  call sieve(kia,kib,irecnn,kkmx,1,klim)
      do jf=3,norma(2)+2
      karr(jf-2)=norma(jf)
      end do
      ilen=norma(2)
      if (ipr(iz).lt.10000)goto 5002
      ilen2=2
      kbarr(1)=int(ipr(iz)/10000)
      kbarr(2)=ipr(iz)-kbarr(1)*10000
      goto 5003
5002  ilen2=1
      kbarr(1)=ipr(iz)
5003  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.ne.0)goto 500
      kr(iz,ikk)=kia
      ikk=ikk+1
      if (ikk.eq.11)goto 510
500   end do
      print *,'ikk',ikk
510   end do
!     write solutions to file
      write (2)kr
12345 close (unit=2)
!      open (unit=1,file='rargkr',access='direct',form=&
!     'formatted',recl=27000,status='old')
      do iz2=1,10
      write(1,12347,rec=iz2)(kr(iz,iz2),iz=1,klim)
      end do
12347 format(12000i7)
      close (unit=1)
!     change parameter
12349 do i=2,klim
      rel=ipr(i)
      lgpr(i)=log10(rel)
      end do
      
      print *,'kr13',kr(13,1)
      print *,'kr3',kr(3,1)
      print *,'kr19',kr(19,1)
      
      
      print *,'ipr',(ipr(jf),jf=1,200)
      
      call sieve(1,1,irecnn,kkmx,10,klim)
      
      goto 460
      print *,'leadc con',leadc(1)
      if (leadc(1).eq.0)goto 460
      print *,'leadc',(leadc(jf),jf=2,leadc(1)+1)
      call sieve(1,1,irecnn,kkmx,7,klim)
      do jf=1,ijfak(1)+2
      ifak(jf)=ijfak(jf)
      ifakn(jf)=ijfakn(jf)
      end do
      
      call sieve(1,1,irecnn,kkmx,8,klim)
      do jf=1,ijfak(1)+2
      jfak(jf)=ijfak(jf)
      jfakn(jf)=ijfakn(jf)
      end do
      print *,'ifak',(ifak(jf),jf=1,ifak(1)+2)
      print *,'ifakn',(ifakn(jf),jf=1,ifakn(1)+2)
      print *,'jfak',(jfak(jf),jf=1,jfak(1)+2)
      print *,'jfakn',(jfakn(jf),jf=1,jfakn(1)+2)
      irecnn=irecnn+1
      do jf=1,20
      littr(jf)=0
      normar(jf)=0
      end do
      do jf=1,ifak(1)
      iabp(jf)=ifak(jf+2)
      iabpn(jf)=ifakn(jf+2)
      ity(jf)=1
      end do
      nin=ifak(1)
      ity(nin)=ity(nin)+2*ifak(2)
      do jf=1,jfak(1)
      iabp(nin+jf)=jfak(jf+2)
      iabpn(nin+jf)=jfakn(jf+2)
      ity(nin+jf)=2
      end do
      icur=nin+jfak(1)
      ity(icur)=ity(icur)+2*jfak(2)
      print *,'icur',icur,'iabp',(iabp(jf),jf=1,icur),'iabpn',&
      (iabpn(jf),jf=1,icur),'ity',(ity(jf),jf=1,icur)
      kia=0
      kib=0
      littr(1)=m1(1)
!      write(4,*)irecnn,kia,kib,icur,(littr(j1),j1=1,20),(normar(j2)&
!      ,j2=1,20)
      
!      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      do kkx=2,ipr(klim)
      kib=kkx
      
      
      do jf=2,klim
      if (irecnn.eq.3000)goto 51
      if (kib.eq.ipr(jf))goto 14
      if (kib.gt.ipr(jf))goto 50
      goto 15
50    end do
      goto 19
      

      
15    if (kib.lt.10000)goto 10
      norma(1)=0
      norma(2)=2
      norma(3)=kib/10000
      norma(4)=kib-norma(3)*10000
      goto 11
10    norma(1)=0      
      norma(2)=1
      norma(3)=kib
11    call sieve(0,kib,irecnn,kkmx,9,klim)
      print *,'ijfak',(ijfak(jk),jk=1,ijfak(1)+2),'kib',kib
      print *,'ijfakn',(ijfakn(jk),jk=1,ijfakn(1)+2)
      
      goto 21
14    if (kib.gt.ipr(2000))goto 19
      do jk=2,leadc(1)+1      
      if (kib.eq.leadc(jk))goto 20
      end do
      if (kr(jf,1).eq.9999999)goto 19
20    a=a      
      ijfak(1)=1
      ijfak(2)=0
      ijfak(3)=kib
      ijfakn(1)=1
      ijfakn(2)=0
      ijfakn(3)=1
21    ind1=3
      ind2=3
      ind3=1
      print *,'ijfak',(ijfak(jk),jk=1,ijfak(1)+2)
      
      if (ijfak(2).eq.1)goto 19
      
31    if (ifak(ind1).lt.ijfak(ind2))goto 30
      if (ifak(ind1).gt.ijfak(ind2))goto 35
      iabp(ind3)=ifak(ind1)
      iabpn(ind3)=ifakn(ind1)+ijfakn(ind2)
      ity(ind3)=1
      ind1=ind1+1
      ind2=ind2+1
      ind3=ind3+1
      if (ind1.gt.ifak(1)+2-ifak(2))goto 36
      if (ind2.gt.ijfak(1)+2)goto 37
      goto 31
30    iabp(ind3)=ifak(ind1)
      iabpn(ind3)=ifakn(ind1)
      ity(ind3)=1
      ind1=ind1+1
      ind3=ind3+1
      if (ind1.gt.ifak(1)+2-ifak(2))goto 36
      goto 31 
35    iabp(ind3)=ijfak(ind2)
      iabpn(ind3)=ijfakn(ind2)
      ity(ind3)=1
      ind2=ind2+1
      ind3=ind3+1
      if (ind2.gt.ijfak(1)+2)goto 37
      goto 31
36    if (ind2.gt.ijfak(1)+2)goto 38
      do jf=ind2,ijfak(1)+2
      iabp(ind3)=ijfak(jf)
      iabpn(ind3)=ijfakn(jf)
      ity(ind3)=1
      ind3=ind3+1
      end do
      goto 38
37    do jf=ind1,ifak(1)+2-ifak(2)
      iabp(ind3)=ifak(jf)
      iabpn(ind3)=ifakn(jf)
      ity(ind3)=1
      ind3=ind3+1
      end do
38    ind1=3
      ind2=3
      if (ifak(2).eq.0)goto 41
      iabp(ind3)=1
      iabpn(ind3)=1
      ity(ind3)=3
      ind3=ind3+1


!     end of ijfak loop
      
41    if (jfak(ind1).lt.ijfak(ind2))goto 40
      if (jfak(ind1).gt.ijfak(ind2))goto 45
      iabp(ind3)=jfak(ind1)
      iabpn(ind3)=jfakn(ind1)+5*ijfakn(ind2)
      ity(ind3)=2
      ind1=ind1+1
      ind2=ind2+1
      ind3=ind3+1
      if (ind1.gt.jfak(1)+2-jfak(2))goto 46
      if (ind2.gt.ijfak(1)+2)goto 47
      goto 41
40    iabp(ind3)=jfak(ind1)
      iabpn(ind3)=jfakn(ind1)
      ity(ind3)=2
      ind1=ind1+1
      ind3=ind3+1
      if (ind1.gt.jfak(1)+2-jfak(2))goto 46
      goto 41 
45    iabp(ind3)=ijfak(ind2)
      iabpn(ind3)=ijfakn(ind2)*5
      ity(ind3)=2
      ind2=ind2+1
      ind3=ind3+1
      if (ind2.gt.ijfak(1)+2)goto 47
      goto 41
46    if (ind2.gt.ijfak(1)+2)goto 48
      do jf=ind2,ijfak(1)+2
      iabp(ind3)=ijfak(jf)
      iabpn(ind3)=ijfakn(jf)*5
      ity(ind3)=2
      ind3=ind3+1
      end do
      goto 48
47    do jf=ind1,jfak(1)+2-jfak(2)
      iabp(ind3)=jfak(jf)
      iabpn(ind3)=jfakn(jf)
      ity(ind3)=2
      ind3=ind3+1
      end do
48    if (jfak(2).eq.0)goto 49
      iabp(ind3)=1
      iabpn(ind3)=1
      ity(ind3)=4
      ind3=ind3+1
49    icur=ind3-1      
      irecnn=irecnn+1
!      write(4,*)irecnn,kia,kib,icur,(littr(j1),j1=1,20),(normar(j2)&
!      ,j2=1,20)
      
!      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
!      print *,'hit no',irecnn,'kib',kib

19    end do
      stop

51    a=a



460   if(ia.eq.1)goto 300      
      if (ib.eq.1)goto 300
202   iabia=abs(ia)      
      iabib=abs(ib)
      if (iabib.gt.iabia)goto 220
      ibig =iabia
      little =iabib
      call subgcd(ibig,little,igcd2)
      if (igcd2.eq.1)goto 300
      ia=ia+1
      if (ia.gt.isivl)goto 3404
      goto 460
220   ibig=iabib
      little=iabia
      call subgcd(ibig,little,igcd2)
      if (igcd2.eq.1)goto 300
      ia =ia+1
      if (ia.gt.isivl)goto 3404
      goto 460
      


300   kia=ia
      kib=ib
      
      kndx(ia*2-1)=1
      

3007  kia=ia *(-1)
      
      kndx(ia*2)=1
3004  ia=ia+1
      if (ia.gt.isivl)goto 3404
      goto 460

3404  a=a     
      
!      print *,'phase 1 recno=',ib
!      if (ib.gt.10)goto 420
     if (ib.gt.isivl)goto 420
3401  format(4000(i1,f9.5))      
!      goto 16159
      call fscomp(kib,vlgls,1,klim)
!      do iz=1,klim/3
      do iz=2,klim
!      print *,'iz',iz
!      do iz=2,4500
      do jf=3,m2(2)+2
      karr(jf-2)=m2(jf)
      end do
      ilen=m2(2)
      if (ipr(iz).lt.10000)goto 1800
      iprar(1)=0
      iprar(2)=2
      iprar(3)=ipr(iz)/10000
      iprar(4)=ipr(iz)-iprar(3)*10000
      goto 1801
1800  iprar(1)=0
      iprar(2)=1
      iprar(3)=ipr(iz)
1801  do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      do jf=1,irlen
      karr(jf+2)=irrr(jf)
      end do
      karr(2)=irlen
      karr(1)=m2(1)
      if (m2(1).eq.0)goto 1802
      do jf=1,iprar(2)+2
      kbarr(jf)=iprar(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
1802  do jf=1,karr(2)+2
      karb(jf)=karr(jf)
      end do
      do jf=1,iprar(2)+2
      kara(jf)=iprar(jf)
      end do
      
      call mpgcd
      
      do jf=3,karv(2)+2
      karr(jf-2)=karv(jf)
      end do
      ilen=karv(2)
      if (kib.lt.10000)goto 1804 
      kbarr(1)=kib/10000
      kbarr(2)=kib-kbarr(1)*10000
      ilen2=2
      goto 1805
1804  kbarr(1)=kib
      ilen2=1
1805  call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,m1(2)+2
      kbarr(jf-2)=m1(jf)
      end do
      ilen2=m1(2)
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
      if (irlen.eq.0)goto 1602
      do jf=1,irlen
      karr(jf+2)=irrr(jf)
      end do
      karr(1)=m1(1)
      karr(2)=irlen
      if (karr(1).eq.0)goto 1810
      do jf=1,iprar(2)+2
      kbarr(jf)=iprar(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
1810  if (karr(2).eq.1)goto 1812
      if (karr(2).eq.0)goto 1602
      markb=karr(3)*10000+karr(4)
      goto 1603
1812  markb=karr(3)
      goto 1603
1602  markb=0
      
      goto 1604

1603  markbn=ipr(iz)-markb
1604  markm=markb
      markbn=ipr(iz)-markb
16041 if ((markb.gt.200000).and.(markbn.gt.200000))goto 1615
      if ((markb.gt.200000).or.(markb.eq.0))goto 1606
      if (kndx(markb*2-1).ne.1)goto 1606
      krecarr(markb*2-1)=krecarr(markb*2-1)-lgpr(iz)
1606  if ((markbn.gt.200000).or.(markbn.eq.0))goto 16062
      if (kndx(markbn*2).ne.1)goto 16062
      krecarr(markbn*2)=krecarr(markbn*2)-lgpr(iz)
16062  markb=markb+ipr(iz)
      markbn=markbn+ipr(iz)
      goto 16041

1615  end do
16159 call fscomp(kib,vlgls,2,klim)
      

        
        





!     change parameter      
!      do iz=2,4500
      do iz=1,klim
      if (ipr(iz).lt.10000) goto 8012
      iprar(2)=2
      iprar(3)=int(ipr(iz)/10000)
      iprar(4)=ipr(iz)-iprar(3)*10000
      goto 8014
8012  iprar(2)=1      
      iprar(3)=ipr(iz)
8014  do iz2=1,5
      if (kr(iz,iz2).eq.9999999)goto 615
      if (kr(iz,iz2).eq.0)goto 602
      if (kib.lt.10000)goto 8004
      karr(1)=int(kib/10000)
      karr(2)=kib-karr(1)*10000
      ilen=2
      goto 8006
      
      
8004  karr(1)=kib
      ilen=1
8006  if (kr(iz,iz2).lt.10000)goto 8008      
      kbarr(1)=int(kr(iz,iz2)/10000)
      kbarr(2)=kr(iz,iz2)-kbarr(1)*10000
      ilen2=2
      goto 8010


8008  kbarr(1)=kr(iz,iz2)
      ilen2=1
8010  call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen2=iprar(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 602
      if (irlen.eq.1)goto 603
      markb=irrr(1)*10000+irrr(2)
      goto 604
602   markb=0
      goto 604
603   markb=irrr(1)
      markbn=ipr(iz)-markb
604   markm=markb
      markbn=ipr(iz)-markb
6041  if ((markb.gt.200000).and.(markbn.gt.200000))goto 615
      if ((markb.gt.200000).or.(markb.eq.0))goto 606
      if (kndx(markb*2-1).ne.1)goto 606
      krecarr(markb*2-1)=krecarr(markb*2-1)-lgpr(iz)
606   if ((markbn.gt.200000).or.(markbn.eq.0))goto 6062
      if (kndx(markbn*2).ne.1)goto 6062
      krecarr(markbn*2)=krecarr(markbn*2)-lgpr(iz)
6062  markb=markb+ipr(iz)
      markbn=markbn+ipr(iz)
      goto 6041
615   end do
      end do
      
      
6151  do kz=1,400000
      if (kndx(kz).ne.1)goto 710
!     was 116.0       
      if (krecarr(kz).gt.vlgls)goto 710
      isgn=mod(kz,2)
      kia=int(kz/2)
      if (isgn.eq.1)goto 699
      kia =kia*(-1)
      goto 702
               
699   kia=kia+1
      

702   kpoin=0
      call sieve(kia,kib,irecnn,kkmx,kpoin,klim)
!      if (irecnn.ge.klim*2+klim/4)goto 420
      if (irecnn.ge.4200)goto 420
!      if (irecnn.ge.4200)goto 420
      if (kpoin.ne.3)goto 710
      kndx(kz)=3
710   end do
      ia=1
      ib=ib+1
! temp inst.
      
      if (ib.gt.isivl)goto 420
      if (ib.gt.8243)goto 420
      do jf=1,400000
      kndx(jf)=0
      krecarr(jf)=127.0
      end do
      
      
      goto 460
      

      






420   print *,'max no of primes=',kkmx
      
      print *,'irecnn=',irecnn
      write (2,4001,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),&
      (ia4(jf),jf=1,10),(ia3(jf),jf=1,10),(ia2(jf),jf=1,10),&
      (ia1(jf),jf=1,10),(ia0(jf),jf=1,10),(m1(jf),jf=1,10),&
      (m2(jf),jf=1,10),irecnn,izz2,klim,izz4,izz5,izz6
      close(unit=2)
      
      
      end
      
      subroutine sieve(kia,kib,irecnn,kkmx,kpoin,klim)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(20),litd(20),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(20),normar(20)
      real krecarr
!     change parameter      
!      lprx1=ipr(klim/3)                                      
      lprx1=ipr(klim)                                      
!      lprx1=ipr(2000)
      nnsw=1
      icur=0
      larp(1)=int(lprx1/10000)
      larp(2)=lprx1 -larp(1) *10000
      
      
      do i=1,50
      ity(i)=0
      iabp(i)=0
      iabpn(i)=0
      end do
      if (kpoin.eq.10)goto 10
      if (kpoin.eq.7)goto 30
      if (kpoin.eq.8)goto 40
      if (kpoin.eq.9)goto 1300
!      print *,'leadc',(leadc(jf),jf=1,20)
      
      itemp=abs(kia)
      
      
      
      if(itemp.lt.10000)goto 1002
      ilen=2
      ilen2=2
      karr(1)=int(itemp/10000)
      ktemp(1,1)=karr(1)
      karr(2)=itemp-karr(1)*10000
      ktemp(1,2)=karr(2)
      ktempl(1)=2
      goto 1004
1002  ilen=1 
      karr(1)=itemp
      ktemp(1,1)=itemp
      ilen2=ilen
      ktempl(1)=ilen
      
1004  do i=1,ilen
      kbarr(i)=karr(i)
      end do
      do j=2,5
      call mpmul(ilen,ilen2,ilen3)
      do i=1,ilen3
      kbarr(i)=kcarr(i)
      ktemp(j,i)=kcarr(i)
      end do
      ktempl(j)=ilen3
      
      ilen2 =ilen3
      end do
      
      itemp=abs(kib)
      if(itemp.lt.10000)goto 1006
      karr(1)=int(itemp/10000)
      ktemp2(1,1)=karr(1)
      karr(2)=itemp-karr(1)*10000
      ktemp2(1,2)=karr(2)
      ilen=2
      ktempl2(1)=ilen
      goto 1008
1006  karr(1)=itemp
      ilen=1
      ktemp2(1,1)=itemp
      ktempl2(1)=ilen
1008  do i=1,ilen
      kbarr(i)=karr(i)
      end do
      do i=1,20
      norma(i)=0
      end do
      
      ilen2 =ilen
      do j=2,5
      call mpmul(ilen,ilen2,ilen3)
      
      do i=1,ilen3
      kbarr(i)=kcarr(i)
      ktemp2(j,i)=kbarr(i)
      end do
      ktempl2(j)=ilen3
      ilen2 =ilen3
      end do
      
      if(kia.lt.0)goto 1009
      ksgn=0
      goto 1012
1009  ksgn=1
1012  if(ia5(1).eq.1)goto 1014
      isgn=0
      goto 1016
1014  isgn=1
1016  ilen=ia5(2)
      if (ia5(2).eq.0)goto 2001
      do i=1,ilen
      karr(i)=ia5(i+2)
      end do
      
      do i=1,ktempl(5)
      kbarr(i)=ktemp(5,i)
      
      end do
      
      ilen2=ktempl(5)
      call mpmul(ilen,ilen2,ilen3)
      norma(2)=ilen3
      do i=1,ilen3
      norma(i+2)=kcarr(i)
      end do
      norma(1)=mod((ksgn*5 +isgn),2)
      do jf=1,norma(2)+2
      ipol1(1,jf)=norma(jf)
      end do
      goto 20011
2001  ipol1(1,1)=0    
      ipol1(1,2)=0
20011 if (ia4(2).eq.0)goto 2002
      do i=1,ia4(2)
      karr(i)=ia4(i+2)
      end do
      
      ilen=ia4(2)
      do i=1,ktempl(4)
      kbarr(i)=ktemp(4,i)
      end do
      ilen2=ktempl(4)
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do i=1,ilen3
      karr(i)=kcarr(i)
      end do
      ilen2=ktempl2(1)
      do i=1,ktempl2(1)
      kbarr(i)=ktemp2(1,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      if(kia.lt.0)goto 1018
      ksgn2=0
      goto 1020
1018  ksgn2=0
1020  kbarr(1)=mod(ia4(1) +ksgn*4 ,2)
      do jf=1,kbarr(2)+2
      ipol1(2,jf)=kbarr(jf)
      end do
      
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do

      
      call mpadd(0)
      

      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      goto 20021
2002  ipol1(2,1)=0      
      ipol1(2,2)=0
20021 ilen=ia3(2)
      if (ia3(2).eq.0)goto 2003
      do i=1,ilen
      karr(i)=ia3(i+2)
      end do
      ilen2=ktempl(3)
      do i=1,ilen2
      kbarr(i)=ktemp(3,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do i=1,ilen3
      karr(i)=kcarr(i)
      end do
      ilen2=ktempl2(2)
      do i=1,ilen2
      kbarr(i)=ktemp2(2,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia3(1)+ksgn*3,2)
      do jf=1,kbarr(2)+2
      ipol1(3,jf)=kbarr(jf)
      end do
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      goto 20031
2003  ipol1(3,1)=0
      ipol1(3,2)=0
      
20031 ilen=ia2(2)
      if (ia2(2).eq.0)goto 2004
      do i=1,ilen
      karr(i)=ia2(i+2)
      end do
      ilen2 =ktempl(2)
      do i=1,ilen2
      kbarr(i)=ktemp(2,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do i=1,ilen
      karr(i)=kcarr(i)
      end do
      ilen2 =ktempl2(3)
      do i=1,ilen2
      kbarr(i)=ktemp2(3,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia2(1) +ksgn*2 ,2)
      do jf=1,kbarr(2)+2
      ipol1(4,jf)=kbarr(jf)
      end do




      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      goto 20041
2004  ipol1(4,1)=0
      ipol1(4,2)=0
20041 ilen=ia1(2)
      if (ia1(2).eq.0)goto 2005
      do i=1,ilen
      karr(i)=ia1(i+2)
      end do
      ilen2=ktempl(1)
      do i=1,ilen2
      kbarr(i)=ktemp(1,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do i=1,ilen
      karr(i)=kcarr(i)
      end do
      ilen2=ktempl2(4)
      do i=1,ilen2
      kbarr(i)=ktemp2(4,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia1(1) +ksgn ,2)
      do jf=1,kbarr(2)+2
      ipol1(5,jf)=kbarr(jf)
      end do
      
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      goto 20051
2005  ipol1(5,1)=0
      ipol1(5,2)=0
      
20051 ilen=ia0(2)
      if (ia0(2).eq.0)goto 2006
      do i=1,ilen
      karr(i)=ia0(i+2)
      end do
      ilen2=ktempl2(5)
      do i=1,ilen2
      kbarr(i)=ktemp2(5,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia0(1) ,2)
      do jf=1,kbarr(2)+2
      ipol1(6,jf)=kbarr(jf)
      end do
      
      do i=1,norma(2) +2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      goto 20061
2006  ipol1(6,1)=0
      ipol1(6,2)=0
20061 if (kpoin.eq.1)goto 1500 
      

      ilen=ktempl(1)
      do i=1,ilen
      karr(i)=ktemp(1,i)
      end do
      ilen2=m2(2)
      do i=1,ilen2
      kbarr(i)=m2(i+2)
      end do
      call mpmul(ilen,ilen2,ilen3)
      litt(2)=ilen3
      do i=1,ilen3
      litt(i+2)=kcarr(i)
      end do
      litt(1)=mod(ksgn+m2(1),2)
      
      
      ilen=ktempl2(1)
      do i=1,ilen
      karr(i)=ktemp2(1,i)
      end do
      ilen2 =m1(2)
      do i=1,ilen2
      kbarr(i)=m1(i+2)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(m1(1),2)
      do i=1,litt(2) +2
      karr(i)=litt(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2)  +2
      litt(i)=kcarr(i)
      end do
      if (litt(2).eq.0)goto 1500
      if((norma(2).eq.1).and.(norma(3).lt.2))goto 1500
      
      

!     finish of calculation of number to be sieved
      do i=1,kcarr(2)+2
      litz(i)=litt(i)
      litd(i)=litt(i)
      end do
      goto 11
10    do jf=1,ia5(2)+2      
      litz(jf)=ia5(jf)
      litd(jf)=ia5(jf)
      end do
      goto 11
30    do jf=1,m1(2)+2
      litz(jf)=m1(jf)
      litd(jf)=m1(jf)
      end do
      goto 11
40    do jf=1,ia0(2)+2
      norma(jf)=ia0(jf)
      
      end do
      goto 1300
11    a=a
      litd(1)=0
      
      iflim=10000
      iab=1
      icur=0
      
      do i=1,50
      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
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
1300  if (kpoin.eq.10)goto 1500
       
      nnsw=0
      if (kpoin.eq.7)goto 1500
      nnsw=1
      do i=1,norma(2)+2
      litt(i)=norma(i)
      litd(i)=norma(i)
      end do
      print *,'bigstop'
      mprx1=ipr(klim)
      larp(1)=mprx1/10000
      larp(2)=mprx1-larp(1)*10000
      
      
      
      
      litd(1)=0
      print *,(litd(i),i=1,litd(2) +2),'stpge','ia=',kia,'kib=',kib,&
      irecnn
      iab=1
      do i=1,50
      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1306  iab=iab+1
      iprx1= ipr(iab)
      if (leadc(1).eq.0)goto 22
      do jjf=2,leadc(1)+1
      if (leadc(jjf).eq.iprx1)goto 1307
      end do
!      if (ipr(iab).eq.2)goto 1307
!      if (ipr(iab).eq.5)goto 1307
!      if (ipr(iab).eq.7)goto 1307
!      if (ipr(iab).eq.17)goto 1307
!      if (ipr(iab).eq.317)goto 1307
22     if (kr(iab,1).eq.9999999)goto 1306
1307  a=a
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
      ilen =litd(2)
      if (iprx1.lt.10000)goto 1368
      kbarr(1)=iarp(1)
      kbarr(2)=iarp(2)
      ilen2=2
      goto 1370
1368  kbarr(1)=iprx1
      ilen2=1
1370  call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      
!      if ((icont.eq.0).and.(ipqt(1).eq.1))goto 1386
      if(irlen.gt.0)goto 1306
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
      if ((icont.eq.1).and.(ipqt(1).eq.1))goto 1386
      goto 1366
1384  iabpn(icur)=iabpn(icur)+1
      if ((icont.eq.1).and.(ipqt(1).eq.1))goto 1386
      goto 1366
      
1386  a=a
      
      nnsw=0
      if (kpoin.eq.8)goto 1500
      if (kpoin.eq.9)goto 1500
      irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kib,(litt(i),i=1,litt(2) +2)
      do i=1,20
      normar(i)=0
      littr(i)=0
      end do
      do i=1,litz(2) +2
      littr(i)=litz(i)
      end do
      do i=1,norma(2)+2
      normar(i)=norma(i)
      end do
! note write insts. removed      
!      write(4,*)irecnn,kia,kib,icur,(littr(j1),j1=1,20),(normar(j2)&
!      ,j2=1,20)
      
!      write(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      kpoin=3
      

      if (irecnn.lt.160000)goto 1500
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
      
      goto 1500
1400  if (litd(2).ne.iprex(2))goto 1402
      do i=3,litd(2)+2
      if (litd(i).ne.iprex(i))goto 1402
      end do
      goto 1404
1402  if (litd(2).gt.2)goto 3020
      if (litd(2).eq.1)goto 1366
      if (litd(3).lt.larp(1))goto 1366
      if (litd(3).gt.larp(1))goto 3020
      if (litd(4).le.larp(2))goto 1366


3020  do i=3,litd(2)+2
      mnum(i-2)=litd(i)
      end do
      inlen=litd(2)
      
      call mpprime(icorp,inlen)
      if(icorp.eq.1)goto 1500
      do i=1,litd(2)+2
      iprex(i)=litd(i)
      end do
1404  ipfar(1)=ipfar(1)+2
      goto 1366
      call mpprime(icorp,inlen)
1500  if (kpoin.eq.10)goto 14
      if (kpoin.eq.7)goto 31
      if (kpoin.eq.8)goto 31
      if (kpoin.eq.9)goto 31
      goto 15
14    if (icur.eq.0)goto 16
      leadc(1)=icur
      do jf=1,icur
      leadc(jf+1)=iabp(jf)
      end do
      goto 15
31    ijfak(1)=icur
      ijfakn(1)=icur
      if (icur.eq.0)goto 32
      do jf=1,icur
      ijfak(jf+2)=iabp(jf)
      ijfakn(jf+2)=iabpn(jf)
      end do
32    if (nnsw.eq.0)goto 33
      ijfak(1)=ijfak(1)+1
      ijfak(2)=1
      ijfak(icur+3)=1
      ijfakn(icur+3)=1
      goto 15
33    ijfak(2)=0
      ijfakn(2)=0
      goto 15






16    leadc(1)=0
15    a=a
!      print *,'nnsw',nnsw
      return
      end
      
      subroutine fscomp(kib,vlgls,kjk,klim)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      dimension ixold(6,50),ixnew(6,50),ipol2(6,50),ispol(6,50)
      dimension ipowr(5,6,50),itempb(20)
      real krecarr
      if (kjk.eq.2)goto 60
      if (kib.lt.10000)goto 62
      
      karr(1)=kib/10000
      karr(2)=kib-karr(1)*10000
      ilen=2
      goto 63
62    karr(1)=kib
      ilen=1
63    do jf=3,m1(2)+2
      kbarr(jf-2)=m1(jf)
      end do
      ilen2=m1(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      itempb(jf+2)=kcarr(jf)
      end do
      itempb(2)=ilen3
      itempb(1)=mod(m1(1)+1,2)
      do jf=1,itempb(2)+2
      karr(jf)=itempb(jf)
      end do
      do jf=1,m2(2)+2
      kbarr(jf)=m2(jf)
      end do
      iconz=-1
      do i=1,200000
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      iconz=iconz+2
      if (kndx(iconz).ne.1)goto 64
      rel=zlog(karr(3))+4*karr(2)-4
      krecarr(iconz)=rel
64    end do
      do jf=1,itempb(2)+2
      karr(jf)=itempb(jf)
      end do
      iconz=0
      do i=1,200000
      call mpadd(1)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      iconz=iconz+2
      if (kndx(iconz).ne.1)goto 66
      rel=zlog(karr(3))+4*karr(2)-4
      krecarr(iconz)=rel
66    end do
      goto 100
















60    do i=1,5
      do j=1,6
      ipowr(i,j,1)=0
      ipowr(i,j,2)=0
      end do
      end do
      do i=1,5
      do j=i,6
      ipowr(i,j,1)=0
      ipowr(i,j,2)=1
      end do
      end do
      ipowr(1,1,3)=1
      ipowr(1,2,3)=5
      ipowr(1,3,3)=10
      ipowr(1,4,3)=10
      ipowr(1,5,3)=5
      ipowr(1,6,3)=1
      ipowr(2,2,3)=1
      ipowr(2,3,3)=4
      ipowr(2,4,3)=6
      ipowr(2,5,3)=4
      ipowr(2,6,3)=1
      ipowr(3,3,3)=1
      ipowr(3,4,3)=3
      ipowr(3,5,3)=3
      ipowr(3,6,3)=1
      ipowr(4,4,3)=1
      ipowr(4,5,3)=2
      ipowr(4,6,3)=1
      ipowr(5,5,3)=1
      ipowr(5,6,3)=1











      do ibig=1,2
      if (ibig.eq.2)goto 50
      call sieve(1,kib,irecnn,kkmax,1,klim)
!      print *,'ipol12',(ipol1(2,jf),jf=1,ipol1(2,2)+2)
      
      do i=1,6
      do jf=1,ipol1(i,2)+2
      ispol(i,jf)=ipol1(i,jf)
      end do
      end do
      goto 52
50    ipowr(1,2,1)=1
      ipowr(1,4,1)=1
      ipowr(1,6,1)=1
      ipowr(2,3,1)=1
      ipowr(2,5,1)=1
      ipowr(3,4,1)=1
      ipowr(3,6,1)=1
      ipowr(4,5,1)=1
      ipowr(5,6,1)=1
      do i=1,6
      do jf=1,ispol(i,2)+2
      ipol1(i,jf)=ispol(i,jf)
      end do
      end do


52   do jf=1,ipol1(6,2)+2 
     ixold(6,jf)=ipol1(6,jf)
     end do
     do kk1=1,5
     do jf=1,ipol1(6,2)+2
     ipol2(6,jf)=ipol1(6,jf)
     end do
     do i=1,5
     ipol2(i,1)=0
     ipol2(i,2)=0
     end do
     do ii=kk1,5
     do jj=kk1,6
     if (ipol1(ii,2).eq.0)goto 30
     do jf=3,ipol1(ii,2)+2
     karr(jf-2)=ipol1(ii,jf)
     end do
     ilen=ipol1(ii,2)
     do jf=3,ipowr(ii,jj,2)+2
     kbarr(jf-2)=ipowr(ii,jj,jf)
     end do
     ilen2=ipowr(ii,jj,2)
     call mpmul(ilen,ilen2,ilen3)
     do jf=1,ilen3
     karr(jf+2)=kcarr(jf)
     end do
     karr(1)=mod(ipol1(ii,1)+ipowr(ii,jj,1),2)
     karr(2)=ilen3
     do jf=1,ipol2(jj,2)+2
     kbarr(jf)=ipol2(jj,jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     ipol2(jj,jf)=kcarr(jf)
     end do
30   end do
!     print *,'iipol21',(ipol2(ii,jf),jf=1,ipol2(ii,2)+2),'ii',ii,'ibig',ibig
     end do
     
     do jj=kk1+1,6
     do jf=1,ipol2(jj,2)+2
     karr(jf)=ipol2(jj,jf)
     end do
     do jf=1,ipol1(jj,2)+2
     kbarr(jf)=ipol1(jj,jf)
     end do
     call mpadd(1)
     do jf=1,kcarr(2)+2
     ipol1(jj,jf)=kcarr(jf)
     end do
     end do
     do jf=1,ipol1(6,2)+2
     ixold(6-kk1,jf)=ipol1(6,jf)
     end do
!    end kk1 loop
!     print *,'kk1',kk1,'ibig',ibig
     end do
     
     iconz=-2+ibig
     do jz=1,200000
     do kz=6,2,-1
     do jf=1,ixold(kz,2)+2
     karr(jf)=ixold(kz,jf)
     end do
     do jf=1,ixold(kz-1,2)+2
     kbarr(jf)=ixold(kz-1,jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     ixnew(kz,jf)=kcarr(jf)
     end do
     end do
     do kz=2,6
     if (ixnew(kz,2).eq.0)goto 15
     do jf=1,ixnew(kz,2)+2
     ixold(kz,jf)=ixnew(kz,jf)
     end do
     goto 16
15   ixold(kz,1)=0
     ixold(kz,2)=0
16   end do
     iconz=iconz+2
     if (kndx(iconz).ne.1)goto 14
     if (krecarr(iconz).gt.vlgls-0.5)goto 19
     if (ixnew(6,2).eq.0)goto 17
     rel=zlog(ixnew(6,3))+4*ixnew(6,2)-4
     krecarr(iconz)=rel
     goto 14
19   kndx(iconz)=0
     goto 14
17   krecarr(iconz)=0
     stop
14   a=a
     end do
     end do
     
     
100  return
     end















      
      
      
      
      
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      real krecarr
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      real krecarr
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      dimension kdum(50),isub(50)
      real krecarr
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      real krecarr
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      dimension ie(50),itwo(100),nnum(50)
      real krecarr
      
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

100   a=a
      if(ileng.ne.1)goto 110
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      real krecarr
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      
      dimension karu(50),karv1(50),karv3(50),karqq(50)
      dimension kart3(50),kart1(50)
      real krecarr
      
      
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
!      print *,'kart3',(kart3(i),i=1,kart3(2)+2)
      
      
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
!      print *,'kart1',(kart1(i),i=1,kart1(2)+2)
      
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000),norma(20),leadc(20),ijfak(20),ijfakn(20)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common kr(10001,10),kndx(400000),krecarr(400000),zlog(10000)
      common ipol1(6,50)
      dimension karae(50),karde(50),karpe(50),karh(50),ktemp(50)
      real krecarr
      
      
      
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
