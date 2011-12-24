     program bobmodj7 
! comes from bobmodj5 :uses different field for extracting roots     
! was  program bern1 
!     first stage multi-precisioning of "bern",computes multi-precision      
!     complex cube and square roots mod p

!     first in generalised GNFS suite    
!     smaller sieving interval
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(55),iansar(55),narc(55)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      dimension ireclen(100),irecarr(200),memfil1(8,30),memfil2(10,30)
      dimension iptest(100),nstack(30,100),idstack(30),ipemp1(100)
      dimension indstak(20),itemp1(100),jdis(100),jcc(100),mpre(100)
      dimension ipre(100),memlis(30),icurr(100),nfact2(100),limp(3)
      dimension memfil3(4,60),memfil4(15,100),iprar(100),idar(60),idsol(60)
      
      kpf(1)=0
      kpf(2)=1
      kpf(3)=2
      kqf(1)=0
      kqf(2)=1
      kqf(3)=3
      krf(1)=0
      krf(2)=1
      krf(3)=130
      ip(1)=0
      ip(2)=1
      ip(3)=131
!      call subbern(3,isok)
!      stop

!       call bobecm
!      stop
      open (unit=3,file='modjlen3',access='direct',form=&
      'formatted',recl=600,status='old')
      open (unit=2,file='modjpol3',access='sequential')
      read (3,1005,rec=1)(ireclen(jf),jf=1,79)
      open (unit=1,file='recl.dat',access='direct',form=&
      'formatted',recl=390000,status='old')
      read (1,1006,rec=1)(ipr(jf),jf=1,65000)
1006  format (65000i6)      
      print *,'ipr10833',ipr(10833),'ipr10834',ipr(10834)
      print *,'ipr63951',ipr(63951),'ipr63952',ipr(63952)
      
      
      iccc=0
      icd1=2
      icd2=0
      icd3=0
      icd4=0
      imax=0
1005  format (1000i4)
200   iccc=iccc+1      
      limt=ireclen(iccc)
      read (2,*,end=400)(irecarr(jf),jf=1,limt)
      if (irecarr(3).gt.3)goto 201
      if (limt.le.imax)goto 199
      imax=limt
199   if (irecarr(3).eq.1)goto 197
      if (irecarr(3).eq.3)goto 300
      icd2=icd2+1
      do jf=1,limt
      memfil2(icd2,jf)=irecarr(jf)
      end do
      goto 198
300   icd3=icd3+1
      do jf=1,limt
      memfil3(icd3,jf)=irecarr(jf)
      end do
!      print *,'dis 3',(irecarr(jk),jk=1,limt)
      goto 200
197   icd1=icd1+1
      do jf=1,limt
      memfil1(icd1,jf)=irecarr(jf)
      end do
!      print *,'dis1',(irecarr(jk),jk=1,limt)
      goto 200
198   a=a
!      print *,'dis2',(irecarr(jf),jf=1,limt)
      goto 200
201   a=a
      icd4=icd4+1
      do jf=1,limt
      memfil4(icd4,jf)=irecarr(jf)
      end do
!      if (irecarr(2).ne.172)goto 200
      print *,'dis4',(irecarr(jf),jf=1,limt)
      goto 200      
400   close (unit=2)      
      close (unit=3)
      
      print *,'maxlen dis <4',imax,'icd1',icd1,'icd2',icd2,'icd3',icd3
      print *,'icd4',icd4
      
      print *,'memfil1',(memfil1(3,jf),jf=1,memfil1(3,4))
      print *,'memfil3 2',(memfil3(2,jf),jf=1,memfil3(2,4))
      print *,'memfil3 3',(memfil3(3,jf),jf=1,memfil3(3,4))
      
      
      memfil3(1,1)=1
      memfil3(1,2)=379
      memfil3(1,3)=3
      memfil3(1,4)=39
      memfil3(1,5)=7
      memfil3(1,6)=9
      memfil3(1,7)=10
      memfil3(1,8)=0
      memfil3(1,9)=7
      memfil3(1,10)=364
     memfil3(1,11)=3954
      memfil3(1,12)=410
       memfil3(1,13)=4624
       memfil3(1,14)=2390
       memfil3(1,15)=1824
       memfil3(1,16)=6144
       memfil3(1,17)=1
       memfil3(1,18)=9
       memfil3(1,19)=1
       memfil3(1,20)=2156
       memfil3(1,21)=7791
       memfil3(1,22)=98
       memfil3(1,23)=8087
       memfil3(1,24)=6719
       memfil3(1,25)=5385
       memfil3(1,26)=2832
       memfil3(1,27)=1536
       memfil3(1,28)=0
        memfil3(1,29)=10
        memfil3(1,30)=15
        memfil3(1,31)=4436
        memfil3(1,32)=4
        memfil3(1,33)=7689
        memfil3(1,34)=119
        memfil3(1,35)=4802
        memfil3(1,36)=4601
        memfil3(1,37)=8074
        memfil3(1,38)=1514
        memfil3(1,39)=8544
      mpre(1)=0
      mpre(2)=0
      memfil1(1,2)=3
      memfil1(2,2)=4
      nstcon=1
      ndiscon=1
      indstak(1)=1
      memlis(1)=1
      idstack(nstcon)=memfil1(ndiscon,2)
!     input iptest here      
      iptest(1)=0
      iptest(2)=11
      iptest(3)=4
      iptest(4)=5994
      iptest(5)=8113
      iptest(6)=4788
      iptest(7)=6846
      iptest(8)=3102
      iptest(9)=2172
      iptest(10)=8895
      iptest(11)=2230
      iptest(12)=3430
      iptest(13)=1839
      print *,'number length radix 10000?'
      read *,iptest(2)
      print *,'number to be tested?'
      read *,(iptest(jf),jf=3,iptest(2)+2)
      if (iptest(2).ne.1)goto 520
      ntest(1)=0 
      ntest(2)=iptest(2)
      ntest(3)=iptest(3)
!      ntest(4)=iptest(4)
!      ntest(5)=iptest(5)
      nstack(1,1)=0
      nstack(1,2)=iptest(2)
      nstack(1,3)=iptest(3)
!      nstack(1,4)=iptest(4)
!      nstack(1,5)=iptest(5)
      goto 522
520   a=a
      do jf=1,iptest(2)+2
      ip(jf)=iptest(jf)
      ipd(jf)=ip(jf)
      nstack(nstcon,jf)=iptest(jf)
      end do
      
      
      
      call rmill(indp)
      if (indp.eq.1)goto 10
523   print *,'number is composite at outset',(iptest(jf),jf=1,iptest(2)+2)
      stop
10    idar(1)=3
      idar(2)=4
      idar(3)=7
      idar(4)=11
      idar(5)=19
      idar(6)=43
      idar(7)=67
      idar(8)=163
      idar(9)=8
      idar(10)=15
      idar(11)=51
      idar(12)=52
      idar(13)=88
      idar(14)=115
      idar(15)=123
      idar(16)=148
      idar(17)=232
      idar(18)=235
      idar(19)=267
      idar(20)=427
      idar(21)=20
      idar(22)=24
      idar(23)=40
      idar(24)=331
      idar(25)=1411
      idar(26)=120
      idar(27)=132
      idar(28)=168
      idar(29)=228
      idar(30)=280
      idar(31)=312
      idar(32)=340
      idar(33)=372
      idar(34)=408
      idar(35)=520
      idar(36)=708
      idar(37)=715
      idar(38)=760
      idar(39)=795
      idar(40)=315
      idar(41)=379
      idar(42)=643
      idar(43)=291
      idar(44)=259
      idar(45)=203
      idar(46)=1027
      idar(47)=23
      iddcon=0
      do i=1,47
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=idar(i)
!      call nacchia(inach)
      if (inach.ne.0)goto 530
      iddcon=iddcon+1
      idsol(iddcon)=idar(i)
530   end do
!      print *,'iddcon',iddcon,'idsols',(idsol(jf),jf=1,iddcon)
!      stop
      memind1=1
      memind2=0
      memind3=0
      memind4=0
1     kdcorn(1)=1
      kdcorn(2)=1
!      kdcorn(3)=88
!      kdcorn(3)=379
     kdcorn(3)=memfil1(memind1,2)
!      print *,'kdcorn3',kdcorn(3)
      
4     call nacchia(inach)
      if (memind2.ne.3)goto 101
!      print *,'nach',inach,'kdcorn',(kdcorn(jf),jf=1,kdcorn(2)+2)
!      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      
101   a=a
      
      jok=2
      if (inach.eq.0)goto 2
      memind1=memind1+1
      if (memind1.gt.8)goto 3
      memlis(nstcon)=memind1
      goto 1
3     memind2=memind2+1
      indstak(nstcon)=2
      if (memind2.gt.10)goto 5
      memlis(nstcon)=memind2
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=memfil2(memind2,2)
      goto 4
5     memind3=memind3+1
      indstak(nstcon)=3
      if (memind3.gt.1)goto 301
      memlis(nstcon)=memind3
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=memfil3(memind3,2)
      goto 4
301   memind4=memind4+1
      indstak(nstcon)=4
      if (memind4.gt.15)goto 600
      memlis(nstcon)=memind4
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=memfil4(memind4,2)
      goto 4

600   if (nstcon.ne.1)goto 70

59    print *,'failure to get started',' memind2=',memind2
      stop
2     a=a
      
      if (memind1.eq.1)goto 401
      if (memind1.eq.2)goto 410
      if (memind1.le.8)goto 80
      if (memind2.gt.10)goto 310
      indt=memfil2(memind2,5)
      do jf=7,8+indt
      marr(jf-6)=memfil2(memind2,jf)
      mbarr(jf-6)=memfil2(memind2,jf)
      end do
!      print *,'marr',(marr(jf),jf=1,marr(2)+2)
      
      call menmul
      do jf=1,mcarr(2)+2
      itemp1(jf)=mcarr(jf)
      end do
      indt2=memfil2(memind2,6)
      do jf=9+indt,10+indt+indt2
      print *,'jf',jf,'file',memfil2(memind2,jf)
      marr(jf-8-indt)=memfil2(memind2,jf)
      end do
!      print *,'conmarr',(marr(jf),jf=1,marr(2)+2)
      
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=4
      call menmul
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,itemp1(2)+2
      karr(jf)=itemp1(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      if (mcarr(1).eq.0)goto 6
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      jdis(jf)=kcarr(jf)
      kard(jf)=kcarr(jf)
      end do
7     do jf=1,ip(2)+2
      karp(jf)=ip(jf)
      end do
      call mpkron(k)
      if (k.eq.1)goto 8
11    memind2=memind2+1
      if (memind2.gt.10)goto 5
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=memfil2(memind2,2)
      goto 4
6     do jf=1,mcarr(2)+2
      jdis(jf)=mcarr(jf)
      kard(jf)=mcarr(jf)
      end do
      goto 7
8     do jf=1,jdis(2)+2
      ncom(1,jf)=jdis(jf)
      end do
      ncom(2,1)=0
      ncom(2,2)=0
      call cornsq
      if (nsq.eq.0)goto 11
      if (isqurar(1,2,2).ne.0)goto 11
      do jf=7,8+indt
      marr(jf-6)=memfil2(memind2,jf)
      print *,'marrj-6',marr(jf-6)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
!   situation  computing minus b where b might be negative     
      if (mcarr(1).eq.1)goto 12
      do jf=1,mcarr(2)+2
      kbarr(jf)=mcarr(jf)
      end do
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      call mpadd(1)
!      print *,'firkcarr',(kcarr(jf),jf=1,kcarr(2)+2)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      goto 13
12    do jf=1,mcarr(2)+2 
      karr(jf)=mcarr(jf)
      end do
      karr(1)=0
13    do jf=1,isqurar(1,1,2)+2
      kbarr(jf)=isqurar(1,1,jf)
      end do
      call mpadd(0)
!      print *,'kcarr',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      itemp1(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      marr(jf)=ip(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+2
      karr(jf)=mdarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,itemp1(2)+2
      mbarr(jf)=itemp1(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      itemp1(jf)=mcarr(jf)
      karr(jf)=mcarr(jf)
      end do
!      print *,'jdis',(jdis(jf),jf=1,jdis(2)+2)
!      print *,'isqr',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2)
      
!      print *,'j',(itemp1(jf),jf=1,itemp1(2)+2)
      
      
86    kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1728
      call mpadd(1)
      
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      if (mcarr(1).eq.0)goto 14
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karb(jf)=kcarr(jf)
      end do
      goto 15
14    do jf=1,mcarr(2)+2       
      karb(jf)=mcarr(jf)
      end do
      
15    do jf=1,ip(2)+2
      kara(jf)=ip(jf)
      end do
      call mpgcd
      do jf=1,karv(2)+2
      marr(jf)=karv(jf)
      end do
      do jf=1,itemp1(2)+2
      mbarr(jf)=itemp1(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      if (mcarr(1).eq.0)goto 16
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      jcc(jf)=kcarr(jf)
      end do
      goto 17
16    do jf=1,mcarr(2)+2
      jcc(jf)=mcarr(jf)
      end do
17    a=a
      print *,'jcc',(jcc(jf),jf=1,jcc(2)+2)
      
      do jf=1,jcc(2)+2
      marr(jf)=jcc(jf)
      end do
      mbarr(1)=1
      mbarr(2)=1
      mbarr(3)=3
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      
      
      
      
      
      if (mcarr(1).eq.0)goto 18
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      jaa(jf)=kcarr(jf)
      end do
      goto 19
18    do jf=1,mcarr(2)+2 
      jaa(jf)=mcarr(jf)
      end do
19    a=a
      do jf=1,jcc(2)+2
      karr(jf)=jcc(jf)
      kbarr(jf)=jcc(jf)
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
      jbb(jf)=mcarr(jf)
      end do
420   print *,'jaa',(jaa(jf),jf=1,jaa(2)+2)
      print *,'jbb',(jbb(jf),jf=1,jbb(2)+2)
!      print *,'kdcorn',(kdcorn(jf),jf=1,kdcorn(2)+2)
      
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do jf=1,kbcorn(2)+2
      kbarr(jf)=kbcorn(jf)
      end do
      call mpadd(0) 
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ntest(jf)=kcarr(jf)
      end do
      do jf=1,ntest(2)+2
      ip(jf)=ntest(jf)
      ipd(jf)=ip(jf)
      mmbig(jf)=ntest(jf)
      end do
      ithcon=0
      iffcon=0
      imsw=0
      imsw2=0
      imsw3=0
35    call rmill(indp)
!      print *,'book1',' imsw',imsw,'indp',indp
      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      
      if (indp.eq.0)goto 33
      
      
      do jf=2,ip(2)+2
      if (ip(jf).lt.mmbig(jf))goto 50
      end do
      
      if (imsw2.eq.1)goto 55
      goto 40



33    call subrho(iok)
      print *,'iok',iok,'nfact1',(nfact1(jf),jf=1,nfact1(2)+2),'d',kdcorn(3)
      print *,'nstcon',nstcon,'kycorn',(kycorn(jf),jf=1,kycorn(2)+2),&
      'ithcon',ithcon
      if (iok.eq.0)goto 20
      
      imsw3=1
      if (imsw.eq.1)goto 40
      goto 202
20    do jf=1,ntest(2)+2
      marr(jf)=ntest(jf)
      end do
      do jf=1,nfact1(2)+2
      mbarr(jf)=nfact1(jf)
      end do
      call mendiv
      do jf=1,mdarr(2)+2
      nfact2(jf)=mdarr(jf)
      end do
      imsw=1
!      print *,'nfact1',(nfact1(jf),jf=1,nfact1(2)+2)
!      print *,'nfact2',(nfact2(jf),jf=1,nfact2(2)+2)
202     do jf=1,mpre(2)+2 
       if (mpre(jf).ne.iptest(jf))goto 53
      end do
!      if (imsw.eq.0)goto 53
!      do jf=1,mpre(2)+2
!      if (mpre(jf).ne.iptest(jf))goto 53
!      end do
      if (iok.eq.1)goto 40
      goto 54
53    do jf=1,lcorn(2)+2
      karr(jf)=lcorn(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=2
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
! 53    do jf=1,lcorn(2)+2
!      marr(jf)=lcorn(jf)
!      end do
      inj=marr(2)+2
!      inj=lcorn(2)+2
      marr(inj+1)=0
      marr(inj+2)=0
      marr(2)=marr(2)+2
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+2
      itemp1(jf)=mdarr(jf)
      end do
      do jf=1,iptest(2)+2
      mpre(jf)=iptest(jf)
      end do
      imsw=1








!      itemp1(2)=itemp1(2)+2
      
      do jf=1,itemp1(2)+2
      ipre(jf)=itemp1(jf)
      end do
61    do jf=1,itemp1(2)+2
      marr(jf)=itemp1(jf)
      end do
      do jf=1,ipre(2)+2
      mbarr(jf)=ipre(jf)
      end do
      call mendiv
      do jf=1,mdarr(2)+2
      karr(jf)=mdarr(jf)
      end do
      do jf=1,ipre(2)+2
      kbarr(jf)=ipre(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+2
      icurr(jf)=mdarr(jf)
      end do
      do jf=2,icurr(2)+2
      if (icurr(jf).lt.ipre(jf))goto 62
      if (icurr(jf).gt.ipre(jf))goto 63
      end do
      goto 63
62   do jf=1,icurr(2)+2
     ipre(jf)=icurr(jf)
     end do
     goto 61
63   do jf=1,ipre(2)+2
     icurr(jf)=ipre(jf)
     karr(jf)=ipre(jf)
     end do
     print *,'4th root=',(icurr(jf),jf=1,icurr(2)+2)
     
     kbarr(1)=0
     kbarr(2)=2
     kbarr(3)=1
     kbarr(4)=5000
     kbarr(5)=0
     
     call mpadd(0)
     do jf=1,kcarr(2)+2
     marr(jf)=kcarr(jf)
     end do
     marr(2)=marr(2)-1
     do jf=1,marr(2)+2
     mbarr(jf)=marr(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     ipemp1(jf)=mcarr(jf)
     end do
!     print *,'ipemp1',(ipemp1(jf),jf=1,ipemp1(2)+2)
  
      
54    do jf=2,nfact1(2)+2
      if (nfact1(jf).lt.ipemp1(jf))goto 30
      if (nfact1(jf).gt.ipemp1(jf))goto 32
      end do
      goto 32       
30    do jf=1,nfact2(2)+2
      if (nfact2(jf).lt.ipemp1(jf))goto 40
      if (nfact2(jf).gt.ipemp1(jf))goto 42
      end do
      goto 42
40    if (imsw2.eq.1)goto 55
      do jf=1,iptest(2)+2
      karr(jf)=iptest(jf)
      end do
      do jf=1,kbcorn(2)+2
      kbarr(jf)=kbcorn(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ip(jf)=kcarr(jf)
      ipd(jf)=ip(jf)
      mmbig(jf)=kcarr(jf)
      ntest(jf)=kcarr(jf)
      end do
      imsw2=1
!      print *,'secip',(ip(jf),jf=1,ip(2)+2)
      
      
      goto 35
32    do jf=1,nfact1(2)+2
      ntest(jf)=nfact1(jf)
      
      ip(jf)=nfact1(jf)
      ipd(jf)=ip(jf)
      end do
      
      goto 35
42    do jf=1,nfact2(2)+2
      ntest(jf)=nfact2(jf)
      
      ip(jf)=nfact2(jf)
      ipd(jf)=ip(jf)
      end do

      goto 35
50    nn=1
51    nn=nn*607
!      print *,'stop we are expecting'
      
      nn=mod(nn,5003)
      marr(1)=0
      marr(2)=1
      marr(3)=nn
      do jf=1,iptest(2)+2
      mbarr(jf)=iptest(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kard(jf)=mcarr(jf)
      end do
      do jf=1,iptest(2)+2
      karp(jf)=iptest(jf)
      end do
      call mpkron(k)
!      print *,'k',k,'kard',(kard(jf),jf=1,kard(2)+2)
!      print *,'iptest',(iptest(jf),jf=1,iptest(2)+2)
      if (k.eq.1)goto 51
      if (memind1.ne.1)goto 422
      iaas(1,1)=0
      iaas(1,2)=1
      iaas(1,3)=nn
      iaas(2,1)=0
      iaas(2,2)=0
      do jf=1,iptest(2)+2
      karr(jf)=iptest(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=3
      call mendiv
      do jf=1,mdarr(2)+2
      ipn(jf)=mdarr(jf)
      end do
      do jf=1,iptest(2)+2
      ip(jf)=iptest(jf)
      ipd(jf)=ip(jf)
      end do
      
      call sub516
      if ((ibprod(1,2).eq.1).and.(ibprod(1,3).eq.1))goto 51
      
   



422   jgg(1)=0
      jgg(2)=1
      jgg(3)=nn
      do jf=1,iptest(2)+2
      ip(jf)=iptest(jf)
      ipd(jf)=ip(jf)
      end do
      print *,'stop right'
      print *,'ip',(ip(jf),jf=1,ip(2)+2),'kdcorn',kdcorn(3)
      print *,'mmbig',(mmbig(jf),jf=1,mmbig(2)+2),'jgg',jgg(3)
      print *,'ntest',(ntest(jf),jf=1,ntest(2)+2)
      print *,'jgg',jgg(3),'jbb',jbb(3),'kbcorn',kbcorn(3)
      print *,'jaa',(jaa(jf),jf=1,jaa(2)+2),'jbb',(jbb(jf),jf=1,jbb(2)+2)
      print *,'kbcorn',(kbcorn(jf),jf=1,kbcorn(2)+2)
      print *,'kycorn',(kycorn(jf),jf=1,kycorn(2)+2)
      
      
      
      
      call bobecm(jok)
!  jok=0 means go forward,jok=1 means backtrack
      print *,'jokstop jok',jok,'memind1',memind1,'memind2',memind2
      
      if (jok.eq.1)goto 70
      if (jok.eq.2)goto 551
      if (ntest(2).gt.3)goto 76
      if ((ntest(3).gt.5000).and.(ntest(2).eq.3))goto 76
522   limp(1)=101
      limp(2)=10007
      limp(3)=799999
      ippx=ntest(2)
      llimp=limp(ippx)
!      print *,'ipr63951',ipr(63951),'ipr63952',ipr(63952)
83    do i=2,63951
      if (ipr(i).gt.llimp)goto 1000
      if (ipr(i).lt.10000)goto 77
      iprar(1)=0
      iprar(2)=2
      iprar(3)=ipr(i)/10000
      iprar(4)=ipr(i)-10000*iprar(3)
      goto 78
77    iprar(1)=0
      iprar(2)=1
      iprar(3)=ipr(i)
78    do jf=1,ntest(2)+2
      marr(jf)=ntest(jf)
      end do
      do jf=1,iprar(2)+2
      mbarr(jf)=iprar(jf)
      end do
      call mendiv
      if (mcarr(2).eq.0)goto 524
      end do
      goto 1000
524   if ((mdarr(2).eq.1).and.(mdarr(3).eq.1))goto 1000
      goto 70
76    nstcon=nstcon+1
      if (nstcon.gt.30)goto 1002
      do jf=1,ntest(2)+2
      nstack(nstcon,jf)=ntest(jf)
      iptest(jf)=ntest(jf)
      ip(jf)=ntest(jf)
      ipd(jf)=ip(jf)
      end do
      idstack(nstcon)=memfil1(1,2)
     indstak(nstcon)=1
     goto 10
70   nstcon=nstcon-1
     if (ntest(2).eq.1)goto 523
     if ((nstcon.eq.0).and.(jok.eq.1))goto 99
     if (nstcon.eq.0)goto 59
     if ((indstak(nstcon).eq.1).and.(memlis(nstcon).eq.8))goto 73
     if (indstak(nstcon).eq.2)goto 74
     if (indstak(nstcon).eq.3)goto 306
     if (indstak(nstcon).eq.4)goto 640
71   memlis(nstcon)=memlis(nstcon)+1
     indstak(nstcon)=1
     memind1=memlis(nstcon)
     memind2=0
     memind3=0
     memind4=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil1(memind1,2)
715  do jf=1,nstack(nstcon,2)+2
     iptest(jf)=nstack(nstcon,jf)
     ip(jf)=nstack(nstcon,jf)
     ipd(jf)=ip(jf)
     end do
     goto 4
!306  if (memlis(nstcon).eq.1)goto 70
306  if (memlis(nstcon).eq.1)goto 642     
     memlis(nstcon)=memlis(nstcon)+1
     indstak(nstcon)=3
     memind3=memlis(nstcon)
     memind1=9
     memind2=11
     memind4=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(memind3,2)
     goto 715
308  indstak(nstcon)=3
     memlis(nstcon)=1
     memind1=9
     memind2=11
     memind3=1
     memind4=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(1,2)
     goto  715

73   memind1=9
     memind2=1
     memind3=0
     memind4=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil2(1,2)
     indstak(nstcon)=2
     memlis(nstcon)=1
     goto 715
74   if (memlis(nstcon).eq.10)goto 308
     indstak(nstcon)=2
     memlis(nstcon)=memlis(nstcon)+1
     memind2=memlis(nstcon)
     memind1=9
     memind3=0
     memind4=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil2(memind2,2)
     goto 715
640  a=a
! backtracking on 4th degee polynomial
     if (memlis(nstcon).eq.15)goto 70
     memlis(nstcon)=memlis(nstcon)+1
     indstak(nstcon)=4
     memind4=memlis(nstcon)
     memind1=9
     memind2=11
     memind3=2
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil4(memind4,2)
     goto 715
642  indstak(nstcon)=4
     memlis(nstcon)=1
     memind1=9
     memind2=11
     memind3=2
     memind4=1
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil4(memind4,2)
     goto 715

55   a=a
     if ((kdcorn(3).eq.3).and.(ithcon.eq.0))goto 450
     if ((kdcorn(3).eq.3).and.(ithcon.eq.1))goto 470
     if ((kdcorn(3).eq.4).and.(iffcon.eq.0))goto 460
551  do jf=1,iptest(2)+2
     ip(jf)=iptest(jf)
     ipd(jf)=ip(jf)
     end do
     memind1=memind1+1
     jok=2
     if (memind1.gt.8)goto 56 
     memlis(nstcon)=memind1
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil1(memind1,2)
     
     goto 4
56   indstak(nstcon)=2
     memind2=memind2+1
     if (memind2.gt.10)goto 304
     memlis(nstcon)=memind2
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil2(memind2,2)
     
     goto 4
450  do jf=1,kycorn(2)+2
     marr(jf)=kycorn(jf)
     end do
     mbarr(1)=0
     mbarr(2)=1
     mbarr(3)=3
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,kbcorn(2)+2
     kbarr(jf)=kbcorn(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     marr(jf)=kcarr(jf)
     end do
     mbarr(1)=0 
     mbarr(2)=1
     mbarr(3)=2
     call mendiv
     do jf=1,mdarr(2)+2
     kbcorn(jf)=mdarr(jf)
     kbarr(jf)=mdarr(jf)
     end do
     
     do jf=1,iptest(2)+2
     karr(jf)=iptest(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     karr(jf)=kcarr(jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(0)
     do jf=1,kcarr(2)+2
     mmbig(jf)=kcarr(jf)
     ip(jf)=kcarr(jf)
     ipd(jf)=ip(jf)
     ntest(jf)=kcarr(jf)
     end do
     ithcon=1
     imsw=0
     imsw2=0
     imsw3=0
     goto 35
460  do jf=1,kycorn(2)+2
     karr(jf)=kycorn(jf)
     kbarr(jf)=kycorn(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     kbcorn(jf)=kcarr(jf)
     karr(jf)=kcarr(jf)
     end do
     do jf=1,iptest(2)+2
     kbarr(jf)=iptest(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     karr(jf)=kcarr(jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(0)
     do jf=1,kcarr(2)+2
     mmbig(jf)=kcarr(jf)
     ip(jf)=kcarr(jf)
     ipd(jf)=ip(jf)
     ntest(jf)=kcarr(jf)
     end do
     iffcon=1
     imsw=0
     imsw2=0
     imsw3=0
     goto 35
470  do jf=1,kycorn(2)+2
     marr(jf)=kycorn(jf)
     end do
     mbarr(1)=0
     mbarr(2)=1
     mbarr(3)=3
     call menmul
     do jf=1,mcarr(2)+2
     kbarr(jf)=mcarr(jf)
     end do
     do jf=1,kbcorn(2)+2
     karr(jf)=kbcorn(jf)
     end do
     call mpadd(1)
     
     
     do jf=1,kcarr(2)+2
     kbcorn(jf)=kcarr(jf)
     kbarr(jf)=kcarr(jf)
     end do
     do jf=1,iptest(2)+2
     karr(jf)=iptest(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     karr(jf)=kcarr(jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(0)
     do jf=1,kcarr(2)+2
     mmbig(jf)=kcarr(jf)
     ip(jf)=kcarr(jf)
     ipd(jf)=ip(jf)
     ntest(jf)=kcarr(jf)
     end do
     ithcon=2
     imsw=0
     imsw2=0
     imsw3=0
     goto 35

304  indstak(nstcon)=3      
     memind3=memind3+1
!     if (memind3.gt.1)goto 70
     if (memind3.gt.1)goto 630
     memlis(nstcon)=memind3
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(memind3,2)
     goto 4
630  indstak(nstcon)=4
! polynomial of 4th degree involved here     
     memind4=memind4+1
     if (memind4.gt.15)goto 70
     memlis(nstcon)=memind4
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil4(memind4,2)
     goto 4

     
80   a=a
!  linear minimum polynomial used here
     indt3=memfil1(memind1,5) 
     do jf=6,7+indt3
     marr(jf-5)=memfil1(memind1,jf)
     end do
     marr(1)=mod(marr(1)+1,2)
     print *,'marr',(marr(jf),jf=1,marr(2)+2)
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.0)goto 81
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     karr(jf)=kcarr(jf)
     itemp1(jf)=kcarr(jf)
     end do
!     print *,'karr',(karr(jf),jf=1,karr(2)+2)
     goto 86
81   do jf=1,mcarr(2)+2 
     karr(jf)=mcarr(jf)
     itemp1(jf)=mcarr(jf)
     end do
!     print *,'karr',(karr(jf),jf=1,karr(2)+2)
     
     goto 86
310  if (memind3.gt.1)goto 610
!       polynomial of third degree involved here 
     indt=memfil3(memind3,5)
     do jf=8,9+indt
     marr(jf-7)=memfil3(memind3,jf)
     end do
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 311
     do jf=1,mcarr(2)+2
     jpol(2,jf)=mcarr(jf)
     end do
     goto 312
311  do jf=1,mcarr(2)+2      
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(2,jf)=kcarr(jf)
     end do






     
312  indt2=memfil3(memind3,6)
     print *,'did reach here'
     
     do jf=10+indt,11+indt+indt2
     marr(jf-9-indt)=memfil3(memind3,jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 313
     do jf=1,mcarr(2)+2
     jpol(3,jf)=mcarr(jf)
     end do
     goto 314
313  do jf=1,mcarr(2)+2     
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(3,jf)=kcarr(jf)
     end do








314  indt3=memfil3(memind3,7)
     do jf=12+indt+indt2,13+indt+indt2+indt3
     marr(jf-11-indt-indt2)=memfil3(memind3,jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 315
     do jf=1,mcarr(2)+2
     jpol(4,jf)=mcarr(jf)
     end do
     goto 316
315  do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(4,jf)=kcarr(jf)
     end do
316  jpold=3
     do jf=1,ip(2)+2
     ipd(jf)=ip(jf)
     end do

     call  bobfacp2(jpold)
     print *,'memind3',memind3
     print *,'nfsol',(nfsol(1,jf),jf=1,nfsol(1,2)+2)
     print *,'memfil3',(memfil3(memind3,jk),jk=1,memfil3(memind3,4))
     
!     if (isok.eq.1)goto 551
     do jf=1,nfsol(1,2)+2
     karr(jf)=nfsol(1,jf)
     itemp1(jf)=nfsol(1,jf)
     end do
     
     goto 86
610  a=a
! polynomial of 4th degree involved here
     
     indt=memfil4(memind4,5)
     indt2=memfil4(memind4,6)
     indt3=memfil4(memind4,7)
     indt4=memfil4(memind4,8)
     do jf=9,10+indt
     marr(jf-8)=memfil4(memind4,jf)
     end do
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 611
     do jf=1,mcarr(2)+2
     jpol(2,jf)=mcarr(jf)
     end do
     goto 612
611  do jf=1,mcarr(2)+2 
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(2,jf)=kcarr(jf)
     end do
612  do jf=11+indt,12+indt+indt2
     marr(jf-10-indt)=memfil4 (memind4,jf)
     end do
     
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 613
     do jf=1,mcarr(2)+2
     jpol(3,jf)=mcarr(jf)
     end do
     goto 614
613  do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(3,jf)=kcarr(jf)
     end do
614  do jf=13+indt+indt2,14+indt+indt2+indt3
     marr(jf-12-indt-indt2)=memfil4(memind4,jf)
     end do
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 615 
     do jf=1,mcarr(2)+2
     jpol(4,jf)=mcarr(jf)
     end do
     goto 616
615  do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(4,jf)=kcarr(jf)
     end do
616  do jf=15+indt+indt2+indt3,16+indt+indt2+indt3+indt4
     marr(jf-14-indt-indt2-indt3)=memfil4(memind4,jf)
     end do
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 617
     do jf=1,mcarr(2)+2
     jpol(5,jf)=mcarr(jf)
     end do
     goto 618
617  do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     kbarr(jf)=ip(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jpol(5,jf)=kcarr(jf)
     end do
618  jpold=4
     do jf=1,ip(2)+2
     ipd(jf)=ip(jf)
     end do
     call bobfacp2(jpold)
     do jf=1,nfsol(1,2)+2
     karr(jf)=nfsol(1,jf)
     itemp1(jf)=nfsol(1,jf)
     end do
     goto 86
     


401  jaa(1)=0
     jaa(2)=0
     do jf=1,ip(2)+2
     karr(jf)=ip(jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(1)
     do jf=1,kcarr(2)+2
     jbb(jf)=kcarr(jf)
     end do
     goto 420
410  jbb(1)=0
     jbb(2)=0
     do jf=1,ip(2)+2
     karr(jf)=ip(jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(1)
     do jf=1,kcarr(2)+2
     jaa(jf)=kcarr(jf)
     end do
     goto 420




99   print *,'number is composite'
     stop
1002 print *,'memory too small'
     stop
1000 print *,'path length=',nstcon
     print *,'number is definitely prime',(nstack(1,jf),jf=1,nstack(1,2)+2)
     end





      
      subroutine nacchia(inach)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
 
      dimension ipre(200),icurr(200)
      dimension isqrtd(200),kacorn(200)
      dimension krcorn(200),ktemp1(200),kchek(200),ktest(200)
      
      goto 778
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=43
      ncom(2,1)=0
      ncom(2,2)=0
      ip(1)=0
      ip(2)=1
      ip(3)=101
      call cornsq
      stop

!     call rmill(indp)
778   a=a
      do i=ip(2)+3,ip(2)+6
      ip(i)=0
      end do
      ip(2)=ip(2)+4
      do jf=1,ip(2)+2
      ipre(jf)=ip(jf)
      end do
1     do jf=1,ip(2)+2
      marr(jf)=ip(jf)
      end do
      do jf=1,ipre(2)+2
      mbarr(jf)=ipre(jf)
      end do
      call mendiv
      do jf=1,mdarr(2)+2
      karr(jf)=mdarr(jf)
      end do
      do jf=1,ipre(2)+2
      kbarr(jf)=ipre(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+ 2
      icurr(jf)=mdarr(jf)
      end do
      do jf=2,icurr(2)+2
      if (icurr(jf).lt.ipre(jf))goto 2
      if (icurr(jf).gt.ipre(jf))goto 3
      end do
      goto 3
2     do jf=1,icurr(2)+2
      ipre(jf)=icurr(jf)
      end do
      goto 1
3     do jf=1,ipre(2)+2
      icurr(jf)=ipre(jf)
      marr(jf)=ipre(jf)
      end do
      print *,'sq. root=',(icurr(jf),jf=1,icurr(2)+2)
!      print *,'marr',(marr(jf),jf=1,marr(2)+2)
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
!      print *,'karr',(karr(jf),jf=1,karr(2)+2)
!      kbarr(1)=0
!      kbarr(2)=2
!      kbarr(3)=5000
!      kbarr(4)=0
!      kbarr(5)=0
!      call mpadd(0)
      
      do jf=1,karr(2)+2
      lcorn(jf)=karr(jf)
      end do
      lcorn(2)=lcorn(2)-2
!      print *,'lcorn',(lcorn(jf),jf=1,lcorn(2)+2)
!  must be negative discriminant here      
      ip(2)=ip(2)-4
      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      
!      kdcorn(1)=1
!      kdcorn(2)=1
!      kdcorn(3)=19
      do jf=1,kdcorn(2)+2
      marr(jf)=kdcorn(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kard(jf)=kcarr(jf)
      end do
      do jf=1,ip(2)+2
      karp(jf)=ip(jf)
      end do
      print *,'oknach1'
      call mpkron(k)
      print *,'oknach2'
      if (k.ne.1)goto 99
      do jf=1,kard(2)+2
      ncom(1,jf)=kard(jf)
      end do
      ncom(2,1)=0
      ncom(2,2)=0
      call cornsq
      if (nsq.eq.0)goto 99
      if (isqurar(1,2,2).ne.0)goto 98
      ivan=isqurar(1,1,2)+2
      jvan=kdcorn(2)+2
!      print *,'ivan',ivan,'jvan',jvan
!      print *,'isq',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2)
!      print *,'kdcorn',(kdcorn(jf),jf=1,kdcorn(2)+2)
      
      
      if (mod(isqurar(1,1,ivan),2).ne.mod(kdcorn(jvan),2))goto 11
      do jf=1,isqurar(1,1,2)+2
      isqrtd(jf)=isqurar(1,1,jf)
      end do
!      print *,'okstop'
      
      goto 12
11    do jf=1,isqurar(1,1,2)+2
      kbarr(jf)=isqurar(1,1,jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      isqrtd(jf)=kcarr(jf)
      end do
12    do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      kbarr(jf)=ip(jf)
      end do
!  put double prime in kacorn      
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kacorn(jf)=kcarr(jf)
      end do
      do jf=1,isqrtd(2)+2
      kbcorn(jf)=isqrtd(jf)
      end do
15    print *,'kacorn',(kacorn(jf),jf=1,kacorn(2)+2)
!      print *,'kbcorn',(kbcorn(jf),jf=1,kbcorn(2)+2)
      

      do jf=2,kbcorn(2)+2
      if (kbcorn(jf).gt.lcorn(jf))goto 14
      if (kbcorn(jf).lt.lcorn(jf))goto 20
      end do
      goto 20

14    do jf=1,kacorn(2)+2
      marr(jf)=kacorn(jf)
      end do
      do jf=1,kbcorn(2)+2
      mbarr(jf)=kbcorn(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      krcorn(jf)=mcarr(jf)
      end do
      do jf=1,kbcorn(2)+2
      kacorn(jf)=kbcorn(jf)
      end do
      do jf=1,krcorn(2)+2
      kbcorn(jf)=krcorn(jf)
      end do
      
      goto 15
!  consider absolute value of discriminant
20    do jf=1,kbcorn(2)+2
      marr(jf)=kbcorn(jf)
      mbarr(jf)=kbcorn(jf)
      end do
!      print *,'2kbcorn',(kbcorn(jf),jf=1,kbcorn(2)+2)
       
      
      call menmul
      do jf=1,mcarr(2)+2
      ktemp1(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      marr(jf)=ip(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=4
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ktemp1(2)+2
      kbarr(jf)=ktemp1(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
!      print *,'2marr',(marr(jf),jf=1,marr(2)+2)
      
      do jf=1,kdcorn(2)+2
      mbarr(jf)=kdcorn(jf)
      end do
      mbarr(1)=0
      call mendiv
      if (mcarr(2).ne.0)goto 97
      do jf=1,mdarr(2)+2
      ktest(jf)=mdarr(jf)
      end do
      
      do jf=1,ktest(2)+2
      ipre(jf)=ktest(jf)
      end do
31    do jf=1,ktest(2)+2
      marr(jf)=ktest(jf)
      end do
      do jf=1,ipre(2)+2
      mbarr(jf)=ipre(jf)
      end do
      call mendiv
      do jf=1,mdarr(2)+2
      karr(jf)=mdarr(jf)
      end do
      do jf=1,ipre(2)+2
      kbarr(jf)=ipre(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+ 2
      icurr(jf)=mdarr(jf)
      end do
      do jf=2,icurr(2)+2
      if (icurr(jf).lt.ipre(jf))goto 32
      if (icurr(jf).gt.ipre(jf))goto 33
      end do
      goto 33
32    do jf=1,icurr(2)+2
      ipre(jf)=icurr(jf)
      end do
      goto 31
33    do jf=1,ipre(2)+2
      icurr(jf)=ipre(jf)
      kycorn(jf)=ipre(jf)
      end do
!      print *,'sq. root=',(icurr(jf),jf=1,icurr(2)+2)
!      print *,'marr',(marr(jf),jf=1,marr(2)+2)
      
      do jf=1,icurr(2)+2
      marr(jf)=icurr(jf)
      mbarr(jf)=icurr(jf)
      end do
      call menmul
!      print *,'okb'
      do jf=1,mcarr(2)+2
      kchek(jf)=mcarr(jf)
      end do
      do jf=2,kchek(2)+2
      print *,'kchek',kchek(jf),'ktest',ktest(jf),'jf',jf
      if (kchek(jf).ne.ktest(jf))goto 96
      end do
      print *,'solution x=',(kbcorn(jf),jf=1,kbcorn(2)+2)
      print *,'solution y=',(icurr(jf),jf=1,icurr(2)+2)
      inach=0
      goto 100
96    print *,'no solution: square test fails'
      inach=1
      goto 100
97    print *,'no solution: division test fails'
      inach=1
      goto 100
98    print *,'complex square root'
      inach=1
      goto 100
99    print *,'kronecker value=',k
      inach=1
       


       
100  return
     end
     
     
      subroutine addstar(in,out)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)

      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      character*(*) in,out
      out =in
      i =len(out)
      do while(out(i:i) ==' ')
      out(i:i)='*'
      i = i-1
      end do 
      return 
      end
      subroutine rmill(indp)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(55),iansar(55),narc(55)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      
      dimension ipsm(40),ixr(200),ntp(200),ibb(200)
!      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      
      
      
      idb=1
      ipsm(1)=2
      ipsm(2)=3
      ipsm(3)=5
      ipsm(4)=7
      ipsm(5)=11
      ipsm(6)=13
      ipsm(7)=17
      ipsm(8)=19
      ipsm(9)=23
      ipsm(10)=29
      ipsm(11)=31
      ipsm(12)=37
      ipsm(13)=41
      ipsm(14)=43
      ipsm(15)=47
      ipsm(16)=53
      ipsm(17)=59
      ipsm(18)=61
      ipsm(19)=67
      ipsm(20)=71
      ipsm(21)=73
      ipsm(22)=79
      ipsm(23)=83
      ipsm(24)=89
      ipsm(25)=97
      



!      ip(1)=0
!      ip(2)=4
!      ip(3)=1000
!      ip(4)=0
!      ip(5)=0
!      ip(6)=919
      
      do jf=1,ip(2)+2
      marr(jf)=ip(jf)
      end do
      
      do i=1,25
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=ipsm(i)
      call mendiv
      
      if (mcarr(2).eq.0)goto 660
      end do
!      print *,'ok1'
      
      
      
      
      
      
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      
      ntp(jf)=kcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      ind1=0
      
      do jf=1,ntp(2)+2
      marr(jf)=ntp(jf)
      end do
246   call mendiv
      if (mcarr(2).ne.0)goto 260
      do jf=1,mdarr(2)+2
      marr(jf)=mdarr(jf)
      end do
     ind1=ind1+1
     goto 246
260  itt=ind1
     do jf=1,marr(2)+2
     ipn(jf)=marr(jf)
     end do
     ixr(1)=0
     ixr(2)=1 
     ixr(3)=1 
264  do jf=1,ixr(2)+2 
     marr(jf)=ixr(jf) 
     end do 
     mbarr(1)=0 
     mbarr(2)=2 
     mbarr(3)=4035
     mbarr(4)=3607
     call menmul
!     print *,'ok2'
     do jf=1,mcarr(2)+2
     marr(jf)=mcarr(jf)
     end do
     do jf=1,ip(2)+2
     mbarr(jf)=ip(jf)
     end do
     call mendiv
     do jf=1,mcarr(2)+2
     ixr(jf)=mcarr(jf)
     end do
     do jf=1,ixr(2)+2
     iaas(1,jf)=ixr(jf)
     end do
     iaas(2,1)=0
     iaas(2,2)=0
     call sub516
     if ((icprod(1,2).eq.1).and.(icprod(1,3).eq.1))goto 700
     do jf=1,ntp(2)+2
     if (icprod(1,jf).ne.ntp(jf))goto 614
     end do
     goto 700
614  do jf=1,icprod(1,2)+2
     ibb(jf)=icprod(1,jf)
     end do
     iee=1
630  if (iee.gt.itt-1)goto 660
     do jf=1,ibb(2)+2
     marr(jf)=ibb(jf)
     mbarr(jf)=ibb(jf)
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
     ibb(jf)=mcarr(jf)
     end do
     do jf=1,ntp(2)+2
     if (ibb(jf).ne.ntp(jf))goto 646
     end do
     goto 700
646  iee=iee+1 
     goto 630
700  if (idb.eq.20)goto 750      
     idb=idb+1
     goto 264 
660  print *,(ip(jf),jf=1,ip(2)+2),'is composite'      
     indp=0
     goto 1000
      
750  print *,(ip(jf),jf=1,ip(2)+2),'is very probably prime'     
     indp=1
1000 return
     end 

     
      
      
      
      subroutine subgcd2
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
     
      do jf=3,kara(2)+2
      karr(jf-2)=kara(jf)
      end do
      ilen=kara(2)
      do jf=3,karb(2)+2
      kbarr(jf-2)=karb(jf)
      end do
      ilen2=karb(2)
1     call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 2
      do jf=1,ilen2
      karr(jf)=kbarr(jf)
      end do
      ilen=ilen2
      do jf=1,irlen
      kbarr(jf)=irrr(jf)
      end do
      ilen2=irlen
      goto 1
2     kard(1)=0
      kard(2)=ilen2
      do jf=1,ilen2
      kard(jf+2)=kbarr(jf)
      end do
      return 
      end



      
      
      
      subroutine subrho(iok)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100) 
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
     
     dimension ipre(100),nsub(100),itest(100),ijgcd(100),ijgpre(100)
     ijgpre(1)=0
     ijgpre(2)=1
     ijgpre(3)=1
     ijcon=0
     if (ntest(2).gt.30)goto 30
!     if (ntest(2).gt.25)goto 30
     jdtest=131271
     ijtexx=100000
     goto 31
30   jdtest=400000     
     ijtexx=200000
31   ipolc=1 
!     goto 999
     idiff=2 
     icon=4 
     ipre(1)=0 
     ipre(2)=1 
     ipre(3)=677 
1    do jf=1,ipre(2)+2
     nsub(jf)=ipre(jf)
     end do
     do i=1,icon
     do jf=3,ipre(2)+2
     karr(jf-2)=ipre(jf)
     kbarr(jf-2)=ipre(jf)
     end do
     ilen=ipre(2)
     ilen2=ipre(2)
     call mpmul(ilen,ilen2,ilen3)
     do jf=1,ilen3
     karr(jf+2)=kcarr(jf)
     end do
     karr(1)=0
     karr(2)=ilen3
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=ipolc
     call mpadd(0)
     do jf=3,kcarr(2)+2
     karr(jf-2)=kcarr(jf)
     end do
     ilen=kcarr(2)
     do jf=3,ntest(2)+2
     kbarr(jf-2)=ntest(jf)
     end do
     ilen2=ntest(2)
     call mpdiv(ilen,ilen2,irlen,icont,iswq)
     if (irlen.eq.0)goto 99
     do jf=1,irlen
     ipre(jf+2)=irrr(jf)
     end do
     ipre(1)=0
     ipre(2)=irlen
     if (i.le.idiff)goto 3
     do jf=1,ipre(2)+2
     karr(jf)=ipre(jf)
     end do
     do jf=1,nsub(2)+2
     kbarr(jf)=nsub(jf)
     end do
     call mpadd(1)
     if (kcarr(2).eq.0)goto 99
     isgn=kcarr(1)
     do jf=3,kcarr(2)+2
     karr(jf-2)=kcarr(jf)
     end do
     ilen=kcarr(2)
     do jf=3,ntest(2)+2
     kbarr(jf-2)=ntest(jf)
     end do
     ilen2=ntest(2)
     call mpdiv(ilen,ilen2,irlen,icont,iswq)
     if (irlen.eq.0)goto 99
     if (isgn.eq.0)goto 2
     do jf=1,irlen
     karr(jf+2)=irrr(jf)
     end do






     karr(1)=isgn
     karr(2)=irlen
     do jf=1,ntest(2)+2
     kbarr(jf)=ntest(jf)
     end do
     call mpadd(1)
     do jf=1,kcarr(2)+2
     itest(jf)=kcarr(jf)
     end do
     goto 25
2    do jf=1,irlen
     itest(jf+2)=irrr(jf)
     end do
     itest(1)=0
     itest(2)=irlen
25   do jf=1,itest(2)+2
     marr(jf)=itest(jf)
     end do
     do jf=1,ijgpre(2)+2
     mbarr(jf)=ijgpre(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     marr(jf)=mcarr(jf)
     end do
     do jf=1,ntest(2)+2
     mbarr(jf)=ntest(jf)
     end do
     call mendiv
     if (mcarr(2).eq.0)goto 9999
     do jf=1,mcarr(2)+2
     
     ijgpre(jf)=mcarr(jf)
     end do
     ijcon=ijcon+1
     if (ijcon.eq.200)goto 9999
     goto 201
9999 do jf=1,ijgpre(2)+2
     itest(jf)=ijgpre(jf)
     end do

     ijcon=0






     do jf=1,ntest(2)+2
     kara(jf)=ntest(jf)
     end do
     do jf=1,itest(2)+2
     karb(jf)=itest(jf)
     end do
     call subgcd2
     if (kard(2).gt.1)goto 90
     if (kard(3).gt.1)goto 90
201  if (i+icon.gt.jdtest)goto 999
3    end do
     idiff=icon
     icon=icon+icon
     print *,'icon',icon
!     if (icon.gt.100000000)goto 999
!     if (icon.gt.jdtest)goto 999
     goto 1
90   do jf=1,kard(2)+2
     nfact1(jf)=kard(jf)
     end do
     iok=0
     goto 100
99   ipolc=ipolc+1
     if (ipolc.gt.9000)goto 999
     idiff=2
     icon=4
     ipre(1)=0
     ipre(2)=1
     ipre(3)=676+ipolc
     goto 1
999  iok=1
! p-1 test applied after rho has failed
     if (ntest(2).le.30)goto 100
     icnt=0
     
     l1=ijtexx
     iprc=1
     ix1(1,1)=0
     ix1(1,2)=1
     ix1(1,3)=2
     ix2(1,1)=0
     ix2(1,2)=1
     ix2(1,3)=2
114  iprc=iprc+1
     if (ipr(iprc).gt.l1)goto 302
     iprx1=ipr(iprc)

1141 iq=iprx1
     iq1=iq
     l=l1/iq
128  if (iq1.gt.l)goto 131
     iq1=iq1*iq
     goto 128
131  a=a
!131  print *,'iq1',iq1,'l',l,'iprc',iprc,'ipr',ipr(iprc)
     iee=1
138  if (iee.gt.iq1)goto 146
     iee=iee*2
     goto 138
146  iee=iee/2
     iq1=iq1-iee
147  if (iee.eq.1)goto 300
     iee=iee/2
     do jf=1,ix2(1,2)+2
     marr(jf)=ix2(1,jf)
     mbarr(jf)=ix2(1,jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     marr(jf)=mcarr(jf)
     end do
     do jf=1,ntest(2)+2
     mbarr(jf)=ntest(jf)
     end do
     call mendiv
     if (mcarr(2).eq.0)goto 360
     do jf=1,mcarr(2)+2
     ix2(1,jf)=mcarr(jf)
     end do
     if (iq1.lt.iee)goto 147
     iq1=iq1-iee
1601 do jf=1,ix1(1,2)+2
     marr(jf)=ix1(1,jf)
     end do
     do jf=1,ix2(1,2)+2
     mbarr(jf)=ix2(1,jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     marr(jf)=mcarr(jf)
     end do
     do jf=1,ntest(2)+2
     mbarr(jf)=ntest(jf)
     end do
     call mendiv
     if (mcarr(2).eq.0)goto 360
     do jf=1,mcarr(2)+2
     ix2(1,jf)=mcarr(jf)
     end do
     goto 147
300  if ((ix2(1,2).eq.1).and.(ix2(1,3).eq.1))goto 1400

     do jf=1,ix2(1,2)+2
     ix1(1,jf)=ix2(1,jf)
     ijgpre(jf)=ix2(1,jf)
     end do
     icnt=icnt+1
     if (icnt.eq.1000)goto 302
     goto 114
302  icnt=0
     
     do jf=1,ix2(1,2)+2
     karr(jf)=ix2(1,jf)
     end do
     kbarr(1)=0
     kbarr(2)=1
     kbarr(3)=1
     call mpadd(1)
     do jf=1,kcarr(2)+2
     karb(jf)=kcarr(jf)
     end do
     do jf=1,ntest(2)+2
     kara(jf)=ntest(jf)
     end do
     call subgcd2
     if (kard(2).gt.1)goto 1500
     if (kard(3).gt.1)goto 1500
     if (ipr(iprc).gt.l1)goto 800
     goto 114
360  print *,'problem to be investigated'
3601 iok=1
     goto 100
1400 print *,'solution trivial'
     if (iprc.eq.2)goto 3601
     do jf=1,ijgpre(2)+2
     ix1(1,jf)=ijgpre(jf)
     ix2(1,jf)=ijgpre(jf)
     end do
     goto 147
800  print *,'no solution from p-1 algorithm'
     iok=1
     goto 100
1500 do jf=1,kard(2)+2
     nfact1(jf)=kard(jf)
     end do
     print *,'solution from p-1 algorithm',(nfact1(jf),jf=1,nfact1(2)+2)
     iok=0
     
100  return
     
     end




      
      
      
      

      

      subroutine bobecm(jok)
!     elliptic curve method multiple parallel inverses phase 1 only
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
     
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      dimension ixp(1,100),iyp(1,100),ibb(1,100),idd(1,100)
      dimension iprod2(400),iprod(400),nn(10)
     dimension iprecod(100),idy(400),idx(400),iee(100),iquar(100),itemp1(100)
!  dummy variable iaa not used in multi-precisioning      
      jjcon=0
      iaa=1
      do jf=1,ip(2)+2
      n(jf)=ip(jf)
      end do
      itag=0
      ithcon=0
      iffcon=0
      ifsw=0
      nn(1)=0
      nn(2)=1
      nn(3)=1
401   do jf=1,nn(2)+2
      marr(jf)=nn(jf)
      end do
      if (jjcon.ge.1000)goto 27
      mbarr(1)=0
      mbarr(2)=2
      mbarr(3)=4035
      mbarr(4)=3607
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      mbarr(1)=0
      mbarr(2)=2
      mbarr(3)=1000
      mbarr(4)=0
      call mendiv
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      nn(jf)=mcarr(jf)
      end do





!      print *,'ifsw',ifsw
!      if (jjcon.ge.20)goto 28
      do jf=1,nn(2)+2
      marr(jf)=nn(jf)
      mbarr(jf)=nn(jf)
      end do
      
      
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      ncom(1,jf)=mcarr(jf)
      end do
      do jf=1,jaa(2)+2
      marr(jf)=jaa(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ncom(1,2)+2
      kbarr(jf)=ncom(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      
      kbarr(jf)=kcarr(jf)
      end do
      do jf=1,jbb(2)+2
      karr(jf)=jbb(jf)
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
      ncom(1,jf)=mcarr(jf)
      kard(jf)=mcarr(jf)
      end do
      
      do jf=1,ip(2)+2
      karp(jf)=ip(jf)
      end do
      call mpkron(k)
      if (k.eq.1)goto 402
      goto 401
402   ncom(2,1)=0      
      ncom(2,2)=0
      call cornsq
      print *,'stop wrong','ifsw',ifsw,'jjcon',jjcon
!      if (ifsw.ne.2)goto 105
!      print *,'mmbig',(mmbig(jf),jf=1,mmbig(2)+2) 
!      print *,'n',(n(jf),jf=1,n(2)+2)
!      print *,'jgg',(jgg(jf),jf=1,jgg(2)+2)
!      print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
!      print *,'squrar',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2)
!      print *,'ntest',(ntest(jf),jf=1,ntest(2)+2)
      print *,'nn',(nn(jf),jf=1,nn(2)+2)
!      print *,'jaa sec',(jaa(jf),jf=1,jaa(2)+2)
!      print *,'jbb sec',(jbb(jf),jf=1,jbb(2)+2)
!      print *,'secstop'
      
105   a=a
      if (isqurar(1,2,2).ne.0)goto 401
      if (nsq.eq.0)goto 291
      do jf=1,nn(2)+2
      ix1(1,jf)=nn(jf)
      ix2(1,jf)=nn(jf)
      end do
!      ix1(1,1)=0
!      ix1(1,2)=1
!      ix1(1,3)=nn
!      ix2(1,1)=0
!      ix2(1,2)=1
!      ix2(1,3)=nn

      
      
      do jf=1,isqurar(1,1,2)+2
      iy1(1,jf)=isqurar(1,1,jf)
      iy2(1,jf)=isqurar(1,1,jf)
      end do
! temporary only      
!      ix1(1,1)=0
!      ix1(1,2)=1
!      ix1(1,3)=2
!      ix2(1,1)=0
!      ix2(1,2)=1
!      ix2(1,3)=2
!      iy1(1,1)=0
!      iy1(1,2)=1
!      iy1(1,3)=124
!      iy2(1,1)=0
!      iy2(1,2)=1
!      iy2(1,3)=124
!      jaa(1)=0
!      jaa(2)=1
!      jaa(3)=6703











      do jf=1,mmbig(2)+2
      iquar(jf)=mmbig(jf)
      end do
      iesw=1
      marr(1)=0
      marr(2)=1
      marr(3)=2
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
1     call menmul
      do jf=2,mcarr(2)+2
      if (mcarr(jf).lt.iquar(jf))goto 2
      if (mcarr(jf).gt.iquar(jf))goto 3
      end do
      goto 2
3     do jf=1,mcarr(2)+2
      iee(jf)=mcarr(jf)
      end do
      goto 4
2     do jf=1,mcarr(2)+2
      mbarr(jf)=mcarr(jf)
      end do
      goto 1
4     do jf=1,iee(2)+2      
      marr(jf)=iee(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+2
      iee(jf)=mdarr(jf)
      end do
!      print *,'iee',(iee(jf),jf=1,iee(2)+2)
      
      
      if ((iesw.eq.2).or.(ifsw.eq.2))goto 114
      goto 114

10    iesw=2
      do jf=1,mmbig(2)+2
      marr(jf)=mmbig(jf)
      end do
      do jf=1,ntest(2)+2
      mbarr(jf)=ntest(jf)
      end do
!      print *,'point e3 ok'
      call mendiv
      do jf=1,mdarr(2)+2
      iquar(jf)=mdarr(jf)
      end do
      
      marr(1)=0
      marr(2)=1
      marr(3)=2
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      goto 1
430   ithcon=ithcon+1   
      if (ithcon.gt.6)goto 29
      do jf=1,jbb(2)+2
      marr(jf)=jbb(jf)
      end do
      do jf=1,jgg(2)+2
      mbarr(jf)=jgg(jf)
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
      jbb(jf)=mcarr(jf)
      end do
      goto 401
440   iffcon=iffcon+1
      if (iffcon.gt.4)goto 29
      do jf=1,jaa(2)+2
      marr(jf)=jaa(jf)
      end do
      do jf=1,jgg(2)+2
      mbarr(jf)=jgg(jf)
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
      jaa(jf)=mcarr(jf)
      end do
      goto 401

       


11    a=a
      if (kdcorn(3).eq.3)goto 430
      if (kdcorn(3).eq.4)goto 440
      if (ifsw.eq.2)goto 29
      
!      print *,'point e4'
!      print *,'jaa1',(jaa(jf),jf=1,jaa(2)+2)
!      print *,'jbb1',(jbb(jf),jf=1,jbb(2)+2)
!      print *,'jgg',(jgg(jf),jf=1,jgg(2)+2)
!      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      do jf=1,jgg(2)+2
      marr(jf)=jgg(jf)
      mbarr(jf)=jgg(jf)
      end do
      call menmul
!      do jf=1,mcarr(2)+2
!      itemp1(jf)=mcarr(jf)
!      mbarr(jf)=mcarr(jf)
!      end do
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      itemp1(jf)=mcarr(jf)
      mbarr(jf)=mcarr(jf)
      end do



!      print *,'mbarr mid',(mbarr(jf),jf=1,mbarr(2)+2)
!      print *,'jaa again',(jaa(jf),jf=1,jaa(2)+2)
!      do jf=1,jaa(2)+2
!      mcarr(jf)=jaa(jf)
!      end do
!      do jf=1,mcarr(2)+2
!      marr(jf)=mcarr(jf)
!      end do
      do jf=1,jaa(2)+2
      marr(jf)=jaa(jf)
      end do
!      print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
!      print *,'jaa 3',(jaa(jf),jf=1,jaa(2)+2)
!      print *,'marr midmid',(marr(jf),jf=1,marr(2)+2)
!      print *,'mbarr midmid',(mbarr(jf),jf=1,mbarr(2)+2)
      
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
!      print *,'marr mid',(marr(jf),jf=1,marr(2)+2)
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      jaa(jf)=mcarr(jf)
      end do
      do jf=1,itemp1(2)+2
      marr(jf)=itemp1(jf)
      end do
      do jf=1,jgg(2)+2
      mbarr(jf)=jgg(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,jbb(2)+2
      mbarr(jf)=jbb(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
!      print *,'point e4 ok'
      do jf=1,mcarr(2)+2
      jbb(jf)=mcarr(jf)
      end do
!      print *,'jaa2',(jaa(jf),jf=1,jaa(2)+2)
!      print *,'jbb2',(jbb(jf),jf=1,jbb(2)+2)
      
      ifsw=2
      goto 401
27    jok=2
      goto 1600


28    jok=0
      goto 1600
29    jok=1
!      print *,'ip',(ip(jf),jf=1,ip(2)+2)
!      print *,'mmbig',(mmbig(jf),jf=1,mmbig(2)+2)
!      print *,'ntest',(ntest(jf),jf=1,ntest(2)+2)
      print *,'kdcorn',kdcorn(3),'ithcon',ithcon,'iffcon',iffcon
      goto 1600
291   print *,'square root expected but not present kdcorn=',kdcorn(3)
      jok=1
      goto 1600






      
      
      irecnn=0
      kkmx=0
      

      


      n(1)=0
      n(2)=1
      n(3)=5839
      do jf=1,n(2)+2
      iprecod(jf)=n(jf)
      end do
104    do i=1,200
      
      
      
      
      do jf=3,iprecod(2)+2
      karr(jf-2)=iprecod(jf)
      end do
      ilen=iprecod(2)
      kbarr(1)=2
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (icont.eq.0)goto 700
      iprecod(1)=0
      iprecod(2)=icont
      do jf=1,icont
      iprecod(jf+2)=ipqt(jf)
      end do
      if (irlen.gt.0)goto 40
      end do
      print *,'errorstart'
      stop
40    do jf=1,iprecod(2)+2
      karr(jf)=iprecod(jf)
      end do
      do jf=1,iprecod(2)+2
      kbarr(jf)=iprecod(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      karr(jf)=kcarr(jf)
      end do
      do jf=1,irlen
      kbarr(jf+2)=irrr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=irlen
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iprecod(jf)=kcarr(jf)
      end do
      do i=1,200
      do jf=3,iprecod(2)+2
      karr(jf-2)=iprecod(jf)
      end do
      ilen=iprecod(2)
      kbarr(1)=3
      ilen2=1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (icont.eq.0)goto 700
      iprecod(1)=0
      iprecod(2)=icont
      do jf=1,icont
      iprecod(jf+2)=ipqt(jf)
      end do
      if (irlen.gt.0)goto 53
      end do
      print *,'errorstart2'
      stop
53    do jf=3,iprecod(2)+2
      karr(jf-2)=iprecod(jf)
      
      end do
      ilen=iprecod(2)
      kbarr(1)=3
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      
      
      
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
      do jf=1,irlen
      kbarr(jf+2)=irrr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=irlen
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iprecod(jf)=kcarr(jf)
      end do
!      print *,'iprecod',(iprecod(jf),jf=1,iprecod(2)+2)
      goto 114
64    do i=1,1
! 64    do i=1,30
      ix1(i,1)=0
      ix1(i,2)=1
      ix1(i,3)=2
      ix2(i,1)=0
      ix2(i,2)=1
      ix2(i,3)=2
      iy1(i,1)=0
      iy1(i,2)=1
      iy1(i,3)=5440
      iy2(i,1)=0
      iy2(i,2)=1
      iy2(i,3)=5440
      end do
      
       
114   a=a
!      iquar(1)=0
!      iquar(2)=1
!      iquar(3)=8068 
!      iee(1)=0
!      iee(2)=1
!      iee(3)=4096







       do i=1,1
!      do i=1,30
      do jf=1,ix2(i,2)+2
      ixp(i,jf)=ix2(i,jf)
      end do
      do jf=1,iy2(i,2)+2
      iyp(i,jf)=iy2(i,jf)
      end do
      end do
      iq=12436
      iq1=iq
     ieb=8192
!     print *,'iee fir',(iee(jf),jf=1,iee(2)+2) 
      
       
! 128   if (iq1.gt.l)goto 131
!      iq1=iq1*iq
!      goto 128
! 131   print *,'iaaiq1',iaa,iq1
!      iee=1
! 138   if (iee.gt.iq1)goto 146
!      iee=iee*2
!      goto 138
!146   iee=iee/2
      do jf=1,iquar(2)+2
      karr(jf)=iquar(jf)
      end do
      do jf=1,iee(2)+2
      kbarr(jf)=iee(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      iquar(jf)=kcarr(jf)
      end do
!      print *,'iquar 1',(iquar(jf),jf=1,iquar(2)+2)
      iq1=iq1-ieb
! 147   if (ieb.eq.1)goto 300      
 147   if ((iee(2).eq.1).and.(iee(3).eq.1))goto 300
      do jf=1,iee(2)+2
      marr(jf)=iee(jf)
      end do
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call mendiv
      do jf=1,mdarr(2)+2
      iee(jf)=mdarr(jf)
      end do
!      print *,'iquar',(iquar(jf),jf=1,iquar(2)+2),'iq1',iq1
      ieb=ieb/2
!      print *,'iee',(iee(jf),jf=1,iee(2)+2),'ieb',ieb 
!      print *,'ix2',(ix2(1,jf),jf=1,ix2(1,2)+2)
!      print *,'iy2',(iy2(1,jf),jf=1,iy2(1,2)+2)
      
      
      

!     iq1=iq1-iee
! 147   if (iee.eq.1)goto 300
!      iee=iee/2
      do jf=1,iy2(1,2)+2
      karr(jf)=iy2(1,jf)
      kbarr(jf)=karr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      idy(jf)=kcarr(jf)
      end do
      
      do jf=3,idy(2)+2
      karr(jf-2)=idy(jf)
      end do
      ilen=idy(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      indz=1
      if (irlen.eq.0)goto 360
      
      do jf=1,irlen
      idy(jf+2)=irrr(jf)
      end do
      idy(2)=irlen
      if (idy(1).eq.0)goto 148
      do jf=1,idy(2)+2
      karr(jf)=idy(jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      idy(jf)=kcarr(jf)
      end do

148   do jf=1,idy(2)+2
      ibb(1,jf)=idy(jf)
      iprod(jf)=idy(jf)
      idd(1,jf)=idy(jf)
      end do

      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      do jf=1,ibb(1,2)+2
      karb(jf)=ibb(1,jf)
      end do

!      do jf=1,ibb(30,2)+2
!      karb(jf)=ibb(30,jf)
!      end do
!     print *,'hello1'
      call mpgcd
!     print *,'hello2' 
      if (kard(2).gt.1) goto 1500
      if (kard(3).gt.1)goto 1500
      do jf=1,karv(2)+2
      icc(1,jf)=karv(jf)
      end do
!      print *,'icc',(icc(1,jf),jf=1,icc(1,2)+2)


            
      call subeq(iaa,itag)
!      print *,'ix2',(ix2(1,jf),jf=1,ix2(1,2)+2),'hello3'
!      print *,'iy2',(iy2(1,jf),jf=1,iy2(1,2)+2)
      
      if (itag.eq.1)goto 360
      do jf=2,iquar(2)+2
      if (iquar(jf).lt.iee(jf))goto 147
      if (iquar(jf).gt.iee(jf))goto 149
      end do
!      if (iq1.lt.ieb)goto 147
      iq1=iq1-ieb
149   do jf=1,iquar(2)+2
      karr(jf)=iquar(jf)
      end do
      do jf=1,iee(2)+2
      kbarr(jf)=iee(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      iquar(jf)=kcarr(jf)
      end do


      
!      if (iq1.lt.iee)goto 147
!      iq1=iq1-iee
      do jf=1,ix1(1,2)+2
      karr(jf)=ix1(1,jf)
      end do
      do jf=1,ix2(1,2)+2
      kbarr(jf)=ix2(1,jf)
      end do
      call mpadd(1)
      
      do jf=1,kcarr(2)+2
      idx(jf)=kcarr(jf)
      end do
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      indz=7
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      ibb(1,jf+2)=irrr(jf)
      end do
      ibb(1,2)=irlen
      ibb(1,1)=idx(1)
      
      if (ibb(1,1).eq.0)goto 160
      do jf=1,ibb(1,2)+2
      karr(jf)=ibb(1,jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ibb(1,jf)=kcarr(jf)
      end do
      
160   do jf=1,ibb(1,2)+2      
      
      idd(1,jf)=ibb(1,jf)
      iprod(jf)=ibb(1,jf)
      end do
      if (ibb(1,2).eq.0)goto 360

      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
!      do jf=1,ibb(30,2)+2
!      karb(jf)=ibb(30,jf)
!      end do
      do jf=1,ibb(1,2)+2
      karb(jf)=ibb(1,jf)
      end do
!      print *,'hello4'
      call mpgcd
!      print *,'hello5'
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      do jf=1,karv(2)+2
      icc(1,jf)=karv(jf)
      end do


            
      call subneq(iaa,itag)
!      print *,'hello7'
      
      
      
      if (itag.eq.1)goto 360
      
      
      
      
      
      goto 147
300   do jf=1,ix2(1,2)+2
      ix1(1,jf)=ix2(1,jf)
      end do
      do jf=1,iy2(1,2)+2
      iy1(1,jf)=iy2(1,jf)
      end do
      print *,'normstop'
      
      if (iesw.eq.1)goto 11
      goto 28
      stop

! 300   do i=1,30
!      do jf=1,ix2(i,2)+2
!      ix1(i,jf)=ix2(i,jf)
!      end do
!      do jf=1,iy2(i,2)+2
!      iy1(i,jf)=iy2(i,jf)
!      end do
      

!      end do
!      goto 114
360   iaa=iaa +i+1
      itag=0
      jjcon=jjcon+1
!      print *,'ix2',(ix2(1,jf),jf=1,ix2(1,2)+2)
!      print *,'iy2',(iy2(1,jf),jf=1,iy2(1,2)+2)
!      print *,'indz',indz,'primestop'
       print *,'iesw360',iesw
       
!      if (ifsw.eq.2)goto 601
      if (iesw.eq.1)goto 10
      if (iesw.gt.1)goto 401
601   stop
      goto 64
600   iaa=iaa+30
      
      goto 64
700   print *,'all 2s and 3s'
      stop
1500  print *,'endwell prime=',(iquar(jf),jf=1,iquar(2)+2)
      print *,'prime divisor found on curve',(jaa(jf),jf=1,jaa(2)+2)
      print *,'prime divisor=',(kard(jf),jf=1,kard(2)+2)
      do jf=3,n(2)+2
      karr(jf-2)=n(jf)
      end do
      ilen=n(2)
      do jf=3,kard(2)+2
      kbarr(jf-2)=kard(jf)
      end do
      ilen2=kard(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      print *,'cofactor=0',icont,(ipqt(jf),jf=1,icont)
      goto 29
      stop






      
      
5     format(65000i6)      
1600  return
      end
      subroutine subeq(iaa,itag)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),mxxx(5),myyy(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      dimension m1(400),ixt(400),iyt(400),ixtemp(400),ms1(2000)
!      print *,'subeq ix2',(ix2(1,jf),jf=1,ix2(1,2)+2)
!      print *,'subeq iy2',(iy2(1,jf),jf=1,iy2(1,2)+2)
      
      do i=1,1
!      do i=1,30
      do jf=1,jaa(2)+2
      m1(jf)=jaa(jf)
      end do
      

      
!      m1(3)=iaa+i-1
!      m1(1)=0
!      m1(2)=1

      if (ix2(i,2).eq.0)goto 22
      do jf=3,ix2(i,2)+2
      karr(jf-2)=ix2(i,jf)
      kbarr(jf-2)=ix2(i,jf)
      end do
      ilen=ix2(i,2)
      ilen2=ix2(i,2)
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen',ilen,'ilen2',ilen2,'ilen3',ilen3
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      kbarr(1)=3
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen 3 2',ilen3
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      
      karr(2)=ilen3
      karr(1)=0
!      kbarr(1)=0
!      kbarr(2)=1
!      kbarr(3)=m1(3)
     do jf=1,m1(2)+2
     kbarr(jf)=m1(jf)
     end do

      
      call mpadd(0)
      
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      print *,'irlen 1',irlen
      if (irlen.eq.0)goto 20
      do jf=1,irlen
      m1(jf+2)=irrr(jf)
      end do
      m1(2)=irlen
      goto 22
20    m1(1)=0
      m1(2)=0
      m1(3)=0
      goto 24



      
22    do jf=3,m1(2)+2
      karr(jf-2)=m1(jf)
      end do
      ilen=m1(2)
      do jf=3,icc(i,2)+2
      kbarr(jf-2)=icc(i,jf)
      end do
      ilen2=icc(i,2)
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen3',ilen3
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      print *,'irlen 2',irlen
      if (irlen.eq.0)goto 71
      do jf=1,irlen
      m1(jf+2)=irrr(jf)
      end do
      m1(2)=irlen
      do jf=3,m1(2)+2
      karr(jf-2)=m1(jf)
      kbarr(jf-2)=m1(jf)
      end do
      ilen=m1(2)
      ilen2=m1(2)
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen3 3',ilen3
      do jf=1,ilen3
      ms1(jf+2)=kcarr(jf)
      end do
      ms1(1)=0
      ms1(2)=ilen3

24    if (ix2(i,2).eq.0)goto 30
      do jf=1,ix2(i,2)+2
      karr(jf)=ix2(i,jf)
      kbarr(jf)=karr(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ixt(jf)=kcarr(jf)
      end do
      ixt(1)=mod(kcarr(1)+1,2)
      goto 32
30    ixt(1)=0
      ixt(2)=0
      ixt(3)=0
32    do jf=1,ixt(2)+2
      karr(jf)=ixt(jf)
      end do
      do jf=1,ms1(2)+2
      kbarr(jf)=ms1(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ixt(jf)=kcarr(jf)
      end do
      do jf=1,ix2(i,2)+2
      karr(jf)=ix2(i,jf)
      end do
      do jf=1,ixt(2)+2
      kbarr(jf)=ixt(jf)
      end do
      call mpadd(1)
      
      do jf=1,kcarr(2)+2
      ixtemp(jf)=kcarr(jf)
      end do
      if (ixtemp(2).eq.0)goto 50
      do jf=3,ixtemp(2)+2
      karr(jf-2)=ixtemp(jf)
      end do
      ilen=ixtemp(2)
      if (m1(2).eq.0)goto 60
      do jf=3,m1(2)+2
      kbarr(jf-2)=m1(jf)
      end do
      ilen2=m1(2)
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen3 4',ilen3
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(2)=ilen3
      kbarr(1)=mod(m1(1)+ixtemp(1),2)
      goto 52
50    kbarr(1)=0
      kbarr(2)=0
      kbarr(3)=0
52    do jf=1,iy2(i,2)+2
      karr(jf)=iy2(i,jf)
      end do
      karr(1)=mod(iy2(i,1)+1,2)
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      iyt(jf)=kcarr(jf)
      end do
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      print *,'irlen3 3',irlen
      if (irlen.eq.0)goto 60
      do jf=1,irlen
      iy2(i,jf+2)=irrr(jf)
      end do
      iy2(i,1)=iyt(1)
      iy2(i,2)=irlen
      
      
      goto 62
60    iy2(i,1)=0
      iy2(i,2)=0
      iy2(i,3)=0
      print *,'neq'
      itag=1
      goto 80
62    do jf=3,ixt(2)+2
      karr(jf-2)=ixt(jf)
      end do
      ilen=ixt(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      if (irlen.eq.0)goto 70
      do jf=1,irlen
      ix2(i,jf+2)=irrr(jf)
      end do
      ix2(i,2)=irlen
      ix2(i,1)=ixt(1)
      if (ix2(i,1).eq.0)goto 72
      do jf=1,ix2(i,2)+2
      karr(jf)=ix2(i,jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ix2(i,jf)=kcarr(jf)
      end do
      goto 72
70    ix2(i,1)=0
      ix2(i,2)=0
      ix2(i,3)=0
      itag=1
      goto 80
      
71    print *,'factor found m=',(m1(jf),jf=1,m1(2)+2)
      do kk=1,20
      a=a
      end do
      
      
      
      
      itag=1
      goto 80
      
      
72    if (iy2(i,1).eq.0)goto 74
      do jf=1,iy2(i,2)+2
      karr(jf)=iy2(i,jf)
      end do
      
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      iy2(i,jf)=kcarr(jf)
      end do
74    a=a
      
      
      
      end do
80    return
      print *,'endeq'
      end
      subroutine subneq(iaa,itag)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),mxxx(5),myyy(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      dimension m1(400),ixt(400),iyt(400),ixtemp(400)
      
      
      inz=0
      do i=1,1
!      do i=1,30
      do jf=1,iy1(i,2)+2
      karr(jf)=iy1(i,jf)
      end do
      do jf=1,iy2(i,2)+2
      kbarr(jf)=iy2(i,jf)
      end do
      call mpadd(1)
      
      isgn=kcarr(1)
      
      if (kcarr(2).eq.0)goto 10
      do jf=3,kcarr(2) +2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,icc(i,2)+2
      kbarr(jf-2)=icc(i,jf)
      end do
      ilen2=icc(i,2)
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
      
      inz=1
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      m1(jf+2)=irrr(jf)
      end do
      m1(1)=mod(isgn+icc(i,1),2)
      m1(2)=irlen
      goto 12
10    m1(1)=0
      m1(2)=0
      m1(3)=0
12    do jf=2,ix1(i,2)+2  
      karr(jf)=ix1(i,jf)
      end do
      
      karr(1)=mod(ix1(i,1)+1,2)
      do jf=2,ix2(i,2)+2
      kbarr(jf)=ix2(i,jf)
      end do
      kbarr(1)=mod(ix2(i,1)+1,2)
      
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ixtemp(jf)=kcarr(jf)
      end do
      if (m1(2).eq.0)goto 20
      do jf=3,m1(2)+2
      karr(jf-2)=m1(jf)
      kbarr(jf-2)=m1(jf)
      end do
      ilen=m1(2)
      ilen2=m1(2)
      call mpmul(ilen,ilen2,ilen3)
      
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(2)=ilen3
      kbarr(1)=0
      do jf=1,ixtemp(2)+2
      karr(jf)=ixtemp(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ixt(jf)=kcarr(jf)
      end do
      goto 22
20    do jf=1,ixtemp(2)+2
      ixt(jf)=ixtemp(jf)
      end do
22    do jf=3,ixt(2)+2 
      karr(jf-2)=ixt(jf)
      end do
      
      ilen=ixt(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      inz=2
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      ix2(i,jf+2)=irrr(jf)
      end do
      ix2(i,2)=irlen
      ix2(i,1)=ixt(1)
      if (ix2(i,1).eq.0)goto 32
      do jf=1,ix2(i,2)+2
      karr(jf)=ix2(i,jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      ix2(i,jf)=kcarr(jf)
      end do
      
32    do jf=1,ix2(i,2)+2 
      kbarr(jf)=ix2(i,jf)
      end do
      
      do jf=1,ix1(i,2)+2
      karr(jf)=ix1(i,jf)
      end do
      call mpadd(1)
      
      isgn=kcarr(1)
      inz=3
      if (m1(2).eq.0)goto 360
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do
      do jf=3,m1(2)+2
      karr(jf-2)=m1(jf)
      end do
      ilen=m1(2)
      ilen2=kcarr(2)
      
      call mpmul(ilen,ilen2,ilen3)
      
      kbarr(1)=mod(isgn+m1(1),2)
      kbarr(2)=ilen3
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      do jf=2,iy1(i,2)+2
      karr(jf)=iy1(i,jf)
      end do
      karr(1)=mod(iy1(i,1)+1,2)
      call mpadd(0)
      
      isgn=kcarr(1)
      ilen=kcarr(2)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      inz=4
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      iy2(i,jf+2)=irrr(jf)
      end do
      iy2(i,2)=irlen
      iy2(i,1)=isgn
      if (isgn.eq.0)goto 40
      do jf=1,iy2(i,2)+2
      karr(jf)=iy2(i,jf)
      end do
      
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      iy2(i,jf)=kcarr(jf)
      end do
      goto 40
360   iy2(i,1)=0
      iy2(i,2)=0
      iy2(i,3)=0
      print *,'sneq','inz',inz
      print *,'icc',(icc(i,jf),jf=1,icc(i,2)+2),'i',i
      itag=1
      goto 80
      stop
40    a=a      
      
      end do
80    return
      end


      

      
      
      
                         
      







       
       subroutine bobfacp2(jpold)
!      program for trying out different methods of poly. factorization
!      involves powering to compute  GCD's of     
!      polynomials       
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
     print *,'did reach 2 jpold ',jpold  
       
       
       
       
       issz=0
       ncony=0
       nfcon=0
       ibbsw=0
       jpol(1,1)=0
       jpol(1,2)=1
       jpol(1,3)=1
       
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       nnh(jf)=mdarr(jf)
       end do

       ihalf(1)=0
       ihalf(2)=1
       
       ihalf(3) =1
37     do jf=2,ihalf(2)+2
       if (ihalf(jf).lt.nnh(jf))goto 79
       if (ihalf(jf).gt.nnh(jf))goto 38
       end do
79     do jf=1,ihalf(2)+2 
       karr(jf)=ihalf(jf)
       kbarr(jf)=ihalf(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ihalf(jf)=kcarr(jf)
       end do
       goto 37
38     do jf=1,ihalf(2)+2
       marr(jf)=ihalf(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ihalf(jf)=mdarr(jf)
       end do
       do jf=1,nnh(2)+2
       karr(jf)=nnh(jf)
       end do
       
       
       do jf=1,ihalf(2)+2
       
       kbarr(jf)=ihalf(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nnh(jf)=kcarr(jf)
       end do

       
       
       

       if (jpold.eq.1)goto 100
       if (jpold.eq.2)goto 200
       if (jpold.eq.3)goto 300
       if (jpold.eq.4)goto 400
       if (jpold.eq.5)goto 500
       stop
100    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       karb(jf)=igarray(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,igarray(2,2)+2
       mbarr(jf)=igarray(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
       goto 1000








200    call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
!       print *,'nfcon',nfcon,'sols',isol(1),isol(2)
       goto 1000
300    a=a
       
       call split(jpold,jpold2)
       print *,'jpold2',jpold2
       
!       if (issz.eq.1)goto 309
       if (jpold2.eq.2)goto 350
301    if (jpold2.eq.3)goto 350
303    if (jpold2.eq.4)goto 350       
       do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       karb(jf)=igarray(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,igarray(2,2)+2
       mbarr(jf)=igarray(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
       
       
       if (ibbsw.eq.1)goto 410
       if (ibbsw.eq.2)goto 510
       do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       idegg=jpold-jpold2
       call subbw2(idegd,idegr)
324    do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do

       
       nfcon=nfcon+2
       goto 1000
309    stop       
350    a=a
       
       
       do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,iqt(1,2)+2
       karb(jf)=iqt(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,iqt(2,2)+2
       mbarr(jf)=iqt(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
!       if (issz.eq.1)goto 309
       
       
       
!       print *,'nfcon',nfcon,'nfsol',nfsol(nfcon),'ibbsw',ibbsw
       
       if (ibbsw.eq.1)goto 430
       if (ibbsw.eq.2)goto 530
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       igarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw2(idegd,idegr)
       goto 324
400    call split(jpold,jpold2)
       print *,'400 jpold2',jpold2
       
       if (jpold2.eq.2)goto 420
       ibbsw=1
       goto 301
410    ibbsw=0
       do jf=1,4
       do jk=1,iqt(jf,2)+2
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       jpold=3
       goto 300

420    do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       igarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
       do jf=1,3
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       idegg=2
       call subbw2(idegd,idegr)
       goto 324
430    do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       jpold=jpold2
       ibbsw=0
       issz=1
      
       goto 300
500    call split(jpold,jpold2)
       print *,'jpold2',jpold2
       
       if (jpold2.eq.2)goto 520
       if (jpold2.eq.3)goto 522
       ibbsw=2
       goto 301
510    ibbsw=0
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       jpold=jpold2
       goto 400
520    call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
       do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       jpold=3
       ibbsw=0
       goto 300
522    do jf=1,3
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       print *,((isol(jf,jk),jk=1,isol(jf,2)+2),jf=1,2)
       
       nfcon=nfcon+2
       do jf=1,4
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       jpold=jpold2
       ibbsw=0
       goto 300
530    ibbsw=0
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       jpold=jpold2
       goto 400
1000   a=a
      print *,'nfcon',nfcon,'sols',((nfsol(jf,jk),jk=1,nfsol(jf,2)+2),&       
       jf=1,nfcon)
       
              
       end
       subroutine split(jpold,jpold2)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)

       
       
       
       dimension igg1(10,55),igg2(10,55),iapar(55)
       dimension ie(55),nn(55)
       print *,'did reach 3'
       
       idegp=1
       iapar(1)=0
       iapar(2)=1
       iapar(3)=1
61     igg1(1,1)=0
       igg1(1,2)=1
       igg1(1,3)=1
       do jf=1,iapar(2)+2
       igg1(2,jf)=iapar(jf)
       igg2(2,jf)=iapar(jf)
       end do
       
       igg2(1,1)=0
       igg2(1,2)=1
       igg2(1,3)=1
       
       idegn=1
       do jf=1,ihalf(2)+2
       ie(jf)=ihalf(jf)
       end do
       do jf=1,nnh(2)+2
       nn(jf)=nnh(jf)
       end do
       
       
       print *,'ihalf',(ihalf(jf),jf=1,ihalf(2)+2)
       print *,'nnh',(nnh(jf),jf=1,nnh(2)+2)
       print *,'ie',(ie(jf),jf=1,ie(2)+2)
!       if (iapar(3).ne.2)goto 4
       
       
4      if ((ie(2).eq.1).and.(ie(3).eq.1))goto 30
       print *,'ie',ie(3),'nn',nn(3)
       do jf=1,ie(2)+2
       marr(jf)=ie(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ie(jf)=mdarr(jf)
       end do

       
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult1(jf,jk)=igg2(jf,jk)
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       idegm=idegn
       call multy(idegm,idegn)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 10
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       isarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg)
       print *,'big'
       
       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
6      do jf=2,nn(2)+2
       if (nn(jf).lt.ie(jf))goto 4
       if (nn(jf).gt.ie(jf))goto 12
       end do
       goto 12

10     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 6
12     do jf=1,nn(2)+2
       karr(jf)=nn(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nn(jf)=kcarr(jf)
       end do





       print *,'coming thru'
       do jf=1,idegp+1
       do jk=1,igg1(jf,2)+2
       mult1(jf,jk)=igg1(jf,jk)
       end do
       end do
       idegm=idegp
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       call multy(idegm,idegn)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 20
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       isarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg)
       print *,'small'
       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
       goto 4
20     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 4
30     print *,'igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
       do jf=1,igg2(idegn+1,2)+2
       karr(jf)=igg2(idegn+1,jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       if (kcarr(1).eq.0)goto 32
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
32     do jf=1,kcarr(2)+2        
       igg2(idegn+1,jf)=kcarr(jf)
       end do







       
       do i=1,idegn+1
       
       if (igg2(i,2).eq.0)goto 33
       goto 34
33     end do
!       print *,'split gcd',(jpol(jf),jf=1,jpold+1)
!       iapar=iapar+1
!       if (iapar.gt.ipd)goto 71
       do jf=1,iapar(2)+2
       karr(jf)=iapar(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iapar(jf)=kcarr(jf)
       end do
       
       goto 61
71     print *,'problem in split'
       stop
34     idegs=idegn+1-i       
       print *,'idegs',idegs
       do jf=1,idegs+1
       do jk=1,igg2(jf+idegn-idegs,2)+2
       isarray(jf,jk)=igg2(jf+idegn-idegs,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       ibarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=jpold
       call subgcd(idegs,idegb,idegg)
       print *,'split gcd',((igarray(jf,jk),jk=1,igarray(jf,2)+2),&
       jf=1,idegg+1),'idegg',idegg,'idegb',idegb
       
       iapar=iapar+1
!       if (iapar.gt.ipd)goto 41
       do jf=1,iapar(2)+2
       karr(jf)=iapar(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iapar(jf)=kcarr(jf)
       end do
       
        
42     if ((idegg.eq.0).or.(idegg.eq.jpold))goto 61       
       do jf=1,idegg+1
       do jk=1,igarray(jf,2)+2
       isarray(jf,jk)=igarray(jf,jk)
       jpol2(jf,jk)=igarray(jf,jk)
       end do
       end do
       jpold2=idegg
       idegs=idegg
       call subbw4(idegs,idegb,idegg)
       
!       print *,'iquo',(iqt(jf),jf=1,idegb-idegs+1)
       
       
       return
       end
       subroutine menmul
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)

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
9      mcarr(1)=0
       mcarr(2)=0
10     return
       end
       subroutine mendiv
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
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
6      if (icont.eq.0)goto 12
       do jf=1,icont
       mdarr(jf+2)=ipqt(jf)
       end do
       mdarr(2)=icont
       mdarr(1)=mod(marr(1)+mbarr(1),2)
       goto 15
11     print *,'halted : attempted division by zero'
       stop
10     mcarr(1)=0
       mcarr(2)=0
       goto 6
9      mcarr(1)=0 
       mcarr(2)=0
12     mdarr(1)=0
       mdarr(2)=0
15     return
       end



       subroutine subbw2(idegd,idegr)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
       
       dimension irowp(20),irown(20),ic(20)
       dimension iran(20),nnumb1(55),nnumb2(55),nnumb3(55)
       dimension nnumb4(110),nnumb5(110)
       dimension ninv(55),iprar(55),nsol(2,55),ixold(6,25),ixnew(6,25)
       
       dimension itemp(55)
       


19     a=a
       
       do jf=1,ipd(2)+2
       iprar(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       nnumb1(jf)=igarray(1,jf)
       end do
       do jf=1,igarray(2,2)+2
       nnumb2(jf)=igarray(2,jf)
       end do
       do jf=1,igarray(3,2)+2
       nnumb3(jf)=igarray(3,jf)
       end do






312    do jf=1,iprar(2)+2
       kara(jf)=iprar(jf)
       end do
       do jf=1,igarray(1,2)+2
       karr(jf)=igarray(1,jf)
       kbarr(jf)=igarray(1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       karb(jf)=mcarr(jf)
       end do






       
       
!       itemp=2*igarray(1)
316    call mpgcd
       do jf=1,karv(2)+2
       ninv(jf)=karv(jf)
       end do
       do jf=3,nnumb2(2)+2
       karr(jf-2)=nnumb2(jf)
       kbarr(jf-2)=nnumb2(jf)
       end do
       ilen=nnumb2(2)
       ilen2=nnumb2(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       nnumb4(jf+2)=kcarr(jf)
       end do
       nnumb4(2)=ilen3
       nnumb4(1)=0
       print *,'nnumb4',(nnumb4(jk),jk=1,nnumb4(2)+2)
       do jf=3,nnumb1(2)+2
       karr(jf-2)=nnumb1(jf)
       end do
       ilen=nnumb1(2)
       do jf=3,nnumb3(2)+2
       kbarr(jf-2)=nnumb3(jf)
       end do
       ilen2=nnumb3(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf)=kcarr(jf)
       end do
       ilen2=ilen3
       karr(1)=4
       ilen=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=0

       
       do jf=1,nnumb4(2)+2
       karr(jf)=nnumb4(jf)
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


       
       
       do jf=1,irlen
       ncom(1,jf+2)=irrr(jf)
       end do
       ncom(1,1)=0
       ncom(1,2)=irlen
       if (isgn.eq.0)goto 318
       do jf=2,ncom(1,2)+2
       karr(jf)=ncom(1,jf)
       end do
       karr(1)=1


       do jf=1,iprar(2)+2
       kbarr(jf)=iprar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ncom(1,jf)=kcarr(jf)
       end do
       ncom(2,1)=0
       ncom(2,2)=0
318    print *,'ncom',(ncom(1,jk),jk=1,ncom(1,2)+2)
       
       call cornsq
       print *,'nsq',nsq
       

       
       do jf=1,igarray(2,2)+2
       kbarr(jf)=igarray(2,jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,isqurar(1,1,2)+2
       kbarr(jf)=isqurar(1,1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       if (mcarr(1).eq.0)goto 100
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       goto 102
100    do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
102    do jf=1,ninv(2)+2
       mbarr(jf)=ninv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       isol(1,jf)=mcarr(jf)
       end do
       print *,'isol1',(isol(1,jf),jf=1,isol(1,2)+2)
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,igarray(2,2)+2
       kbarr(jf)=igarray(2,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,isqurar(1,1,2)+2
       kbarr(jf)=isqurar(1,1,jf)
       end do
       call mpadd(1)
       
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       
       if (mcarr(1).eq.0)goto 150
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       goto 152
150    do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
152    do jf=1,ninv(2)+2
       mbarr(jf)=ninv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       isol(2,jf)=mcarr(jf)
       end do
       print *,'isol2',(isol(2,jf),jf=1,isol(2,2)+2)
       
       
       

       
       
400    a=a              
       return
       end
       subroutine subgcd(idegs,idegb,idegg)
       
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       dimension iws(20,55),iwb(20,55),itempb(20,55),mul(55),ib(55)
       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
       jf=1,idegb+1)
       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
       jf=1,idegs+1)

       
       do i= 1,idegs +1
       do jf=1,isarray(i,2)+2
       
       iws(i,jf) = isarray(i,jf)
       end do
       end do
       do i=1,idegb+1
       do jf=1,ibarray(i,2)+2
       iwb(i,jf) = ibarray(i,jf)
       end do
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
!       print *,'iwbd',iwbd,'iwsd',iwsd
       do jf=1,iws(1,2)+2
       ib(jf)=iws(1,jf)

       end do
       if (ib(2).ne.0)goto 101
       print *,'worries'
       stop
101    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,ib(2)+2
       karb(jf)=ib(jf)
       end do
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       
       

501    do i =1,loopl
       
       if (iwb(i,2).eq.0)goto 60
       if (mul(2).ne.0)goto 5019
       print *,'worries'
       stop
5019   do jf=3,iwb(i,2)+2
       karr(jf-2)=iwb(i,jf)
       end do
       ilen=iwb(i,2)
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ipd(2)+2
       kbarr(jf-2)=ipd(jf)
       end do
       ilen2=ipd(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       iqt(i,jf+2)=irrr(jf)
       end do
       iqt(i,1)=0
       iqt(i,2)=irlen
       
       
509    iwb(i,1)=0
       iwb(i,2)=0
       
       
       
       

       do j =2,iwsd +1
       do jf=1,iws(j,2)+2
       marr(jf)=iws(j,jf)
       end do
       do jf=1,iqt(i,2)+2
       mbarr(jf)=iqt(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,iwb(i+j-1,2)+2
       karr(jf)=iwb(i+j-1,jf)
       end do
       call mpadd(1)
       if (kcarr(1).eq.0)goto 20
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
       
20     do jf=1,kcarr(2)+2
       iwb(i+j-1,jf)=kcarr(jf)
       end do
       end do
60     end do
       do i=1,iwbd+1
       if (iwb(i,2).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       do jj=1,itempbd
       do jf=1,iwb(i+jj-1,2)+2
       itempb(jj,jf)=iwb(i+jj-1,jf)
       end do
       end do
       do i=1,iwsd+1
       do jf=1,iws(i,2)+2
       iwb(i,jf)=iws(i,jf)
       end do
       end do
       iwbd=iwsd
       do jj=1,itempbd
       do jf=1,itempb(jj,2)+2
       iws(jj,jf)=itempb(jj,jf)
       end do
       end do
       iwsd=itempbd-1
       goto 10
100    idegg=iwsd
       do i=1,iwsd+1
       do jf=1,iws(i,2)+2
       igarray(i,jf)=iws(i,jf)
       end do
       end do
       return
       end

       subroutine multy(idegm,idegn)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)

       
       do i=1,idegm+idegn+1
       do jf=1,50
       mult3(i,jf) =0
       end do
       end do
        
       do i =1,idegm +1
       do j =1,idegn+1
       
       do jf=1,mult1(i,2)+2
       marr(jf)=mult1(i,jf)
       end do
       do jf=1,mult2(j,2)+2
       mbarr(jf)=mult2(j,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,mult3(i+j-1,2)+2
       karr(jf)=mult3(i+j-1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       
       call mendiv
       do jf=1,mcarr(2)+2
       mult3(i+j-1,jf)=mcarr(jf)
       end do
       end do
       end do


!       mult3(i+j-1)=mult3(i+j-1) +mult1(i) *mult2(j)
!       mult3(i+j-1) =mod(mult3(i+j-1),ipd)
!       print *,'ok','idegm',idegm,'idegn',idegn
       return
       end
       subroutine subbw6(ia,ib,iv)
       
       
       ib =mod(ib,ia)
       
       
       
       if(ib.ge.0)goto 10
       ib =ia + ib
10     iu =1
       
       id = ia
       if(ib.eq.0)goto 888
       iv1=0
       iv3 =ib
815    if(iv3.eq.0)goto 830
       iqq = int(id/iv3)
       it3 =id -iqq*iv3
       it1 =iu -iqq*iv1
       iu =iv1
       id = iv3
       iv1 = it1
       iv3 = it3
       goto 815
830    iv =(id -ia *iu)/ib
       if(iu.le.0)goto 870
       iv =iv *(1-ia)
       iv =mod(iv,ia)
       if(iv.ge.0)goto 870
       iv=iv+ia
870    a=a
       if (id.gt.1)goto 890
       
       goto 890
888    iv =ib
890    return
       end
       subroutine subbw4(idegs,idegb,idegg)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
       dimension iws(20,55),iwb(20,55),itempb(20,55),mul(55),ib(55)
!       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
!       jf=1,idegb+1)
!       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
!       jf=1,idegs+1)
       
       
       
       do i= 1,idegs +1
       do jf=1,isarray(i,2)+2
       
       iws(i,jf) = isarray(i,jf)
       end do
       end do
       do i=1,idegb+1
       do jf=1,ibarray(i,2)+2
       iwb(i,jf) = ibarray(i,jf)
       end do
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       do jf=1,iws(1,2)+2
       ib(jf)=iws(1,jf)

       end do
       if (ib(2).ne.0)goto 101
       print *,'worries'
       stop
101    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,ib(2)+2
       karb(jf)=ib(jf)
       end do
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       
       

501    do i =1,loopl
       
       if (iwb(i,2).eq.0)goto 60
       if (mul(2).ne.0)goto 5019
       print *,'worries'
       stop
5019   do jf=3,iwb(i,2)+2
       karr(jf-2)=iwb(i,jf)
       end do
       ilen=iwb(i,2)
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ipd(2)+2
       kbarr(jf-2)=ipd(jf)
       end do
       ilen2=ipd(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       iqt(i,jf+2)=irrr(jf)
       end do
       iqt(i,1)=0
       iqt(i,2)=irlen
!       print *,'iqti',i,(iqt(i,jf),jf=1,iqt(i,2)+2)
       
509    iwb(i,1)=0
       iwb(i,2)=0
       
       
       
       

       do j =2,iwsd +1
       do jf=1,iws(j,2)+2
       marr(jf)=iws(j,jf)
       end do
       do jf=1,iqt(i,2)+2
       mbarr(jf)=iqt(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,iwb(i+j-1,2)+2
       karr(jf)=iwb(i+j-1,jf)
       end do
       call mpadd(1)
       if (kcarr(1).eq.0)goto 20
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
       
20     do jf=1,kcarr(2)+2
       iwb(i+j-1,jf)=kcarr(jf)
       end do
       end do
60     end do
       do i=1,iwbd+1
       if (iwb(i,2).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       
      do jj =1,itempbd
      do jf=1,iwb(jj+i-1,2)+2
      irarray(jj,jf) =iwb(jj+i-1,jf)
      end do
      end do
      idegr = itempbd -1
      goto 110

100   idegr =0
110   a=a
!     print *,'quos',(iqt(jf),jf=1,loopl)
!      print *,'rem',((irarray(jf,jk),jk=1,irarray(jf,2)+2),jf=1,idegr+1)
      idegg=idegr
      return
      end
      
       

       
       
      
      


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      
      
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
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      
      dimension kdum(600),isub(600)
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
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      
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

      
      

      subroutine mpgcd
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(55),iansar(55),narc(55)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
      
      
      dimension karu(190),karv1(190),karv3(190),karqq(190)
      dimension kart3(190),kart1(190)
      
      
      
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
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      
      
      
      dimension karae(190),karde(190),karpe(190),karh(190),ktemp(190)
      
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
       
      
          
       

       subroutine cornsq
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       
     
       
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       do jf=1,ipd(2)+2
       ip(jf)=ipd(jf)
       end do
       
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       
       iconz=1
       nn=2
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
       end do
       
       
       
!       do jf=3,ip(2)+2
!       karr(jf-2)=ip(jf)
!       kbarr(jf-2)=ip(jf)
!       end do
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
!       print *,'firprecod',(iprecod(jf),jf=1,iprecod(2)+2)
       
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
!       print *,'precod',(iprecod(jf),jf=1,iprecod(2)+2),'irem',irem
       
       
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
!       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
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
132    nn=nn+1        
       nn=mod(nn,5003)
       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
       kard(1)=0
       kard(2)=1
       kard(3)=nn
       do jf=1,ip(2)+2
       karp(jf)=ip(jf)
       end do
       call mpkron(k)
       if (k.ne.-1)goto 132
       
!       nn=nn*607
!       nn=mod(nn,1000)
!       iaas(2,1)=0
!       iaas(2,2)=1
!       iaas(2,3)=nn
       iaas(2,1)=0
       iaas(2,2)=0
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       if (ibsw.eq.1)goto 1620
112    call sub516
!       print *,'ok4'
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
!       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
!       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
!       ),jf=1,ibcn(2,2)+2)
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
       nsq=1
       
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
       
       
       
       
       
       
       
                   
       
       
       
          
       
       
       
       subroutine sub516
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
       
       
       
       
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)

       
       

       
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
!       print *,'ie',(ie(jf),jf=1,ie(2)+2)
       
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
!       print *,'ok2 516'
       
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
620    a=a
! 620    print *,'ok3 516'
       return
       end


       subroutine sub1100
       common ibarray(20,55),isarray(20,55),igarray(20,55),inv(25)
       common mult1(20,55),mult2(20,55),mult3(20,55),mult4(20,55)
       common mult5(20,55)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,55)
       common ipp(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(55)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20)
       common kara(190),karb(190),kard(190),karp(190),karv(190)
       common iarq(2),ncom(2,100),irarray(20,55),jpol(10,55),jpol2(10,55)
       common iqt(10,55),ipd(55)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(55),ihalf(55)
       
       
       
       
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120),khsq(120)
      common nfsol(5,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       do jf=1,ipd(2)+2
       ip(jf)=ipd(jf)
       end do
       
       if ((iacn(1,2).eq.0).or.(ibcn(1,2).eq.0))goto 4
       do jf=3,iacn(1,2)+2
       karr(jf-2)=iacn(1,jf)
       end do
       ilen=iacn(1,2)
       do jf=3,ibcn(1,2)+2
       kbarr(jf-2)=ibcn(1,jf)
       end do
       ilen2=ibcn(1,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       icprod(1,jf+2)=kcarr(jf)
       end do
       icprod(1,1)=0
       icprod(1,2)=ilen3
       goto 41
4      icprod(1,1)=0 
       icprod(1,2)=0
41     if ((iacn(2,2).eq.0).or.(ibcn(2,2).eq.0))goto 42

       do jf=3,iacn(2,2)+2
       karr(jf-2)=iacn(2,jf)
       end do
       ilen=iacn(2,2)
       do jf=3,ibcn(2,2)+2
       kbarr(jf-2)=ibcn(2,jf)
       end do
       ilen2=ibcn(2,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=ilen3
       do jf=1,icprod(1,2)+2
       karr(jf)=icprod(1,jf)
       end do
       goto 51
42     do jf=1,icprod(1,2)+2       
       kcarr(jf)=icprod(1,jf)
       end do
       goto 52

51     call mpadd(1)
52     if (kcarr(2).eq.0)goto 1
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
       if (irlen.eq.0)goto 1
       if (isgn.eq.0)goto 2
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
       icprod(1,jf)=kcarr(jf)
       end do
       goto 3
1      icprod(1,1)=0
       icprod(1,2)=0
       goto 3
2      do jf=1,irlen
       icprod(1,jf+2)=irrr(jf)
       end do
       icprod(1,1)=isgn
       icprod(1,2)=irlen
       
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
       
       
       
       



      
      
