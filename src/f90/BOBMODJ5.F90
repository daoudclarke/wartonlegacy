     program bobmodj5 
! was  program bern1 
!     first stage multi-precisioning of "bern",computes multi-precision      
!     complex cube and square roots mod p

!     first in generalised GNFS suite    
!     smaller sieving interval
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
      dimension ireclen(100),irecarr(200),memfil1(8,30),memfil2(10,30)
      dimension iptest(100),nstack(25,100),idstack(25),ipemp1(100)
      dimension indstak(20),itemp1(100),jdis(100),jcc(100),mpre(100)
      dimension ipre(100),memlis(25),icurr(100),nfact2(100),limp(3)
      dimension memfil3(4,60)
      
!      call subbern(3,isok)
!      stop

!       call bobecm
!      stop
      open (unit=3,file='modjlen2',access='direct',form=&
      'formatted',recl=600,status='old')
      open (unit=2,file='modjpol2',access='sequential')
      read (3,1005,rec=1)(ireclen(jf),jf=1,79)
      open (unit=1,file='recl.dat',access='direct',form=&
      'formatted',recl=390000,status='old')
      read (1,1006,rec=1)(ipr(jf),jf=1,65000)
1006  format (65000i6)      
      iccc=0
      icd1=2
      icd2=0
      icd3=0
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
      goto 200
198   a=a
!      print *,'dis2',(irecarr(jf),jf=1,limt)
      goto 200
201   a=a
!      if (irecarr(2).gt.714)goto 200
      print *,'dis4',(irecarr(jf),jf=1,limt)
      goto 200      
400   close (unit=2)      
      close (unit=3)
      
      print *,'maxlen dis <4',imax,'icd1',icd1,'icd2',icd2,'icd3',icd3
      print *,'memfil1',(memfil1(3,jf),jf=1,memfil1(3,4))
      print *,'memfil3 2',(memfil3(2,jf),jf=1,memfil3(2,4))
      print *,'memfil3 3',(memfil3(3,jf),jf=1,memfil3(3,4))
      
      
!      memfil3(1,10)=301
!      memfil3(1,11)=9767
!      memfil3(1,12)=9720
!      memfil3(1,15)=91
!      memfil3(1,16)=2937
!      memfil3(1,17)=9192
!      memfil3(1,18)=4024
!      memfil3(1,21)=36
!      memfil3(1,22)=7589
!      memfil3(1,23)=7833
!      memfil3(1,24)=7742
!      memfil3(1,25)=4823
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







      do jf=1,iptest(2)+2
      ip(jf)=iptest(jf)
      nstack(nstcon,jf)=iptest(jf)
      end do
      
      
      call rmill(indp)
      if (indp.eq.1)goto 10
      print *,'number is composite at outset',(iptest(jf),jf=1,iptest(2)+2)
      stop






      
      
10    memind1=1
      memind2=0
      memind3=0
1     kdcorn(1)=1
      kdcorn(2)=1
!      kdcorn(3)=88
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
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=memfil3(memind3,2)
      goto 4
301   if (nstcon.ne.1)goto 70

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
      mmbig(jf)=kcarr(jf)
      ntest(jf)=kcarr(jf)
      end do
      imsw2=1
!      print *,'secip',(ip(jf),jf=1,ip(2)+2)
      
      
      goto 35
32    do jf=1,nfact1(2)+2
      ntest(jf)=nfact1(jf)
      
      ip(jf)=nfact1(jf)
      end do
      
      goto 35
42    do jf=1,nfact2(2)+2
      ntest(jf)=nfact2(jf)
      
      ip(jf)=nfact2(jf)
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
      end do
      
      call sub516
      if ((ibprod(1,2).eq.1).and.(ibprod(1,3).eq.1))goto 51
      
   



422   jgg(1)=0
      jgg(2)=1
      jgg(3)=nn
      do jf=1,iptest(2)+2
      ip(jf)=iptest(jf)
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
      limp(1)=101
      limp(2)=10007
      limp(3)=ipr(65000)
      ippx=ntest(2)
      llimp=limp(ippx)
83    do i=2,65000
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
      if (mcarr(2).eq.0)goto 70
      end do
      goto 1000


76    nstcon=nstcon+1
      if (nstcon.gt.25)goto 1002
      do jf=1,ntest(2)+2
      nstack(nstcon,jf)=ntest(jf)
      iptest(jf)=ntest(jf)
      ip(jf)=ntest(jf)
      end do
      idstack(nstcon)=memfil1(1,2)
     indstak(nstcon)=1
     goto 10
70   nstcon=nstcon-1
     if ((nstcon.eq.0).and.(jok.eq.1))goto 99
     if (nstcon.eq.0)goto 59
     if ((indstak(nstcon).eq.1).and.(memlis(nstcon).eq.8))goto 73
     if (indstak(nstcon).eq.2)goto 74
     if (indstak(nstcon).eq.3)goto 306
71   memlis(nstcon)=memlis(nstcon)+1
     indstak(nstcon)=1
     memind1=memlis(nstcon)
     memind2=0
     memind3=0
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil1(memind1,2)
715  do jf=1,nstack(nstcon,2)+2
     iptest(jf)=nstack(nstcon,jf)
     ip(jf)=nstack(nstcon,jf)
     end do
     goto 4
306  if (memlis(nstcon).eq.1)goto 70
     memlis(nstcon)=memlis(nstcon)+1
     indstak(nstcon)=3
     memind3=memlis(nstcon)
     memind1=9
     memind2=11
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(memind3,2)
     goto 4
308  indstak(nstcon)=3
     memlis(nstcon)=1
     memind1=9
     memind2=11
     memind3=1
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(1,2)
     goto 4

73   memind1=9
     memind2=1
     memind3=0
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
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil2(memind2,2)
     goto 715
55   a=a
     if ((kdcorn(3).eq.3).and.(ithcon.eq.0))goto 450
     if ((kdcorn(3).eq.3).and.(ithcon.eq.1))goto 470
     if ((kdcorn(3).eq.4).and.(iffcon.eq.0))goto 460
551  do jf=1,iptest(2)+2
     ip(jf)=iptest(jf)
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
     ntest(jf)=kcarr(jf)
     end do
     ithcon=2
     imsw=0
     imsw2=0
     imsw3=0
     goto 35

304  indstak(nstcon)=3      
     memind3=memind3+1
     if (memind3.gt.1)goto 70
     memlis(nstcon)=memind3
     kdcorn(1)=1
     kdcorn(2)=1
     kdcorn(3)=memfil3(memind3,2)
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
310  a=a
!  polynomial of third degree involved here 
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
     kpf(jf)=mcarr(jf)
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
     kpf(jf)=kcarr(jf)
     end do






     
312  indt2=memfil3(memind3,6)
     do jf=10+indt,11+indt+indt2
     marr(jf-9-indt)=memfil3(memind3,jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 313
     do jf=1,mcarr(2)+2
     kqf(jf)=mcarr(jf)
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
     kqf(jf)=kcarr(jf)
     end do








314  indt3=memfil3(memind3,7)
     do jf=12+indt+indt2,13+indt+indt2+indt3
     marr(jf-11-indt-indt2)=memfil3(memind3,jf)
     end do
     call mendiv
     if (mcarr(1).eq.1)goto 315
     do jf=1,mcarr(2)+2
     krf(jf)=mcarr(jf)
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
     krf(jf)=kcarr(jf)
     end do





316  call subbern(3,isok)
     print *,'isok',isok,'memind3',memind3
     print *,'iroota',(iroota(1,jf),jf=1,iroota(1,2)+2)
     print *,'memfil3',(memfil3(memind3,jk),jk=1,memfil3(memind3,4))
     
     if (isok.eq.1)goto 551
     do jf=1,iroota(1,2)+2
     karr(jf)=iroota(1,jf)
     itemp1(jf)=iroota(1,jf)
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
 
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
      kbarr(1)=0
      kbarr(2)=2
      kbarr(3)=5000
      kbarr(4)=0
      kbarr(5)=0
      call mpadd(0)
      
      do jf=1,kcarr(2)+2
      lcorn(jf)=kcarr(jf)
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

      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
     
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100) 
     common kpf(120),kqf(120),krf(120),iroota(2,120)
     
     dimension ipre(100),nsub(100),itest(100)
     if (ntest(2).gt.35)goto 30
     jdtest=100000
     goto 31
30   jdtest=200000     
31   ipolc=1 
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
25   do jf=1,ntest(2)+2
     kara(jf)=ntest(jf)
     end do
     do jf=1,itest(2)+2
     karb(jf)=itest(jf)
     end do
     call subgcd2
     if (kard(2).gt.1)goto 90
     if (kard(3).gt.1)goto 90
3    end do
     idiff=icon
     icon=icon+icon
     print *,'icon',icon
!     if (icon.gt.100000000)goto 999
     if (icon.gt.jdtest)goto 999
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
100  return
     end




      
      
      subroutine menmul
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
      dimension kdum(400),isub(400)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)

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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
      dimension karu(100),karv1(100),karv3(100),karqq(100)
      dimension kart3(100),kart1(100)
      
      
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)

      dimension karae(120),karde(120),karpe(120),karh(120),ktemp(120)
      
      
      
      
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
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       
       
       common ipr(65000),norma(50)
       common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
       
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
       nn=mod(nn,1000)
       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
       nn=nn*607
       nn=mod(nn,1000)
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

       subroutine cornsq
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       
       
       common ipr(65000),norma(50)
       common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
       
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
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
       
       
       
       
       
       
       
       subroutine cub5
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       
       common ipr(65000),norma(50)
       common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
       
       dimension itempz(2),iaa(2,100),inv2(100),minusb(2,100)
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
       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
       ),jf=1,ibcn(2,2)+2)
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
       print *,'nsq',nsq,'should not be'
       ncubr=0
       stop
       
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
1510   print *,'major error'
       stop











       
       
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
       print *,'icubar',icubar(1,1,2),icubar(1,1,3)
       
       call bwq5
       print *,'icubar13',icubar(1,1,2),icubar(1,1,3)
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
       print *,'icubar1',(icubar(1,1,jf),jf=1,icubar(1,1,2)+2)
       print *,'icubarcom',(icubar(1,1,jf),jf=1,icubar(1,1,2)+2)
       
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
       
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
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
       
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       do ibig=1,2
       if (icprod(ibig,2).eq.0)goto 1518
       do jf=1,icprod(ibig,2)+2
       kbarr(jf)=icprod(ibig,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       icprod(ibig,jf)=kcarr(jf)
       end do
1518   end do 
       goto 1521
1520  print *,'2icubar',(icubar(2,1,jf),jf=1,icubar(2,1,2)+2) 
      print *,'2icubarcomp',(icubar(2,2,jf),jf=1,icubar(2,2,2)+2)
      print *,'3icubar',(icubar(3,1,jf),jf=1,icubar(3,1,2)+2)
      print *,'3icubar',(icubar(3,2,jf),jf=1,icubar(3,2,2)+2)
      print *,'1icubar',(icubar(1,1,jf),jf=1,icubar(1,1,2)+2)
      print *,'1icubarc',(icubar(1,2,jf),jf=1,icubar(1,2,2)+2)
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
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       
       
       common ipr(65000),norma(50)
       common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
       
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
!       print *,'ok2'
       
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
! 620    print *,'ok3'
       return
       end


       subroutine sub1100
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       
       
       common ipr(65000),norma(50)
       common kara(90),karb(90),kard(90),karp(90),karv(90)

       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
       
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
       
       
       
       








      

      subroutine bobecm(jok)
!     elliptic curve method multiple parallel inverses phase 1 only
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
     
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
!      do i=2,30
!      do jf=1,iy2(i,2)+2
!      karr(jf)=iy2(i,jf)
!      kbarr(jf)=karr(jf)
!      end do
!      call mpadd(0)
!      do jf=1,kcarr(2)+2
!      idd(i,jf)=kcarr(jf)
!      end do
!      do jf=3,idd(i,2)+2
!      karr(jf-2)=idd(i,jf)
!      end do
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen=idd(i,2)
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=2
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      idd(i,jf+2)=irrr(jf)
!      end do
!      idd(i,2)=irlen
!      if (idd(i,1).eq.0)goto 150
!      do jf=1,idd(i,2)+2
!      karr(jf)=idd(i,jf)
!      end do
!      do jf=1,n(2)+2
!      kbarr(jf)=n(jf)
!      end do
!      call mpadd(0)
!      do jf=1,kcarr(2)+2
!      idd(i,jf)=kcarr(jf)
!      end do
! 150   do jf=3,idd(i,2)+2
!      karr(jf-2)=idd(i,jf)
!      end do
!      ilen=idd(i,2)
      
!      do jf=3,iprod(2)+2
!      kbarr(jf-2)=iprod(jf)
!      end do
!      ilen2=iprod(2)
!      
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=3
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      iprod(jf+2)=irrr(jf)
!      end do
!      iprod(2)=irlen
      
!      do jf=1,iprod(2)+2
!      ibb(i,jf)=iprod(jf)
!      end do
      
      
!      end do
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

!      do jf=3,karv(2)+2
!      karr(jf-2)=karv(jf)
!      iprod2(jf)=karv(jf)
!      end do
!      ilen=karv(2)
!      iprod2(2)=ilen
!      iprod2(1)=0
!     do jf=3,ibb(29,2)+2
!      kbarr(jf-2)=ibb(29,jf)
!      end do
!      ilen2=ibb(29,2)
!      call mpmul(ilen,ilen2,ilen3) 
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=4
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      icc(30,jf+2)=irrr(jf)
!      end do
!      icc(30,2)=irlen
!      icc(30,1)=0
!      do i=1,28
!      do jf=3,iprod2(2)+2
!      karr(jf-2)=iprod2(jf)
!      end do
!      ilen=iprod2(2)
!      do jf=3,idd(31-i,2)+2
!      kbarr(jf-2)=idd(31-i,jf)
!      end do
!      ilen2=idd(31-i,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
      
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      iprod2(jf+2)=irrr(jf)
!      end do
!      iprod2(1)=0
!      iprod2(2)=irlen
!      do jf=3,iprod2(2)+2
!      karr(jf-2)=iprod2(jf)
!      end do
!      ilen=iprod2(2)






!      do jf=3,ibb(29-i,2)+2
!      kbarr(jf-2)=ibb(29-i,jf)
!      end do
!      ilen2=ibb(29-i,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
      
      
      
      
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=5
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      icc(30-i,jf+2)=irrr(jf)
!      end do
!      icc(30-i,2)=irlen
!      icc(30-i,1)=0
!      
!      end do
!      do jf=3,iprod2(2)+2
!      karr(jf-2)=iprod2(jf)
!      end do
!      ilen=iprod2(2)
!      do jf=3,idd(2,2)+2
!      kbarr(jf-2)=idd(2,jf)
!      end do
!      ilen2=idd(2,2)
      
      
      
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=6
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      icc(1,jf+2)=irrr(jf)
!      end do
!      icc(1,1)=0
!      icc(1,2)=irlen
            
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
!      do i=2,30
!      do jf=1,ix1(i,2)+2
!      karr(jf)=ix1(i,jf)
!      end do
!      do jf=1,ix2(i,2)+2
!      kbarr(jf)=ix2(i,jf)
!      end do
!      call mpadd(1)
!      do jf=3,kcarr(2)+2
!      karr(jf-2)=kcarr(jf)
!      end do
!      ilen=kcarr(2)
!      idd(i,1)=kcarr(1)
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=8
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      idd(i,jf+2)=irrr(jf)
!      karr(jf)=irrr(jf)
!      end do
!      ilen=irlen
!      idd(i,2)=irlen
!      if (idd(i,1).eq.0)goto 191
!      do jf=1,idd(i,2)+2
!      karr(jf)=idd(i,jf)
!      end do
!      do jf=1,n(2)+2
!      kbarr(jf)=n(jf)
!      end do
!      call mpadd(0)
!      do jf=1,kcarr(2)+2
!      idd(i,jf)=kcarr(jf)
      
!      end do
!      do jf=3,kcarr(2)+2
!      karr(jf-2)=kcarr(jf)
!      end do
!      ilen=kcarr(2)


! 191   a=a      
      
      
!      do jf=3,iprod(2)+2
!      kbarr(jf-2)=iprod(jf)
!      end do
!      ilen2=iprod(2)
!      call mpmul(ilen,ilen2,ilen3)
      
      
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=9
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      iprod(jf+2)=irrr(jf)
!      end do
!      iprod(2)=irlen
!      do jf=1,iprod(2)+2
!      ibb(i,jf)=iprod(jf)
!      end do
      
!      end do
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

!      do jf=3,karv(2)+2
!      karr(jf-2)=karv(jf)
!      iprod2(jf)=karv(jf)
!      end do
!      ilen=karv(2)
!      iprod2(2)=ilen
!      iprod2(1)=0
!      do jf=3,ibb(29,2)+2
!      kbarr(jf-2)=ibb(29,jf)
!      end do
!      ilen2=ibb(29,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=10
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      icc(30,jf+2)=irrr(jf)
!      end do
!      icc(30,2)=irlen
!      icc(30,1)=0
      
!      do i=1,28
!      do jf=3,iprod2(2)+2
!      karr(jf-2)=iprod2(jf)
!      end do
!      ilen=iprod2(2)
!      do jf=3,idd(31-i,2)+2
!      kbarr(jf-2)=idd(31-i,jf)
!      end do
!      ilen2=idd(31-i,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=11
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      iprod2(jf+2)=irrr(jf)
!      karr(jf)=irrr(jf)
!      end do
!      iprod2(2)=irlen
!      iprod2(1)=0
!      ilen=irlen
!      do jf=3,ibb(29-i,2)+2
!      kbarr(jf-2)=ibb(29-i,jf)
!      end do
!      ilen2=ibb(29-i,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=12
!      if (irlen.eq.0)goto 360

!      do jf=1,irlen
       
!      icc(30-i,jf+2)=irrr(jf)
!      end do
!      icc(30-i,2)=irlen
!      icc(30-i,1)=0
      
      
      
!      end do
!      do jf=3,iprod2(2)+2
!      karr(jf-2)=iprod2(jf)
!      end do
!      ilen=iprod2(2)
!      do jf=3,idd(2,2)+2
!      kbarr(jf-2)=idd(2,jf)
!      end do
!      ilen2=idd(2,2)
!      call mpmul(ilen,ilen2,ilen3)
!      do jf=1,ilen3
!      karr(jf)=kcarr(jf)
!      end do
!      ilen=ilen3
!      do jf=3,n(2)+2
!      kbarr(jf-2)=n(jf)
!      end do
!      ilen2=n(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      indz=13
!      if (irlen.eq.0)goto 360
!      do jf=1,irlen
!      icc(1,jf+2)=irrr(jf)
!      end do
!      icc(1,1)=0
!      icc(1,2)=irlen
!      itag=0
!     print *,'hello6'
            
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
      if (iesw.eq.2)goto 401
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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
      stop
      goto 72
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      
      common ipr(65000),norma(50)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
      
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
      
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


      

      
      
      subroutine subbern(idegind,isok) 
    
!     smaller sieving interval
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common ipr(65000),norma(50)
      
      
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
      common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
      common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
     common ntest(100),nfact1(100),kdcorn(200),lcorn(200),mmbig(100),jgg(100)
     common kpf(120),kqf(120),krf(120),iroota(2,120)
     
      dimension krr(120),iden(120),invcom(120),jnvcom(2,120)
      dimension kq(120),ka(120),kb(120)
      dimension kp(120),inv2(100),inv3(100),inv4(100),inv8(100)
      dimension inv16(100),inv27(100),kcub(6) 
      dimension khsq(120),itempr(120),kicr(120),ians2(2,120)
      dimension khcr(120),icuba(6,2,100)
      dimension kcapa(2,50)
      idegind=3
!      ip(1)=0
!      ip(2)=1
!      ip(3)=131
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do jf=1,kpf(2)+2
      kbarr(jf)=kpf(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kp(jf)=kcarr(jf)
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
      
      print *,'inverses 4 8 27 ',inv4(3),inv8(3),inv27(3)
      print *,'inv16',inv16(3),'inv3',inv3(3),'inv2',inv2(3)
      
      

      
! to use this routine for solving cubic equations set 2nd coeff to kpf      
! 3rd coeff to kqf , const to krf AND minus 2nd to kp
!      kp(1)=0
!      kp(2)=1
!      kp(3)=129
!      kpf(1)=0
!      kpf(2)=1
!      kpf(3)=2
!      kqf(1)=0
!      kqf(2)=1
!      kqf(3)=3
!      krf(1)=0
!      krf(2)=1
!      krf(3)=77
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
      if (khsq(2).eq.0)goto 200
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
200   do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,kb(2)+2
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
      marr(jf)=mcarr(jf)
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
      kicr(jf)=mcarr(jf)
      end do
      if (kicr(2).eq.0)goto 999
      ians(2,1)=0
      ians(2,2)=0
      ncub=0
      if (kicr(1).eq.1)goto 201
      do jf=1,kicr(2)+2
      karr(jf)=kicr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)  
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kicr(jf)=kcarr(jf)
      end do
      goto 201
60    do jf=1,khsq(2)+2 
      ncom(1,jf)=khsq(jf)
      end do
      ncom(2,1)=0
      ncom(2,2)=0
      call bwq5
      if (nsq.eq.0)goto 999
      do ibig=1,2
      do jf=1,isqurar(1,ibig,2)+2
      ians(ibig,jf)=isqurar(1,ibig,jf)
      
      end do
      print *,'ians',(ians(ibig,jk),jk=1,ians(ibig,2)+2)
      end do
      
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      
      do ibig=1,2
      if (ians(ibig,2).eq.0)goto 100
      do jf=1,ians(ibig,2)+2
      kbarr(jf)=ians(ibig,jf)
      end do
      
      call mpadd(1)
      do jf=1,kcarr(2)+2
      ians2(ibig,jf)=kcarr(jf)
      end do
      goto 101
100   ians2(ibig,1)=0 
      ians2(ibig,2)=0
101   end do
      print *,'ians2re',(ians2(1,jk),jk=1,ians2(1,2)+2)
      print *,'ians2im',(ians2(2,jk),jk=1,ians2(2,2)+2)
      
      do jf=1,kb(2)+2
      kbarr(jf)=kb(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ians(1,2)+2
      kbarr(jf)=ians(1,jf)
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
      khcr(jf)=mcarr(jf)
      end do
      
      do jf=1,khcr(2)+2
      ncom(1,jf)=khcr(jf)
      end do
      do jf=1,ians(2,2)+2
      ncom(2,jf)=ians(2,jf)
      end do
      print *,'ncom2',(ncom(2,jk),jk=1,ncom(2,2)+2)
      print *,'khcr',(khcr(jf),jf=1,khcr(2)+2)
      
      call cub5
      print *,'khcr',(khcr(jf),jf=1,khcr(2)+2)
      print *,'ians2',(ians2(2,jf),jf=1,ians2(2,2)+2),'ncubr',ncubr
      
      ncub=ncubr
      if (ncub.eq.0)goto 202
      do jbig=1,ncub
      do ibig=1,2
      do jf=1,icubar(jbig,ibig,2)+2
      icuba(jbig,ibig,jf)=icubar(jbig,ibig,jf)
      end do
      end do
      end do
202   do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      print *,'ok1'
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,kb(2)+2
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
      print *,'ok2'
      do jf=1,mcarr(2)+2
      marr(jf)=mcarr(jf)
      end do
      do jf=1,inv2(2)+2
      mbarr(jf)=inv2(jf)
      end do
      call menmul
      print *,'ok25'
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      do jf=1,ians2(1,2)+2
      kbarr(jf)=ians2(1,jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      print *,'ok3'
      
      do jf=1,ip(2)+2
      mbarr(jf)=ip(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      kicr(jf)=mcarr(jf)
      end do
      if (kicr(1).eq.0)goto 112
      do jf=1,kicr(2)+2
      karr(jf)=kicr(jf)
      end do
      do jf=1,ip(2)+2
      kbarr(jf)=ip(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kicr(jf)=kcarr(jf)
      end do
112   a=a
      print *,'khcr',(khcr(jf),jf=1,khcr(2)+2)
      print *,'kicr',(kicr(jf),jf=1,kicr(2)+2)
      print *,'ianscom',(ians(2,jf),jf=1,ians(2,2)+2)
      print *,'khsq',(khsq(jf),jf=1,khsq(2)+2)
      print *,'icuba',(icuba(1,1,jf),jf=1,icuba(1,1,2)+2)
      print *,'icubacom',(icuba(1,2,jf),jf=1,icuba(1,2,2)+2)
201   do jf=1,kicr(2)+2
      ncom(1,jf)=kicr(jf)
      end do
      do jf=1,ians2(2,2)+2
      ncom(2,jf)=ians2(2,jf)
      end do
      
      
      call cub5
      if (ncubr.eq.0)goto 110
      goto 111
110   print *,'ncubr',ncubr,'ncub',ncub      
      if (ncub+ncubr.eq.0)goto 999
      stop
111   do jbig=1,ncubr
      do ibig=1,2
      do jf=1,(icubar(jbig,ibig,2)+2)
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
      enddo
      inxyz=1




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
      kbarr(jf)=mcarr(jf)
      end do
      do jf=1,iden(2)+2
      karr(jf)=iden(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      call mendiv
      do jf=1,mcarr(2)+2
      iden(jf)=mcarr(jf)
      end do
      do jf=1,iden(2)+2
      karb(jf)=iden(jf)
      end do
      call mpgcd
      do jf=1,karv(2)+2
      invcom(jf)=karv(jf)
      end do
      do jf=1,ip(2)+2
      karr(jf)=ip(jf)
      end do
      do jf=1,kcapa(2,2)+2
      kbarr(jf)=kcapa(2,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      marr(jf)=kcarr(jf)
      end do
      do jf=1,invcom(2)+2
      mbarr(jf)=invcom(jf)
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
      jnvcom(2,jf)=mcarr(jf)
      end do
      do jf=1,jnvcom(2,2)+2
      marr(jf)=jnvcom(2,jf)
      end do
      do jf=1,kcapa(2,2)+2
      mbarr(jf)=kcapa(2,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(0)
      do jf=1,kcarr(2)+2
      itempr(jf)=kcarr(jf)
      end do
      do jf=1,kcapa(1,2)+2
      karb(jf)=kcapa(1,jf)
      end do
      call mpgcd
      do jf=1,karv(2)+2
      marr(jf)=karv(jf)
      end do
      do jf=1,itempr(2)+2
      mbarr(jf)=itempr(jf)
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
      if (inxyz.eq.1)goto 121
!      if (inxyz.eq.2)goto 122
      
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
      print *,'kkbig',kkbig,'kicr',(kicr(jf),jf=1,kicr(2)+2)
      print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
      print *,'icuba',(icuba(kkbig,1,jf),jf=1,icuba(kkbig,1,2)+2)
      print *,'icubac',(icuba(kkbig,2,jf),jf=1,icuba(kkbig,2,2)+2)
      
      if (kkbig.eq.ncub)goto 312
      if (iroota(2,2).eq.0)goto 312
      kkbig=kkbig+1
      goto 160
312   isok=0
      goto 1000
    





203   stop
999   print *,'no solutions mod',(ip(jf),jf=1,ip(2)+2)
      isok=1
1000  return
      end
      
                        
                         
      







