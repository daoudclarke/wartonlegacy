      
      program bern1 
!     first stage multi-precisioning of "bern",computes multi-precision      
!     complex cube and square roots mod p

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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
!      dimension ka3(50),ka2(50),ka1(50),kc(50),ks1(50),ks2(50)
!      dimension ks3(50),ks4(50),ks12(50)
!      dimension kp(50),inv2(50),inv3(50),inv4(50),inv8(50)
!      dimension inv16(50),inv27(50),inv256(50),kcub(6),kzans(30),ksol(4) 
      dimension ipsm(30),ntp(200),ixr(200),ibb(200),ipre(200),icurr(200)
      dimension lcorn(200),kdcorn(200),isqrtd(200),kacorn(200)
      dimension kbcorn(200),krcorn(200),ktemp1(200),kchek(200),ktest(200)
      goto 100
!      goto 778
      iaas(1,1)=0
      iaas(1,2)=1
      iaas(1,3)=2
      iaas(2,1)=0
      iaas(2,2)=0
      ipn(1)=0
      ipn(2)=1
      ipn(3)=15
      ipn(4)=3666
      ip(1)=0
      ip(2)=2
      ip(3)=100
      ip(4)=999
      call sub516
      print *,'ibprod',(ibprod(1,jf),jf=1,ibprod(1,2)+2)
      print *,'2ibprod',(ibprod(2,jf),jf=1,ibprod(2,2)+2)
      stop

      
      ip(1)=0
      ip(2)=1
      ip(3)=2029
      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=8
      ncom(2,1)=0
      ncom(2,2)=0
      call cub5
      stop
!      goto 777
      goto 778
100   ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=53
      ncom(2,1)=0
      ncom(2,2)=0
      ip(1)=0
      ip(2)=1
      ip(3)=449
      ip(4)=5678
      ip(5)=9012
      ip(6)=3456
      ip(7)=7890
      ip(8)=1234
      ip(9)=5678
      ip(10)=9097


      call cornsq
      stop

!      goto 777
778   ip(1)=0
      ip(2)=5
      ip(3)=877
      ip(4)=0
      do i=5,8
      ip(i)=0
      end do
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
      print *,'marr',(marr(jf),jf=1,marr(2)+2)
      mbarr(1)=0
      mbarr(2)=1
      mbarr(3)=2
      call menmul
      do jf=1,mcarr(2)+2
      karr(jf)=mcarr(jf)
      end do
      print *,'karr',(karr(jf),jf=1,karr(2)+2)
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
      print *,'lcorn',(lcorn(jf),jf=1,lcorn(2)+2)
!  must be negative discriminant here      
      ip(2)=ip(2)-4
      print *,'ip',(ip(jf),jf=1,ip(2)+2)
      
      kdcorn(1)=1
      kdcorn(2)=1
      kdcorn(3)=708
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
      call mpkron(k)
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
      print *,'ivan',ivan,'jvan',jvan
      print *,'isq',(isqurar(1,1,jf),jf=1,isqurar(1,1,2)+2)
      print *,'kdcorn',(kdcorn(jf),jf=1,kdcorn(2)+2)
      
      
      if (mod(isqurar(1,1,ivan),2).ne.mod(kdcorn(jvan),2))goto 11
      do jf=1,isqurar(1,1,2)+2
      isqrtd(jf)=isqurar(1,1,jf)
      end do
      print *,'okstop'
      
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
      print *,'kbcorn',(kbcorn(jf),jf=1,kbcorn(2)+2)
      

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
      print *,'2kbcorn',(kbcorn(jf),jf=1,kbcorn(2)+2)
       
      
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
      print *,'2marr',(marr(jf),jf=1,marr(2)+2)
      
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
      marr(jf)=ipre(jf)
      end do
      print *,'sq. root=',(icurr(jf),jf=1,icurr(2)+2)
      print *,'marr',(marr(jf),jf=1,marr(2)+2)
      
      do jf=1,icurr(2)+2
      marr(jf)=icurr(jf)
      mbarr(jf)=icurr(jf)
      end do
      call menmul
      print *,'okb'
      do jf=1,mcarr(2)+2
      kchek(jf)=mcarr(jf)
      end do
      do jf=2,kchek(2)+2
      print *,'kchek',kchek(jf),'ktest',ktest(jf),'jf',jf
      if (kchek(jf).ne.ktest(jf))goto 96
      end do
      print *,'solution x=',(kbcorn(jf),jf=1,kbcorn(2)+2)
      print *,'solution y=',(icurr(jf),jf=1,icurr(2)+2)
      stop
96    print *,'no solution: square test fails'
      stop
97    print *,'no solution: division test fails'
      stop
98    print *,'complex square root'
      stop
99    print *,'kronecker value=',k
      stop
       


      stop


777   idb=1
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
      



      ip(1)=0
      ip(2)=1
      ip(3)=8209
      ip(4)=0
      ip(5)=0
      ip(6)=919
      
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
      print *,'ok1'
      
      
      
      
      
      
      
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
     print *,'ok2'
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
     stop 
750  print *,(ip(jf),jf=1,ip(2)+2),'is very probably prime'     
     stop

      ncom(1,1)=0
      ncom(1,2)=1
      ncom(1,3)=264
      ncom(2,1)=1
      ncom(2,2)=1
      ncom(2,3)=259
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


1000  end
      
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
      common ncom(2,100),iprar(50),iansar(50),narc(50)
      common ip(100),ians(2,100)
      common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
      common ibprod(2,100),icprod(2,100)
      common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
      common marr(200),mbarr(200),mcarr(400),mdarr(200)


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
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(1,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100)                               
       dimension ibbb(100),ittt(100),iyyy(100)
       dimension ipow(100),ixxx(100)
       
       
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
       end do
       
       
       

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
       ibig=ibig-1
       jrrr=ibig
       
       print *,'ibig',ibig
       
       
       do jf=1,iprecod(2)+2
       karr(jf)=iprecod(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2-irem
!       call mpadd(0)
       call mpadd(1)
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
       
       do jf=1,iaa(1,2)+2
       marr(jf)=iaa(1,jf)
       end do
       do jf=1,ix(1,2)+2
       mbarr(jf)=ix(1,jf)
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
       do jf=1,ix(1,2)+2
       mbarr(jf)=ix(1,jf)
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
       ibbb(jf)=mcarr(jf)
       end do
       do jf=1,iaa(1,2)+2
       marr(jf)=iaa(1,jf)
       end do
       do jf=1,ix(1,2)+2
       mbarr(jf)=ix(1,jf)
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
       ixxx(jf)=mcarr(jf)
       end do
       



        
132    nn=nn*607        
       nn=mod(nn,5003)
       if (ip(2).gt.1)goto 403
       nn=mod(nn,ip(3))

403    kard(1)=0
       kard(2)=1
       kard(3)=nn
       
       do jf=1,ip(2)+2
       karp(jf)=ip(jf)
       end do
       call mpkron(k)
       if (k.ne.-1)goto 132
       


       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
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
!       if (ibsw.eq.1)goto 1620
112    call sub516
       print *,'ok4'
       do ibig=1,1
       do jf=1,ibprod(ibig,2)+2
       iz(ibig,jf)=ibprod(ibig,jf)
       iyyy(jf)=ibprod(ibig,jf)
       end do
       end do
       print *,'iyyy',(iyyy(jf),jf=1,iyyy(2)+2) ,'nn',nn,'ipn',(ipn(jf)&
       ,jf=1,ipn(2)+2)
       print *,'ibbb',(ibbb(jf),jf=1,ibbb(2)+2)
       print *,'ixxx',(ixxx(jf),jf=1,ixxx(2)+2)
       
402    if ((ibbb(2).eq.1).and.(ibbb(3).eq.1))goto 228
       do jf=1,ibbb(2)+2
       ipow(jf)=ibbb(jf)
       end do
       
       do ibig=1,jrrr
       do jf=1,ipow(2)+2
       marr(jf)=ipow(jf)
       mbarr(jf)=ipow(jf)
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
       ipow(jf)=mcarr(jf)
       end do
       if ((ipow(2).eq.1).and.(ipow(3).eq.1))goto 401
       end do
       print *,'error stop z'
       stop
401    if (ibig.eq.jrrr)goto 200
       print *,'ibig2',ibig
       
       ipn(1)=0
       ipn(2)=1
       ipn(3)=jrrr-ibig-1
       if (ipn(3).eq.0)goto 404
       iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=2
       iaas(2,1)=0
       iaas(2,2)=0
       call sub516
       do jf=1,ibprod(1,2)+2
       ipn(jf)=ibprod(1,jf)
       end do
       do jf=1,iyyy(2)+2
       iaas(1,jf)=iyyy(jf)
       end do
       iaas(2,1)=0
       iaas(2,2)=0
       call sub516
       goto 405
404    do jf=1,iyyy(2)+2
       ibprod(1,jf)=iyyy(jf)
       end do
             
405    do jf=1,ibprod(1,2)+2
       ittt(jf)=ibprod(1,jf)
       marr(jf)=ibprod(1,jf)
       mbarr(jf)=ibprod(1,jf)
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
       iyyy(jf)=mcarr(jf)
       end do
       print *,'secitt',(ittt(jf),jf=1,ittt(2)+2)
       print *,'seciyyy',(iyyy(jf),jf=1,iyyy(2)+2)
       
       jrrr=ibig
       do jf=1,ixxx(2)+2
       marr(jf)=ixxx(jf)
       end do
       do jf=1,ittt(2)+2
       mbarr(jf)=ittt(jf)
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
       ixxx(jf)=mcarr(jf)
       end do
       do jf=1,ibbb(2)+2
       marr(jf)=ibbb(jf)
       end do
       do jf=1,iyyy(2)+2
       mbarr(jf)=iyyy(jf)
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
       ibbb(jf)=mcarr(jf)
       end do
       goto 402


       
       
200    print *,'no square root exists'
       stop

228    print *,'square root=',(ixxx(jf),jf=1,ixxx(2)+2),'ipn',ipn(3)
!       print *,'square root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       print *,'jrrr',jrrr
       do ibig=1,1
       do jf=1,ixxx(2)+2
       isqurar(1,ibig,jf)=ixxx(jf)
       end do
       end do
       nsq=1
       
       goto 300
280    print *,'no square root exists',' iconz',iconz
       stop
       
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
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       dimension itempz(2),iaa(2,100),inv2(100),minusb(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100),nn(5)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
       iconz=1
       nn(1)=0
       nn(2)=1
       nn(3)=1
!       nn=1
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
       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
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
132    do jf=1,nn(2)+2
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
       iaas(1,jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
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
       iaas(2,jf)=mcarr(jf)
       end do



!132    nn=nn*607        
!       nn=mod(nn,1000)
!       iaas(1,1)=0
!       iaas(1,2)=1
!       iaas(1,3)=nn
!       nn=nn*607
!       nn=mod(nn,1000)
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
       print *,'nsq',nsq
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
       stop
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

       
       stop
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
232    stop
300    return
       end

       
       
          
       
       
       
       subroutine sub516
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
       print *,'ok2'
       
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
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)

       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
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
!200    print *,'icprod1',(icprod(1,jf),jf=1,icprod(1,2)+2)
!       print *,'icprod2',(icprod(2,jf),jf=1,icprod(2,2)+2)
!       print *,'iacn1',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',&
!       (iacn(2,jf),jf=1,iacn(2,2)+2)
!       print *,'ibcn1',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',&
!       (ibcn(2,jf),jf=1,ibcn(2,2)+2)
       
       return
       end
       
       
       
       subroutine bigb5(noyes)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100),ians(2,100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),ncubr,icubar(6,2,100)
       common marr(200),mbarr(200),mcarr(400),mdarr(200)
       dimension ibase(50),ipow(50),ibst(50),nnum(50),ie(50)
       
       ibigsw=0
       do jf=3,ncom(1,2)+2
       karr(jf-2)=ncom(1,jf)
       end do
       ilen=ncom(1,2)
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











