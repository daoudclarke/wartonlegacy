     

      program xobmpbak 
!     first of two stages of MPQS double large prime variant procedures   
!     8000 stage 1 primes and 10000 stage2 primes
!     contains small prime variant (no sieving on primes < 100      
!     also large prime variant extension      
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
      dimension kr(30000,2),mm(20),mres(20),n(50),msqr(20)
      dimension mmsq(20),mstem(20),zlog(10000)
      dimension iaeq(50),ibeq(50),iceq(50),ninv(50),ixx(5)
      dimension ibbeq(50),kkm(30000,2)
      dimension narr(2,50),nbarr(2,50),ncarr(2,50),iback(50)
      
      real krecarr(10000),lglm
!      open(unit=2,file='hquat',access='direct',form=&
!      'formatted',recl=410,status='new')
!      do ii=869,872
!      read (2,72222,rec=ii)irecnn,iconq,(narr(ii,jf),jf=1,50),&
!      (nbarr(ii,jk),jk=1,50)
!      print *,'irecnn=',irecnn,'iconq',iconq,'a=',(narr(ii,jf),jf=1,&
!       narr(ii,2)+2),'b=',(nbarr(ii,jk),jk=1,nbarr(ii,2)+2)
!       end do
!      close(unit=2)
!      stop
      ntest(1)=0
      ntest(2)=6
      ntest(3)=100
      ntest(4)=99
      ntest(5)=3000
      ntest(6)=92
      ntest(7)=1091
      ntest(8)=4553
      ntest(9)=8883
      ntest(10)=1243
      
      iok=0
      call subrho(iok)
      print *,'subrho',(nfact1(jf),jf=1,nfact1(2)+2)
      stop


      do jf=1,10000
      ibbp(jf)=999999999
      nbbp(jf)=0
      end do
      imatc=0
      iccz=0
      ijdpr=0
      iprerec=0
      iprec2=0
      irecone=0
      iconq=1
      iconq2=1
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
      end do
      iaeq(1)=0
      iaeq(2)=7
      iaeq(3)=745
      iaeq(4)=3559
      iaeq(5)=9249
      iaeq(6)=9929
      iaeq(7)=8988
      iaeq(8)=305
      iaeq(9)=7889
      
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
      n(2)=18
      n(3)=111
      n(4)=1111
      n(5)=1111
      n(6)=1111
      n(7)=1111
      n(8)=1111
      n(9)=1111
      n(10)=1111
      n(11)=1111
      n(12)=1111
      n(13)=1111
      n(14)=1111
      n(15)=1111
      n(16)=1111
      n(17)=1111
      n(18)=1111
      n(19)=1111
      n(20)=1111
      
      do jf=1,50
      ncom(jf)=n(jf)
      end do
      do jf=1,iaeq(2)+2
      iprar(jf)=iaeq(jf)
      end do
      call bigb5(noyes)
      print *,'bigans',(iansar(jf),jf=1,iansar(2)+2)
      print *,'noyes',noyes
      
      
      
      
      
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
      open(unit=4,file='tnx',access='sequential')
      open(unit=7,file='kaze',access='direct',form=&
      'formatted',recl=270000,status='old')
       open (unit=8,file='nbbp',access='direct',form=&
       'formatted',recl=30000,status='old')
81111 format (10000i3)
!      read (7,6111,rec=1)(kr(jf,1),jf=1,30000)
!      read (7,6111,rec=2)(kr(jf,2),jf=1,30000)
!      print *,'kr1',(kr(jf,1),jf=1,11)
!      print *,'kr2',(kr(jf,2),jf=1,11)
!      close (unit=3)
!      close(unit=4)
!      close(unit=7)
!      stop
      read(3,5,rec=1) (ipr(i),i=1,65000)
5     format(65000i6)      
      
      print *,'no of primes<800001=',ipr(1),ipr(2)
      
      
      rel=ipr(20000)
!     instruction changes for small prime variant      
      lglm=log10(rel)+9.0 
      key =1
      print *,'number?'
      read *,istring
      call addstar(istring,ostring)
      
      print *,ostring
      write (gstring,1)ostring
1     format(a200)
      print 1,gstring
      
      read (gstring,2) (loc(i),i=1,200)
      print *,loc(9)
      
      print *,loc(15)


2     format(200a1)      
      itemp =ichar(loc(9)) -48
      jtemp =ichar(loc(2)) -48
      ktemp =ichar(loc(10))-48
      ltemp =ichar(loc(15))-48

3     format(4i10)
      print 3,itemp,jtemp,ktemp,ltemp
      do i =1,200
      ii= 201-i
      itemp = ichar(loc(ii))
      if (itemp.le.47)goto 10
      if (itemp.ge.58)goto 10
      
      goto 12
10    end do      
      goto 100
12    jtemp =int(ii/4)      
      irem1 =ii -jtemp *4
      do i =1,50
      karr(i) =0
      end do
      do j=1,jtemp
      print *,'ii',ii
      do k =1,4
      kemp =ichar(loc(ii+1-k))-48
      karr(jtemp+2-j) =karr(jtemp+2-j) +kemp *10**(k-1)
      end do
      ii=ii-4
      end do
      IF (iREM1.GT.0)GOTO 24
      DO J =1,JTEMP
      karr(j)=karr(j+1)
      end do
      jvr=jtemp
      goto 26
24    do k =1,irem1
      kemp =ichar(loc(ii+1-k))-48
      karr(1) =karr(1) +kemp *10**(k -1)
      end do
      jvr =jtemp+1
26    print *,'karrs',(karr(i),i=1,jvr)
      if (irem1.gt.0)goto 28
      ilen =jtemp
      goto 30









28    ilen =jtemp+1  
      
      


30    print *,'ilen',ilen
      inlen =ilen
      do i=1,inlen
      mnum(i)=karr(i)
      end do
      
      ilen2 =3
      kbarr(1)=8926
      kbarr(2)=4071
      kbarr(3)=6193
      kbarr(4)=8987
      kbarr(5) =6573
      
      
      call mpmul(ilen,ilen2,ilen3)
      
      
32    print *,(kcarr(i),i=1,ilen3)
      
      
      karr(1)=7
      karr(2)=5480
      karr(3)=1247
      karr(4) =7879
      karr(5)=132
      karr(6)=4834
      karr(7)=5243
      karr(8)=3610
      karr(9)=8280
      karr(10)=5162
      karr(11)=5180
      karr(12)=5673

      kbarr(1)=5
      kbarr(2)=1234
      kbarr(3)=5678
      kbarr(4)=9123
      kbarr(5)=4631
      
      ilen =12
      ilen2=5
      do iv =1,1
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      end do
      
      if (icont.eq.0)goto 101
      print *,'quotient',(ipqt(i),i=1,icont)
101   if(irlen.eq.0)goto 102      
      print *,'remainder',(irrr(i),i=1,irlen)
      goto 102
100   print *,'input error'
102   karr(1)=1
      

      

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
      
      call mpprime(icorp,inlen)
      
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
      kard(2)=1
      kard(3)=543
      kard(4)=3
      karp(1)=0
      karp(2)=1
      karp(3)=17
      karp(4)=7
      do jf=1,n(2)+2
      kard(jf)=n(jf)
      end do
      call mpkron(k)
      print *,'k',k
      






      
      m1(1)=0
      m1(2)=1
      m1(3)=120
      m1(4)=5172
      m2(1)=0
      m2(2)=1
      m2(3)=113
      m2(4)=5169
      ia0(1)=0
      ia0(2)=1
      ia0(3)=11
      ia0(4)=2489
      ia1(1)=0
      ia1(2)=1
      ia1(3)=105
      ia1(4)=282
      ia2(1)=0
      ia2(2)=1
      ia2(3)=85
      ia2(4)=392
      ia3(1)=0
      ia3(2)=1
      ia3(3)=39
      ia3(4)=1210
      ia4(1)=0
      ia4(2)=1
      ia4(3)=118
      ia4(4)=3543
      ia5(1)=1
      ia5(2)=1
      ia5(3)=191
      ia5(4)=8312
      
      iblim=0
      iulim=2000
      lprx1=100000
      icon=0
      isivl=2000
      
      
      
      
400   do i=1,30000
      
      kr(i,1)=999999999
      kr(i,2)=999999999
      end do
      
      
      goto 4004
4002  print *,'factor found early=2'
      stop


4004  iva=n(2)+2
      kr(2,1)=mod(n(iva),2)
      if (kr(2,1).eq.0)goto 4002
      do jf=1,n(2)+2
      kard(jf)=n(jf)
      end do
      
      iccn=0
      karp(1)=0
      do i=3,20000
      do jf=1,n(2)+2
      kard(jf)=n(jf)
      end do


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
      if (ip.ne.17)goto 5901
      
5901  call bwq5(ip,ians)
      kr(i,1)=ians
      kr(i,2)=ip-ians
      iccn=iccn+1
      print *,'kr=',kr(i,1),'ipr=',ipr(i),'iccn=',iccn
      if (iccn.eq.8000)goto 512
590   end do
      print *,'loop too short'
      stop

      
      
      
      
      
      
      

512   print *,'square root phase finished'
!     change parameter      
      do i=2,20000
      if (kr(i,1).eq.999999999)goto 592
      rel=ipr(i)
      gpr(i)=log10(rel)
      goto 594
592   gpr(i)=99.0      
594   end do      
      rel=ipr(20000)
!     changes for large prime variant     
      lglm=log10(rel)+9.0 
      
      
      print *,'i',i
      print *,'krs',(kr(jf,1),jf=1,50)
      print *,'krs2',(kr(jf,2),jf=1,50)
      write(7,6111,rec=1)(kr(jf,1),jf=1,30000)
      write(7,6111,rec=2)(kr(jf,2),jf=1,30000)
6111  format(30000i9)


      print *,'end y stage'
      open(unit=2,file='iquat',access='direct',form=&
      'formatted',recl=410,status='old')
      print *,'okk'
      do jf=3,iaeq(2)+2
      mnum(jf-2)=iaeq(jf)
      end do
      inlen=iaeq(2)
55555 a=a
      
      
5122  itemp=mod(mnum(inlen),8)
      if (itemp.ne.1)goto 5121
5123  karr(1)=0
      karr(2)=inlen
      do jf=1,inlen
      karr(jf+2)=mnum(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=2
      call mpadd(0)
      inlen=kcarr(2)
      do jf=3,kcarr(2)+2
      mnum(jf-2)=kcarr(jf)
      end do
      goto 5122
5121  call mpprime(icorp,inlen)
      if (icorp.eq.0)goto 5123
      do jf=1,n(2)+2
      kard(jf)=n(jf)
      end do
      do jf=1,inlen
      karp(jf+2)=mnum(jf)
      end do
      karp(1)=0
      karp(2)=inlen
      call mpkron(k)
      if (k.lt.0)goto 5123
      do jf=1,inlen
      iprar(jf+2)=mnum(jf)
      end do
      iprar(1)=0
      iprar(2)=inlen
      
      call bigb5(noyes)
      inlen=iprar(2)
      do jf=3,iprar(2)+2
      mnum(jf-2)=iprar(jf)
      end do

      
      if (noyes.eq.1)goto 5123
      do jf=1,iansar(2)+2
      nbarr(1,jf)=iansar(jf)
      kbarr(jf)=iansar(jf)
      end do
      goto 11167
      do jf=1,iprar(2)+2
      karr(jf)=iprar(jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      nbarr(2,jf)=kcarr(jf)
      end do
11167 do jf=1,iprar(2)+2
      narr(1,jf)=iprar(jf)
      narr(2,jf)=iprar(jf)
      end do
      do ik=1,2
      do jf=3,nbarr(ik,2)+2
      karr(jf-2)=nbarr(ik,jf)
      kbarr(jf-2)=nbarr(ik,jf)
      end do
      ilen=nbarr(ik,2)
      ilen2=ilen
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
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      isgn=kcarr(1)
      do jf=3,narr(1,2)+2
      kbarr(jf-2)=narr(1,jf)
      end do
      ilen2=narr(1,2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      ncarr(ik,jf+2)=ipqt(jf)
      end do
      ncarr(ik,2)=icont
      ncarr(ik,1)=isgn
      end do
      
      


!     start of big loop      
      do ixyz=1,1
      do jf=1,narr(1,2)+2
      iaeq(jf)=narr(1,jf)
      end do
      do jf=1,nbarr(ixyz,2)+2
      ibeq(jf)=nbarr(ixyz,jf)
      end do
      do jf=1,ncarr(ixyz,2)+2
      iceq(jf)=ncarr(ixyz,jf)
      end do
      do jf=1,ibeq(2)+2
      karr(jf)=ibeq(jf)
      kbarr(jf)=ibeq(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      ibbeq(jf)=kcarr(jf)
      end do
      read(7,6111,rec=1)(kr(jf,1),jf=1,20000)
      read(7,6111,rec=2)(kr(jf,2),jf=1,20000)



!     finding solutions to particular equation for all mods     
      do iz=2,20000
      if (kr(iz,1).eq.999999999)goto 899
!     compute inverse of 2*a
      
      
      do jf=3,iaeq(2)+2
      karr(jf-2)=iaeq(jf)
      end do
      ilen=iaeq(2)
      
      if (ipr(iz).lt.10000)goto 800
      kbarr(1)=ipr(iz)/10000
      kbarr(2)=ipr(iz)-kbarr(1)*10000
      ilen2=2
      
      goto 802
800   kbarr(1)=ipr(iz)
      ilen2=1
802   call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,irlen
      karb(jf+2)=irrr(jf)
      end do
      karb(1)=0
      karb(2)=irlen
      do jf=1,ilen2
      kara(jf+2)=kbarr(jf)
      iprar(jf+2)=kbarr(jf)
      end do
      kara(1)=0
      kara(2)=ilen2
      iprar(1)=0
      iprar(2)=ilen2
      
      call mpgcd
      do jf=1,karv(2)+2
      ninv(jf)=karv(jf)
      end do
!     determine numerator
      do jf=3,ibeq(2)+2
      karr(jf-2)=ibeq(jf)
      end do
      ilen=ibeq(2)
      do jf=3,iprar(2)+2
      kbarr(jf-2)=iprar(jf)
      end do
      ilen2=iprar(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      isum=0
      do jf=irlen,1,-1
      isum=isum+irrr(jf)*10000**(irlen-jf)
      end do
      item=ipr(iz)-isum+kr(iz,1)
      ibx=isum
      do ind1=1,2
      if (kr(iz,ind1).eq.999999999)goto 816
      if (item.lt.10000)goto 804
      karr(1)=item/10000
      karr(2)=item-karr(1)*10000
      ilen=2
      goto 806
804   karr(1)=item
      ilen=1
806   do jf=3,ninv(2)+2
      kbarr(jf-2)=ninv(jf)
      end do
      ilen2=ninv(2)
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
      if (irlen.eq.0)goto 810
      isum=0
      do jf=irlen,1,-1
      isum=isum+irrr(jf)*10000**(irlen-jf)
      end do
      
      kr(iz,ind1)=isum
      
      goto 812
810   kr(iz,ind1)=0
812   item=2*ipr(iz)-ibx+kr(iz,2)
      if (item.lt.0)goto 814
      goto 816
814   item=item+ipr(iz)
816   end do
899   end do
      
     goto 820










      
      

820   do i=3,20000
      if (kr(i,1).eq.999999999)goto 5160
      if (kr(i,1).lt.kr(i,2))goto 5160
      itemp=kr(i,1)
      kr(i,1)=kr(i,2)
      kr(i,2)=itemp
5160  end do
      









      print *,'end of inversions' 

      icdn=-20000
      



      
      
!     instruction changes for small prime variant      
599   do jz=2,20000
      
      if (gpr(jz).gt.90.0)goto 650
      
      if (icdn.gt.-20000)goto 5991
      markl=mod(200000000+kr(jz,1),ipr(jz))
      if (markl.gt.10000)goto 6501
      goto 606
      
5991  markl=kkm(jz,1)
      if (markl.gt.10000)goto 610
      

      
      
606   krecarr(markl)=krecarr(markl)-gpr(jz)
608   markl=markl+ipr(jz)
      
      if (markl.gt.10000)goto 610
      krecarr(markl)=krecarr(markl)-gpr(jz)
      goto 608
610   kkm(jz,1)=markl-10000
6101  if (ipr(jz).eq.2)goto 650
      if (kr(jz,2).eq.999999)goto 650
      if (icdn.gt.-20000)goto 6141
      markr=mod(200000000+kr(jz,2),ipr(jz))
      if (markr.gt.10000)goto 6502
      goto 6142
      
6141  markr=kkm(jz,2)
612   if (markr.gt.10000)goto 6502
6142  krecarr(markr)=krecarr(markr)-gpr(jz)
      markr=markr+ipr(jz)
      
      goto 612
614   markr=kr(jz,2)
      goto 612
6501  kkm(jz,1)=markl-10000
      goto 6101
6502  kkm(jz,2)=markr-10000
650   end do
      
!     start of last phase
      
      
668   do kz=1,10000
      
!     instruction changes for small prime variant  limit incr. by 4    
      if (krecarr(kz).gt.99.0)goto 690
!     compute accurately      
      if (icdn.lt.0)goto 7000
      ixx(1)=0
      if (kz.eq.10000)goto 910
      iddn=icdn
      kzx=kz
      if (iddn.eq.0)goto 900
901   if (iddn.lt.10000)goto 902
      ixx(2)=3
      ixx(3)=iddn/10000
      ixx(4)=iddn-ixx(3)*10000 
      ixx(5)=kzx
      goto 920
902   ixx(2)=2
      ixx(3)=iddn 
      ixx(4)=kzx
      goto 920
900   ixx(2)=1
      ixx(3)=kzx
      goto 920
910   kzx=0
      iddn=icdn+1
      goto 901
920   do jf=3,iaeq(2)+2
      karr(jf-2)=iaeq(jf)
      end do
      ilen=iaeq(2)
      do jf=3,ixx(2)+2
      kbarr(jf-2)=ixx(jf)
      end do
      ilen2=ixx(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
      do jf=1,ibbeq(2)+2
      kbarr(jf)=ibbeq(jf)
      end do
      call mpadd(0)
      
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,ixx(2)+2
      kbarr(jf-2)=ixx(jf)
      end do
      ilen2=ixx(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=0
      karr(2)=ilen3
      do jf=1,iceq(2)+2
      kbarr(jf)=iceq(jf)
      end do
      call mpadd(0)
9777  do jf=1,kcarr(2)+2
      norma(jf)=kcarr(jf)
      end do
      rel=zlog(norma(3))+4*norma(2)-4
      rel2=rel-127.0 +krecarr(kz)
      if (rel2.gt.lglm)goto 690
      kkb=0
      kia=kz
      call sieve(kia,kkb,irecnn,kkmx,icdn,iconq,iccz,imatc,&
      ijdpr,idpr1,idpr2,jdpr1,jdpr2)
      if (kkb.eq.0)goto 690
      do jf=3,ixx(2)+2
      karr(jf-2)=ixx(jf)
      end do
      ilen=ixx(2)
      
      
      do jf=3,iaeq(2)+2
      kbarr(jf-2)=iaeq(jf)
      end do
      ilen2=iaeq(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      if (icdn.ge.0)goto 61111
      karr(1)=1
      goto 61112
      
61111 karr(1)=0
61112 karr(2)=ilen3
      do jf=1,ibeq(2)+2
      kbarr(jf)=ibeq(jf)
      end do
      call mpadd(0)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      iback(jf)=kcarr(jf)
      end do
      iback(1)=0
      iback(2)=kcarr(2)
      ilen=kcarr(2)
      ilen2=2
      kbarr(1)=idpr1
      kbarr(2)=idpr2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      isum=0
      do jf=irlen,1,-1
      isum=isum+irrr(jf)*10000**(irlen-jf)
      end do
      
      kr(20000+iccz-kkb+1,1)=isum
      kr(20000+iccz-kkb+1,2)=idpr1*10000+idpr2-isum
      if (kkb.eq.1)goto 690
      do jf=3,iback(2)+2
      karr(jf-2)=iback(jf)
      end do
      ilen=iback(2)
      kbarr(1)=jdpr1
      kbarr(2)=jdpr2
      ilen2=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      isum=0
      do jf=irlen,1,-1
      isum=isum+irrr(jf)*10000**(irlen-jf)
      end do
      kr(20000+iccz,1)=isum
      kr(20000+iccz,2)=jdpr1*10000+jdpr2-isum
      goto 690
      

7000  ixx(1)=0
      if (kz.eq.10000)goto 7010      
      iddn=abs(icdn+1)
      kzx=10000-kz
7001  if (iddn.eq.0)goto 7090
      if (iddn.lt.10000)goto 7002
      ixx(2)=3
      ixx(3)=iddn/10000
      ixx(4)=iddn-ixx(3)*10000
      ixx(5)=kzx
      goto 7020
7002  ixx(2)=2
      ixx(3)=iddn
      ixx(4)=kzx
      goto 7020
7010  iddn=abs(icdn+1)
      kzx=0
      goto 7001
7090  ixx(2)=1
      ixx(3)=kzx

7020  do jf=3,iaeq(2)+2
      karr(jf-2)=iaeq(jf)
      end do
      ilen=iaeq(2)
      do jf=3,ixx(2)+2
      kbarr(jf-2)=ixx(jf)
      end do
      ilen2=ixx(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=1
      karr(2)=ilen3
      do jf=1,ibbeq(2)+2
      kbarr(jf)=ibbeq(jf)
      end do
      call mpadd(0)
      isgn=kcarr(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,ixx(2)+2
      kbarr(jf-2)=ixx(jf)
      end do
      ilen2=ixx(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      karr(1)=mod(isgn+1,2)
      karr(2)=ilen3
      do jf=1,iceq(2)+2
      kbarr(jf)=iceq(jf)
      end do
      call mpadd(0)
      goto 9777
690   end do
      icdn=icdn+1
      
      do jf=1,10000
      krecarr(jf)=127.0
      end do
      if (icdn.eq.20001)goto 720
      goto 599










720   print *,'maximum no of primes=',kkmx,'ixyz=',ixyz
      print *,'irecnn=',irecnn,'irecpl=',irecpl
      if (iprerec.eq.irecnn)goto 69001
      irecdif=irecnn-iprerec
      if (irecdif.ne.1)goto 69002
      irecone=irecone+1
69002 iprerec=irecnn
      print *,'irecone=',irecone
      iconq=iconq+1
      write (2,72222,rec=iconq)irecnn,iconq,(narr(ixyz,jf),jf=1,50),(nbarr&
      (ixyz,jk),jk=1,50)
      if (iccz.ge.3000)goto 730
72222 format(i6,i4,50i4,50i4)      
!     end of big loop      
69001 a=a       
      end do
      if (iprec2.eq.irecnn)goto 7303
      iprec2=irecnn 
      iconq2=iconq2+1


      
7303  do jf=1,narr(1,2)+2
      karr(jf)=narr(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=2
      call mpadd(0)
      inlen=kcarr(2)
      do jf=1,inlen
      mnum(jf)=kcarr(jf+2)
      end do
      goto 55555
730   print *,'run completed',' no of records=',irecnn,'irecone=',irecone
      read (2,72222,rec=1)irecnn,ispac,(narr(1,jf),jf=1,50),&
      (nbarr(1,jk),jk=1,50)
      write (2,72222,rec=1)irecnn,iconq,(narr(1,jf),jf=1,50),&
      (nbarr(1,jk),jk=1,50)
      read (7,6111,rec=1)(kr(jf,1),jf=1,20000)
      read (7,6111,rec=2)(kr(jf,2),jf=1,20000)
      do i=20001,30000
      if (kr(i,1).eq.999999999)goto 7111
      if (kr(i,1).lt.kr(i,2))goto 7111
      itemp=kr(i,1)
      kr(i,1)=kr(i,2)
      kr(i,2)=itemp
7111  end do      
      print *,'iccz',iccz,'imatc',imatc,'ijdpr',ijdpr
      print *,'iaeq',(iaeq(jf),jf=1,iaeq(2)+2)
      write (7,6111,rec=1)(kr(jf,1),jf=1,30000)
      write (7,6111,rec=2)(kr(jf,2),jf=1,30000)
      write (8,81111,rec=1)(nbbp(jf),jf=1,10000)
      close (unit=8)
      close (unit=2)
      close(unit=3)
      
      
      
      
      
      
      close(unit=4)
      close(unit=7)
      end
      
      subroutine sieve(kia,kkb,irecnn,kkmx,icdn,iconq,iccz,imatc,&
      ijdpr,idpr1,idpr2,jdpr1,jdpr2)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(50),litd(50),ity(50),iabp(50),iabpn(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(50)
      dimension littr(50),normar(50),nfact2(50)
      
      
!     change parameter      
11    lprx1=ipr(20000)                                      
      
      
      
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
      ,'iconq',iconq
      do i=1,50

      iprex(i)=0
      end do
      iprex(2)=litd(2)
      ipfar(1)=1
      ipfar(2)=0
1010  iab=iab+1
      if (iab.eq.20000)goto 1490
      if (gpr(iab).gt.90.0)goto 1010
       
      iprx1=ipr(iab)
      
      iarp(1)=int(iprx1/10000)
      iarp(2)=iprx1-iarp(1) *10000
      ibsw=0
      
      
      
      if (iarp(1).gt.larp(1))goto 1490
      if (iarp(1).lt.larp(1))goto 1162
      if (iarp(2).gt.larp(2))goto 1490
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
      if (icorp.eq.1)goto 1490
      do i=1,litd(2)+2
      iprex(i)=litd(i)
      end do
1204  ipfar(1)=ipfar(1) +2
      goto 1166


      

1300  irecnn=irecnn+1
      print *,'hit no',irecnn,'a=',kia,'b=',kkb,(litt(i),i=1,litt(2) +2)&
      ,'iconq',iconq ,'iccz',iccz,'ijdpr',ijdpr,'imatc',imatc
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
      if (irecnn.lt.40000)goto 1500
13861 close(unit=4)
      print *,'max no of primes=',kkmx
      stop
1490  if (litd(2).gt.4)goto 1500      
      if ((litd(2).eq.2).and.(litd(3).le.120))goto 2000
      if (litd(2).lt.4)goto 2002
      if ((litd(3).eq.1).and.(litd(4).le.4400))goto 2002
      goto 1500
2002  inlen=litd(2)
      do jf=3,litd(2)+2
      mnum(jf-2)=litd(jf)
      end do
      call mpprime(icorp,inlen)
      if (icorp.eq.1)goto 1500
      do jf=1,litd(2)+2
      ntest(jf)=litd(jf)
      end do
      call subrho(iok)
      if (iok.eq.1)goto 1500
      do jf=3,ntest(2)+2
      karr(jf-2)=ntest(jf)
      end do
      ilen=ntest(2)
      do jf=3,nfact1(2)+2
      kbarr(jf-2)=nfact1(jf)
      end do
      ilen2=nfact1(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      do jf=1,icont
      nfact2(jf+2)=ipqt(jf)
      end do
      nfact2(1)=0
      nfact2(2)=icont
      print *,'nfact1',(nfact1(jf),jf=1,nfact1(2)+2)
      print *,'nfact2',(nfact2(jf),jf=1,nfact2(2)+2)

      
      
      
      if ((nfact1(2).gt.2).or.(nfact2(2).gt.2))goto 1500
      if ((nfact1(3).gt.120).or.(nfact2(3).gt.120))goto 1500
      icur=icur+1
      iabp(icur)=nfact1(3)*10000+nfact1(4)
      iabpn(icur)=1
      ity(icur)=1
      itz=iabp(icur)
      idpr1=nfact1(3)
      idpr2=nfact1(4)
      jdpr1=nfact2(3)
      jdpr2=nfact2(4)
      
      
      itz2=nfact2(3)*10000+nfact2(4)
      if (itz.eq.itz2)goto 2006
      
      icur=icur+1
      iabp(icur)=itz2
      iabpn(icur)=1
      ity(icur)=1
      
      ijdpr=ijdpr+1
      

2006  if (iccz.eq.0)goto 20201
      do jf=1,iccz
      if (itz.ne.ibbp(jf))goto 2007
      nbbp(jf)=nbbp(jf)+1
      imatc=imatc+1
      kkb=0
      goto 2008
2007  end do
20201 iccz=iccz+1
      ibbp(iccz)=itz
      nbbp(iccz)=nbbp(iccz)+1
      kkb=1
      print *,'ibbp',itz
      if (itz.eq.itz2)goto 2011
      goto 2012
2008  if (itz.eq.itz2)goto 2011
      goto 2012

      

      
2011  nbbp(iccz)=nbbp(iccz)+1      
      iabpn(icur)=iabpn(icur)+1
      goto 1300
2012  do jf=1,iccz
      if (itz2.ne.ibbp(jf))goto 2019
      nbbp(iccz)=nbbp(iccz)+1
      imatc=imatc+1
      goto 1300
2019  end do
      iccz=iccz+1
      ibbp(iccz)=itz2
      nbbp(iccz)=nbbp(iccz)+1
      kkb=kkb+1
      print *,'ibbp',itz2
      
      
      
      
      goto 1300
2000  icur=icur+1
      iabp(icur)=litd(3)*10000+litd(4)
      iabpn(icur)=iabpn(icur)+1
      ity(icur)=1
      if (iccz.eq.0)goto 1499
      do jf=1,iccz
      if (iabp(icur).ne.ibbp(jf))goto 1498
      nbbp(jf)=nbbp(jf)+1
      imatc=imatc+1
      goto 1300
1498  end do
1499  iccz=iccz+1
      ibbp(iccz)=iabp(icur)
      nbbp(iccz)=nbbp(iccz)+1
      kkb=1
      idpr1=litd(3)
      idpr2=litd(4)
      print *,'ibbp',ibbp(iccz)
      goto 1300


1500  return
      

      
      end
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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

      subroutine subgcd2
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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


      subroutine mpgcd
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
      common ipr(65000),norma(50),gpr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common ncom(50),iprar(50),iansar(50),narc(50)
      common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
       common ipr(65000),norma(50),gpr(65000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
       common ipr(65000),norma(50),gpr(65000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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
       common ipr(65000),norma(50),gpr(65000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
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


       subroutine subrho(iok)
       common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(50),gpr(65000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common ncom(50),iprar(50),iansar(50),narc(50)
       common ibbp(10000),nbbp(10000),ntest(50),nfact1(50)
       dimension ipre(50),nsub(50),itest(50)
       print *,'ntest',(ntest(jf),jf=1,ntest(2)+2)
       ipolc=1
       idiff=2
       icon=4
       ipre(1)=0
       ipre(2)=1
       ipre(3)=677
1      do jf=1,ipre(2)+2 
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
       do jf =1,irlen
       ipre(jf+2)=irrr(jf)
       end do
       ipre(1)=0
       ipre(2)=irlen
       if (i.le.idiff)goto 3
       do jf=1,ipre(2)+2
       karr(jf)=ipre(jf)
       end do
       do jf=1,nsub(2) +2
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
2      do jf=1,irlen
       itest(jf+2)=irrr(jf)
       end do
       itest(1)=0
       itest(2)=irlen
25     do jf=1,ntest(2)+2
       kara(jf)=ntest(jf)
       end do
       do jf=1,itest(2)+2
       karb(jf)=itest(jf)
       end do
       call subgcd2
       if (kard(2).gt.1)goto 90
       if (kard(3).gt.1)goto 90
       
3      end do       
       idiff=icon
       icon=icon+icon
       print *,'icon',icon
       if (icon.gt.10000000)goto 999
       goto 1
90     do jf=1,kard(2)+2
       nfact1(jf)=kard(jf)
       end do
       iok=0
       goto 100
99     ipolc=ipolc+1
       if (ipolc.gt.9000)goto 999
       idiff=2
       icon=4
       ipre(1)=0
       ipre(2)=1
       ipre(3)=676+ipolc
       goto 1
999    iok=1
100    return
       end









