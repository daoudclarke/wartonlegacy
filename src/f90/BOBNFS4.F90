      

      program bobnfs4
! generalisation of bobmid
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
      dimension kron(100)
      irecnn=0
      ipoin=0
      open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
!      open(unit=4,file='nfsf2',access='sequential')




      
      read(3,5,rec=1) (ipr(i),i=1,65000)
5     format(65000i6)      
      
      print *,'no of primes<800001=',ipr(1)
      close(3)
      open (unit=5,file='nfspar',access='direct',form=&
      'formatted',recl=488,status='old')
      read (5,501,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),(ia4(jf),jf=1,10),&
      (ia3(jf),jf=1,10),(ia2(jf),jf=1,10),(ia1(jf),jf=1,10),&
      (ia0(jf),jf=1,10),(m1(jf),jf=1,10),(m2(jf),jf=1,10),irecper,klim,kkll,&
      izz4,kmx1,izz6
      print *,'klim',klim,'irecper',irecper,'kkll',kkll
!      klim=2000
      write (5,501,rec=1)(n(jf),jf=1,30),(ia5(jf),jf=1,10),(ia4(jf),jf=1,10),&
      (ia3(jf),jf=1,10),(ia2(jf),jf=1,10),(ia1(jf),jf=1,10),&
      (ia0(jf),jf=1,10),(m1(jf),jf=1,10),(m2(jf),jf=1,10),irecper,klim,kkll,&
      izz4,kmx1,izz6
      
      kkll2=13+kkll
501   format (30i4,80i4,6i8)      
      open(unit=1,file='nfsf4',access='direct',form&
      ='formatted',recl=kkll2,status='old')
      read(1,1,rec=4200)match(i*2-1),match(i*2),isgn,(nar(jf),jf=1,kkll)
      
      print *,(nar(jf),jf=kkll-51,kkll) 
      
      do i=1,irecper
      read(1,1,rec=i)match(i*2-1),match(i*2),isgn,(nar(jf),jf=1,kkll)
1     format (i6,i6,i1,11000i1)      
      
      end do
!      read(2,2,rec=1)(match(jf),jf=1,400)
      kkll4=12*irecper
      open (unit=7,file='nfsf5',access='direct',form=&
      'formatted',recl=kkll4,status='old')
701   format (10000(i6,i6))
      write (7,701,rec=1)(match(jf),jf=1,2*irecper)
      print *,'match8399',match(8399),'match8400',match(8400)
      

2     format(400i5)      
      kara(1)=0
      kara(2)=3
      kara(3)=1
      kara(4)=2345
      kara(5)=6791
      karb(1)=0
      karb(2)=1
      karb(3)=1234
      call mpgcd
      print *,(karv(i),i=1,karv(2)+2)
      print *,'gcd=',(kard(i),i=1,kard(2) +2)
      kia =199
      kib=70
      kbarr(1)=79
      kbarr(2)=117
      ilen2=2
      karp(1)=0
      karp(2)=2
      karp(3)=79
      
      karp(4)=6139
      karr(1)=kib
      ilen=1
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=0
      karr(1)=0
      karr(2)=1
      karr(3)=kia

      call mpadd(1)
      do i=1,kcarr(2)+2
      kard(i)=kcarr(i)
      end do
      call mpkron(k)
      kard(3)=2
      kard(4)=8818
      kard(1)=0
      kard(2)=2
      karp(1)=0
      karp(2)=1
      karp(3)=617
      call mpkron(k)
      print *,'ok up to here'
      ilen=ia5(2)
      do jf=3,ia5(2)+2
      karr(jf-2)=ia5(jf)
      end do
      ilen2=1
      kbarr(1)=5
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      idrv5(jf+2)=kcarr(jf)
      end do
      idrv5(2)=ilen3
      idrv5(1)=ia5(1)
      if (ia4(2).eq.0)goto 101
      do jf=3,ia4(2)+2
      karr(jf-2)=ia4(jf)
      end do
      ilen=ia4(2)
      ilen2=1
      kbarr(1)=4
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      idrv4(jf+2)=kcarr(jf)
      end do
      idrv4(2)=ilen3
      idrv4(1)=ia4(1)
      goto 102
101   idrv4(1)=0
      idrv4(2)=0
102   if (ia3(2).eq.0)goto 103
      do jf=3,ia3(2)+2
      karr(jf-2)=ia3(jf)
      end do
      ilen=ia3(2)
      ilen2=1
      kbarr(1)=3
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      idrv3(jf+2)=kcarr(jf)
      end do
      idrv3(2)=ilen3
      idrv3(1)=ia3(1)
      goto 104
103   idrv3(1)=0 
      idrv3(2)=0
104   if (ia2(2).eq.0)goto 105
      do jf=3,ia2(2)+2
      karr(jf-2)=ia2(jf)
      end do
      ilen=ia2(2)
      ilen2=1
      kbarr(1)=2
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      idrv2(jf+2)=kcarr(jf)
      end do
      idrv2(2)=ilen3
      idrv2(1)=ia2(1)
      goto 106
105   idrv2(1)=0
      idrv2(2)=0
106   do jf=1,ia1(2)+2
      idrv1(jf)=ia1(jf)
      end do
     print *,'ok up to here 2'
      
      
      call kieve(1,1,789003,0,klim,kkll,irecper)
      close(unit=1)
      
      close(unit=8)
      close (unit=5)
      close (unit=7)
      
      end

      
      subroutine kieve(kia,kib,iqrx1,ipoin,klim,kkll,irecper)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
      dimension ktemp(5,50),ktempl(5),ktemp2(5,50),ktempl2(5)
      dimension litt(20),litd(20),norma(50)
      dimension iprex(50),larp(2),iarp(2),ipfar(4),litz(20)
      dimension littr(20),normar(50),imark(720)
      dimension kron(50),khit(1000),izzy(100)
      dimension kr(10000,10),kktema(5),kktemb(5)
      lprx1=789000                                      
!      kmx1=70
!      kmx2=119
      
      ihit =0
      
      

      
      print *,'idrv',(idrv1(i),i=1,idrv1(2)+2)
      
      
!      larp(1)=int(lprx1/10000)
!      larp(2)=lprx1 -larp(1) *10000
      
      
      
      
      do i=1,720
      imark(i)=0
      end do
! start here      
      kklim=klim
      kkll2=(kklim+300)*7
      open (unit=8,file='nfsf1',access='direct',form=&
      'formatted',recl=kkll2,status='old')
      do i=1,10
      read (8,601,rec=i)(kr(jf,i),jf=1,klim+300)
      end do
601   format (11000i7)
      print *,'ok1'
      lmain=klim+1 
       print *,'krlmains',(kr(lmain,jt),jt=1,10)
      do lmain=klim+1,klim+300
      do ijj=1,8
      if (kr(lmain,ijj).eq.9999999)goto 700
      

      
      
      
      
      
!      do kia=1,1000
       kia=kr(lmain,ijj)
       print *,'krlmains',(kr(lmain,jt),jt=1,9)
       print *,'lmain',lmain,'kia',kia  ,'iprlmain',ipr(lmain)
       
      
      
      
      
      
      
      
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
      
      
      do i=1,50
      norma(i)=0
      normar(i)=0
      end do
      goto 610
      
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
      
      
2001  if (ia4(2).eq.0)goto 2002
      do i=1,ia4(2)
      karr(i)=ia4(i+2)
      end do
      
      ilen=ia4(2)
      do i=1,ktempl(4)
      kbarr(i)=ktemp(4,i)
      end do
      ilen2=ktempl(4)
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia4(1) +ksgn*4 ,2)
      
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do

      
      call mpadd(0)
      

      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      
      
2002  ilen=ia3(2)
      if (ia3(2).eq.0)goto 2003
      do i=1,ilen
      karr(i)=ia3(i+2)
      end do
      ilen2=ktempl(3)
      do i=1,ilen2
      kbarr(i)=ktemp(3,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      
      kbarr(1)=mod(ia3(1)+ksgn*3,2)
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      
      
2003  ilen=ia2(2)
      if (ia2(2).eq.0)goto 2004
      do i=1,ilen
      karr(i)=ia2(i+2)
      end do
      ilen2 =ktempl(2)
      do i=1,ilen2
      kbarr(i)=ktemp(2,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(ia2(1) +ksgn*2 ,2)
      






      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)  +2
      norma(i)=kcarr(i)
      end do
      
      
2004  ilen=ia1(2)
      if (ia1(2).eq.0)goto 2005
      do i=1,ilen
      karr(i)=ia1(i+2)
      end do
      ilen2=ktempl(1)
      do i=1,ilen2
      kbarr(i)=ktemp(1,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      
      kbarr(1)=mod(ia1(1) +ksgn ,2)
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      
      
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      
2005  do i=1,ia0(2)+2      
      kbarr(i)=ia0(i)
      end do
      
      kbarr(1)=mod(ia0(1) ,2)
      do i=1,norma(2) +2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      normar(i)=kcarr(i)
      end do
      
       
! middle starts here      
610  a=a






      
      do i=3,idrv5(2)+2
      karr(i-2)=idrv5(i)
      end do
      ilen =idrv5(2)
      do i=1,ktempl(4)
      kbarr(i)=ktemp(4,i)
      end do
      ilen2=ktempl(4)
      call mpmul(ilen,ilen2,ilen3)
      norma(2)=ilen3
      do i=1,ilen3
      norma(i+2)=kcarr(i)
      end do
      
      
      norma(1)=mod(idrv5(1) +ksgn*4,2)
      
      do i=3,idrv4(2)+2
      karr(i-2)=idrv4(i)
      end do
      ilen=idrv4(2)
      ilen2=ktempl(3)
      do i=1,ilen2
      kbarr(i)=ktemp(3,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(idrv4(1) +ksgn*3,2)
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      
      do i=3,idrv3(2)+2
      karr(i-2)=idrv3(i)
      end do
      ilen=idrv3(2)
      do i=1,ktempl(2)
      kbarr(i)=ktemp(2,i)
      end do
      ilen2=ktempl(2)
      call mpmul(ilen,ilen2,ilen3)
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      kbarr(1)=mod(idrv3(1)+ksgn*2,2)
      kbarr(2)=ilen3
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      
      do i=3,idrv2(2)+2
      karr(i-2)=idrv2(i)
      end do
      ilen=idrv2(2)
      ilen2=ktempl(1)
      do i=1,ktempl(1)
      kbarr(i)=ktemp(1,i)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(1)=mod(idrv2(1) +ksgn,2)
      kbarr(2)=ilen3
      do i=1,ilen3
      kbarr(i+2)=kcarr(i)
      end do
      do i=1,norma(2) +2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2)+2
      norma(i)=kcarr(i)
      end do
      
      do i=1,idrv1(2)+2
      kbarr(i)=idrv1(i)
      
      end do
      
      
      
      do i=1,norma(2)+2
      karr(i)=norma(i)
      end do
      call mpadd(0)
      do i=1,kcarr(2) +2
      norma(i)=kcarr(i)
      end do
      
      ilen2=2
      
!      do kpr=1,300
!      if (imark(kpr).eq.1)goto 3502
!      do i=3,normar(2)+2
!      karr(i-2)=normar(i)
!      end do
!      ilen=normar(2)
!      iqrx1=ipr(klim+kpr)
      iqrx1=ipr(lmain)
      if (iqrx1.lt.10000)goto 151
      iarq(1)=0
      iarq(2)=2
      iarq(3)=iqrx1/10000
      iarq(4)=iqrx1-iarq(3)*10000
!      iarq(1)=int(iqrx1/10000)
!      iarq(2)=iqrx1-iarq(1)*10000
      goto 152
151   iarq(1)=0
      iarq(2)=1
      iarq(3)=iqrx1
152  do i=3,iarq(2)+2
     kbarr(i-2)=iarq(i)
     end do

      ilen2=iarq(2)
!      kbarr(1)=iqrx1
!      kbarr(2)=iarq(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
!      if (irlen.ne.0)goto 3502
      do i=3,norma(2)+2
      karr(i-2)=norma(i)
      end do
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 7001
      
      do jj=1,irecper
      kkia=match(2*jj-1)
      kkib=match(2*jj)
      if (abs(kkia).lt.10000)goto 161
      kktema(1)=0
      kktema(2)=2
      kktema(3)=abs(kkia)/10000
      kktema(4)=abs(kkia)-kktema(3)*10000
      goto 162
161   kktema(1)=0      
      kktema(2)=1
      kktema(3)=abs(kkia)
162   if (abs(kkib).lt.10000)goto 163
      kktemb(1)=0
      kktemb(2)=2
      kktemb(3)=abs(kkib)/10000
      kktemb(4)=abs(kkib)-kktemb(3)*10000
      goto 164
163   kktemb(1)=0
      kktemb(2)=1
      kktemb(3)=abs(kkib)
164   if (kkia.lt.0)goto  165
      ksgna=0
      goto 166
165   ksgna=1
166   if (kkib.lt.0)goto 167
      ksgnb=0
      goto 168
167   ksgnb=1
168   do jf=1,ktempl(1)
      karr(jf)=ktemp(1,jf)
      end do
      ilen=ktempl(1)
      do jf=3,kktemb(2)+2
      kbarr(jf-2)=kktemb(jf)
      end do
      ilen2=kktemb(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      kbarr(2)=ilen3
      kbarr(1)=mod(ksgnb+ksgn,2)
      do jf=1,kktema(2)+2
      karr(jf)=kktema(jf)
      end do
      karr(1)=ksgna
      call mpadd(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,iarq(2)+2
      kbarr(jf-2)=iarq(jf)
      end do
      ilen2=iarq(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 7001



       end do




!      mtes=kkia-kkib*kia
!      if (mod(mtes,iqrx1).eq.0)goto 3502
      
      ihit=ihit+1
      print *,'hitn no.=',ihit,'kia=',kia,iqrx1
      print *,'normar',(normar(jf),jf=1,normar(2)+2)
      print *,'norma',(norma(jf),jf=1,norma(2)+2) 
      imark(ihit)=ipr(lmain)
      khit(ihit)=kia
      if (ihit.eq.51)goto 701
      goto 700
7001  end do
700   end do
701   a=a      
      print *,'khits',(khit(jf),jf=1,50)
      print *,'imarks',(imark(jf),jf=1,50)
      
      
      do jj=1,irecper

      iconty=0
      do ll=1,51
!      if (imark(ll).eq.0)goto 3510
      iconty=iconty+1
!      karp(3)=imark(ll)
       if (imark(ll).lt.10000)goto 3507
       karp(1)=0
       karp(2)=2
       karp(3)=imark(ll)/10000
       karp(4)=imark(ll)-karp(3)*10000
       goto 3508
3507   karp(1)=0
       karp(2)=1
       karp(3)=imark(ll)
3508   if (match(2*jj).lt.10000)goto 3509
       ilen=2
       karr(1)=match(2*jj)/10000
       karr(2)=match(2*jj)-karr(1)*10000
       goto 3510
3509   karr(1)=match(2*jj)        
       ilen=1
3510   if (khit(ll).lt.10000)goto 3511
       kbarr(1)=khit(ll)/10000
       kbarr(2)=khit(ll)-kbarr(1)*10000
       ilen2=2
       goto 3512
3511   kbarr(1)=khit(ll)       
       ilen2=1
3512   call mpmul(ilen,ilen2,ilen3)       
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=ilen3
       if (abs(match(2*jj-1)).lt.10000)goto 3513
       karr(2)=2
       karr(3)=abs(match(2*jj-1))/10000
       karr(4)=abs(match(2*jj-1))-karr(3)*10000
       goto 3514
3513   karr(2)=1
       karr(3)=abs(match(2*jj-1))
3514   if (match(2*jj-1).lt.0)goto 3515
       karr(1)=0
       goto 3516
3515   karr(1)=1
3516   call mpadd(1)
       isgn=kcarr(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,karp(2)+2
       kbarr(jf-2)=karp(jf)
       end do
       ilen2=karp(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 3522
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(2)=irlen
       kbarr(1)=isgn
       if (isgn.eq.0)goto 3520
       do jf=1,karp(2)+2
       karr(jf)=karp(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       kard(jf)=kcarr(jf)
       end do
       goto 3521
3520   do jf=1,kbarr(2)+2
       kard(jf)=kbarr(jf)
       end do
       goto 3521
       


3522   print *,'error slipped through'
       stop
3521   a=a






       
       
       
       


!      karp(3)=ipr(klim+ll)
!      karp(2)=1
!      karp(1)=0
!      ktia=khit(ll)
!      itemp=match(2*jj-1)-match(2*jj)*ktia
!      iatemp=abs(itemp)
!      if (iatemp.lt.10000)goto 35021
!      kard(2)=2
!      kard(3)=int(iatemp/10000)
!      kard(4)=iatemp-kard(3)*10000
!      goto 35022
!35021 kard(2)=1
!      kard(3)=iatemp

!35022 if (itemp.lt.0)goto 3504
!      kard(1)=0
!      goto 3505
!3504  kard(1)=1
!      do jf=3,karp(2)+2
!      kbarr(jf-2)=karp(jf)
!      end do
!      ilen2=karp(2)
!      do jf=3,kard(2)+2
!      karr(jf-2)=kard(jf)
!      end do
!      ilen=kard(2)
!      call mpdiv(ilen,ilen2,irlen,icont,iswq)
!      if (irlen.eq.0)goto 3507
!      do jf=1,irlen
!      kard(jf+2)=irlen(jf)
!      end do
!      kard(2)=irlen













3505  a=a
!3505  print *,'kard',kard(2),kard(3),karp(3),jj,iconty
!      print *,'matches',match(2*jj-1),match(2*jj),khit(ll)
      
      call mpkron(kx)
      
      kron(iconty)=mod(kx+3,3)-1
!      if (iconty.eq.51)goto 3506
      end do
      read(1,1,rec=jj)kkia,kkib,isgn,(nar(jf),jf=1,kkll)
      
!      do jl=1,51
!      nar(kkll-51+jl)=kron(jl)
!      end do
1     format(i6,i6,i1,11000i1)
      write(1,1,rec=jj)kkia,kkib,isgn,(nar(jf),jf=1,kkll-51),&
      (kron(jl),jl=1,51)
      print *,'kronjl',(kron(jl),jl=1,51),'jj',jj
      
3506  end do
      
      
      
      print *,'imarks',(imark(jf),jf=1,51)
      print *,'iprs 111 112',ipr(111),ipr(112)
      
      do jf =1,irecper
!      print *,match(2*jf-1),match(2*jf),(kron(jf,jk),jk=1,10)
      end do
      print *,'kaze',(khit(jf),jf=1,51)
      
      return
      end


   
      
      
      subroutine addstar(in,out)

      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      common mnum(50)
      common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
      
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
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      common ipr(65000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(5)
      common match(20000)
      common idrv5(10),idrv4(10),idrv3(10),idrv2(10),idrv1(10)
      common nar(10000),n(100)
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
      kx=k
      return
      end
