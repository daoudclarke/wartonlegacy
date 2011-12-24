      program bobecm6
!     elliptic curve method multiple parallel inverses phase 1 and 2
!     uses polynomials in continuation phase 20 degree
!     fast multiplication version (twice)     
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension ixp(30,50),iyp(30,50),ibb(30,50),idd(30,50),iprod(400)
      dimension iprod2(400)
      dimension iprecod(50),idy(400),idx(400)
      dimension ixdf(21,30,50),iydf(21,30,50)
      dimension kdoub2(21,300)
      dimension ixdf2(21,30,50),iydf2(21,30,50),iqarr(300),ie(300)
      dimension ixpow(21,30,50),iypow(21,30,50)
      iaa=31
      ilen=16
      ilen2=16
      do jf=1,ilen
      karr(jf)=jf+2
      kbarr(jf)=jf+3
      end do
      do i=1,10000
      call mpmul(ilen,ilen2,ilen3)
      end do
      print *,'ilen3',ilen3,(kcarr(jf),jf=1,ilen3)
      stop
      ilen=298
      ilen2=298
      
      do i=1,99
      karr(i)=1234
      kbarr(i)=5678
      end do
      do i=100,298
      karr(i)=1111
      kbarr(i)=2222
      end do
      do i=1,10000
      call mpmul(ilen,ilen2,ilen3)
      end do
      
!      print *,'ilen3',ilen3
!      print *,(kcarr(jf),jf=1,ilen3)
      print *,'ilen3',ilen3
      stop
      karr(1)=1
      karr(2)=2
      karr(3)=3
      karr(4)=0
      karr(5)=0
      karr(6)=6
      karr(7)=7
      karr(8)=8
      karr(9)=9
      karr(10)=10
      ilen=10
      kbarr(1)=10
      kbarr(2)=0
      ilen2=2
      call mpmul(ilen,ilen2,ilen3)
      print *,'ilen3',ilen3
      print *,'prod',(kcarr(jf),jf=1,ilen3)
      stop
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      



      
      irecnn=0
      iprdl=799999
      kkmx=0
      open(unit=1,file='polynom2',access='direct',form=&
      'formatted',recl=1200,status='old')
      
      
      
      goto 4993
      
      
      
      inlen=2
      mnum(1)=80
      mnum(2)=1
      do k=1,15
      icon=1
499   call mpprime(icorp,inlen)
      if (icorp.eq.0)goto 500
      ipr(icon)=mnum(1)*10000+mnum(2)
      if (icon.eq.65000)goto 502
      icon=icon+1
500   karr(1)=0
      karr(2)=2
      karr(3)=mnum(1)
      karr(4)=mnum(2)
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=2
      call mpadd(0)
      mnum(1)=kcarr(3)
      mnum(2)=kcarr(4)
      inlen=2
      goto 499
502   write(2,4992,rec=k)(ipr(jf),jf=1,65000)
      end do
4992  format(1500i4,1500i4)
      close(unit=2)
      stop
4993  open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
      read(3,5,rec=1)(ipr(jf),jf=1,65000)
      
      print *,'code length radix 10000?'
      read *,n(2)
      print *,'enter code'
      read *,(n(jf),jf=3,n(2)+2)
      
      print *,'expected max. length of min. factor in powers of 10?'
      read *,nr
      print *,'phase 2 multiple? of 50'
      read *,iph2
      n(1)=0
      print *,'n',(n(jf),jf=1,n(2)+2)
      econ=0.4342945
      rel=nr/econ
      pow=rel*log(rel)

      
      print *,'pow',pow
      pow=pow**0.5
      print *,'pow',pow
      rel2=exp(pow)
      pow2=1/(2**0.5)
      l1=rel2**pow2
      print *,'l1=',l1
      
      l2=l1*iph2
      lten=(l1/10)*10
      n(1)=0
      do jf=1,n(2)+2
      iprecod(jf)=n(jf)
      end do
      do i=1,200
      
      
      
      
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
      print *,'iprecod',(iprecod(jf),jf=1,iprecod(2)+2)
      

64    do i=1,30
      ix1(i,1)=0
      ix1(i,2)=0
      ix1(i,3)=0
      ix2(i,1)=0
      ix2(i,2)=0
      ix2(i,3)=0
      iy1(i,1)=0
      iy1(i,2)=1
      iy1(i,3)=1
      iy2(i,1)=0
      iy2(i,2)=1
      iy2(i,3)=1
      end do
      ibbsw=0
      ibbsw2=0
      iprc=1 
114   iprc=iprc+1
      iprx1=ipr(iprc)
      
      if (iprx1.gt.l1)goto 800 
      
      iq=iprx1
      iq1=iq
      
      l=l1/iq
       
128   if (iq1.gt.l)goto 131
      iq1=iq1*iq
      goto 128
131   print *,'iaaiq1',iaa,iq1
      iee=1
138   if (iee.gt.iq1)goto 146
      iee=iee*2
      goto 138
146   iee=iee/2
      iq1=iq1-iee
147   if (iee.eq.1)goto 300
      iee=iee/2
1471  do jf=1,iy2(1,2)+2
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
      do i=2,30
      do jf=1,iy2(i,2)+2
      karr(jf)=iy2(i,jf)
      kbarr(jf)=karr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      idd(i,jf)=kcarr(jf)
      end do
      do jf=3,idd(i,2)+2
      karr(jf-2)=idd(i,jf)
      end do
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen=idd(i,2)
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      indz=2
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      idd(i,jf+2)=irrr(jf)
      end do
      idd(i,2)=irlen
      if (idd(i,1).eq.0)goto 150
      do jf=1,idd(i,2)+2
      karr(jf)=idd(i,jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      idd(i,jf)=kcarr(jf)
      end do
150   do jf=3,idd(i,2)+2
      karr(jf-2)=idd(i,jf)
      end do
      ilen=idd(i,2)
      
      do jf=3,iprod(2)+2
      kbarr(jf-2)=iprod(jf)
      end do
      ilen2=iprod(2)
      
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
      indz=3
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      iprod(jf+2)=irrr(jf)
      end do
      iprod(2)=irlen
      
      do jf=1,iprod(2)+2
      ibb(i,jf)=iprod(jf)
      end do
      
      
      end do
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      do jf=1,ibb(30,2)+2
      karb(jf)=ibb(30,jf)
      end do
!     print *,'hello1'
      call mpgcd
      
      if (kard(2).gt.1) goto 1500
      if (kard(3).gt.1)goto 1500
      do jf=3,karv(2)+2
      karr(jf-2)=karv(jf)
      iprod2(jf)=karv(jf)
      end do
      ilen=karv(2)
      iprod2(2)=ilen
      iprod2(1)=0
      do jf=3,ibb(29,2)+2
      kbarr(jf-2)=ibb(29,jf)
      end do
      ilen2=ibb(29,2)
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
      indz=4
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      icc(30,jf+2)=irrr(jf)
      end do
      icc(30,2)=irlen
      icc(30,1)=0
      do i=1,28
      do jf=3,iprod2(2)+2
      karr(jf-2)=iprod2(jf)
      end do
      ilen=iprod2(2)
      do jf=3,idd(31-i,2)+2
      kbarr(jf-2)=idd(31-i,jf)
      end do
      ilen2=idd(31-i,2)
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
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      iprod2(jf+2)=irrr(jf)
      end do
      iprod2(1)=0
      iprod2(2)=irlen
      do jf=3,iprod2(2)+2
      karr(jf-2)=iprod2(jf)
      end do
      ilen=iprod2(2)






      do jf=3,ibb(29-i,2)+2
      kbarr(jf-2)=ibb(29-i,jf)
      end do
      ilen2=ibb(29-i,2)
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
      indz=5
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      icc(30-i,jf+2)=irrr(jf)
      end do
      icc(30-i,2)=irlen
      icc(30-i,1)=0
      
      end do
      do jf=3,iprod2(2)+2
      karr(jf-2)=iprod2(jf)
      end do
      ilen=iprod2(2)
      do jf=3,idd(2,2)+2
      kbarr(jf-2)=idd(2,jf)
      end do
      ilen2=idd(2,2)
      
      
      
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
      indz=6
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      icc(1,jf+2)=irrr(jf)
      end do
      icc(1,1)=0
      icc(1,2)=irlen
      
      

      call subeq(iaa)
      if (ibbsw.gt.0)goto 777
      
      if (iq1.lt.iee)goto 147
      iq1=iq1-iee
1601  do jf=1,ix1(1,2)+2
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
      do i=2,30
      do jf=1,ix1(i,2)+2
      karr(jf)=ix1(i,jf)
      end do
      do jf=1,ix2(i,2)+2
      kbarr(jf)=ix2(i,jf)
      end do
      call mpadd(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      idd(i,1)=kcarr(1)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      indz=8
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      idd(i,jf+2)=irrr(jf)
      karr(jf)=irrr(jf)
      end do
      ilen=irlen
      idd(i,2)=irlen
      if (idd(i,1).eq.0)goto 191
      do jf=1,idd(i,2)+2
      karr(jf)=idd(i,jf)
      end do
      do jf=1,n(2)+2
      kbarr(jf)=n(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      idd(i,jf)=kcarr(jf)
      
      end do
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)


191   a=a      
      
      
      do jf=3,iprod(2)+2
      kbarr(jf-2)=iprod(jf)
      end do
      ilen2=iprod(2)
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
      indz=9
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      iprod(jf+2)=irrr(jf)
      end do
      iprod(2)=irlen
      do jf=1,iprod(2)+2
      ibb(i,jf)=iprod(jf)
      end do
      
      end do
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      do jf=1,ibb(30,2)+2
      karb(jf)=ibb(30,jf)
      end do
      call mpgcd
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      do jf=3,karv(2)+2
      karr(jf-2)=karv(jf)
      iprod2(jf)=karv(jf)
      end do
      ilen=karv(2)
      iprod2(2)=ilen
      iprod2(1)=0
      do jf=3,ibb(29,2)+2
      kbarr(jf-2)=ibb(29,jf)
      end do
      ilen2=ibb(29,2)
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
      indz=10
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      icc(30,jf+2)=irrr(jf)
      end do
      icc(30,2)=irlen
      icc(30,1)=0
      
      do i=1,28
      do jf=3,iprod2(2)+2
      karr(jf-2)=iprod2(jf)
      end do
      ilen=iprod2(2)
      do jf=3,idd(31-i,2)+2
      kbarr(jf-2)=idd(31-i,jf)
      end do
      ilen2=idd(31-i,2)
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
      indz=11
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      iprod2(jf+2)=irrr(jf)
      karr(jf)=irrr(jf)
      end do
      iprod2(2)=irlen
      iprod2(1)=0
      ilen=irlen
      do jf=3,ibb(29-i,2)+2
      kbarr(jf-2)=ibb(29-i,jf)
      end do
      ilen2=ibb(29-i,2)
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
      indz=12
      if (irlen.eq.0)goto 360

      do jf=1,irlen















      icc(30-i,jf+2)=irrr(jf)
      end do
      icc(30-i,2)=irlen
      icc(30-i,1)=0
      
      
      
      end do
      do jf=3,iprod2(2)+2
      karr(jf-2)=iprod2(jf)
      end do
      ilen=iprod2(2)
      do jf=3,idd(2,2)+2
      kbarr(jf-2)=idd(2,jf)
      end do
      ilen2=idd(2,2)
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
      indz=13
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      icc(1,jf+2)=irrr(jf)
      end do
      icc(1,1)=0
      icc(1,2)=irlen
      itag=0
!     print *,'hello'
             
      call subneq(iaa,itag)
      if (ibbsw.gt.0)goto 777
      if (itag.eq.1)goto 360
      
      
      
      
      
      goto 147
300   if (ibbsw2.eq.1)goto 802
      
      do i=1,30
      do jf=1,ix2(i,2)+2
      ix1(i,jf)=ix2(i,jf)
      end do
      do jf=1,iy2(i,2)+2
      iy1(i,jf)=iy2(i,jf)
      
      
      end do




      end do
      
      goto 114
360   iaa=iaa+i+1
      itag=0
      print *,'i',i,'iprx1',iprx1
      print *,'ix2',(ix2(1,jf),jf=1,ix2(1,2)+2)
      print *,'iy2',(iy2(1,jf),jf=1,iy2(1,2)+2)
      print *,'indz',indz
      
      goto 64
600   iaa=iaa+30
      goto 64
700   print *,'all 2s and 3s'
      stop
777   if (ibbsw.eq.1)goto 2040
      if (ibbsw.eq.2)goto 2044
      if (ibbsw.eq.3)goto 822
      if (ibbsw.eq.4)goto 902
      











800   if (iph2.eq.1)goto 600
      icph=1
      do i=1,30
      do jf=1,ix2(i,2)+2
      ixp(i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iyp(i,jf)=iy2(i,jf)
      end do
      end do
!     setting up ibr1       
      ibr1=1
801   iq1=lten
      ibbsw=0
      ibbsw2=1
      goto 131

      
802   do i=1,30
      do jf=1,ix2(i,2)+2
      ix1(i,jf)=ix2(i,jf)
      ixpow(ibr1,i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iypow(ibr1,i,jf)=iy2(i,jf)
      iy1(i,jf)=iy2(i,jf)
      end do
      end do
      if (ibr1.eq.20)goto 803
      ibr1=ibr1+1
      goto 801
      
!     setting up ibig      
803   ibig=1
809   irs=(ibig*(45-ibig))/2 -21
      print *,'ibig',ibig
      ibr2=22-ibig
      irv=irs
      do i3=1,ibr2
      read (1,1111,rec=irv)(kdoub2(i3,jf),jf=1,300)
1111  format(300i4)      
      
      irv=irv+1
      end do
!     setting up ibig2      
      
      ibig2=1
807   do jf=1,300
      iqarr(jf)=kdoub2(ibr2+1-ibig2,jf)
      end do
      if (ibig2.ne.1)goto 810
      do i=1,30
      do jf=1,ixp(i,2)+2
      ix1(i,jf)=ixp(i,jf)
      ix2(i,jf)=ixp(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iyp(i,2)+2
      iy1(i,jf)=iyp(i,jf)
      iy2(i,jf)=iyp(i,jf)
      end do
      end do
      
      
      ibbsw3=1
      
      goto 2000
806   do i=1,30
      do jf=1,ix2(i,2)+2
      ixdf(ibig,i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iydf(ibig,i,jf)=iy2(i,jf)
      end do
      end do




824   if (ibig2.eq.ibr2)goto 808   
      ibig2=ibig2+1
      print *,'ibig2',ibig2,'ibig',ibig
      goto 807
808   if (ibig.eq.21)goto 900      
      ibig=ibig+1
      goto 809
810   ivvr=ibig2-1
      do i=1,30
      do jf=1,ixpow(ivvr,i,2)+2
      ix2(i,jf)=ixpow(ivvr,i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iypow(ivvr,i,2)+2
      iy2(i,jf)=iypow(ivvr,i,jf)
      end do
      end do
      
      if (iqarr(2).ne.1)goto 8104
      if (iqarr(3).ne.1)goto 8104
      do i=1,30
      do jf=1,ixdf(ibig,i,2)+2
      ix1(i,jf)=ixdf(ibig,i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iydf(ibig,i,2)+2
      iy1(i,jf)=iydf(ibig,i,jf)
      end do
      end do
      ibbsw=3
      goto 1601
8104  ibbsw3=2
      do i=1,30
      do jf=1,ix2(i,2)+2
      ix1(i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iy1(i,jf)=iy2(i,jf)
      end do
      end do
      
      goto 2000
820   do i=1,30
      do jf=1,ixdf(ibig,i,2)+2
      ix1(i,jf)=ixdf(ibig,i,jf)
      end do
      end do
      ibbsw=3
      do i=1,30
      do jf=1,iydf(ibig,i,2)+2
      iy1(i,jf)=iydf(ibig,i,jf)
      end do
      end do
      
      goto 1601
822   do i=1,30
      do jf=1,ix2(i,2)+2
      ixdf(ibig,i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iydf(ibig,i,jf)=iy2(i,jf)
      end do
      end do
      
      goto 824
!     setting up igg
900   igg=1
      ibbsw=4
901   do i=1,30
      do jf=1,ixdf(igg,i,2) +2
      ix1(i,jf)=ixdf(igg,i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iydf(igg,i,2)+2
      iy1(i,jf)=iydf(igg,i,jf)
      end do
      end do
      do i=1,30
      do jf=1,ixdf(igg+1,i,2)+2
      ix2(i,jf)=ixdf(igg+1,i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iydf(igg+1,i,2)+2
      iy2(i,jf)=iydf(igg+1,i,jf)
      end do
      end do
      goto 1601
902   do i=1,30
      do jf=1,ix2(i,2)+2
      ixdf2(igg,i,jf)=ix2(i,jf)
      end do
      end do
      do i=1,30
      do jf=1,iy2(i,2)+2
      iydf2(igg,i,jf)=iy2(i,jf)
      end do
      end do
      if (igg.eq.20)goto 910
      igg=igg+1
      goto 901
910   if (icph.eq.iph2)goto 600
      print *,'icph',icph
      icph=icph+1
      do igg=1,20
      do i=1,30
      do jf=1,ixdf2(igg,i,2)+2
      ixdf(igg,i,jf)=ixdf2(igg,i,jf)
      end do
      end do
      end do
      do igg=1,20
      do i=1,30
      do jf=1,iydf2(igg,i,2)+2
      iydf(igg,i,jf)=iydf2(igg,i,jf)
      end do
      end do
      end do
      goto 900





2000  karr(1)=2
      kbarr(1)=2
      ilen=1
      ilen2=1
      inlen=iqarr(2)
2004  call mpmul(ilen,ilen2,ilen3)
      if (ilen3.lt.inlen)goto 2022
      if (ilen3.gt.inlen)goto 2020
      do i=1,ilen3
      if (iqarr(i+2).lt.kcarr(i))goto 2020
      if (iqarr(i+2).gt.kcarr(i))goto 2022
      end do
      goto 2022
2020  ie(2)=ilen3
      do i=1,ilen3
      ie(i+2)=kcarr(i)
      end do
      goto 2030
2022  do i=1,ilen3
      kbarr(i)=kcarr(i)
      end do
      ilen2=ilen3
      goto 2004
2030  do i=1,ie(2)
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
      ie(1)=0
      do i=1,iqarr(2)+2
      karr(i)=iqarr(i)
      end do
      do i=1,ie(2)+2
      kbarr(i)=ie(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2)+2
      iqarr(i)=kcarr(i)
      end do
2031  if (ie(2).gt.1)goto 2032
      if (ie(3).eq.1)goto 2100
2032  do i=1,ie(2)
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
      ibbsw=1
      goto 1471
2040  if (iqarr(2).lt.ie(2))goto 2031
      if (iqarr(2).gt.ie(2))goto 2042
      do i=3,iqarr(2)+2
      if (iqarr(i).lt.ie(i))goto 2031
      if (iqarr(i).gt.ie(i))goto 2042
      end do
2042  do i=1,iqarr(2)+2
      karr(i)=iqarr(i)
      end do
      do i=1,ie(2)+2
      kbarr(i)=ie(i)
      end do
      call mpadd(1)
      do i=1,kcarr(2)+2
      iqarr(i)=kcarr(i)
      end do
      ibbsw=2
      goto 1601
2044  goto 2031

2100  if (ibbsw3.eq.1)goto 806
      if (ibbsw3.eq.2)goto 820








!     exiting to different points
1500  print *,'endwell prime=',ipr(iprc)
      print *,'prime divisor found on curve',iaa
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







      
      
5     format(65000i6)      
      
      end
      subroutine subeq(iaa)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension m1(400),ixt(400),iyt(400),ixtemp(400),ms1(2000)
      
      
      
      
      
      
      
      
      
      
      do i=1,30
      m1(3)=iaa+i-1
      m1(1)=0
      m1(2)=1

      if (ix2(i,2).eq.0)goto 22
      do jf=3,ix2(i,2)+2
      karr(jf-2)=ix2(i,jf)
      kbarr(jf-2)=ix2(i,jf)
      end do
      ilen=ix2(i,2)
      ilen2=ix2(i,2)
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      kbarr(1)=3
      ilen2=1
      call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      end do
      
      karr(2)=ilen3
      karr(1)=0
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=m1(3)

      
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
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
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
      if (iyt(2).lt.200)goto 11111
      print *,'iytlen',iyt(2)
      stop
      
11111 do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      end do
      ilen=kcarr(2)
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
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
      stop
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
      
      
      
      
      
      
      
      stop
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
      return
      end
      subroutine subneq(iaa,itag)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension m1(400),ixt(400),iyt(400),ixtemp(400)
      
      
      inz=0
      do i=1,30
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
      stop
40    a=a      
      
      end do
      return
      end



























                    
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension karrp(800),kbarrp(800),karra(2,800),kbarra(2,800)
      dimension ilena(2),ilena2(2),ibdf(800)
      dimension iarr(800),iadf(800)
      ilenp=ilen
      ilenp2=ilen2
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      if (ilenp.gt.ilenp2)goto 1
      ilenas=ilenp2/2
      if (mod(ilenp2,2).eq.0)goto 2
      ilenas=ilenas+1
      goto 2
1     ilenas=ilenp/2
      if (mod(ilenp,2).eq.0)goto 2
      ilenas=ilenas+1

2    a=a      
     if (ilenp.lt.ilenas+1)goto 210
      ilena(1)=ilenp-ilenas
      
      do jf=1,ilena(1)
      karra(1,jf)=karr(jf)
      end do
      
      
      
211   ired=ilena(1)+1 
      ijk=1
      do jf=ired,ired+ilenas-1
      if (karr(jf).eq.0)goto 204
      ilena(2)=ilenas+1-ijk
      do jk=1,ilena(2)
      karra(2,jk)=karr(jf-1+jk)
      end do
      goto 212
204   ijk=ijk+1
      end do
      ilena(2)=0

      goto 212
210   ilena(1)=0
      if (ilenp.eq.ilenas)goto 211
      ilena(2)=ilenp
      do jf=1,ilenp
      karra(2,jf)=karr(jf)
      end do
          
212   a=a
!      print *,ilena(1),ilena(2)
      do i=1,2
      if (ilena(i).eq.0)goto 213
!      print *,'i',i,'karra',(karra(i,jf),jf=1,ilena(i))
213   end do
      
      
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
      goto 312
304   ijk=ijk+1 
      end do
      ilena2(2)=0
      goto 312
310   ilena2(1)=0
      if (ilenp2.eq.ilenas)goto 311
      ilena2(2)=ilenp2
      do jf=1,ilenp2
      kbarra(2,jf)=kbarr(jf)
      end do



      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      end do
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      end do
      ilen2=ilena2(1)
      isub=2
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'ilen3',ilen3
      goto 540
535   iarr(1)=0
      iarr(2)=0
      goto 545
540   do jf=1,ilen3
      iarr(jf+2)=kcarr(jf)
      end do
      

!      print *,'point  1'
      iarr(1)=0
      iarr(2)=ilen3
      karr(1)=iarr(1)
      karr(2)=iarr(2)+ilenas
      do jf=1,ilen3
      karr(jf+2)=iarr(jf+2)
      end do
      do jf=ilen3+3,ilen3+2 +ilenas 
      karr(jf)=0
      end do
      kbarr(1)=0 
      kbarr(2)=karr(2)+ilenas
      do jf=3,karr(2)+2
      kbarr(jf)=karr(jf)
      end do
      do jf=karr(2)+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'iarr1',iarr(2)
      

545   if (ilena(1).eq.0)goto 546
      do jf=1,ilena(1)
      karr(jf+2)=karra(1,jf)
      end do
      karr(1)=0
      karr(2)=ilena(1)
      goto 547
546   karr(1)=0
      karr(2)=0
      
547   if (ilena(2).eq.0)goto 548
      do jf=1,ilena(2)
      kbarr(jf+2)=karra(2,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena(2)
      goto 549
548   kbarr(1)=0
      kbarr(2)=0
549   call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
      if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5495  karr(1)=0       
      karr(2)=0
5497  if (ilena2(1).eq.0)goto 5498
      
      do jf=1,ilena2(1)
      kbarr(jf+2)=kbarra(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena2(1)
      goto 5499
5498  kbarr(1)=0      
      kbarr(2)=0
5499  call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do 
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      end do
      ilen=iadf(2)
      isub=3
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'2ilen3',ilen3
      
550   do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      
      
      kbarr(1)=mod(iadf(1)+isgn,2)
!      print *,'kbarr1',kbarr(1)
      kbarr(2)=ilen3+ilenas
      
      
      do jf=ilen3+3,ilen3+2+ilenas
      kbarr(jf)=0
      end do
!      print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
      
!      print *,'iarr2',iarr(2)
      
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'second iarr',iarr(2)
      
552   if (ilena(2).eq.0)goto 553      
      do jf=1,ilena(2)
      karr(jf)=karra(2,jf)
      end do
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      end do
      ilen2=ilena2(2)
      isub=4 
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'3ilen3',ilen3
      
560   do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      kbarr(jf+2)=kcarr(jf)
      end do
      
      karr(1)=0
      kbarr(1)=0
      karr(2)=ilen3
      
      kbarr(2)=ilen3+ilenas
      do jf=ilen3+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kbarr(jf)=kcarr(jf)
      end do
!      print *,'last',(kcarr(jf),jf=1,kcarr(2)+2)
!      print *,'iarrlen',iarr(2)
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      ilen3=kcarr(2)
!      print *,'fir ilen3',ilen3
      
      do jf=3,ilen3 +2
      kcarr(jf-2)=kcarr(jf)
      
      end do
!      print *,'kcarrml',(kcarr(jf),jf=1,ilen3)
      
      goto 110
553   ilen3=iarr(2)
      do jf=1,ilen3
      kcarr(jf)=iarr(jf+2)
      end do
      goto 110
      stop
      
      jdf=ilenp-(i-1)*ilenas
      

!      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     a=a
10    do i=1,ilen+ilen2
      kcarr(i) =0
      end do
!      print *,'two',ilen,ilen2
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
!      print *,'four ilen3',ilen3
      if(kcarr(1).ne.0)goto 100
      ilen3 =ilen3-1
!      print *,'three ilen3',ilen3
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

100   if (isub.eq.2)goto 540
      if (isub.eq.3)goto 550
      if (isub.eq.4)goto 560
110   ilen=ilenp
      ilen2=ilenp2
      do jf=1,ilen
      karr(jf)=karrp(jf)
      end do
      do jf=1,ilen2
      kbarr(jf)=kbarrp(jf)
      end do
      return


      end

      subroutine mpdiv(ilen,ilen2,irlen,icont,iswq)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension kdum(800),isub(800)
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
      subroutine mpmul2(ilen,ilen2,ilen3)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension karrp(800),kbarrp(800),karra(2,800),kbarra(2,800)
      dimension ilena(2),ilena2(2),ibdf(800)
      dimension iarr(800),iadf(800)
      ilenp=ilen
      ilenp2=ilen2
      
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      if (ilenp.gt.ilenp2)goto 1
      ilenas=ilenp2/2
      if (mod(ilenp2,2).eq.0)goto 2
      ilenas=ilenas+1
      goto 2
1     ilenas=ilenp/2
      if (mod(ilenp,2).eq.0)goto 2
      ilenas=ilenas+1

2    a=a      
     if (ilenp.lt.ilenas+1)goto 210
      ilena(1)=ilenp-ilenas
      
      do jf=1,ilena(1)
      karra(1,jf)=karr(jf)
      end do
      
      
      
211   ired=ilena(1)+1 
      ijk=1
      do jf=ired,ired+ilenas-1
      if (karr(jf).eq.0)goto 204
      ilena(2)=ilenas+1-ijk
      do jk=1,ilena(2)
      karra(2,jk)=karr(jf-1+jk)
      end do
      goto 212
204   ijk=ijk+1
      end do
      ilena(2)=0

      goto 212
210   ilena(1)=0
      if (ilenp.eq.ilenas)goto 211
      ilena(2)=ilenp
      do jf=1,ilenp
      karra(2,jf)=karr(jf)
      end do

      
      
      
212   a=a
!      print *,ilena(1),ilena(2)
      do i=1,2
      if (ilena(i).eq.0)goto 213
!      print *,'i',i,'karra',(karra(i,jf),jf=1,ilena(i))
213   end do
      
      
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
      goto 312
304   ijk=ijk+1 
      end do
      ilena2(2)=0
      goto 312
310   ilena2(1)=0
      if (ilenp2.eq.ilenas)goto 311
      ilena2(2)=ilenp2
      do jf=1,ilenp2
      kbarra(2,jf)=kbarr(jf)
      end do
      
      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      end do
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      end do
      ilen2=ilena2(1)
      isub=2
      goto 10
535   iarr(1)=0
      iarr(2)=0
      goto 545
540   do jf=1,ilen3
      iarr(jf+2)=kcarr(jf)
      end do
      

!      print *,'point  1'
      iarr(1)=0
      iarr(2)=ilen3
      karr(1)=iarr(1)
      karr(2)=iarr(2)+ilenas
      do jf=1,ilen3
      karr(jf+2)=iarr(jf+2)
      end do
      do jf=ilen3+3,ilen3+2 +ilenas 
      karr(jf)=0
      end do
      kbarr(1)=0 
      kbarr(2)=karr(2)+ilenas
      do jf=3,karr(2)+2
      kbarr(jf)=karr(jf)
      end do
      do jf=karr(2)+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'iarr1',iarr(2)
      

545   if (ilena(1).eq.0)goto 546
      do jf=1,ilena(1)
      karr(jf+2)=karra(1,jf)
      end do
      karr(1)=0
      karr(2)=ilena(1)
      goto 547
546   karr(1)=0
      karr(2)=0
      
547   if (ilena(2).eq.0)goto 548
      do jf=1,ilena(2)
      kbarr(jf+2)=karra(2,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena(2)
      goto 549
548   kbarr(1)=0
      kbarr(2)=0
549   call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
      if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5495  karr(1)=0       
      karr(2)=0
5497  if (ilena2(1).eq.0)goto 5498
      
      do jf=1,ilena2(1)
      kbarr(jf+2)=kbarra(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena2(1)
      goto 5499
5498  kbarr(1)=0      
      kbarr(2)=0
5499  call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do 
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      end do
      ilen=iadf(2)
      isub=3
      goto 10
550   do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      
      
      kbarr(1)=mod(iadf(1)+isgn,2)
!      print *,'kbarr1',kbarr(1)
      kbarr(2)=ilen3+ilenas
      
      
      do jf=ilen3+3,ilen3+2+ilenas
      kbarr(jf)=0
      end do
!      print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
      
!      print *,'iarr2',iarr(2)
      
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'second iarr',iarr(2)
      
552   if (ilena(2).eq.0)goto 553      
      do jf=1,ilena(2)
      karr(jf)=karra(2,jf)
      end do
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      end do
      ilen2=ilena2(2)
      isub=4 
      goto 10
560   do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      kbarr(jf+2)=kcarr(jf)
      end do
      
      karr(1)=0
      kbarr(1)=0
      karr(2)=ilen3
      
      kbarr(2)=ilen3+ilenas
      do jf=ilen3+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kbarr(jf)=kcarr(jf)
      end do
!      print *,'last',(kcarr(jf),jf=1,kcarr(2)+2)
!      print *,'iarrlen',iarr(2)
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      ilen3=kcarr(2)
!      print *,'fir ilen3',ilen3
      
      do jf=3,ilen3 +2
      kcarr(jf-2)=kcarr(jf)
      
      end do
!      print *,'kcarrml',(kcarr(jf),jf=1,ilen3)
      
      goto 110
553   ilen3=iarr(2)
      do jf=1,ilen3
      kcarr(jf)=iarr(jf+2)
      end do
      goto 110
      stop





      
      
      jdf=ilenp-(i-1)*ilenas




!      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     a=a
10    do i=1,ilen+ilen2
      kcarr(i) =0
      end do
!      print *,'two',ilen,ilen2
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
!      print *,'four ilen3',ilen3
      if(kcarr(1).ne.0)goto 100
      ilen3 =ilen3-1
!      print *,'three ilen3',ilen3
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

100   if (isub.eq.2)goto 540
      if (isub.eq.3)goto 550
      if (isub.eq.4)goto 560
110   ilen=ilenp
      ilen2=ilenp2
      do jf=1,ilen
      karr(jf)=karrp(jf)
      end do
      do jf=1,ilen2
      kbarr(jf)=kbarrp(jf)
      end do
      return


      end

      
       



      subroutine mpadd(isora)
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
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
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
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
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
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
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension karu(90),karv1(90),karv3(90),karqq(90)
      dimension kart3(90),kart1(90)
      
      
      
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
      common karr(800),kbarr(800),kcarr(1600),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      common ix1(30,50),iy1(30,50),ix2(30,50),iy2(30,50),icc(30,50)
      common n(50)
      dimension karae(90),karde(90),karpe(90),karh(90),ktemp(90)
      
      
      
      
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
