      

      program bobpq1
!     phases 1 and 2 of p-1 method
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
      common n(50)
      dimension ixp(50),iyp(50),ix1(50),ix2(50)
      dimension iprod2(400)
      dimension iprecod(50),idy(400),idx(400)
      dimension ixl(50),iyl(50),ixdf(500,50),iydf(1000,50),igcd(50)
      iaa=1
      
      karr(1)=23
      karr(2)=2274
      ilen=2
      kbarr(1)=10
      kbarr(2)=0
      ilen2=2
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      



      
      irecnn=0
      jrecnn=0
      iprdl=799999
      kkmx=0
      open(unit=2,file='ecmprim',access='direct',form=&
      'formatted',recl=520000,status='old')
      open(unit=4,file='ecmprim2',access='direct',form=&
      'formatted',recl=585000,status='old')
5992  format (65000i9)
!     read (2,4992,rec=14)(ipr(jf),jf=1,65000)
!     print *,ipr(65000)
!     read (2,4992,rec=15)(ipr(jf),jf=1,65000)
!     print *,ipr(1)
      

      goto 4993
!     read (2,4992,rec=15)(ipr(jf),jf=1,65000)
!     print *,(ipr(jf),jf=1,65000)
      
      
      inlen=2
      mnum(1)=80
      mnum(2)=1
      do k=1,15
      icon=1
499   call mpprime(icorp,inlen)
      if (icorp.eq.0)goto 500
      ipr(icon)=mnum(1)*10000+mnum(2)
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
      
      if (icon.eq.65001)goto 502
      
      goto 499
502   write(2,4992,rec=k)(ipr(jf),jf=1,65000)
      end do
4992  format(65000i8)
      close(unit=2)
      stop
4993  open(unit=3,file ='recl.dat',access='direct',form =&
      'formatted',recl=390000,status='old')
      read(3,5,rec=1)(ipr(jf),jf=1,65000)
      
      print *,'code length radix 10000?'
      read *,n(2)
      print *,'enter code'
      read *,(n(jf),jf=3,n(2)+2)
      
      print *,'first stage limit'
      read *,l1
      print *,'phase 2 multiple?'
      read *,iph2
      l2=l1*iph2
      n(1)=0
      print *,'n',(n(jf),jf=1,n(2)+2)
      icnt=0
      iprc=1
      ibbsw2=0
      ix1(1)=0
      ix1(2)=1
      ix1(3)=2
      ix2(1)=0
      ix2(2)=1
      ix2(3)=2
114   iprc=iprc+1 
      if (iprc.eq.65001)goto 330
      
      if (ipr(iprc).gt.l1)goto 800
      
      if (ipr(iprc).gt.l1)goto 800
!      if (ipr(iprc).eq.iprdl)goto 330
      if ((ipr(iprc).eq.0).and.(irecnn.eq.0))goto 330
      
      iprx1=ipr(iprc)
1141  do jf=1,ix2(2)+2
      ixp(jf)=ix2(jf)
      end do
      
      
      
      
      iq=iprx1
      iq1=iq
      
      l=l1/iq
       
128   if (iq1.gt.l)goto 131
      
      iq1=iq1*iq
      goto 128
131   if (icnt.gt.0)goto 132
      print *,'iq1',iq1
132   iee=1
138   if (iee.gt.iq1)goto 146
      iee=iee*2
      goto 138
146   iee=iee/2
      iq1=iq1-iee
147   if (iee.eq.1)goto 300
      iee=iee/2
1471  do jf=3,ix2(2)+2
      karr(jf-2)=ix2(jf)
      kbarr(jf-2)=ix2(jf)
      end do
      ilen=ix2(2)
      ilen2=ilen
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
      if (irlen.eq.0)goto 370
      do jf=1,irlen
      ix2(jf+2)=irrr(jf)
      end do
      ix2(1)=0
      ix2(2)=irlen
      
      if (iq1.lt.iee)goto 147
      iq1=iq1-iee
1601  do jf=3,ix1(2)+2
      karr(jf-2)=ix1(jf)
      end do
      do jf=3,ix2(2)+2
      kbarr(jf-2)=ix2(jf)
      end do
      ilen=ix1(2)
      ilen2=ix2(2)
      call mpmul(ilen,ilen2,ilen3)
      
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen=ilen3
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 370
      do jf=1,irlen
      ix2(jf+2)=irrr(jf)
      end do  
      ix2(1)=0

      ix2(2)=irlen
      
      
      
      
      
      
      
      goto 147
300   do jf=1,ix2(2)+2
      ix1(jf)=ix2(jf)
      
      end do
      icnt=icnt+1
      if (n(2).lt.3)goto 302
      if (icnt.eq.1000)goto 302
      
      if (ibbsw2.eq.1)goto 780
      if (ibbsw2.eq.2)goto 785
!      if (iprx1.ne.iprdl)goto 3021
!      iprdl=ipr(65000)
!      iprc=0
3021  goto 114
302   icnt=0
      
      if ((ix2(2).eq.1).and.(ix2(3).eq.1))goto 1400
      do jf=1,ix2(2)+2
      karr(jf)=ix2(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      karb(jf)=kcarr(jf)
      end do
      
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      
      call subgcd2
      
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      if (ibbsw2.eq.1)goto 780
      if (ibbsw2.eq.2)goto 785
      goto 114
330   irecnn=irecnn+1
      if (irecnn.eq.88)goto 710

      read(2,4992,rec=irecnn)(ipr(jf),jf=1,65000)
      
      iprc=1
      iprx1=ipr(iprc)
      iprdl=ipr(65000)
      goto 1141
360   print *,'prime pass'
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do

      do ii=1,icnt
      do jf=1,iydf(ii,2)+2
      karb(jf)=iydf(ii,jf)
      end do
      call subgcd2
      do jj=2,kard(2)+2
      if (kard(jj).lt.n(jj))goto 3603
      goto 3604

3603  if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
3604  end do
      end do
      print *,'number probable prime'
370   print *,'problem to be investigated'
      stop
700   jrecnn=jrecnn+1
      if (jrecnn.eq.319)goto 710
      read(4,5992,rec=jrecnn)(ipr(jf),jf=1,65000)
      goto 902

710   print *,'limit reached'
      do jf=1,igcd(2)+2
      karb(jf)=igcd(jf)
      end do
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      call subgcd2
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      stop
   





800   a=a
      do jf=1,ix2(2)+2
      ixl(jf)=ix2(jf)
      karr(jf)=ix2(jf)
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=1,kcarr(2)+2
      igcd(jf)=kcarr(jf)
      end do
      if (iph2.eq.1)goto 710


      
      icnt=0
      do jf=3,ixp(2)+2
      karr(jf-2)=ixp(jf)
      kbarr(jf-2)=ixp(jf)
      end do
      ilen=ixp(2)
      ilen2=ilen
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
      ixdf(1,jf+2)=irrr(jf)
      end do
      ixdf(1,1)=0
      ixdf(1,2)=irlen
!      do i=1,49
      do i=1,499
      do jf=3,ixdf(i,2)+2
      karr(jf-2)=ixdf(i,jf)
      end do
      do jf=3,ixdf(1,2)+2
      kbarr(jf-2)=ixdf(1,jf)
      end do
      ilen=ixdf(i,2)
      ilen2=ixdf(1,2)
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
      ixdf(i+1,jf+2)=irrr(jf)
      end do
      ixdf(i+1,1)=0
      ixdf(i+1,2)=irlen
      end do
      
      
      

805   iq1=ipr(iprc)-ipr(iprc-1)
      if (iq1.gt.1000)goto 820
      
      idcn=iq1/2
      

      do jf=3,ixl(2)+2
      karr(jf-2)=ixl(jf)
      end do
      do jf=3,ixdf(idcn,2)+2
      kbarr(jf-2)=ixdf(idcn,jf)
      end do
      ilen=ixl(2)
      ilen2=ixdf(idcn,2)
806   call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen=ilen3
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      
      
      if (irlen.eq.0)goto 360
      do jf=1,irlen
      ixl(jf+2)=irrr(jf)
      end do
      ixl(1)=0
      ixl(2)=irlen


      icnt=icnt+1
      goto 830
!     if (icnt.eq.100)goto 830
8201  iprc=iprc+1
!      print *,'ipr2',ipr(iprc)
      if (iprc.eq.65001)goto 900
      if (ipr(iprc).gt.l2)goto 710

      
      if ((ipr(iprc).eq.0).and.(irecnn.eq.0))goto 900
      
      goto 805
820   ibbsw2=1
      do jf=1,ixp(2)+2
      ix1(jf)=ixp(jf)
      ix2(jf)=ixp(jf)
      end do
      
      goto 131

830   do jf=1,ixl(2)+2
      karr(jf)=ixl(jf)
      
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      iydf(icnt,jf)=kcarr(jf)
      end do
      ilen=kcarr(2)
      iydf(icnt,1)=kcarr(1)
      iydf(icnt,2)=kcarr(2)
      do jf=3,igcd(2)+2
      kbarr(jf-2)=igcd(jf)
      end do
      ilen2=igcd(2)
      call mpmul(ilen,ilen2,ilen3)
      do jf =1,ilen3
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
      igcd(jf+2)=irrr(jf)
      end do
      igcd(1)=0
      igcd(2)=irlen





      if (icnt.ne.1000)goto 8201
      print *,'ipr2',ipr(iprc)
      icnt=0
      do jf=1,igcd(2)+2
      karb(jf)=igcd(jf)
      end do
      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      call subgcd2
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      igcd(1)=0
      igcd(2)=1
      igcd(3)=1
      goto 8201

      
      




780   ibbsw2=0
      
      do jf=3,ixl(2)+2
      karr(jf-2)=ixl(jf)
      end do
      ilen=ixl(2)
      do jf=1,ix2(2)+2
      kbarr(jf-2)=ix2(jf)
      end do
      goto 806
785   ibbsw2=0
      
      do jf=3,ixl(2)+2
      karr(jf-2)=ixl(jf)
      end do
      ilen=ixl(2)
      do jf=3,ix1(2)+2
      kbarr(jf-2)=ix1(jf)
      end do
      ilen2=ix1(2)
      goto 9002
      










900   irecnn=irecnn+1
      icnt=0
!      if (irecnn.eq.20)goto 700
      if (irecnn.ge.88)goto 700
      read(2,4992,rec=irecnn)(ipr(jf),jf=1,65000)
902   iprc=1
      if (ipr(iprc).gt.l2)goto 940
      
      icdn=ipr(1)-iprdl
!      print *,'icdn',icdn,'ipr1',ipr(1),'iprdl',iprdl
      iprdl=ipr(65000)
9001  if (icdn.gt.1000)goto 920

      icdn=icdn/2
      
      do jf=3,ixl(2)+2
      karr(jf-2)=ixl(jf)
      end do
      ilen=ixl(2)
      
      do jf=3,ixdf(icdn,2)+2
      kbarr(jf-2)=ixdf(icdn,jf)
      end do
      ilen2=ixdf(icdn,2)
9002  call mpmul(ilen,ilen2,ilen3)
      do jf=1,ilen3
      karr(jf)=kcarr(jf)
      end do
      ilen=ilen3
      do jf=3,n(2)+2
      kbarr(jf-2)=n(jf)
      end do
      ilen2=n(2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 710
      do jf=1,irlen
      ixl(jf+2)=irrr(jf)
      end do
      ixl(1)=0
      ixl(2)=irlen
      icnt=icnt+1
!     if (icnt.ne.100)goto 9201
      
      do jf=1,ixl(2)+2
      karr(jf)=ixl(jf)
      
      end do
      kbarr(1)=0
      kbarr(2)=1
      kbarr(3)=1
      call mpadd(1)
      do jf=3,kcarr(2)+2
      karr(jf-2)=kcarr(jf)
      iydf(icnt,jf)=kcarr(jf)
      end do
      ilen=kcarr(2)
      iydf(icnt,1)=0
      iydf(icnt,2)=kcarr(2)
      do jf=3,igcd(2)+2
      kbarr(jf-2)=igcd(jf)
      end do
      ilen2=igcd(2)
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
      igcd(jf+2)=irrr(jf)
      end do
      igcd(1)=0
      igcd(2)=irlen
      if (icnt.ne.1000)goto 9201
      print *,ipr(iprc),icdn
      icnt=0
      do jf=1,igcd(2)+2
      karb(jf)=igcd(jf)
      end do




      do jf=1,n(2)+2
      kara(jf)=n(jf)
      end do
      call subgcd2
      if (kard(2).gt.1)goto 1500
      if (kard(3).gt.1)goto 1500
      igcd(1)=0
      igcd(2)=1
      igcd(3)=1


9201  a=a
      iprc=iprc+1
      if (iprc.eq.65001)goto 900
!      print *,ipr(iprc),icdn
      if (ipr(iprc).gt.l2)goto 940
      icdn=ipr(iprc)-ipr(iprc-1)
      goto 9001
      




      
      
920   do jf=1,ixp(2)+2
      ix1(jf)=ixp(jf)
      ix2(jf)=ixp(jf)
      end do
      
      
      ibbsw2=2
      iq1=icdn
      goto 131
930   iprdl=ipr(iprc)
      goto 900

      
940   print *,'limit reached'
      stop
1400  print *,'code trivial or algorithm unsuitable'

1500  print *,'endwell prime=',ipr(iprc)
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
      
      





















                    
      
      
      
      

      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
      common n(50)
      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     print *,'ilenz',ilen,ilen2
10    do i=1,800
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
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
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





      subroutine mpadd(isora)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
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

      subroutine subgcd2
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
      common n(50)
      do jf=3,karb(2)+2
      karr(jf-2)=karb(jf)
      end do
      ilen=karb(2)
      
      
      do jf=3,kara(2)+2
      kbarr(jf-2)=kara(jf)
      end do
      ilen2=kara(2)
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
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
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      
      common ipr(65000)
      common kara(90),karb(90),kard(90),karp(90),karv(90)
      
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
