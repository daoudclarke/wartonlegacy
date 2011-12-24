       program alex
       common ibarray(2000),isarray(1000),igarray(1000),inv(25)
       common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,200),kdoub2(6,200),kdoub3(12,400)
       common ipowar(150)
       dimension kdsol(5,200),jpowar(150),kbper(5,50)
       
       dimension mpol(6),ipol(6),idbarr(12,12),iper(6),ibper(6)
       dimension niv(1000),match(1000),litd(1000),litt(1000),n(100)
       open(unit=1,file='kernel',access='direct',form=&
       'formatted',recl=1000,status='old')
       open(unit=2,file='matches',access='direct',form=&
       'formatted',recl=2000,status='old')
       
       kmx=200
       mmm=kmx
       read(1,1,rec=3)(niv(jk),jk=1,kmx)
1      format(200i5)       
       print *,'niv',(niv(jk),jk=1,kmx)
       
       read(2,2,rec=1)(match(jk),jk=1,kmx*2)
2      format(400i5)       
       n(1)=0
       n(2)=2
       n(3)=3005
       n(4)=3021
       
       m1(1)=0
       m1(2)=1
       m1(3)=13
       m2(1)=0
       m2(2)=1
       m2(3)=18
       litd(1)=0
       litd(2)=1
       litd(3)=1
       iconty=0
       
       do i=1,kmax
       
       if (niv(i).eq.0)goto 700
       iconty=iconty+1
       
       ilen=1
       
       kia=match(2*i-1)
       if (kia.gt.0)goto 651
       isgn=1
       goto 652
651    isgn=0
652    karr(1)=abs(kia)
       ilen=1

       
       ilen2=m2(2)
       do jf=3,m2(2)+2
       kbarr(jf-2)=m2(jf)
       end do
       call mpmul(ilen,ilen2,ilen3)
       litt(2)=ilen3
       do jf=1,ilen3
       litt(jf+2)=kcarr(jf)
       end do
       litt(1)=mod(m2(1)+isgn,2)
       ilen=1
       kib=match(i*2)
       if (kib.gt.0)goto 653
       isgn=1
       goto 654
653    isgn=0
654    karr(1)=abs(kib)
       ilen2=m1(2)
       do jf=1,ilen2
       kbarr(jf)=m1(jf+2)
       end do
       call mpmul(ilen,ilen2,ilen3)
       kbarr(2)=ilen3
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=mod(m1(1)+isgn,2)
       do jf=1,litt(2)+2
       karr(jf)=litt(jf)
       end do
       
       
       call mpadd(1)
       do jf=1,kcarr(2)+2
       litt(jf)=kcarr(jf)
       end do
       ilen=litd(2)
       ilen2=litt(2)
       do jf=1,ilen
       karr(jf)=litd(jf+2)
       end do
       do jf=1,ilen2
       kbarr(jf)=litt(jf+2)
       end do
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       litd(jf+2)=kcarr(jf)
       end do
       litd(2)=ilen3
       litd(1)=mod(litd(1) +litt(1),2)
       print *,'litd',(litd(jf),jf=1,litd(2)+2),'iconty',iconty,i
700    end do
       print *,'litd',(litd(jf),jf=1,litd(2)+2)
       end
      subroutine mpmul(ilen,ilen2,ilen3)
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,200),kdoub2(6,200),kdoub3(12,400)
      common ipowar(150)
      print *,'lenzy',ilen,ilen2
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,200),kdoub2(6,200),kdoub3(12,400)
      common ipowar(150)
      dimension kdum(200),isub(200)
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,200),kdoub2(6,200),kdoub3(6,400)
      common ipowar(150)
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


      
      
      

