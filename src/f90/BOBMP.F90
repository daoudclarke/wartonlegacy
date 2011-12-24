      
      program bobmp
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
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
      goto 26
24    do k =1,irem1
      kemp =ichar(loc(ii+1-k))-48
      karr(1) =karr(1) +kemp *10**(k -1)
      end do
26    print *,karr(1),karr(2),karr(3)
      if (irem1.gt.0)goto 28
      ilen =jtemp
      goto 30









28    ilen =jtemp+1  

30    print *,'ilen',ilen
      ilen2 =5
      kbarr(1)=5
      kbarr(2)=1234
      kbarr(3)=5678
      kbarr(4)=9123
      kbarr(5) =4631

      
      do iv =1,1
      call mpmul(ilen,ilen2)
      
      
      ilen3 =ilen +ilen2
      if (kcarr(1) .ne.0)goto 32
      ilen3 =ilen3 -1
      do i =1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

      
32    print *,(kcarr(i),i=1,ilen3)
      end do
      
      karr(1)=7
      karr(2)=54
      karr(3)=80
      karr(4) =12
      karr(5)=47
      karr(6)=78
      karr(7)=79
      karr(8)=1
      karr(9)=32
      karr(10)=48
      karr(11)=34
      karr(12)=52
      karr(13)=43
      karr(14)=36
      karr(15)=10
      karr(16)=82
      karr(17)=80
      karr(18)=51
      karr(19)=62
      karr(20)=51
      karr(21)=80
      karr(22)=56
      karr(23)=73



      kbarr(1)=5
      kbarr(2)=12
      kbarr(3)=34
      kbarr(4)=56
      kbarr(5)=78
      kbarr(6)=91
      kbarr(7)=23
      kbarr(8)=46
      kbarr(9)=31



      
      ilen =23
      ilen2=9
      do i4 =1,1000
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      end do
      if(iswq.eq.0)goto 102
      if (icont.eq.0)goto 101
      print *,'quotient',(ipqt(i),i=1,icont)
101   if(irlen.eq.0)goto 102      
      print *,'remainder',(irrr(i),i=1,irlen)
      goto 102
100   print *,'input error'
102   end
      subroutine addstar(in,out)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      character*(*) in,out
      out =in
      i =len(out)
      do while(out(i:i) ==' ')
      out(i:i)='*'
      i = i-1
      end do 
      return 
      end

      subroutine mpmul(ilen,ilen2)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      
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
      return
      end

      subroutine mpdiv(ilen,ilen2,irlen,icont,iswq)
      common karr(50),kbarr(50),kcarr(200),ipqt(50),irrr(50)
      dimension kdum(50),isub(50)
      do i =1,ilen
      kdum(i)=karr(i)
      ipqt(i)=0
      end do
      
      icont2 =0
      icont3 =0
      iswm=0
      iswq=1
      ind2 =ilen2
      if (ilen2.gt.ilen)goto 85
      if (ilen2.lt.ilen)goto 10
      if (kbarr(1).gt.karr(1))goto 85
      do i=2,ilen
      if(kbarr(i).gt.karr(i))goto 85
      if(kbarr(i).lt.karr(i))goto 10
      end do

10    kbap=kbarr(1)
      lop =ilen +1 -ilen2
      ll=1
201   do i =1,ilen2
      if (kbarr(i).eq.kdum(ll+i-1))goto 202
      if (kbarr(i).gt.kdum(ll +i-1))goto 203
      
      goto 20
202   end do
      icont =ind2
      iapd =1
      goto 31
203   if(iswm.eq.1)goto 72
      iapd =int((kdum(ll)*1000000 +kdum(ll+1)*10000+kdum(ll+2)*100&
      +kdum(ll+3))/(kbarr(1) *10000 +kbarr(2) *100 +kbarr(3)+1))



      goto 30
20    kbig = kdum(ll) *1000000 +kdum(ll+1)*10000 +kdum(ll+2) *100&
      +kdum(ll+3)
      kbap =kbarr(1) *1000000 + kbarr(2) *10000 + kbarr(3) *100&
      +kbarr(4) +1
      icont =ind2
      
      
      iapd=int(kbig/kbap)
      goto 31
      
      
30    icont =ind2 +1
      


31    ipqt(icont) =ipqt(icont)+iapd
      icont2 =icont2 +1
      
      if (ipqt(icont).lt.100)goto 35
      ipqt(icont)=ipqt(icont)-100
      ipqt(icont-1)=ipqt(icont-1) +1
      if (ipqt(icont-1).lt.100)goto 35
      ipqt(icont-1) =ipqt(icont-1)-100
      ipqt(icont-2)=ipqt(icont-2) +1

35    isub(1)=0
      do i=1,ilen2
      isub(i+1)= kbarr(i)*iapd
      end do
351   do i=1,ilen2 
      ii=ilen2 +2 -i
      if (isub(ii).lt.100)goto 37
      itemp =int(isub(ii)/100)
      irem1 =isub(ii)-itemp *100
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
      kdum(ii+ig) =kdum(ii+ig) +100
      kdum(ii+ig-1) =kdum(ii+ig-1) -1
40    end do
      
      
      
      do i =1,ilens
      if((i.lt.icont+1-ilen2).and.(kdum(i).ne.0))goto 60


      if (kdum(icont +i-ilen2).eq.kbarr(i)) goto 50
      if (kdum(icont+i-ilen2).gt.kbarr(i))goto 60
      goto 70
50    end do
60    if(kdum(ll).ne.0)goto 601
      ll=ll+1      
601   if (int(kdum(ll)/kbarr(1)).lt.3)goto 61
      goto 150
61    iapd =1
      icont3 =icont3 +1
      ipqt(icont) =ipqt(icont)+1
      if (ipqt(icont).lt.100)goto 611
      ipqt(icont)=0
      ipqt(icont-1) =ipqt(icont-1) +1
      if (ipqt(icont-1).lt.100)goto 611
      ipqt(icont-1)=0
      ipqt(icont-2)=ipqt(icont-2) +1
      
      
611   isub(1) =0   
      do i =1,ilen2 
      isub(i+1) =kbarr(i)
      end do
      goto 351
150   kbap = kbarr(1)
      goto 201

70    j=j
      
      
      


      if (icont.eq.ilen)goto 72
      kbap =kbarr(1)
      ll =icont  -ilen2
      icont =icont-1
92    ll =ll+1
      icont =icont +1  
      ind2=icont
      if(kdum(ll).eq.0)goto 92
      if(ilen2.gt.ilen+1-ll)goto 72
      if(ilen2.lt.ilen+1-ll)goto 201
      iswm=1
      
      goto 201








72    do i = 1,ilen
      if (kdum(i).ne.0)goto 80
      end do
      irlen =0
      goto 86
80    do j =1,ilen-i+1
      irrr(j) =kdum(i+j-1)
      end do
      irlen =ilen +1 -i
      goto 86
85    iswq=0
86    print *,'icont2',icont2
      print *,'icont3',icont3
      return
      end




