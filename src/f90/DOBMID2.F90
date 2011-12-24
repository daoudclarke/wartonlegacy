      

      program dobmid2
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      dimension idrv5(5),idrv4(5),idrv3(5),idrv2(5),idrv1(5)
      dimension match(2000),iv(1000)
      dimension ksm(20,1000),ktemp(5,20),ktempl(5)
      irecnn=0
      ipoin=0
      open(unit=1,file='kernel2',access='direct',form=&
      'formatted',recl=2610,status='old')
      
      
      
      
      open(unit=2,file='matches2',access='direct',form=&
      'formatted',recl=5220,status='old')
      open(unit=3,file='square2',access='direct',form=&
      'formatted',recl=2000,status='old')
      open(unit=4,file='minpol2',access='direct',form=&
      'formatted',recl=80,status='old')
      do i=1,6
      do j=1,20
      kmpol(i,j)=0
      end do
      end do
      do i=1,6
      do j=1,500
      kdoub1(i,j)=0
      kdoub2(i,j)=0
      end do 
      end do



      
      m2(4)=9941
      ia0(1)=0
      ia0(2)=1
      ia0(3)=11
      ia0(4)=2813
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
      ia3(4)=843
      ia4(1)=0
      ia4(2)=1
      ia4(3)=118
      ia4(4)=7110
      ia5(1)=1
      ia5(2)=1
      ia5(3)=191
      ia5(4)=1487
      idrv5(1)=1
      idrv5(2)=3
      idrv5(3)=1
      idrv5(4)=7215
      idrv5(5)=7435
      idrv4(1)=0
      idrv4(2)=1
      idrv4(3)=472
      idrv4(4)=8440
      idrv3(1)=0
      idrv3(2)=1
      idrv3(3) =117
      idrv3(4)=2529
      idrv2(1)=0
      idrv2(2)=1
      idrv2(3)=170
      idrv2(4)=784
      idrv1(1)=0
      idrv1(2)=1
      idrv1(3)=105
      idrv1(4)=282
      mm=522
      read(1,1,rec=1)(iv(ii),ii=1,mm)
1     format(522i5)
2     format(1044i5)
      do i=1,mm
      if (iv(i).ne.0)goto 3
      end do
      print *,'file error'
      stop
3     read(2,2,rec=1)(match(ii),ii=1,mm*2)
      ktempl(1)=ia5(2)
      ione =i
      do ii=3,ia5(2)+2
      karr(ii-2)=ia5(ii)
      ktemp(1,ii-2)=ia5(ii)
      end do
      ilen=ia5(2)
      do j=2,4
      ilen2=ia5(2)
      do ii=1,ilen2
      kbarr(ii)=ia5(ii+2)
      end do
      call mpmul(ilen,ilen2,ilen3)
      ilen=ilen3
      do ii=1,ilen3
      karr(ii)=kcarr(ii)
      ktemp(j,ii)=kcarr(ii)
      end do
      ktempl(j)=ilen3
      print *,'ktempl',j,ktempl(j)
      end do
      
      kmpol(1,1)=0
      kmpol(1,2)=1
      kmpol(1,3)=1
      
      if(ia4(2).eq.0)goto 200
      
      kmpol(2,1)=mod(ia4(1),2)
      do ii=2,ia4(2)+2
      kmpol(2,ii)=ia4(ii)
      end do
      goto 202
200   kmpol(2,1)=0
      kmpol(2,2)=0
202   if(ia3(2).eq.0)goto 204
      
      ilen=ia3(2)
      do ii=1,ilen
      karr(ii)=ia3(ii+2)
      end do
      ilen2=ktempl(1)
      do ii=1,ilen2
      kbarr(ii)=ktemp(1,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kmpol(3,2)=ilen3
      do ii=1,ilen3
      kmpol(3,ii+2)=kcarr(ii)
      end do
      kmpol(3,1)=mod(ia3(1)+ia5(1),2)
      goto 206
204   kmpol(3,1)=0
      kmpol(3,2)=0
206   if(ia2(2).eq.0)goto 208
      
      ilen=ia2(2)
      do ii=1,ilen
      karr(ii)=ia2(ii+2)
      end do
      ilen2=ktempl(2)
      do ii=1,ilen2
      kbarr(ii)=ktemp(2,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kmpol(4,2)=ilen3
      do ii=1,ilen3
      kmpol(4,ii+2)=kcarr(ii)
      end do
      kmpol(4,1)=mod(ia2(1)+ia5(1)*2,2)
      goto 210
208   kmpol(4,1)=0
      kmpol(4,2)=0
210   if (ia1(2).eq.0)goto 212
      ilen=ia1(2)
      do ii=1,ilen
      karr(ii)=ia1(ii+2)
      end do
      ilen2=ktempl(3)
      do ii=1,ilen2
      kbarr(ii)=ktemp(3,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kmpol(5,2)=ilen3
      do ii=1,ilen3
      kmpol(5,ii+2)=kcarr(ii)
      end do
      kmpol(5,1)=mod(ia1(1)+ia5(1)*3,2)
      goto 214
212   kmpol(5,1)=0
      kmpol(5,2)=0
214   if(ia0(2).eq.0)goto 216
      ilen=ia0(2)
      do ii=1,ilen
      karr(ii)=ia0(ii+2)
      end do
      ilen2=ktempl(4)
      do ii=1,ilen2
      kbarr(ii)=ktemp(4,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      
      
      
      kmpol(6,2)=ilen3
      do ii=1,ilen3
      kmpol(6,ii+2)=kcarr(ii)
      end do
      kmpol(6,1)=mod(ia0(1)+ia5(1)*4,2)
      print *,'ione',ione
      
      
      goto 218
216   kmpol(6,1)=0 
      kmpol(6,2)=0
      
218   do ii=3,ia5(2)+2
      karr(ii-2)=ia5(ii)
      end do
      ilen=ia5(2)
      kia=match(ione*2-1)
      print *,'kia=',kia
      iabsa=abs(kia)
      if (iabsa.lt.10000)goto 10
      kbarr(1)=int(iabsa/10000)
      kbarr(2)=iabsa-kbarr(1)*10000
      ilen2=2
      goto 12
10    kbarr(1)=iabsa
      ilen2=1
12    call mpmul(ilen,ilen2,ilen3)
      do ii=1,ilen3
      kdoub1(2,ii+2)=kcarr(ii)
      end do
      kdoub1(2,2)=ilen3
      isgn=0
      if (kia.lt.0)goto 20
      goto 22
20    isgn=1
22    kdoub1(2,1)=mod(isgn+ia5(1),2)
      print *,'firkdb2',(kdoub1(2,jf),jf=1,kdoub1(2,2)+2)
      kib=match(2*ione)
      print *,'kib=',kib
      if (kib.lt.10000)goto 24
      kdoub1(1,3)=int(kib/10000)
      kdoub1(1,4)=kib-kdoub1(1,3)*10000
      kdoub1(1,2)=2
      goto 26
24    kdoub1(1,3)=kib
      kdoub1(1,2)=1
26    kdoub1(1,1)=1 
      idegm=1
      print *,'idegm',idegm,kdoub1(1,1),kdoub1(1,2),kdoub1(1,3)
      icc=1
      do i=ione+1,mm
      if (iv(i).eq.0)goto 43
      icc=icc+1
      do ii=3,ia5(2)+2
      karr(ii-2)=ia5(ii)
      end do
      ilen=ia5(2)
      iabsa=abs(match(i*2-1))
      print *,'kia2',match(i*2-1)
      if (iabsa.lt.10000)goto 30
      kbarr(1)=int(iabsa/10000)
      kbarr(2)=iabsa-kbarr(1)*10000
      ilen2=2
      goto 32
30    kbarr(1)=iabsa
      ilen2=1
32    call mpmul(ilen,ilen2,ilen3)
      kdoub2(2,2)=ilen3
      do ii=1,ilen3
      kdoub2(2,ii+2)=kcarr(ii)
      end do
      isgn=0
      if (match(i*2-1).lt.0)goto 34
      goto 36
34    isgn=1
36    kdoub2(2,1)=mod(isgn+ia5(1),2)
      print *,'seckdb2',(kdoub2(2,jf),jf=1,kdoub2(2,2)+2)
      kib=match(i*2)
      print *,'kib2=',kib
      if (kib.lt.10000)goto 40
      kdoub2(1,3)=int(kib/10000)
      kdoub2(1,4)=kib-kdoub2(1,3)*10000
      kdoub2(1,2)=2
      goto 42
40    kdoub2(1,3)=kib
      kdoub2(1,2)=1
42    kdoub2(1,1)=1
      
      idegn=1
      idegp=5
      
      call doubmul(idegm,idegn,idegp,idprod)
      idegm=idprod
      do ii=1,idegm+1
      do jj=1,kdoub3(ii,2)+2
      kdoub1(ii,jj)=kdoub3(ii,jj)
      end do
      end do
      do kf=1,idegm+1
      print *,'k',(kdoub1(kf,jf),jf=1,kdoub1(kf,2)+2)
      end do
      
43    end do
      print *,'icc',icc,'idprod',idprod
      
      idegsm=idprod
      do ii=1,idegm+1
      do jj=1,kdoub3(ii,2)+2
      ksm(ii,jj)=kdoub3(ii,jj)
      end do
      end do
      
      ilen=ia5(2)
      do ii= 1,ilen
      karr(ii)=ia5(ii+2)
      end do

      idegm=4
      idegn=4
      kdoub1(1,1)=0
      kdoub1(1,2)=1
      kdoub1(1,3)=5
      
      do ii=1,idrv4(2)+2
      kdoub1(2,ii)=idrv4(ii)
      end do
      ilen=idrv3(2)
      do ii=1,ilen
      karr(ii)=idrv3(ii+2)
      end do
      ilen2=ktempl(1)
      do ii=1,ilen2
      kbarr(ii)=ktemp(1,ii)
      print *,'ktemp',ktemp(1,ii),ilen2
      end do
      call mpmul(ilen,ilen2,ilen3)
      do ii=1,ilen3
      kdoub1(3,ii+2)=kcarr(ii)
      end do
      kdoub1(3,2)=ilen3
      kdoub1(3,1)=mod(idrv3(1)+ia5(1),2)
      print *,'kdub3',(kdoub1(3,jf),jf=1,kdoub1(3,2)+2)
      ilen=idrv2(2)
      do ii=1,ilen
      karr(ii)=idrv2(ii+2)
      end do
      ilen2=ktempl(2)
      do ii=1,ilen2
      kbarr(ii)=ktemp(2,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do ii=1,ilen3
      kdoub1(4,ii+2)=kcarr(ii)
      end do
      kdoub1(4,2)=ilen3
      kdoub1(4,1)=mod(idrv2(1)+ia5(1)*2,2)
      ilen=idrv1(2)
      do ii=1,ilen
      karr(ii)=idrv1(ii+2)
      end do
      ilen2=ktempl(3)
      do ii=1,ilen2
      kbarr(ii)=ktemp(3,ii)
      end do
      call mpmul(ilen,ilen2,ilen3)
      do ii=1,ilen3
      kdoub1(5,ii+2)=kcarr(ii)
      end do
      kdoub1(5,2)=ilen3
      kdoub1(5,1)=mod(idrv1(1) +ia5(1)*3,2)
      do ii=1,idegm+1
      do jj=1,kdoub1(ii,2)+2
      kdoub2(ii,jj)=kdoub1(ii,jj)
      end do
      print *,'besq',(kdoub2(ii,jf),jf=1,kdoub2(ii,2)+2)
      write(4,1004,rec=ii+1+idprod)(kdoub2(ii,jf),jf=1,20)
      end do
      
      
      call doubmul(idegm,idegn,5,idprod)
      do ii=1,idprod+1
      do jj=1,kdoub3(ii,2)+2
      kdoub2(ii,jj)=kdoub3(ii,jj)
      end do
      print *,'sqdrv',(kdoub2(ii,jf),jf=1,kdoub2(ii,2)+2)
      end do
      kbarr(1)=137
      ilen2=1
      do ii=1,idprod+1
      do jj=3,kdoub2(ii,2)+2
      karr(jj-2)=kdoub2(ii,jj)
      end do
      ilen=kdoub2(ii,2)
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      if (irlen.eq.0)goto 555
      print *,'rems',ii,irrr(1)
555   end do


      
      do ii=1,idegsm+1
      do jj=1,ksm(ii,2)+2
      kdoub1(ii,jj)=ksm(ii,jj)
      end do
      end do
      idegm=idegsm
      idegn=idprod
      call doubmul(idegm,idegn,5,idprod)
      do i=1,idprod+1
      print *,'final',(kdoub3(i,jf),jf=1,kdoub3(i,2)+2)
      write(3,1003,rec=i)(kdoub3(i,jf),jf=1,500)
      end do
      do i=1,idprod+1
      write(4,1004,rec=i)(kmpol(i+1,jk),jk=1,20)
      end do
1003  format(500i4)
1004  format(20i4) 
      print *,'icc',icc
      

500   close(unit=2)
      close(unit=1)
      close(unit=3)
      close(unit=4)
      end
      


   
      
      
      subroutine addstar(in,out)

      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      
      do i=1,2000
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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      dimension kdum(500),isub(500)
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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)

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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
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
      do i=1,irlen
      kart3(i+2)=irrr(i)
      end do
      kart3(2)=irlen
      kart3(1)=0
      print *,'kart3',(kart3(i),i=1,kart3(2)+2)
      
      
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
      print *,'kart1',(kart1(i),i=1,kart1(2)+2)
      
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
      print *,'firkarv',(karv(i),i=1,karv(2)+2)


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
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
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





      subroutine doubmul(idegm,idegn,idegp,idprod)
      common karr(500),kbarr(500),kcarr(2000),ipqt(500),irrr(500)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2)
      
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      
      dimension mmul(1000)
      do i=1,12
      do j=1,1000
      kdoub3(i,j)=0
      end do
      end do
      do i=1,idegm+1
      do j=1,idegn+1
      ilen=kdoub1(i,2)
      ilen2=kdoub2(j,2)
      if ((ilen.eq.0).or.(ilen2.eq.0))goto 10
      do i2=1,ilen
      karr(i2)=kdoub1(i,i2+2)
      end do
      ilen2=kdoub2(j,2)
      do j2=1,ilen2
      kbarr(j2)=kdoub2(j,j2+2)
      end do
      call mpmul(ilen,ilen2,ilen3)
      kbarr(2)=ilen3
      do k2=1,ilen3
      kbarr(k2+2)=kcarr(k2)
      end do
      kbarr(1)=mod(kdoub1(i,1)+kdoub2(j,1),2)
      goto 12
10    kbarr(1)=0      
      kbarr(2)=0
12    do i2=1,kdoub3(i+j-1,2)+2
      karr(i2)=kdoub3(i+j-1,i2)
      end do
      call mpadd(0)
      do k2=1,kcarr(2)+2
      kdoub3(i+j-1,k2)=kcarr(k2)
      end do
      end do
      end do
      
      

      idegg=idegm+idegn
      idd=idegp
      print *,'degs',idegg,idegp
      if (idegg.lt.idegp)goto 100
      
      do j=1,idegg-idegp+1
      do i2=1,kdoub3(j,2) +2
      mmul(i2)=kdoub3(j,i2)
      end do
      do i=1,idd
      do i2=1,kdoub3(j,2)+2
      kdoub3(j,i2)=0
      end do
      ilen=mmul(2)
      if (ilen.eq.0)goto 50
      do i2=3,mmul(2)+2
      karr(i2-2)=mmul(i2)
      end do
      ilen2=kmpol(i+1,2)
      if (ilen2.eq.0)goto 50
      
      do i2=3,kmpol(i+1,2)+2
      kbarr(i2-2)=kmpol(i+1,i2)
      end do
      
      call mpmul(ilen,ilen2,ilen3)
      
      do i2=1,ilen3
      kbarr(i2+2)=kcarr(i2)
      end do
      kbarr(2)=ilen3
      kbarr(1)=mod(mmul(1)+kmpol(i+1,1)+1,2)
      
      
      do i2=1,kdoub3(i+j,2)+2
      karr(i2)=kdoub3(i+j,i2)
      end do
      call mpadd(0)
      do i2=1,kcarr(2)+2
      kdoub3(i+j,i2)=kcarr(i2)
      end do
      
      
      goto 60
50    kdoub3(i+j,1)=0
      kdoub3(i+j,2)=0
60    end do
      end do
      do i=1,idegg+1
      if (kdoub3(i,2).ne.0)goto 70
      end do
      print *,'error algebraic nos'
      stop
70    idprod =idegg+1-i      
      igap=i-1
      print *,'igap',igap
      print *,'kdub3',(kdoub3(mj,2),mj=1,idegg+1-i)
      if (igap.eq.0)goto 102
      
      do i=1,idegg+1-igap
      do j=1,kdoub3(i+igap,2)+2
      kdoub3(i,j)=kdoub3(i+igap,j)
      end do
      end do

      
      goto 102
100   idprod=idegg
102   return 
      end
