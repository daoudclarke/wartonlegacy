       program bobkrz3
!      third generalized exponent recovery program
!      matrix obtained from double quadratic sieve
!      computation of kernel of matrix uses fastest procedure       
! third of 3 kernel-obtaining programs      
      common karr(200),kbarr(200),kcarr(200),ipqt(200),irrr(200)
      common mnum(50)
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      common iarq(2)
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
       
       dimension indic(10000),ipv(3000,20)
       dimension ic(10000),iset(10000)
       dimension nar(3000,20),littr(20)
       dimension n(60),jfak(30),jfreq(30)
       dimension ipd(50),mul(50),isum(50),irno(60)
       dimension ibase(50),jbase(50),itarg(50),jtarg(50)
       dimension iarr1(10000)
       
       do jf=1,60
       irno(jf)=0
       end do
       iig=1
       open (unit=9,file='expparn',access='direct',form=&
       'formatted',recl=240,status='old')
1006   format (60i4)       
       read (9,1006,rec=1)(n(jf),jf=1,60)
       print *,'n',(n(jf),jf=1,n(2)+2)
       close (unit=9)
       
       do jf=1,n(2)+2
       karr(jf)=n(jf)
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
       ipd(jf)=mdarr(jf)
       end do
       print *,'ipd',(ipd(jf),jf=1,ipd(2)+2)
       
       

       open(unit=7,file='exppar',access='direct',form=&
       'formatted',recl=344,status='old')
       read(7,6,rec=1)irecnn,kkll,ijpow,nfak,(ibase(jf),jf=1,20),&
       (jbase(jf),jf=1,20),(itarg(jf),jf=1,20),(jtarg(jf),jf=1,20)
       print *,'kkll',kkll
       close (unit=7)
       
       
      

       
700    nunit=1000000/(kkll*(ipd(2)+2))
       
701    ntem=kkll/nunit 
      nzz=(ntem+2)*nunit 
      nzz3=318+5*nzz
       print *,'nzz3',nzz3
       open (unit=3,file='pindic',access='direct',form=&
       'formatted',recl=nzz3,status='old')
       read (3,1007,rec=1)nzz4,kkll,nzz,ir,(irno(jf),jf=1,60),&
       (indic(jf),jf=1,nzz)
       close (unit=3)
       print *,'nzz4',nzz4,'kkll',kkll,'nzz',nzz,'ir',ir
       print *,'irnos',(irno(jf),jf=1,ir)
       print *,'indic',(indic(jf),jf=1400,nzz)
       
       narc9=nzz*4*(ipd(2)+2)
6      format(i6,i6,i8,i4,20i4,20i4,20i4,20i4)       
       open(unit=1,file='ptran',access='direct',form=&
       'formatted',recl=narc9,status='old')
       
      
      nzz2=narc9
      nzz4=irecnn
      open (unit=2,file='pvec',access='direct',form=&
      'formatted',recl=nzz2,status='old')
      
1007  format (3i5,i3,60i5,20000i5)
       mm=nzz
!       goto 12347
       read(2,2,rec=3)((ipv(j,jf),jf=1,ipd(2)+2),j=1,mm)
       print *,'ipv1468',(ipv(1468,jf),jf=1,ipv(1468,2)+2)
       print *,'ipv1469',(ipv(1469,jf),jf=1,ipv(1469,2)+2)
       print *,'ipv1470',(ipv(1470,jf),jf=1,ipv(1470,2)+2)
       
1008   format (i1,i10,20i4,5000(20i4))
       nzz5=91+kkll*4*(ipd(2)+2)
       open(unit=3,file='zbexp1',access='direct',form=&
       'formatted',recl=nzz5,status='old')
       isum(1)=0
       isum(2)=0
       iggf=15
       do i=1,irno(3)
       print *,'i',i
       read (3,1008,rec=i)kkb,ihitn,(littr(jf),jf=1,20),&
       ((nar(j,jf),jf=1,ipd(2)+2),j=1,kkll)
       if ((nar(iggf,2).eq.0).or.(ipv(i,2).eq.0))goto 12348
       print *,'nar',(nar(iggf,jf),jf=1,nar(iggf,2)+2),'ipv',&
       (ipv(i,jf),jf=1,ipv(i,2)+2)
       
       do jf=1,ipd(2)+2
       marr(jf)=nar(iggf,jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipv(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isum(2)+2
       kbarr(jf)=isum(jf)
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
       isum(jf)=mcarr(jf)
       end do
       print *,'i',i,'isum',(isum(jf),jf=1,isum(2)+2)
12348  end do
       stop
12347  a=a       
       do i=1,mm
       iarr1(i)=0
       end do
       do i=1,mm
       do jf=1,ipd(2)+2
       ipv(i,jf)=0
       end do
       end do
       ir=10
       do i=1,ir
       inj=irno(i)
       iarr1(inj)=1
       end do
       do ibigl=1,ir
       kk1=1
       do jbigl=1,irno(ibigl)-1
       if (iarr1(jbigl).eq.1)goto 100
       iset(kk1)=jbigl
       kk1=kk1+1
100    end do
       iset(kk1)=irno(ibigl)
       ileng=kk1
       print *,'iset1400',(iset(jf),jf=1400,irno(1))
       print *,'indic1462',indic(1462)
       
       print *,'ileng',ileng,'irno1',irno(1)
       
       do kbigl=1,irno(ibigl)
       do jf=1,ipd(2)+2
       ipv(kbigl,jf)=0
       end do
       end do
       inj=irno(ibigl)
       ipv(inj,1)=0
       ipv(inj,2)=1
       ipv(inj,3)=1
       isum(1)=0
       isum(2)=0
       
       do kjay=ileng,2,-1
       inj=iset(kjay-1)
       print *,'inj',inj,'indicinj',indic(inj)
       if (indic(inj).ne.0)goto 12345
       indic(inj)=irno(ibigl)
       print *,'indicinj',inj
12345  inj2=indic(inj)
       read (1,2,rec=inj2)((nar(jk,jf),jf=1,ipd(2)+2),jk=1,mm)
       isum(1)=0
       isum(2)=0
       
       do kkjay=kjay,ileng
       inj=iset(kkjay)
       
       
12346  do jf=1,ipd(2)+2
       marr(jf)=nar(inj,jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipv(inj,jf)
       end do
!       print *,'nar',(nar(inj,jf),jf=1,nar(inj,2)+2)
!       print *,'ipv',(ipv(inj,jf),jf=1,ipv(inj,2)+2)
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isum(2)+2
       kbarr(jf)=isum(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isum(jf)=kcarr(jf)
       end do
       end do
       do jf=1,isum(2)+2
       marr(jf)=isum(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       inj2=iset(kjay-1)
       do jf=1,mcarr(2)+2
       ipv(inj2,jf)=mcarr(jf)
       
       end do
!       print *,'inj2',inj2,'ipv',(ipv(inj2,jf),jf=1,ipv(inj2,2)+2)
       end do
       
       write (2,2,rec=ibigl)((ipv(jk,jf),jf=1,ipd(2)+2),jk=1,mm)
       print *,'ibigl',ibigl
       end do





       
       
       
       
1      format(i1,i10,20i4,5000(20i4))       
       
       










       
       

2      format(5000(20i4))     
       
       
       close(unit=1)
       close(unit=2)
       
       close (unit=9)
       close (unit=7)
       
       end

       subroutine subbw6(ia,ib,ivv)
       ib=mod(ib,ia)
       if (ib.ge.0)goto 10
       ib=ia+ib
10     iu=1
       id=ia
       if (ib.eq.0)goto 888
       iv1=0
       iv3=ib
815    if (iv3.eq.0)goto 830
       iqq=int(id/iv3)
       it3=id-iqq*iv3
       it1=iu-iqq*iv1
       iu=iv1
       id=iv3
       iv1=it1
       iv3=it3
       goto 815
830    iv=(id-ia*iu)/ib
       if (iu.le.0)goto 870
       iv=iv*(1-ia)
       iv=mod(iv,ia)
       if (iv.ge.0)goto 870
       iv=iv+ia
870    a=a
       if (id.gt.1)goto 890
       print *,'inverse',iv
       goto 890
888    iv=ib
890    ivv=iv
       return
       end


      


      
      
      subroutine menmul
      common karr(200),kbarr(200),kcarr(200),ipqt(200),irrr(200)
      common mnum(50)
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      common iarq(2)
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      




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
      common karr(200),kbarr(200),kcarr(200),ipqt(200),irrr(200)
      common mnum(50)
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      
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
11    print *,'halted: attempted division by zero, routine mendiv'
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
      common karr(200),kbarr(200),kcarr(200),ipqt(200)
      common irrr(200)
      common mnum(50)
      
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      
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
      common karr(200),kbarr(200),kcarr(200),ipqt(200)
      common irrr(200)
      common mnum(50)
      
      
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      
      
      
      
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
      common karr(200),kbarr(200),kcarr(200),ipqt(200)
      common irrr(200)
      common mnum(50)
      
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      
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
      common karr(200),kbarr(200),kcarr(200),ipqt(200)
      common irrr(200)
      common mnum(50)
      
      
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      






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
      common karr(200),kbarr(200),kcarr(200),ipqt(200)
      common irrr(200)
      common mnum(50)
      
      
      common kara(200),karb(200),kard(200),karp(200),karv(200)
      
      common iarq(2)
      
      common marr(200),mbarr(200),mcarr(200),mdarr(200)
      
      
      dimension karu(200),karv1(200),karv3(200),karqq(200)
      dimension kart3(200),kart1(200)
      
      
      
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

      


      
      
      
