       program bobnfs12
! for finding good polynomials for number field sieve       
       
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
       dimension kdsol(5,5000),jpowar(5000),kbper(5,5000),iroot(6,8000)
       dimension ivee(5,5000)
       dimension mpol(6),ipol(6),idbarr(12,12),iper(6),ibper(6)
       dimension niv(10000),match(20000),litd(5000),litt(1000),n(100)
       dimension itemp1(1000),itemp2(1000)
       dimension norma(10000),nrem(1000),ipre(10000),icurr(10000)
       dimension ires(1000),icdh(100),mh2(100),mpow1(5,100),mpow2(5,100)
!       dimension normar(1000),kans(2000),kpr2(10000),ipr(65000)
       dimension npow1(5,100),npow2(5,100),iasq(100),jasq(100)
       dimension ja5(10),ja4(10),ja3(10),ja2(10),ja1(10),ja0(10)
       dimension kiat(4),kibt(4),igcd(1000),isquar(5000),iprod(5000)
       dimension ifinv(6,10)
       ibsw1=0
       ibsw2=0
       itcon=0
       ibind=0
       ibsw3=0
       print *,'code length base 10000?'
       read *,n(2)
       print *,'code'
       read *,(n(jf),jf=3,n(2)+2)
       do jf=n(2)+3,n(2)+7
       litd(jf)=0
       litt(jf)=0
       end do
       do jf=1,n(2)+2
       litd(jf)=n(jf)
       litt(jf)=n(jf)
       end do
       litd(2)=litd(2)+5
       litt(2)=litt(2)+5
       do ii=1,4
       ibsw0=0
       do jf=1,litt(2)+2
       marr(jf)=litt(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ii
       call menmul
       do jf=1,mcarr(2)+2
       litd(jf)=mcarr(jf)
       end do






       
       
       inj=litd(2)/2
       do jf=1,inj+2
       ipre(jf)=litd(jf)
       end do
20     do jf=1,ipre(2)+2
       marr(jf)=ipre(jf)
       mbarr(jf)=ipre(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       mbarr(jf)=mcarr(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       mbarr(jf)=mcarr(jf)
       end do
       do jf=1,litd(2)+2
       marr(jf)=litd(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       itemp2(jf)=mdarr(jf)
       end do
       do jf=1,ipre(2)+2
       marr(jf)=ipre(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=4
       call menmul
       do jf=1,mcarr(2)+2
       itemp1(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       do jf=1,itemp2(2)+2
       kbarr(jf)=itemp2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=5
       call mendiv
       do jf=1,mdarr(2)+2
       icurr(jf)=mdarr(jf)
       end do
!       print *,'icurr',(icurr(jf),jf=1,icurr(2)+2)
!       print *,'ipre',(ipre(jf),jf=1,ipre(2)+2)
       if (icurr(2).lt.ipre(2))goto 30
       if (icurr(2).gt.ipre(2))goto 32
       do jf=3,ipre(2)+2
       if (icurr(jf).lt.ipre(jf))goto 30
       if (icurr(jf).gt.ipre(jf))goto 32
       end do
       goto 32
30     do jf=1,icurr(2)+2
       ipre(jf)=icurr(jf)
       end do
       goto 20
32     icurr(2)=icurr(2)-1
33     do jf=1,icurr(2)+2
       m1(jf)=icurr(jf)
       marr(jf)=icurr(jf)
       mbarr(jf)=icurr(jf)
       
       end do
       print *,'m1',(m1(jf),jf=1,m1(2)+2)
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       mbarr(jf)=mcarr(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,icurr(2)+2
       mbarr(jf)=icurr(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,litd(2)+2
       karr(jf)=litd(jf)
       end do
       karr(2)=karr(2)-5
       call mpadd(1)
       do jf=1,kcarr(2)+2
       ires(jf)=kcarr(jf)
       end do
       m2(1)=0
       m2(2)=1
       m2(3)=1
       print *,'ires',(ires(jf),jf=1,ires(2)+2),'ii',ii
       if (ires(2).eq.1)goto 114
!  checking for negative constant coefficient       
       if (ibsw0.eq.1)goto 41
       ibsw0=1
       do jf=1,icurr(2)+2
       karr(jf)=icurr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       icurr(jf)=kcarr(jf)
       end do
       goto 33
41     end do

       do jf=n(2)+3,n(2)+8
       litd(jf)=0
       end do
       do jf=1,n(2)+2
       litd(jf)=n(jf)
       litt(jf)=n(jf)
       end do
       litd(2)=litd(2)+6
       litt(2)=litt(2)+6
       inj=litd(2)/2
       do jf=1,inj+2
       ipre(jf)=litd(jf)
       end do
120     do jf=1,ipre(2)+2
       marr(jf)=ipre(jf)
       mbarr(jf)=ipre(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       mbarr(jf)=mcarr(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       mbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipre(2)+2
       marr(jf)=ipre(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       mbarr(jf)=mcarr(jf)
       end do





       do jf=1,litd(2)+2
       marr(jf)=litd(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       itemp2(jf)=mdarr(jf)
       end do
       do jf=1,ipre(2)+2
       marr(jf)=ipre(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=5
       call menmul
       do jf=1,mcarr(2)+2
       itemp1(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       do jf=1,itemp2(2)+2
       kbarr(jf)=itemp2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=6
       call mendiv
       do jf=1,mdarr(2)+2
       icurr(jf)=mdarr(jf)
       end do
       if (icurr(2).lt.ipre(2))goto 130
       if (icurr(2).gt.ipre(2))goto 132
       do jf=3,ipre(2)+2
       if (icurr(jf).lt.ipre(jf))goto 130
       if (icurr(jf).gt.ipre(jf))goto 132
       end do
       goto 132
130    do jf=1,icurr(2)+2
       ipre(jf)=icurr(jf)
       end do
       goto 120
132     icurr(2)=icurr(2)-1
       do jf=1,icurr(2)+2
       m1(jf)=icurr(jf)
       marr(jf)=icurr(jf)
       mbarr(jf)=icurr(jf)
       m2(jf)=icurr(jf)
       end do
       print *,'m1',(m1(jf),jf=1,m1(2)+2)
      print *,'where are we going m1',(m1(jf),jf=1,m1(2)+2),'m2',&       
      (m2(jf),jf=1,m2(2)+2)



       do jf=1,m1(2)+2
       m2(jf)=m1(jf)
       karr(jf)=m1(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       m1(jf)=kcarr(jf)
       mpow1(1,jf)=kcarr(jf)
       npow1(1,jf)=kcarr(jf)
       marr(jf)=kcarr(jf)
       mbarr(jf)=kcarr(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       m2(jf)=kcarr(jf)
       end do
36     call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       mbarr(jf)=mcarr(jf)
       mpow1(2,jf)=mcarr(jf)
       npow1(2,jf)=mcarr(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       mpow1(4,jf)=mcarr(jf)
       npow1(4,jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       npow1(5,jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       mpow1(5,jf)=mcarr(jf)
       kara(jf)=mcarr(jf)
       end do






       do jf=1,m2(2)+2
       karb(jf)=m2(jf)
       end do
       call subgcd2
       if ((kard(2).eq.1).and.(kard(3).eq.1))goto 40
       do jf=1,m2(2)+2
       karr(jf)=m2(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       m2(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       marr(jf)=m1(jf)
       mbarr(jf)=m1(jf)
       mpow1(1,jf)=m1(jf)
       npow1(1,jf)=m1(jf)
       end do
       print *,'somewhere kard',kard(2),'m1',(m1(jf),jf=1,m1(2)+2),&
       'm2',(m2(jf),jf=1,m2(2)+2)
       goto 36
40     print *,'m1',(m1(jf),jf=1,m1(2)+2)
       print *,'m2',(m2(jf),jf=1,m2(2)+2)
       do jf=1,m1(2)+2
       kara(jf)=m1(jf)
       end do
       do jf=1,m2(2)+2
       karb(jf)=m2(jf)
       mpow2(1,jf)=m2(jf)
       npow2(1,jf)=m2(jf)
       end do
       call subgcd2
       if ((kard(2).eq.1).and.(kard(3).eq.1))goto 38
       goto 52
38     do ii=1,4
       do jf=1,m2(2)+2
       marr(jf)=m2(jf)
       end do
       do jf=1,karb(2)+2
       
       mbarr(jf)=karb(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       npow2(ii+1,jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       karb(jf)=mcarr(jf)
       mpow2(ii+1,jf)=mcarr(jf)
       end do
       call subgcd2
       do jf=1,npow2(ii+1,2)+2
       karb(jf)=npow2(ii+1,jf)
       end do
       if ((kard(2).eq.1).and.(kard(3).eq.1))goto 50
       goto 52
50     end do
       print *,'m1',(m1(jf),jf=1,m1(2)+2)
       print *,'m2',(m2(jf),jf=1,m2(2)+2)
       do jf=1,mpow1(2,2)+2
       marr(jf)=mpow1(2,jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       mpow1(2,jf)=mcarr(jf)
       
       end do
       do jf=1,npow1(2,2)+2
       marr(jf)=npow1(2,jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       npow1(3,jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       mpow1(3,jf)=mcarr(jf)
       end do
       do jf=1,mpow1(4,2)+2
       marr(jf)=mpow1(4,jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       mpow1(4,jf)=mcarr(jf)
       end do
       do jf=1,npow1(4,2)+2
       marr(jf)=npow1(4,jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       npow1(5,jf)=mcarr(jf)
       end do




       do jj=1,5
       print *,'jj',jj,'mpow1',(mpow1(jj,jf),jf=1,mpow1(jj,2)+2),&
       'mpow2',(mpow2(jj,jf),jf=1,mpow2(jj,2)+2)
       print *,'jj',jj,'npow1',(npow1(jj,jf),jf=1,npow1(jj,2)+2),&
       'npow2',(npow2(jj,jf),jf=1,npow2(jj,2)+2)
       
       
       end do
       goto 60
       
52     do jf=1,m2(2)+2
       karr(jf)=m2(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       m2(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       marr(jf)=m1(jf)
       mbarr(jf)=m1(jf)
       end do
       
       goto 36
60     do jf=1,m2(2)+2       
       kara(jf)=m2(jf)
       end do
       do jf=1,mpow1(5,2)+2
       karb(jf)=mpow1(5,jf)
       end do
       call mpgcd
!       print *,'karv5',(karv(jf),jf=1,karv(2)+2)
       do jf=1,n(2)+2
       marr(jf)=n(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ia5(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ia5',(ia5(jf),jf=1,ia5(2)+2)
       print *,'npow1 5',(npow1(5,jf),jf=1,npow1(5,2)+2)
       print *,'n',(n(jf),jf=1,n(2)+2)
       print *,'m2',(m2(jf),jf=1,m2(2)+2)
       do jf=1,npow1(5,2)+2
       mbarr(jf)=npow1(5,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,n(2)+2
       karr(jf)=n(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       
       do jf=1,mpow1(4,2)+2
       karb(jf)=mpow1(4,jf)
       end do
       call mpgcd
!       print *,'karv4',(karv(jf),jf=1,karv(2)+2)
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ia4(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ia4',(ia4(jf),jf=1,ia4(2)+2)
       do jf=1,npow1(4,2)+2
       mbarr(jf)=npow1(4,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow1(3,2)+2
       karb(jf)=mpow1(3,jf)
       end do
       call mpgcd
!       print *,'karv3',(karv(jf),jf=1,karv(2)+2)
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ia3(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ia3',(ia3(jf),jf=1,ia3(2)+2)
       do jf=1,npow1(3,2)+2
       mbarr(jf)=npow1(3,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow1(2,2)+2
       karb(jf)=mpow1(2,jf)
       end do
       call mpgcd
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ia2(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       if (ia2(1).eq.0)goto 85
       if ((ibsw3.eq.1).and.(ibind.eq.1))goto 84
       goto 85
84     ibind=0       
       do jf=1,m2(2)+2
       kbarr(jf)=m2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia2(jf)=kcarr(jf)
       marr(jf)=kcarr(jf)
       end do
85     a=a


       print *,'ia2',(ia2(jf),jf=1,ia2(2)+2)
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       
       do jf=1,npow1(2,2)+2
       mbarr(jf)=npow1(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       print *,'remz',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow1(1,2)+2
       karb(jf)=mpow1(1,jf)
       end do
       call mpgcd
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ia1(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       if (ia1(1).eq.0)goto 87
       if ((ibsw3.eq.1).and.(ibind.eq.1))goto 86
       goto 87
86     ibind=0
       do jf=1,m2(2)+2
       kbarr(jf)=m2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia1(jf)=kcarr(jf)
       marr(jf)=kcarr(jf)
       end do
87     a=a

       print *,'ia1',(ia1(jf),jf=1,ia1(2)+2)
       do jf=1,npow1(1,2)+2
       mbarr(jf)=npow1(1,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m2(2)+2
       mbarr(jf)=m2(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       ia0(jf)=nrem(jf)
       end do
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       if (ibsw3.eq.1)goto 59
! second polynomial        
83     do jf=1,m1(2)+2       
       kara(jf)=m1(jf)
       end do
       do jf=1,mpow2(5,2)+2
       karb(jf)=mpow2(5,jf)
       end do
       call mpgcd
!       print *,'karv5',(karv(jf),jf=1,karv(2)+2)
       do jf=1,n(2)+2
       marr(jf)=n(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ja0(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ja0',(ja0(jf),jf=1,ja0(2)+2)
       print *,'npow2 5',(npow2(5,jf),jf=1,npow2(5,2)+2)
       print *,'n',(n(jf),jf=1,n(2)+2)
       print *,'m1',(m1(jf),jf=1,m1(2)+2)
       do jf=1,npow2(5,2)+2
       mbarr(jf)=npow2(5,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,n(2)+2
       karr(jf)=n(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       
       do jf=1,mpow2(4,2)+2
       karb(jf)=mpow2(4,jf)
       end do
       call mpgcd
!       print *,'karv4',(karv(jf),jf=1,karv(2)+2)
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ja1(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ja1',(ja1(jf),jf=1,ja1(2)+2)
       do jf=1,npow2(4,2)+2
       mbarr(jf)=npow2(4,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow2(3,2)+2
       karb(jf)=mpow2(3,jf)
       end do
       call mpgcd
!       print *,'karv3',(karv(jf),jf=1,karv(2)+2)
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ja2(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       end do
       print *,'ja2',(ja2(jf),jf=1,ja2(2)+2)
       do jf=1,npow2(3,2)+2
       mbarr(jf)=npow2(3,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nremspec',(nrem(jf),jf=1,nrem(2)+2)
       print *,'frac',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow2(2,2)+2
       karb(jf)=mpow2(2,jf)
       end do
       call mpgcd
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ja3(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       if (ja3(1).eq.0)goto 95
       if ((ibsw3.eq.1).and.(ibind.eq.1))goto 94
       goto 95
94     do jf=1,m1(2)+2  
       kbarr(jf)=m1(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja3(jf)=kcarr(jf)
       marr(jf)=kcarr(jf)
       end do
       ibind=0
95     print *,'ja3',(ja3(jf),jf=1,ja3(2)+2)
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       do jf=1,npow2(2,2)+2
       mbarr(jf)=npow2(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       end do
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       print *,'remz',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mpow2(1,2)+2
       karb(jf)=mpow2(1,jf)
       end do
       call mpgcd
       do jf=1,nrem(2)+2
       marr(jf)=nrem(jf)
       end do
       do jf=1,karv(2)+2
       mbarr(jf)=karv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       ja4(jf)=mcarr(jf)
       marr(jf)=mcarr(jf)
       karr(jf)=mcarr(jf)
       end do
       if (ja4(1).eq.0)goto 97
       if ((ibsw3.eq.1).and.(ibind.eq.1))goto 96
       goto 97
96     ibind=0       
       do jf=1,m1(2)+2
       kbarr(jf)=m1(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja4(jf)=kcarr(jf)
       marr(jf)=kcarr(jf)
       end do
97     print *,'ja4',(ja4(jf),jf=1,ja4(2)+2)
       do jf=1,npow2(1,2)+2
       mbarr(jf)=npow2(1,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,nrem(2)+2
       karr(jf)=nrem(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,m1(2)+2
       mbarr(jf)=m1(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       nrem(jf)=mdarr(jf)
       ja5(jf)=nrem(jf)
       end do
       print *,'nrem',(nrem(jf),jf=1,nrem(2)+2)
       
! determining multiples of vectors starts here       
       
59     do jf=1,ja5(2)+2
       marr(jf)=ja5(jf)
       end do
       do jf=1,ia5(2)+2
       mbarr(jf)=ia5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       iprod(jf)=mcarr(jf)
       end do
       do jf=1,ja5(2)+2
       mbarr(jf)=ja5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       isquar(jf)=mcarr(jf)
       end do
       do jf=1,ja4(2)+2
       marr(jf)=ja4(jf)
       end do
! insert 1 remove marks if required to modify routine, repeat       
! for inserts 2 to 8
!       jinj=ja4(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       
       do jf=1,ia4(2)+2
       mbarr(jf)=ia4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ja4(2)+2
       mbarr(jf)=ja4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do




       do jf=1,ja3(2)+2
       marr(jf)=ja3(jf)
       end do
! insert 2       
!       jinj=ja3(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       
       do jf=1,ia3(2)+2
       mbarr(jf)=ia3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ja3(2)+2
       mbarr(jf)=ja3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do


       do jf=1,ja2(2)+2
       marr(jf)=ja2(jf)
       end do
! insert 3       
!       jinj=ja2(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ia2(2)+2
       mbarr(jf)=ia2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ja2(2)+2
       mbarr(jf)=ja2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do

       do jf=1,ja1(2)+2
       marr(jf)=ja1(jf)
       end do
! insert 4       
!       jinj=ja1(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ia1(2)+2
       mbarr(jf)=ia1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ja1(2)+2
       mbarr(jf)=ja1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do

       do jf=1,ja0(2)+2
       marr(jf)=ja0(jf)
       end do

       do jf=1,ia0(2)+2
       mbarr(jf)=ia0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ja0(2)+2
       mbarr(jf)=ja0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do
       do jf=1,iprod(2)+2
       marr(jf)=iprod(jf)
       end do
       marr(2)=marr(2)+1
       inj=marr(2)+2
       marr(inj)=0
       do jf=1,isquar(2)+2
       mbarr(jf)=isquar(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       karr(jf)=mdarr(jf)
       end do
       kbarr(1)=mdarr(1)
       kbarr(2)=1
       kbarr(3)=5000
       call mpadd(0)
       kcarr(2)=kcarr(2)-1
       do jf=2,kcarr(2)+2
       norma(jf)=kcarr(jf)
       end do
       norma(1)=mod(kcarr(1)+1,2)
       print *,'ia5',(ia5(jf),jf=1,ia5(2)+2),'ja5',(ja5(jf),jf=1,ja5(2)+2)
       print *,'ia4',(ia4(jf),jf=1,ia4(2)+2),'ja4',(ja4(jf),jf=1,ja4(2)+2)
       print *,'ia3',(ia3(jf),jf=1,ia3(2)+2),'ja3',(ja3(jf),jf=1,ja3(2)+2)
       print *,'ia2',(ia2(jf),jf=1,ia2(2)+2),'ja2',(ja2(jf),jf=1,ja2(2)+2)
       print *,'ia1',(ia1(jf),jf=1,ia1(2)+2),'ja1',(ja1(jf),jf=1,ja1(2)+2)
       print *,'ia0',(ia0(jf),jf=1,ia0(2)+2),'ja0',(ja0(jf),jf=1,ja0(2)+2)
       print *,'m1',(m1(jf),jf=1,m1(2)+2),'m2',(m2(jf),jf=1,m2(2)+2)
       print *,'iprod',(iprod(jf),jf=1,iprod(2)+2)
       print *,'isquar',(isquar(jf),jf=1,isquar(2)+2)
       print *,'mdarr',(mdarr(jf),jf=1,mdarr(2)+2)
       
       
       
       print *,'multiple',(norma(jf),jf=1,norma(2)+2)
       
       if (norma(2).eq.0)goto 72
       if ((norma(1).eq.1).and.(norma(2).eq.1))goto 61
       goto 62
61    if (norma(3).ne.1)goto 62
       if (ibsw1.ne.0)goto 63
       ibsw1=1
       do jf=1,m1(2)+2
       marr(jf)=m1(jf)
       mbarr(jf)=m1(jf)
       mpow1(1,jf)=m1(jf)
       npow1(1,jf)=m1(jf)
       end do
       do jf=1,m2(2)+2
       karr(jf)=m2(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       m2(jf)=kcarr(jf)
       end do
       if (ibsw3.eq.1)goto 190
       goto 36
63     a=a
! record to be written here with the ja's       
       goto 190
       stop


62     a=a
       do jf=1,norma(2)+2
       marr(jf)=norma(jf)
       end do
       do jf=1,ja5(2)+2
       mbarr(jf)=ja5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia5(2)+2
       kbarr(jf)=ia5(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia5(jf)=kcarr(jf)
       end do

       do jf=1,ja4(2)+2
       mbarr(jf)=ja4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia4(2)+2
       kbarr(jf)=ia4(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia4(jf)=kcarr(jf)
       end do

       do jf=1,ja3(2)+2
       mbarr(jf)=ja3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia3(2)+2
       kbarr(jf)=ia3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia3(jf)=kcarr(jf)
       end do








       do jf=1,ja2(2)+2
       mbarr(jf)=ja2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia2(2)+2
       kbarr(jf)=ia2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia2(jf)=kcarr(jf)
       end do

       do jf=1,ja1(2)+2
       mbarr(jf)=ja1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia1(2)+2
       kbarr(jf)=ia1(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia1(jf)=kcarr(jf)
       end do


       do jf=1,ja0(2)+2
       mbarr(jf)=ja0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ia0(2)+2
       kbarr(jf)=ia0(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ia0(jf)=kcarr(jf)
       end do
       goto 73
       stop
72     ibsw2=1
73     a=a
! swap ia's with ja's
       do jf=1,ia5(2)+2
       marr(jf)=ia5(jf)
       end do
       do jf=1,ja5(2)+2
       mbarr(jf)=ja5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       iprod(jf)=mcarr(jf)
       end do
       do jf=1,ia5(2)+2
       mbarr(jf)=ia5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       isquar(jf)=mcarr(jf)
       end do
       do jf=1,ia4(2)+2
       marr(jf)=ia4(jf)
       end do
! insert 5       
!       jinj=ia4(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ja4(2)+2
       mbarr(jf)=ja4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ia4(2)+2
       mbarr(jf)=ia4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do




       do jf=1,ia3(2)+2
       marr(jf)=ia3(jf)
       end do
! insert 6       
!       jinj=ia3(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ja3(2)+2
       mbarr(jf)=ja3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ia3(2)+2
       mbarr(jf)=ia3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do


       do jf=1,ia2(2)+2
       marr(jf)=ia2(jf)
       end do
! insert 7       
!       jinj=ia2(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ja2(2)+2
       mbarr(jf)=ja2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ia2(2)+2
       mbarr(jf)=ia2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do

       do jf=1,ia1(2)+2
       marr(jf)=ia1(jf)
       end do
! insert 8       
!       jinj=ia1(2)+2
!       marr(jinj+1)=0
!       marr(jinj+2)=0
!       marr(jinj+3)=0
!       marr(2)=marr(2)+3
       do jf=1,ja1(2)+2
       mbarr(jf)=ja1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ia1(2)+2
       mbarr(jf)=ia1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do

       do jf=1,ia0(2)+2
       marr(jf)=ia0(jf)
       end do
       do jf=1,ja0(2)+2
       mbarr(jf)=ja0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,iprod(2)+2
       kbarr(jf)=iprod(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iprod(jf)=kcarr(jf)
       end do
       do jf=1,ia0(2)+2
       mbarr(jf)=ia0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,isquar(2)+2
       kbarr(jf)=isquar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       isquar(jf)=kcarr(jf)
       end do
       do jf=1,iprod(2)+2
       marr(jf)=iprod(jf)
       end do
       marr(2)=marr(2)+1
       inj=marr(2)+2
       marr(inj)=0
       do jf=1,isquar(2)+2
       mbarr(jf)=isquar(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       karr(jf)=mdarr(jf)
       end do
       kbarr(1)=mdarr(1)
       kbarr(2)=1
       kbarr(3)=5000
       call mpadd(0)
       kcarr(2)=kcarr(2)-1
       do jf=2,kcarr(2)+2
       norma(jf)=kcarr(jf)
       end do
       norma(1)=mod(kcarr(1)+1,2)
       print *,'ia5',(ia5(jf),jf=1,ia5(2)+2),'ja5',(ja5(jf),jf=1,ja5(2)+2)
       print *,'ia4',(ia4(jf),jf=1,ia4(2)+2),'ja4',(ja4(jf),jf=1,ja4(2)+2)
       print *,'ia3',(ia3(jf),jf=1,ia3(2)+2),'ja3',(ja3(jf),jf=1,ja3(2)+2)
       print *,'ia2',(ia2(jf),jf=1,ia2(2)+2),'ja2',(ja2(jf),jf=1,ja2(2)+2)
       print *,'ia1',(ia1(jf),jf=1,ia1(2)+2),'ja1',(ja1(jf),jf=1,ja1(2)+2)
       print *,'ia0',(ia0(jf),jf=1,ia0(2)+2),'ja0',(ja0(jf),jf=1,ja0(2)+2)
       print *,'m1',(m1(jf),jf=1,m1(2)+2),'m2',(m2(jf),jf=1,m2(2)+2)
       
       
       
       
       
       print *,'multiple',(norma(jf),jf=1,norma(2)+2)
       if (norma(2).eq.0)goto 172
       if ((norma(1).eq.1).and.(norma(2).eq.1))goto 161
       goto 162
161     if (norma(3).ne.1)goto 162
       if (ibsw1.ne.0)goto 163
       ibsw1=1
       do jf=1,m1(2)+2
       marr(jf)=m1(jf)
       mbarr(jf)=m1(jf)
       mpow1(1,jf)=m1(jf)
       npow1(1,jf)=m1(jf)
       end do
       do jf=1,m2(2)+2
       karr(jf)=m2(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       m2(jf)=kcarr(jf)
       end do
       if (ibsw3.eq.1)goto 190
       goto 36
163     a=a
! record to be written here with the ia's       
       goto 190
       stop


162     a=a
       do jf=1,norma(2)+2
       marr(jf)=norma(jf)
       end do
       do jf=1,ia5(2)+2
       mbarr(jf)=ia5(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja5(2)+2
       kbarr(jf)=ja5(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja5(jf)=kcarr(jf)
       end do

       do jf=1,ia4(2)+2
       mbarr(jf)=ia4(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja4(2)+2
       kbarr(jf)=ja4(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja4(jf)=kcarr(jf)
       end do

       do jf=1,ia3(2)+2
       mbarr(jf)=ia3(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja3(2)+2
       kbarr(jf)=ja3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja3(jf)=kcarr(jf)
       end do








       do jf=1,ia2(2)+2
       mbarr(jf)=ia2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja2(2)+2
       kbarr(jf)=ja2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja2(jf)=kcarr(jf)
       end do

       do jf=1,ia1(2)+2
       mbarr(jf)=ia1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja1(2)+2
       kbarr(jf)=ja1(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja1(jf)=kcarr(jf)
       end do


       do jf=1,ia0(2)+2
       mbarr(jf)=ia0(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ja0(2)+2
       kbarr(jf)=ja0(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ja0(jf)=kcarr(jf)
       end do
       goto 173
       stop
172     ibsw2=ibsw2+1
173     itcon=itcon+1
       if (itcon.eq.4)goto 190
!      if (ibsw2.eq.2)goto 190
      if (ibsw2.eq.1)goto 59

      if (ibsw2.eq.0)goto 59
! polynomial chosen will be one with minimum sum of squares of coefficients


190   a=a
      do jf=1,ia5(2)+2
      marr(jf)=ia5(jf)
      mbarr(jf)=ia5(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      iasq(jf)=mcarr(jf)
      end do
     do jf=1,ia4(2)+2
     marr(jf)=ia4(jf)
     mbarr(jf)=ia4(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,iasq(2)+2
     kbarr(jf)=iasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     iasq(jf)=kcarr(jf)
     end do

     do jf=1,ia3(2)+2
     marr(jf)=ia3(jf)
     mbarr(jf)=ia3(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,iasq(2)+2
     kbarr(jf)=iasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     iasq(jf)=kcarr(jf)
     end do

     do jf=1,ia2(2)+2
     marr(jf)=ia2(jf)
     mbarr(jf)=ia2(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,iasq(2)+2
     kbarr(jf)=iasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     iasq(jf)=kcarr(jf)
     end do

     do jf=1,ia1(2)+2
     marr(jf)=ia1(jf)
     mbarr(jf)=ia1(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,iasq(2)+2
     kbarr(jf)=iasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     iasq(jf)=kcarr(jf)
     end do

     do jf=1,ia0(2)+2
     marr(jf)=ia0(jf)
     mbarr(jf)=ia0(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,iasq(2)+2
     kbarr(jf)=iasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     iasq(jf)=kcarr(jf)
     end do
! j's for i's in summing squares
      do jf=1,ja5(2)+2
      marr(jf)=ja5(jf)
      mbarr(jf)=ja5(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      jasq(jf)=mcarr(jf)
      end do
     do jf=1,ja4(2)+2
     marr(jf)=ja4(jf)
     mbarr(jf)=ja4(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,jasq(2)+2
     kbarr(jf)=jasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jasq(jf)=kcarr(jf)
     end do

     do jf=1,ja3(2)+2
     marr(jf)=ja3(jf)
     mbarr(jf)=ja3(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,jasq(2)+2
     kbarr(jf)=jasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jasq(jf)=kcarr(jf)
     end do

     do jf=1,ja2(2)+2
     marr(jf)=ja2(jf)
     mbarr(jf)=ja2(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,jasq(2)+2
     kbarr(jf)=jasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jasq(jf)=kcarr(jf)
     end do

     do jf=1,ja1(2)+2
     marr(jf)=ja1(jf)
     mbarr(jf)=ja1(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,jasq(2)+2
     kbarr(jf)=jasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jasq(jf)=kcarr(jf)
     end do

     do jf=1,ja0(2)+2
     marr(jf)=ja0(jf)
     mbarr(jf)=ja0(jf)
     end do
     call menmul
     do jf=1,mcarr(2)+2
     karr(jf)=mcarr(jf)
     end do
     do jf=1,jasq(2)+2
     kbarr(jf)=jasq(jf)
     end do
     call mpadd(0)
     do jf=1,kcarr(2)+2
     jasq(jf)=kcarr(jf)
     end do
     print *,'itcon',itcon
     print *,'sum of ia squares',(iasq(jf),jf=1,iasq(2)+2)
     print *,'sum of ja squares',(jasq(jf),jf=1,jasq(2)+2)
        
        ibsw1=0
        ibsw2=0
        ibind=1
        itcon=0
        if (iasq(2).lt.jasq(2))goto 80
        if (iasq(2).gt.jasq(2))goto 81
        do jf=3,iasq(2)+2
        if (iasq(jf).lt.jasq(jf))goto 80
        if (iasq(jf).gt.jasq(jf))goto 81
        end do

80      if (ibsw3.eq.1)goto 110
        ibsw3=1
        goto 83
81    if (ibsw3.eq.1)goto 112       
        ibsw3=1
        goto 60
110     print *,'min vector ia5',(ia5(jf),jf=1,ia5(2)+2),'ia4',& 
        (ia4(jf),jf=1,ia4(2)+2),'ia3',(ia3(jf),jf=1,ia3(2)+2),&
        'ia2',(ia2(jf),jf=1,ia2(2)+2),'ia1',(ia1(jf),jf=1,ia1(2)+2),&
        'ia0',(ia0(jf),jf=1,ia0(2)+2)
        do jf=1,ia5(2)+2
        ifinv(6,jf)=ia5(jf)
        end do
        do jf=1,ia4(2)+2
        ifinv(5,jf)=ia4(jf)
        end do
        do jf=1,ia3(2)+2
        ifinv(4,jf)=ia3(jf)
        end do
        do jf=1,ia2(2)+2
        ifinv(3,jf)=ia2(jf)
        end do
        do jf=1,ia1(2)+2
        ifinv(2,jf)=ia1(jf)
        end do
        do jf=1,ia0(2)+2
        ifinv(1,jf)=ia0(jf)
        end do
1002    do jf=1,ifinv(1,2)+2
        karb(jf)=ifinv(1,jf)
        end do
        karb(1)=0
        do i=2,6
        do jf=1,ifinv(i,2)+2
        kara(jf)=ifinv(i,jf)
        end do
        kara(1)=0
        call subgcd2
        if ((kard(2).eq.1).and.(kard(3).eq.1))goto 1001
        
        do jf=1,kard(2)+2
        karb(jf)=kard(jf)
        end do
        end do
        
        do jf=1,kard(2)+2
        mbarr(jf)=kard(jf)
        end do
        do i=1,6
        do jf=1,ifinv(i,2)+2
        marr(jf)=ifinv(i,jf)
        end do
        call mendiv
        do jf=1,mdarr(2)+2
        ifinv(i,jf)=mdarr(jf)
        end do
        end do
1001    print *,'gcd',(kard(jf),jf=1,kard(2)+2)

        goto 1000
112     print *,'min vector ja5',(ja5(jf),jf=1,ja5(2)+2),'ja4',&
        (ja4(jf),jf=1,ja4(2)+2),'ja3',(ja3(jf),jf=1,ja3(2)+2),&
        'ja2',(ja2(jf),jf=1,ja2(2)+2),'ja1',(ja1(jf),jf=1,ja1(2)+2),&
        'ja0',(ja0(jf),jf=1,ja0(2)+2)
        do jf=1,ja5(2)+2
        ifinv(6,jf)=ja5(jf)
        end do
        do jf=1,ja4(2)+2
        ifinv(5,jf)=ja4(jf)
        end do
        do jf=1,ja3(2)+2
        ifinv(4,jf)=ja3(jf)
        end do
        do jf=1,ja2(2)+2
        ifinv(3,jf)=ja2(jf)
        end do
        do jf=1,ja1(2)+2
        ifinv(2,jf)=ja1(jf)
        end do
        do jf=1,ja0(2)+2
        ifinv(1,jf)=ja0(jf)
        end do
        
        goto 1002
114     print *,'excellent polynomial found constant term=',&
        (ires(jf),jf=1,ires(2)+2),'m1',(m1(jf),jf=1,m1(2)+2),&
        'm2',(m2(jf),jf=1,m2(2)+2)
        do jf=1,10
        ia5(jf)=0
        ia4(jf)=0
        ia3(jf)=0
        ia2(jf)=0
        ia1(jf)=0
        ia0(jf)=0
        end do
        ia5(1)=0
        ia5(2)=1
        ia5(3)=1
        do jf=1,ires(2)+2
        ia0(jf)=ires(jf)
        end do











1000   end
       
       
       
       
       
       subroutine subgcd(idegs,idegb,idegg,ipd)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
       
       dimension iws(20),iwb(20),itempb(20)
       icon=0
       do i= 1,idegs +1
       iws(i) = isarray(i)
       end do
       do i=1,idegb+1
       iwb(i) = ibarray(i)
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       icon = icon+1
       ib=iws(1)
       call subbw6(ipd,ib,iv)

       mul = iv
       do i =1,loopl
       iqt= iwb(i) * mul
       iqt=mod(iqt,ipd)
       iwb(i) =0
       do j =2,iwsd +1
       iwb(i+j-1)=iwb(i+j-1)-iws(j) * iqt
       iwb(i+j-1) =mod(iwb(i+j-1),ipd)
       
       
       if(iwb(i+j-1).ge.0)goto 12
       iwb(i+j-1)=iwb(i+j-1) +ipd
12     end do
       
       
       end do
       
       
       

       do i= 1,iwbd +1

       if(iwb(i).ne.0)goto 20
       end do
       
       goto 100
20     itempbd=iwbd+2 -i
       do jj =1,itempbd
       itempb(jj) = iwb(i+jj-1)
       
       
       end do
       
       
       do i= 1,iwsd + 1
       iwb(i) =iws(i)
       end do
       iwbd = iwsd
       do jj=1,itempbd
       iws(jj) = itempb(jj)
       end do
       iwsd =itempbd-1
       
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8)
       if (icon .ne.4)goto 10
       
       goto 10
100    idegg=iwsd
       do i=1,10
       igarray(i)=0
       end do
       do i = 1,iwsd +1
       igarray(i) = iws(i)
       end do
       return
       end
       subroutine multy(idegm,idegn,ipd)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
       
       
       
       do i=1,idegm+idegn+1
       mult3(i) =0
       mult4(i)=0
       end do
        
       do i =1,idegm +1
       do j =1,idegn+1
       mult3(i+j-1)=mult3(i+j-1) +mult1(i) *mult2(j)
       mult3(i+j-1) =mod(mult3(i+j-1),ipd)
       if(mult3(i+j-1).ge.0)goto 10
       mult3(i+j-1)= mult3(i+j-1) +ipd
10     end do 
       end do
       return
       end
       subroutine subbw6(ia,ib,iv)
       
       
       ib =mod(ib,ia)
       
!       print *,ib
       
       if(ib.ge.0)goto 10
       ib =ia + ib
10     iu =1
!       print *,'iaib',ia,ib
       id = ia
       if(ib.eq.0)goto 888
       iv1=0
       iv3 =ib
815    if(iv3.eq.0)goto 830
       iqq = int(id/iv3)
       it3 =id -iqq*iv3
       it1 =iu -iqq*iv1
       iu =iv1
       id = iv3
       iv1 = it1
       iv3 = it3
       goto 815
830    iv =(id -ia *iu)/ib
       if(iu.le.0)goto 870
       iv =iv *(1-ia)
       iv =mod(iv,ia)
       if(iv.ge.0)goto 870
       iv=iv+ia
870    a=a
! 870    print *,'gcd=',id
       if (id.gt.1)goto 890
!       print *,'inverse=',iv
       goto 890
888    iv =ib
890    return
       end
       subroutine subbw4(idegs,idegb,idegr,ipd)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
       
       
      dimension iqt(100),iwb(100),iws(100)
      
      


      do i =1,idegs+1
      iws(i) =isarray(i)
      end do
      do i=1,idegb+1
      iwb(i) =ibarray(i)
      end do
      iwsd = idegs
      iwbd = idegb
      ib=iws(1)
      loopl = iwbd -iwsd +1
      
      call subbw6(ipd,iws(1),iv)
      
      
      
      
      mul = iv
      do i =1,loopl
      iqt(i) =iwb(i) *mul
      iqt(i)=mod(iqt(i),ipd)
      iwb(i) =0
      do j =2,iwsd + 1
      iwb(i+j-1) =iwb(i+j-1) -iws(j) *iqt(i)
      iwb(i+j-1) =mod(iwb(i+j-1),ipd)
      if(iwb(i+j-1).ge.0)goto 12
      iwb(i+j-1) =iwb(i+j-1) +ipd
12    end do
      end do
      do i =1,iwbd +1
      if(iwb(i).ne.0)goto 20
      end do
      goto 100
20    itempbd =iwbd +2 -i
      do jj =1,itempbd
      igarray(jj) =iwb(i+jj-1)
      end do
      idegr = itempbd -1
      goto 110
100   do i=1,10
      igarray(1) = 0
      end do
      idegr =0
110   a=a 
! 110   print *,'quos',iqt(1),iqt(2),iqt(3)
!      print *,'rem',igarray(1),igarray(2)
      return
      end
      


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
!      print *,'lenzy',ilen,ilen2
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
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      dimension kdum(5000),isub(5000)
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
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
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


      
      
      

      subroutine subgcd2
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
       do jf=3,kara(2)+2
       karr(jf-2)=kara(jf)
       end do
       ilen=kara(2)
       do jf=3,karb(2)+2
       kbarr(jf-2)=karb(jf)
       end do
       ilen2=karb(2)
1      call mpdiv(ilen,ilen2,irlen,icont,iswq)
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
2      kard(1)=0
       kard(2)=ilen2
       if (ilen2.eq.0)goto 3
       do jf=1,ilen2
       kard(jf+2)=kbarr(jf)
       end do
3      return
       end



      
      


      subroutine mpgcd
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      dimension karu(5000),karv1(5000),karv3(5000),karqq(5000)
      dimension kart3(5000),kart1(5000)
      
      
      
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
!      print *,'kart3',(kart3(i),i=1,kart3(2)+2)
      
      
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
!      print *,'firkarv',(karv(i),i=1,karv(2)+2)


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

      
      subroutine doubmul(idegm,idegn,idegp,idprod)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      dimension mmul(5000)
      goto 1
      

      
1     do i=1,12
      do j=1,10000
      kdoub3(i,j)=0
      end do
      end do
      do i=1,idegm+1
      do j=1,idegn+1
      ilen=kdoub1(i,2)
      ilen2=kdoub2(j,2)
      if((ilen.eq.0).or.(ilen2.eq.0))goto 10
      

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
      if (idegg.lt.idegp)goto 102
      
      do j=1,idegg-idegp+1
      do i2=1,kdoub3(j,2) +2
      mmul(i2)=kdoub3(j,i2)
      end do
      do i=1,idd
      do i2=1,kdoub3(j,2)+2
      kdoub3(j,i2)=0
      end do
      ilen=mmul(2)
      if (ilen.eq.0)goto 60
      do i2=3,mmul(2)+2
      karr(i2-2)=mmul(i2)
      end do
      ilen2=kmpol(i,2)
      if (ilen2.eq.0)goto 60
      
      do i2=3,kmpol(i,2)+2
      kbarr(i2-2)=kmpol(i,i2)
      end do
      
      call mpmul(ilen,ilen2,ilen3)
      
      do i2=1,ilen3
      kbarr(i2+2)=kcarr(i2)
      end do
      kbarr(2)=ilen3
      kbarr(1)=mod(mmul(1)+kmpol(i,1)+1,2)
      
      
      do i2=1,kdoub3(i+j,2)+2
      karr(i2)=kdoub3(i+j,i2)
      end do
      call mpadd(0)
      do i2=1,kcarr(2)+2
      kdoub3(i+j,i2)=kcarr(i2)
      end do
     
60    end do
      if (J.ne.50)goto 601
      do jl=1,idegg+1
      print *,'mid',(kdoub3(jl,jf),jf=1,kdoub3(jl,2)+2)
      end do
      stop
      
601   end do
      
      do i=1,idd
      do j=1,kdoub3(i+idprod,2)+2
      kdoub3(i,j)=kdoub3(i+idprod,j)
      end do
      end do
102   return 
      end



      subroutine modmul(idprod)
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      print *,'idprod',idprod
      do i=1,idprod+1
      ilen2=ipowar(2)
      do jf=1,ilen2
      kbarr(jf)=ipowar(jf+2)
      end do
      
      
      ilen=kdoub3(i,2)
      
      do j=3,kdoub3(i,2)+2
      karr(j-2)=kdoub3(i,j)
      end do
      
      call mpdiv(ilen,ilen2,irlen,icont,iswq)
      

      if (irlen.eq.0)goto 90
      
      do jf=1,irlen
      kdoub3(i,jf+2)=irrr(jf)
      end do
      kdoub3(i,2)=irlen
      if (kdoub3(i,1).eq.0)goto 100
      do jf=1,ipowar(2)+2
      karr(jf)=ipowar(jf)
      end do
      kbarr(1)=0
      do jf=2,kdoub3(i,2)+2
      kbarr(jf)=kdoub3(i,jf)
      end do
      call mpadd(1)
      do jf=1,kcarr(2)+2
      kdoub3(i,jf)=kcarr(jf)
      end do
      goto 100
90    kdoub3(i,1)=0    
      kdoub3(i,2)=0
100   end do
      return
      end
      

       
      
      
      subroutine menmul
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      
      




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
       common ibarray(20),isarray(20),igarray(20),inv(25)
       common mult1(20),mult2(20),mult3(20),mult4(20),mult5(20)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(5000),kbarr(5000),kcarr(10000),ipqt(5000),irrr(5000)
       common mnum(50)
       common ia5(10),ia4(10),ia3(10),ia2(10),ia1(10),ia0(10),m1(10),m2(10)
       common kara(5000),karb(5000),kard(5000),karp(5000),karv(5000)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,5000),kdoub2(6,5000),kdoub3(12,10000)
       common ipowar(5000)
       common marr(1000),mbarr(1000),mcarr(1000),mdarr(1000)
      
      
      
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





      

      
