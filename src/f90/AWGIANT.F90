       program awgiant
       common ibarray(2000),isarray(1000),igarray(1000),inv(25)
       common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
       common ipowar(400)
       dimension kdsol(5,400),jpowar(400),kbper(5,500),iroot(12,800)
       
       dimension mpol(6),ipol(6),idbarr(12,12),iper(6),ibper(6)
       dimension niv(1000),match(1000),litd(1000),litt(1000),n(100)
       dimension mm40(100),mm31(100),mm22(100),mm13(100),mm04(100)
       dimension norma(1000),nrem(100),ipre(1000),icurr(1000)
       dimension ires(100),icdh(100),mh2(100),mpow1(5,100),mpow2(5,100)
       dimension normar(100),kans(200),littr(20),iabp(50),iabpn(50)
       dimension ipr(20),modp(20),ity(50)
       
       open(unit=1,file='llrnel2',access='direct',form=&
       'formatted',recl=1000,status='old')
       open(unit=2,file='match4',access='direct',form=&
       'formatted',recl=2000,status='old')
       open(unit=5,file='kits',access='sequential')

       kmx=200
       mab=kmx
       read(1,1,rec=5)(niv(jk),jk=1,kmx)
1      format(200i5)       
       print *,'niv',(niv(jk),jk=1,kmx)
       do i=1,20
       modp(i)=1
       end do
       ipr(1)=109
       ipr(2)=137
       ipr(3)=139
       ipr(4)=149
       ipr(5)=167
       ijsgn=0
       iconty=0
       do i=1,kmx
       read(5, *)irecnn,kia,kib,icur,(littr(j1),j1=1,20),(normar(j2)&
       ,j2=1,20)
       read(5, *)irecnn,(iabp(jf),iabpn(jf),ity(jf),jf=1,icur)
       if (niv(irecnn).ne.1)goto 50
       iconty=iconty+1
       
       print *,'normar',(normar(jz),jz=1,normar(2)+2),'kia',kia,'kib',kib
       do iffy=1,icur
       print *,'iffy',iffy,'iabp',iabp(iffy),'iabpn',iabpn(iffy),'ity',&
       ity(iffy)
       end do
       
       ijsgn=ijsgn+normar(1)
              
       do j2=1,5
       ilen=normar(2)
       do jf=3,normar(2)+2
       karr(jf-2)=normar(jf)
       end do
       
       kbarr(1)=modp(j2)
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       kbarr(1)=ipr(j2)
       ilen2=1
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       modp(j2)=irrr(1)
       end do
       print *,'remz',(modp(jz),jz=1,5)
       if (iconty.ne.3)goto 50
       
50     end do
       print *,'rems',(modp(jf),jf=1,5)
       print *,'ijsgn',ijsgn
       ijsgn=mod(ijsgn,2)
       
       if (ijsgn.eq.0)goto 51
       do i=1,5
       modp(i)=ipr(i)-modp(i)
       end do
       
51     print *,(modp(jf),jf=1,5),'iconty',iconty
!       stop








       
       read(2,2,rec=1)(match(jk),jk=1,kmx*2)
2      format(400i5)       
       n(1)=0
       n(2)=2
       n(3)=3005
       n(4)=3021
       ia5(1)=0
       ia5(2)=1
       ia5(3)=1
       ia4(1)=0
       ia4(2)=1
       ia4(3)=1
       ia3(1)=0
       ia3(2)=1
       ia3(3)=16
       ia2(1)=0
       ia2(2)=1
       ia2(3)=24
       ia1(1)=0
       ia1(2)=1
       ia1(3)=20
       ia0(1)=0
       ia0(2)=1
       ia0(3)=9
       
       
       
       m1(1)=0
       m1(2)=1
       m1(3)=31
       m2(1)=0
       m2(2)=1
       m2(3)=1
       litd(1)=0
       litd(2)=1
       litd(3)=1
       iconty=0
       mab=200
       do i=1,kmx
       
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
       print *,karr(1),ilen,kbarr(1),ilen2
       
       call mpmul(ilen,ilen2,ilen3)
       litt(2)=ilen3
       do jf=1,ilen3
       litt(jf+2)=kcarr(jf)
       end do
       litt(1)=mod(m2(1)+isgn,2)
       ilen=1
       print *,'firlit',(litt(jf),jf=1,litt(2)+2)
       
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
       print *,'kia=',kia,'kib=',kib,'litt',(litt(jf),jf=1,litt(2)+2)
       
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
       
700    end do
       print *,'litd',(litd(jf),jf=1,litd(2)+2)
!      compute rational integer square root       
       
       do j=1,litd(2)+2
       ipre(j)=litd(j)
       end do
720    do j=3,litd(2)+2 
       karr(j-2)=litd(j)
       end do
       ilen=litd(2)
       do j=3,ipre(2)+2
       kbarr(j-2)=ipre(j)
       end do
       ilen2=ipre(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do j=1,icont
       kbarr(j+2)=ipqt(j)
       end do
       kbarr(2)=icont
       kbarr(1)=0
       
       
       do j=1,ipre(2)+2
       karr(j)=ipre(j)
       end do
       call mpadd(0)
       do j=3,kcarr(2)+2
       karr(j-2)=kcarr(j)
       end do
       ilen=kcarr(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do j=1,icont
       icurr(j+2)=ipqt(j)
       end do
       icurr(2)=icont
       icurr(1)=0
       if(icurr(2).lt.ipre(2))goto 730
       do j=3,ipre(2)+2
       if (icurr(j).lt.ipre(j))goto 730
       end do
       goto 732
730    do j=1,icurr(2)+2
       ipre(j)=icurr(j)
       end do
       go to 720
732    print *,'int. root',(icurr(jf),jf=1,icurr(2)+2)
       do j=3,icurr(2)+2
       karr(j-2)=icurr(j)
       kbarr(j-2)=icurr(j)
       end do
       ilen=icurr(2)
       ilen2=ilen
       call mpmul(ilen,ilen2,ilen3)
       print *,'chsq',(kcarr(jf),jf=1,ilen3)
       if(ilen3.ne.litd(2))goto 740
       do j=1,ilen3
       if (kcarr(j).ne.litd(j+2))goto 740
       end do
       goto 742
740    print *,'not perfect square'       
       stop
742    print *,'square ok'      
       do j=3,n(2)+2 
       kbarr(j-2)=n(j)
       end do
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 744
       do j=1,irlen
       ires(j+2)=irrr(j)
       end do
744    ires(1)=0
       ires(2)=irlen
       print *,'ires=',(Ires(jf),jf=1,ires(2)+2)
       close(unit=1)
       close(unit=2)
       
       

       open(unit=3,file='square4',access='direct',form=&
       'formatted',recl=2000,status='old')
       open(unit=4,file='minpol4',access='direct',form&
       ='formatted',recl=80,status='old')
       idprod=4
       do i=1,idprod+1
       read(3,3,rec=i)(kbper(i,jf),jf=1,500)
3      format(500i4)       
       end do
       do i=1,idprod+1
       read(4,4,rec=i)(kmpol(i,jf),jf=1,20)
4      format(20i4)       
       end do
       print *,'kmpol1',(kmpol(1,jf),jf=1,kmpol(1,2)+2)
       ipold=idprod+2                                        
       idegm=idprod
       idegn=idprod
       idd=idprod+1
!       ipd=257
       ipd=167
       call subbw6(ipd,22,iv)
       print *,'iv',iv
       
       
       kbarr(1)=ipd
       ilen2=1
       
       do i=1,idprod+1
       do j=3,kmpol(i,2)+2
       karr(j-2)=kmpol(i,j)
       end do
       ilen=kmpol(i,2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 6
       ipol(i+1)=irrr(1)
       if (kmpol(i,1).eq.0)goto 601
       ipol(i+1)=ipd-ipol(i+1)
       goto 601
6      ipol(i+1)=0
601    end do
       print *,(ipol(jf),jf=2,6)
       
       do i=1,idd
       ilen=kbper(i,2)
       do j=3,kbper(i,2)+2        
       karr(j-2)=kbper(i,j)
       end do
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 5
       mult1(i)=irrr(1)
       if (kbper(i,1).eq.0)goto 501
       mult1(i)=ipd-mult1(i)
       goto 501
5      mult1(i)=0
501    end do
!       goto 55
       print *,'mulz',(mult1(jf),jf=1,idd)
       
       goto 55
!      nnb  these instructions are temporary only
!      This is an application of The Chinese Remainder Theorem  
!      Mod is 167
       mult1(1)=146
       mult1(2)=128
       mult1(3)=62
       mult1(4)=5
       mult1(5)=108
       mult1(1)=74
       mult1(2)=26
       mult1(3)=27
       mult1(4)=87
       mult1(5)=103
!       goto 55
!      Second appliation of CRT       
       kbarr(1)=ipd
       ilen2=1
       do i=1,idd
       read(4,4,rec=i+5)(kans(jf),jf=1,20)
       
       print *,'kans',(kans(jz),jz=1,kans(2)+2)
       do jk=3,kans(2)+2
       karr(jk-2)=kans(jk)
       end do
       ilen=kans(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 53
       mult1(i)=irrr(1)
       if (kans(1).eq.0)goto 54
       mult1(i)=ipd-mult1(i)
       goto 54
53     mult1(i)=0
54     end do
!      These instructions must later be removed
       


       
       
55     ii =1
       
       do i =1,idd
       mat2(i,idd-1) = mult1(i)
       end do
       do i =1,idd
       mult2(i) =mult1(i)
       iper(i) =mult1(i)
       end do
10     call multy(idegm,idegn,ipd)
       idegg=idegm +idegn
       print *,'nmult3',(mult3(jk),jk=1,idegg+1)
       
       ii = ii +1
       if (ii.lt.idd )goto 12
       goto 14

12     do j =1,idegg-idegm
       mul=mult3(j)
       do i =1,idd
       
       mult3(j) =0
       mult3(i+j)=mult3(i+j) +mul *ipol(i+1) *(-1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       if (mult3(i+j).ge.0)goto 1001
       mult3(i+j)=mult3(i+j)+ipd
       
1001   end do
       print *,'try1',(mult3(jk),jk=1,idegg+1)
       end do
       
       do i = 1,idd
       
       mat2(i,idd -ii)=mult3(i+idegm)
       
       
       
       mult1(i)=mult3(i+idegm)
       
       end do
       goto 10
14     do i =1,idd -1
       mat2(i,idd)=0
       end do
       print *,'idegg',idegg,mult3(1),mult3(2),mult3(3),mult3(4),&
       mult3(5)
       
       mat2(i,idd) =1
       do i =1,ipold
       mpol(i) = -1 * ipol(i)
       end do
       do j = 1,idegg-idegm
       mul =mult3(j)
       print *,'mul',mul,mpol(1),mpol(2)
       do i =1,idd
       mult3(j) =0
       mult3(i+j) =mult3(i+j) +mul *mpol(i+1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       end do
       end do
       do i =1,idd
       irhs2(i) =mult3(i+idegg-idegm)
       end do
       do i =1,idd
       do j =1,idd
       mat2(i,j) =mod(mat2(i,j),ipd)
       if (mat2(i,j).ge.0)goto 20
       mat2(i,j) =mat2(i,j) +ipd
20     end do       
       irhs2(i)=mod(irhs2(i),ipd)
       if (irhs2(i).ge.0)goto 22
       irhs2(i) =irhs2(i) +ipd
22     end do
       do i=1,idd
       print *,'firmx',(mat2(i,jk),jk=1,idd),irhs2(i)
       end do
       

       
       call subbw3(idd,ipd)
       print *,'firsl',(isol(jk),jk=1,idd)
!       stop
       do i =1,20
       ip(I) =0
       end do
       
       do i =1,idd
       ip(2*i ) =isol(i) *(-1)
       end do
       idegd = idd*2
       call subbw2(idegd,ipd,idegr)
       print *,'mainidegg',idegr,'idegd',idegd
       print *,igarray(1),igarray(2),igarray(3),igarray(4)
       print *,'isars',isarray(1),isarray(2),isarray(3),isarray(4),&
       isarray(5),isarray(6),isarray(7)
       
       do i =1,10
       do j =1,10
       idbarr(i,j) =0
       end do
       end do
       do i= 1,idd
       mult2(i) =iper(i)
       end do
       do i =1,idegr+1
       idbarr(i,idd) =igarray(i)
       end do
       do lop =1,idegr -1
       do i=1,idd
       mult1(i) = idbarr(lop,i)
       end do
       call multy(idegm,idegn,ipd)
       print *,'mult3',mult3(1),mult3(2),mult3(3),mult3(4),mult3(5)
       
       do j =1,idegg-idegm
       mul=mult3(j)
       do i=1,idd
       mult3(j) =0
       mult3(i+j) =mult3(i+j) +mul *mpol(i+1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       end do
       end do
       do i =1,idd
       mult3(i) =mult3(i +idegm)
       end do
       do i =1,idd
       idbarr(lop +2,i) =idbarr(lop +2,i) +mult3(i)
       idbarr(lop+2,i) =mod(idbarr(lop+2,i),ipd)
       end do
       end do
       do lop =1,idegr +1
       print *,'darrays',(idbarr(lop,jf),jf=1,idd)
       end do
       
       do i=1,idd
       idbarr(idegr+1,i) =ipd -idbarr(idegr+1,i)
       end do
       
       
       do i=1,idd
       if (idbarr(idegr,i).ne.0)goto 30
       end do
       goto 999
       
       
       
30     do k =1,idegm
       do i =1,idd +idegm
       mult3(i) =0
       end do
       
       
       
       do i=1,idd
       mult3(i) =idbarr(idegr,i)
       end do
       do j=1,k
       mul =mult3(j)
       do i =1,idd
       mult3(j) =0
       mult3(i+j)=mult3(i+j) +mul*mpol(i+1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       if (mult3(i+j).ge.0)goto 3001
       mult3(i+j)=mult3(i+j)+ipd
3001   end do
       end do
       do i =1,idd
       mat2(i,idd-k) =mult3(i+k)
       end do
       end do
       do i=1,idd
       mat2(i,idd) =idbarr(idegr,i)
       mat2(i,idd)=mod(mat2(i,idd),ipd)
       if (mat2(i,idd).ge.0)goto 3002
       mat2(i,idd)=mat2(i,idd)+ipd
3002   end do
       do i=1,idd
       irhs2(i) =idbarr(idegr+1,i)
       irhs2(i)=mod(irhs2(i),ipd)
       if (irhs2(i).ge.0)goto 3003
       irhs2(i)=irhs2(i)+ipd
       
3003   end do
       
       do i=1,idd
       print *,(mat2(i,jf),jf=1,idd),irhs2(i)
       end do
       
       call subbw3(idd,ipd)
       
       
       
       print *,'solfir',(isol(jf),jf=1,idd),'ipd',ipd
       

!      check square root       
       print *,'idegm',idegm,'idegn',idegn
       do i=1,idd
       mult1(i)=isol(i)
       mult2(i)=isol(i)
       end do
       call multy(idegm,idegn,ipd)
       idegg =idegm+idegn
       do j=1,idegg-idegm
       mul=mult3(j)
       do i=1,idd
       mult3(j)=0
       mult3(i+j)=mult3(i+j)+mul*ipol(i+1)*(-1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       if (mult3(i+j).ge.0)goto 2991
       mult3(i+j)=mult3(i+j)+ipd
2991   end do
       end do
       do i=1,idd
       mult1(i)=mult3(i+idegm)
       end do
       print *,'ch',(mult1(jk),jk=1,idd)
       
       
       

!      compute inverse of square root mod p      
       do k=1,idegm
       do i=1,idd+idegm
       mult3(i)=0
       end do
       do i=1,idd
       mult3(i)=isol(i)
       end do
       do j=1,k
       mul=mult3(j)
       do i=1,idd
       mult3(j)=0
       mult3(i+j)=mult3(i+j)+mul*mpol(i+1)
       mult3(i+j)=mod(mult3(i+j),ipd)
       if (mult3(i+j).ge.0)goto 1991
       mult3(i+j)=mult3(i+j)+ipd
       
       
1991   end do
       end do
       do i=1,idd
       mat2(i,idd-k)=mult3(i+k)
       end do
       end do
       do i=1,idd
       mat2(i,idd)=isol(i)
       end do
       do i=1,idegm
       irhs2(i)=0
       end do
       irhs2(idd)=1
       do i=1,idd
       print *,(mat2(i,jk),jk=1,idd),irhs2(i)
       end do
       
       call subbw3(idd,ipd)
       do i=1,idd
       print *,'endmat',(mat1(i,j),j=1,idd),irhs1(i)
       end do
       
       print *,'sol2',(isol(jk),jk=1,idd)
       

!      perform lifts
       
199    do i=1,idd
       kdsol(i,3)=isol(i)
       
       kdsol(i,2)=1
       if(isol(i).ge.0)goto 200
       kdsol(i,1)=1
       goto 201
200    kdsol(i,1)=0       
201    end do       
       ipowar(1)=0
       ipowar(2)=1
       ipowar(3)=ipd
       idegm=idprod
       idegn=idprod
       do k=1,6
       
       ilen=ipowar(2)
       ilen2=ilen
       do jf=3,ilen+2
       karr(jf-2)=ipowar(jf)
       kbarr(jf-2)=ipowar(jf)
       end do
       call mpmul(ilen,ilen2,ilen3)
       ipowar(2)=ilen3
       do jf=1,ilen3
       ipowar(jf+2)=kcarr(jf)
       end do
       
       
301    do i=1,idd
       do jf=1,kdsol(i,2)+2
       kdoub1(i,jf)=kdsol(i,jf)
       kdoub2(i,jf)=kdsol(i,jf)
       end do
       end do
       idegp=idprod +1
       
       call doubmul(idegm,idegn,idegp,idprod)
       call modmul(idprod)
       do i=1,idd 
       print *,'now',(kdoub3(i,jk),jk=1,kdoub3(i,2)+2)
       end do
       



       do j=1,idegm+1
       do jf=1,kdoub3(j,2)+2
       kdoub1(j,jf)=kdoub3(j,jf)
       end do
       end do
       
       do i=1,idd
       do jf=1,kbper(i,2)+2
       kdoub2(i,jf)=kbper(i,jf)
       end do
       end do
       
       
       call doubmul(idegm,idegn,idegp,idprod)
       call modmul(idprod)
       
       do i=1,idd
       do jf=1,kdoub3(i,2)+2
       kdoub1(i,jf)=kdoub3(i,jf)
       end do
       
       print *,'ldoub',(kdoub1(i,jk),jk=1,kdoub1(i,2)+2)
       
998    end do
       do jf=1,kdoub1(idd,2)+2
       kbarr(jf)=kdoub1(idd,jf)
       end do
       karr(1)=0
       karr(2)=1
       karr(3)=3
       call mpadd(1)
       do jf=1,kcarr(2)+2
       kdoub1(idd,jf)=kcarr(jf)
       end do
       do i=1,idegm
       kdoub1(i,1)=mod(kdoub1(i,1)+1,2)
       end do
       
       do i=1,idd
       do jf=1,kdsol(i,2)+2
       kdoub2(i,jf)=kdsol(i,jf)
       end do
       end do
       
       
       call doubmul(idegm,idegn,idegp,idprod)
       call modmul(idprod)
       do i=1,idd
       do jf=1,kdoub3(i,2)+2
       kdoub2(i,jf)=kdoub3(i,jf)
       end do
       end do
       ilen=ipowar(2)
       do jf=1,ilen
       karr(jf)=ipowar(jf+2)
       end do
       ilen2=1
       kbarr(1)=2
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
       kdoub1(1,jf)=kcarr(jf)
       end do
       call doubmul(0,idegn,idegp,idprod)
       call modmul(idprod)
       do i=1,idd
       do j=1,kdoub3(i,2)+2
       kdsol(i,j)=kdoub3(i,j)
       
       end do
       print *,'solus',(kdsol(i,jf),jf=1,kdsol(i,2)+2)
       end do
       print *,'power',(ipowar(jf),jf=1,ipowar(2)+2)
       end do
       
!      end of lift procedure
!      start of check square root procedure 
       do i=1,idd
       do j=1,kdsol(i,2)+2
       kdoub1(i,j)=kdsol(i,j)
       end do
       do j=1,kbper(i,2)+2
       kdoub2(i,j)=kbper(i,j)
       end do
       end do
       call doubmul(idegm,idegn,idegp,idprod)
       call modmul(idprod)
       do i=1,idd
       print *,'checks',(kdoub3(i,jf),jf=1,kdoub3(i,2)+2)
       end do
       

       do i=1,idd
       do j=1,kdoub3(i,2)+2
       iroot(i,j)=kdoub3(i,j)
       kdoub1(i,j)=kdoub3(i,j)
       kdoub2(i,j)=kdoub3(i,j)
       end do
       end do
       
       call doubmul(idegm,idegn,idegp,idprod)
       call modmul(idprod)

       do i=1,idd
       print *,'square',(kdoub3(i,jk),jk=1,kdoub3(i,2)+2)
       end do
       
       do i=1,idd
       if (iroot(i,2).lt.100000)goto 1100
       do jf=1,ipowar(2) +2
       kbarr(jf)=ipowar(jf)
       end do
       do jf=1,iroot(i,2)+2
       karr(jf)=iroot(i,jf)
       end do
       call mpadd(1)
       
       do jf=1,kcarr(2)+2
       iroot(i,jf)=kcarr(jf)
       end do
       
       
1100   end do
       do i=1,idd
       print *,'roots',(iroot(i,jk),jk=1,iroot(i,2)+2)
       end do
       do jf=3,iroot(5,2)+2
       karr(jf-2)=iroot(5,jf)
       end do
       ilen=iroot(5,2)
       do jf=3,ipowar(2)+2
       kbarr(jf-2)=ipowar(jf)
       end do
       ilen2=ipowar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       print *,'redroot',(irrr(jf),jf=1,irlen)
       
       
       do i=1,idd
       do jf=1,iroot(i,2)+2
       kdoub1(i,jf)=iroot(i,jf)
       kdoub2(i,jf)=iroot(i,jf)
       end do
       end do
       call doubmul(idegm,idegn,idegp,idprod)
       do i=1,idd
       print *,'raised squares',(kdoub3(i,jf),jf=1,kdoub3(i,2)+2)
       end do
       
       do jf=1,m1(2)+2
       mpow1(1,jf)=m1(jf)
       end do
       do jf=1,m2(2)+2
       mpow2(1,jf)=m2(jf)
       end do
       do jf=3,m1(2)+2
       karr(jf-2)=m1(jf)
       kbarr(jf-2)=m1(jf)
       end do
       ilen=m1(2)
       ilen2=m1(2)
       do i=2,idd-1
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen
       karr(jf)=kcarr(jf)
       end do
       do jf =3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen2=n(2)
       
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1509
       do jf=1,irlen
       karr(jf)=irrr(jf)
       mpow1(i,jf+2)=irrr(jf)
       end do
       do jf=3,m1(2)+2
       kbarr(jf-2)=m1(jf)
       end do
       ilen2=m1(2)
       mpow1(i,2)=irlen
       mpow1(i,1)=mod(m1(1)*i,2)
       ilen=irlen
       end do
       goto 1510
1509   print *,'error in m1'
       stop
1510   do jf=3,m2(2)+2
       karr(jf-2)=m2(jf)
       kbarr(jf-2)=m2(jf)
       end do
       ilen=m2(2)
       ilen2=m2(2)
       do i=2,idd-1
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1519 
       do jf=1,irlen
       karr(jf)=irrr(jf)
       mpow2(i,jf+2)=irrr(jf)
       end do
       ilen=irlen
       mpow2(i,2)=irlen
       mpow2(i,1)=mod(m2(1)*i,2)
       ilen2=m2(2)
       do jf=3,ilen2+2
       kbarr(jf-2)=m2(jf)
       end do
       
       end do
       goto 1520
1519   print *,'error in m2'
       stop
1520   print *,'m14',(mpow1(4,jf),jf=1,mpow1(4,2)+2)
       print *,'m23',(mpow2(3,jf),jf=1,mpow2(3,2)+2)
       do jf=1,100
       normar(jf)=0
       end do
       
       
       if (ia5(2).eq.0)goto 1522
       do jf=3,ia5(2)+2
       karr(jf-2)=ia5(jf)
       end do
       ilen=ia5(2)
       do jf=3,mpow1(4,2)+2
       kbarr(jf-2)=mpow1(4,jf)
       end do
       ilen2=mpow1(4,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       kbarr(1)=5
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       normar(jf+2)=kcarr(jf)
       end do
       normar(2)=ilen3
       normar(1)=mod(ia5(1)+mpow1(4,1),2)
       print *,'normar1',(normar(jf),jf=1,normar(2)+2)

1522   if (ia4(2).eq.0)goto 1524
       do jf=3,ia4(2)+2
       karr(jf-2)=ia4(jf)
       end do
       ilen=ia4(2)
       do jf=3,mpow1(3,2)+2
       kbarr(jf-2)=mpow1(3,jf)
       end do
       ilen2=mpow1(3,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       do jf=3,mpow2(1,2)+2
       kbarr(jf-2)=mpow2(1,jf)
       end do
       ilen2=mpow2(1,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       kbarr(1)=4
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=mod(ia4(1)+mpow1(3,1)+mpow2(1,1),2)
       do jf=1,normar(2)+2
       karr(jf)=normar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       normar(jf)=kcarr(jf)
       end do
       print *,'normar2',(normar(jf),jf=1,normar(2)+2)
1524   if (ia3(2).eq.0)goto 1526
       do jf=3,ia3(2)+2
       karr(jf-2)=ia3(jf)
       end do
       ilen=ia3(2)
       do jf=3,mpow1(2,2)+2
       kbarr(jf-2)=mpow1(2,jf)
       end do
       ilen2=mpow1(2,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       do jf=3,mpow2(2,2)+2
       kbarr(jf-2)=mpow2(2,jf)
       end do
       ilen2=mpow2(2,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3 
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       kbarr(1)=3
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=mod(ia3(1)+mpow1(2,1) +mpow2(2,1),2)
       do jf=1,normar(2)+2
       karr(jf)=normar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       normar(jf)=kcarr(jf)
       end do
       print *,'normar3',(normar(jf),jf=1,normar(2)+2)
1526   if (ia2(2).eq.0)goto 1528
       do jf=3,ia2(2)+2
       karr(jf-2)=ia2(jf)
       end do
       ilen=ia2(2)
       do jf=3,mpow1(1,2)+2
       kbarr(jf-2)=mpow1(1,jf)
       end do
       ilen2=mpow1(1,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       do jf=3,mpow2(3,2)+2
       kbarr(jf-2)=mpow2(3,jf)
       end do
       ilen2=mpow2(3,2)
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       kbarr(1)=2
       ilen2=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=mod(ia2(1)+mpow1(1,1) +mpow2(3,1),2)
       do jf=1,normar(2)+2
       karr(jf)=normar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       normar(jf)=kcarr(jf)
       end do
       print *,'normar4',(normar(jf),jf=1,normar(2)+2)

1528   if (ia1(2).eq.0)goto 1530
       do jf=3,ia1(2)+2
       karr(jf-2)=ia1(jf)
       end do
       ilen=ia1(2)
       do jf=3,mpow2(4,2)+2
       kbarr(jf-2)=mpow2(4,jf)
       end do
       ilen2=mpow2(4,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=mod(ia1(1)+mpow2(4,1),2)
       do jf=1,normar(2)+2
       karr(jf)=normar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       normar(jf)=kcarr(jf)
       end do
       print *,'normar5',(normar(jf),jf=1,normar(2)+2)

1530   do jf=3,normar(2)+2
       karr(jf-2)=normar(jf)
       end do
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen=normar(2)
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1532
       do jf=1,irlen
       normar(jf+2)=irrr(jf)
       end do
       normar(2)=irlen
       goto 1533
1532   print *,'irregular derivative'
       stop
1533   print *,'derivative=',(normar(jf),jf=1,normar(2)+2)
       
       do jf=3,mpow1(3,2)+2
       karr(jf-2)=mpow1(3,jf)
       end do
       ilen=mpow1(3,2)
       do jf=3,mpow2(1,2)+2
       kbarr(jf-2)=mpow2(1,jf)
       end do
       ilen2=mpow2(1,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       mm31(jf+2)=kcarr(jf)
       end do
       mm31(1)=mod(mpow1(3,1)+mpow2(1,1),2)
       mm31(2)=ilen3
       do jf=3,mpow1(2,2)+2
       karr(jf-2)=mpow1(2,jf)
       end do
       ilen=mpow1(2,2)
       do jf=3,mpow2(2,2)+2
       kbarr(jf-2)=mpow2(2,jf)
       end do
       ilen2=mpow2(2,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       mm22(jf+2)=kcarr(jf)
       end do
       mm22(2)=ilen3
       mm22(1)=mod(mpow1(2,1)+mpow2(2,1),2)
       do jf=3,mpow1(1,2)+2
       karr(jf-2)=mpow1(1,jf)
       end do
       ilen=mpow1(1,2)
       do jf=3,mpow2(3,2)+2
       kbarr(jf-2)=mpow2(3,jf)
       end do
       ilen2=mpow2(3,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       mm13(jf+2)=kcarr(jf)
       end do
       mm13(2)=ilen3
       mm13(1)=mod(mpow1(1,1)+mpow2(3,1),2)
       print *,'mm13fir=',(mm13(jf),jf=1,mm13(2)+2)
       
       do jf=3,mpow1(4,2)+2
       karr(jf-2)=mpow1(4,jf)
       end do
       ilen=mpow1(4,2) 
       do jf=3,ia5(2)+2
       kbarr(jf-2)=ia5(jf)
       end do
       ilen2=ia5(2)
       do i=1,idd-1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       end do
       do jf=1,ilen
       mm40(jf+2)=karr(jf)
       end do
       mm40(2)=ilen
       mm40(1)=mod(ia5(1)*4+mpow1(4,1)*4,2)
       print *,'mm40=',(mm40(jf),jf=1,mm40(2)+2)
       

       do jf=3,mm31(2)+2
       karr(jf-2)=mm31(jf)
       end do
       ilen=mm31(2)
       do i=1,idd-2
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       end do
       do jf=1,ilen
       mm31(jf+2)=karr(jf)
       end do
       mm31(2)=ilen
       mm31(1)=mod(ia5(1)*3+mm31(1),2)
       print *,'mm31=',(mm31(jf),jf=1,mm31(2)+2)
       
       do jf=3,mm22(2)+2
       karr(jf-2)=mm22(jf)
       end do
       ilen=mm22(2)
       do i=1,idd-3
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       end do
       do jf=1,ilen
       mm22(jf+2)=karr(jf)
       end do
       mm22(1)=mod(ia5(1)*2+mm22(1),2)
       mm22(2)=ilen
       print *,'mm22',(mm22(jf),jf=1,mm22(2)+2)
       do jf=3,mm13(2)+2
       karr(jf-2)=mm13(jf)
       end do
       ilen=mm13(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       mm13(jf+2)=kcarr(jf)
       end do
       mm13(1)=mod(ia5(1)+mm13(1),2)
       mm13(2)=ilen3
       print *,'mm13=',(mm13(jf),jf=1,mm13(2)+2)
       do jf=1,mpow2(4,2)+2
       mm04(jf)=mpow2(4,jf)
       end do
       print *,'mm04=',(mm04(jf),jf=1,mm04(2)+2)
       
       

       
       do jf=3,iroot(1,2)+2
       karr(jf-2)=iroot(1,jf)
       end do
       ilen=iroot(1,2)
       ilen2=mm40(2)
       do jf=3,mm40(2)+2
       kbarr(jf-2)=mm40(jf)
       end do
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       norma(jf+2)=kcarr(jf)
       end do
       norma(2)=ilen3
       norma(1)=mod(iroot(1,1)+mm40(1),2)
       print *,'norma1',(norma(jf),jf=1,norma(2)+2)
       
       ilen=iroot(2,2)
       do jf=1,ilen
       karr(jf)=iroot(2,jf+2)
       end do
       ilen2=mm31(2)
       do jf=1,ilen2
       kbarr(jf)=mm31(jf+2)
       end do
       call mpmul(ilen,ilen2,ilen3)
       kbarr(2)=ilen3
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=mod(iroot(2,1)+mm31(1),2)
       do jf=1,norma(2)+2
       karr(jf)=norma(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       norma(jf)=kcarr(jf)
       end do
       print *,'norma2',(norma(jf),jf=1,norma(2)+2)
       
       ilen=iroot(3,2)
       do jf=1,ilen
       karr(jf)=iroot(3,jf+2)
       end do
       ilen2=mm22(2)
       do jf=1,ilen2
       kbarr(jf)=mm22(jf+2)
       end do
       call mpmul(ilen,ilen2,ilen3)
       kbarr(2)=ilen3
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=mod(iroot(3,1)+mm22(1),2)
       do jf=1,norma(2)+2
       karr(jf)=norma(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2 
       norma(jf)=kcarr(jf)
       end do
       print *,'norma3',(norma(jf),jf=1,norma(2)+2)
       
       ilen=iroot(4,2)
       do jf=1,ilen
       karr(jf)=iroot(4,jf+2)
       end do
       ilen2=mm13(2)
       do jf=1,ilen2
       kbarr(jf)=mm13(jf+2)
       end do
       call mpmul(ilen,ilen2,ilen3)
       kbarr(2)=ilen3
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=mod(iroot(4,1)+mm13(1),2)
       do jf=1,norma(2)+2
       karr(jf)=norma(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       norma(jf)=kcarr(jf)
       end do
       print *,'norma4',(norma(jf),jf=1,norma(2)+2)
       
       ilen=iroot(5,2)
       do jf=1,ilen
       karr(jf)=iroot(5,jf+2)
       end do
       ilen2=mm04(2) 
       print *,'lens',ilen,ilen2,(mm04(jf),jf=1,mm04(2)+2)
       
       do jf=1,ilen2
       kbarr(jf)=mm04(jf+2)
       end do
       print *,'firkbar',(kbarr(jf),jf=1,ilen2)
       
       call mpmul(ilen,ilen2,ilen3)
       kbarr(2)=ilen3
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=mod(iroot(5,1)+mm04(1),2)
       print *,'seckbar',(kbarr(jf),jf=1,ilen3)
       do jf=1,norma(2)+2
       karr(jf)=norma(jf)
       end do

       call mpadd(0)
       do jf=1,kcarr(2)+2
       norma(jf)=kcarr(jf)
       end do
       print *,'norma5',(norma(jf),jf=1,norma(2)+2)
       print *,'ipowar',(ipowar(jf),jf=1,ipowar(2)+2)
       do jf=3,norma(2)+2
       karr(jf-2)=norma(jf)
       end do
       ilen=norma(2)
       do jf=3,ipowar(2)+2
       kbarr(jf-2)=ipowar(jf)
       end do
       ilen2=ipowar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 6000
       print *,'remhalf',(irrr(jf),jf=1,irlen)
       if (norma(1).eq.0)goto 6002
       do jf=1,irlen
       karr(jf+2)=irrr(jf)
       end do
       karr(1)=1
       karr(2)=irlen
       do jf=1,ipowar(2)+2
       kbarr(jf)=ipowar(jf)
       end do
       call mpadd(0)
       print *,'corrected',(kcarr(jf),jf=1,kcarr(2)+2)
       


6002   a=a
6000   print *,'irregular 1'
       a=a


       
       
       
       
       
       

       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       ilen2=n(2)
       do jf=1,ilen2
       kbarr(jf)=n(jf+2)
       end do
       print *,'karrs',(karr(jf),jf=1,ilen)
       print *,'kbarrs',(kbarr(jf),jf=1,ilen2)
       
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 4000
       do jf=1,irlen
       nrem(jf+2)=irrr(jf)
       end do
       nrem(2)=irlen
       print *,'nremfir',(irrr(jf),jf=1,irlen)
       if (norma(1).eq.0)goto 4002
       
       do jf=1,irlen
       karr(jf+2)=nrem(jf+2)
       end do
       karr(2)=nrem(2)
       karr(1)=1
       kbarr(1)=0
       do jf=2,n(2)+2
       kbarr(jf)=n(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       nrem(jf)=kcarr(jf)
       end do
       goto 4003
4000   print *,'zero warning'
       stop
4002   nrem(1)=0
4003   print *,'n remainders',(nrem(jf),jf=1,nrem(2)+2)
       print *,'ires',(ires(jk),jk=1,ires(2)+2)
       
       icc=iconty/2
       ie=1
5000   if (ie.gt.icc)goto 5001
       ie=ie*2
       go to 5000
5001   ie=ie/2       
       nn=icc-ie
       do jf=3,m2(2)+2
       karr(jf-2)=m2(jf)
       kbarr(jf-2)=m2(jf)
       end do
       ilen=m2(2)
       ilen2=m2(2)
       
5002   if (ie.eq.1)goto 5030       
       ie=ie/2
       call mpmul(ilen,ilen2,ilen3)
       print *,'p5002',nn,ie
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 5012
       do jf=1,irlen
       karr(jf)=irrr(jf)
       kbarr(jf)=irrr(jf)
       end do

       do jf=1,irlen
       mh2(jf+2)=irrr(jf)
       end do
       ilen=irlen
       ilen2=irlen
       
       if (nn.ge.ie)goto 5020
       goto 5002
5020   do jf=3,m2(2)+2       
       kbarr(jf-2)=m2(jf)
       
       end do
       ilen2=m2(2)
       call mpmul(ilen,ilen2,ilen3)
       print *,'p5020',nn,ie
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 5012
       do jf=1,irlen
       karr(jf)=irrr(jf)
       kbarr(jf)=irrr(jf)
       mh2(jf+2)=irrr(jf)
       end do
       ilen=irlen
       ilen2=irlen
       nn=nn-ie
       goto 5002
5012   print *,'error in m2'       
       stop
5030   mh2(2)=ilen
       mh2(1)=mod(m2(1)*nn,2)
       print *,'firmh2',(mh2(jf),jf=1,mh2(2)+2)
       
       
       
       nn=icc+idd-2
       
       ie=1
5100   if (ie.gt.nn)goto 5101
       ie=ie*2
       goto 5100
5101   do jf=3,ia5(2)+2
       karr(jf-2)=ia5(jf)
       kbarr(jf-2)=ia5(jf)
       end do
       ie=ie/2
       print *,ie
       nn=nn-ie
       ilen=ia5(2)
       ilen2=ia5(2)
5102   if (ie.eq.1)goto 5130
       ie=ie/2
       print *,'5102',(karr(jf),jf=1,ilen)
       
       call mpmul(ilen,ilen2,ilen3)
       ilen=ilen3
       do jf=1,ilen
       karr(jf)=kcarr(jf)
       end do
       ilen2=n(2)
       print *,'sec5102',(karr(jf),jf=1,ilen)
       
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 5112
       do jf=1,irlen
       karr(jf)=irrr(jf)
       kbarr(jf)=irrr(jf)
       icdh(jf+2)=irrr(jf)
       end do
       ilen=irlen
       ilen2=irlen
       if(nn.ge.ie)goto 5120
       goto 5102
5120   do jf=3,ia5(2)+2
       kbarr(jf-2)=ia5(jf)
       end do
       ilen2=ia5(2)
       print *,'5120',(karr(jf),jf=1,ilen)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       ilen2=n(2)
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 5112
       do jf=1,irlen
       karr(jf)=irrr(jf)
       kbarr(jf)=irrr(jf)
       icdh(jf+2)=irrr(jf)
       end do
       ilen=irlen
       ilen2=irlen
       nn=nn-ie
       goto 5102
5112   print *,'error in ia5'       
       stop
5130   icdh(2)=ilen       
       icdh(1)=mod(ia5(1)*(icc+3),2)
       print *,'mh2',(mh2(jf),jf=1,mh2(2)+2)
       print *,'icdh',(icdh(jf),jf=1,icdh(2)+2)
       
       
       
       
       do jf=3,icdh(2)+2
       karr(jf-2)=icdh(jf)
       end do
       ilen=icdh(2)
       do jf=3,normar(2)+2
       kbarr(jf-2)=normar(jf)
       end do
       ilen2=normar(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ires(2)+2
       kbarr(jf-2)=ires(jf)
       end do
       ilen2=ires(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kans(jf+2)=kcarr(jf)
       end do
       kans(2)=ilen3
       kans(1)=mod(icdh(1)+normar(1)+ires(1),2)
       do jf=3,mh2(2)+2
       karr(jf-2)=mh2(jf)
       end do
       ilen=mh2(2)
       do jf=3,nrem(2)+2
       kbarr(jf-2)=nrem(jf)
       end do
       ilen2=nrem(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=mod(mh2(1)+nrem(1),2)
       do jf=1,kans(2)+2
       karr(jf)=kans(jf)
       end do
       call mpadd(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,n(2)+2
       kbarr(jf-2)=n(jf)
       end do
       ilen2=n(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       kara(jf+2)=irrr(jf)
       end do
       kara(1)=0
       kara(2)=irlen
       do jf=1,n(2)+2
       karb(jf)=n(jf)
       end do
       call mpgcd
       print *,'gcd=',(kard(jf),jf=1,kard(2)+2)
       stop








       stop
       do i=4,5
       print *,'original number'
       print *,(kbper(i,jk),jk=1,kbper(i,2)+2)
       end do
       
       ipowar(1)=0
       ipowar(2)=2
       ipowar(3)=1
       ipowar(4)=1881
       do i=1,idd 
       do j=1,kbper(i,2)+2
       kdoub3(i,j)=kbper(i,j)
       end do
       end do
       call modmul(idprod)
       do i=1,idd
       print *,'ch2',(kdoub3(i,jk),jk=1,kdoub3(i,2)+2)
       end do
       do i=1,idd
       do j=1,kmpol(i,2)+2
       kdoub3(i,j)=kmpol(i,j)
       end do
       end do
       call modmul(idprod)
       do i=1,idd
       print *,'modreds',(kdoub3(i,jk),jk=1,kdoub3(i,2)+2)
       end do


       
       goto 1000
999    print *,'no square root exists mod',ipd       
1000   end
       
       
       
       
       subroutine subbw2(idegd,ipd,idegr)
        
       
       common  ibarray(2000),isarray(1000),igarray(1000),inv(25)
       common  mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,400),kdoub2(6,400),kdoub(12,800)
       common ipowar(400)
       
       
       
       dimension indic(20),iq(20,20),iqi(20,20),iv(20,20)
       dimension irowp(20),irown(20),ic(20)
       dimension iran(20)
       nmul =1
       ngcd=0
       do i = 1,idegd
       if (ip(i).ge.0)goto 1
       ip(i) = ip(i) + ipd
1      end do
       print *,'ips',(ip(jk),jk=1,idegd)
       
       

       do i =1,20
       indic(i) =0
       end do
       
       nn=(ipd -1)/2
       print *,'nn',nn
       ie =1
2      if (ie.gt.nn)goto 3

       ie = ie *2
       goto 2
3      ie = ie/2
       print *,'ie',ie
       nn = nn-ie
       
       inv(17)=19
       inv(18)=9
       inv(19)=17
       inv(20)=15
       inv(21)=11
       inv(22)=22




       
       
       do i=1,20
       do j=1,20
       iq(i,j) = 0
       end do
       end do
       do i =1,20
       irowp(i) =0
       end do
       irowp(2) =1
       irown(1) =0
       do k = 2,idegd
       do j =1,ipd
       do i =1,idegd
       irown(i+1) =irowp(i) - ip(idegd+1-i) *irowp(idegd+1)
       end do
       do i=1,idegd +1
       irowp(i) =irown(i)
       irowp(i)=mod(irowp(i),ipd) 
       if (irowp(i).ge.0)goto 4
       irowp(i)=irowp(i)+ipd
4      end do
       end do
       
       do i=1,idegd
       iq(k,i) =irowp(i+1)
       end do
       end do
       iq(1,1) =1
       do k = 1,idegd
       print *,iq(k,1),iq(k,2),iq(k,3),iq(k,4),iq(k,5),iq(k,6),iq(k,7),iq(k,8)
       end do
       
       do i=1,idegd
       do j =1,idegd
       iqi(i,j)=iq(i,j)
       end do
       end do
       do i=1,idegd
       iqi(i,i)=iqi(i,i)-1
       if(iqi(i,i).ge.0)goto 12
       iqi(i,i) =iqi(i,i) +ipd 
12     end do       
       ir = 0
       do i=1,20
       
       ic(i)=0
       end do
       do k = 1,idegd
       do j =1,idegd
       if((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
34     ib=iqi(k,j)
       call subbw6(ipd,ib,iggs)
       mul =iggs       
       do i =k,idegd
       iqi(i,j) =(-1) * iqi(i,j)  * mul
       iqi(i,j)=mod(iqi(i,j),ipd) 
       if(iqi(i,j).ge.0)goto 33
       iqi(i,j) =iqi(i,j) +ipd
33     end do
       do jj=1,idegd
       if(jj.eq.j)goto 38
        
       mul =iqi(k,jj)
       do i=k,idegd
       itt = iqi(i,jj) + iqi(i,j) *mul
       iqi(i,jj) =mod(itt,ipd)
       if(iqi(i,jj).ge.0)goto 36
       iqi(i,jj)=iqi(i,jj) +ipd
36     end do
38     end do  
       
       ic(j) =k   
       indic(k) =j
       goto 50

40     ir =ir +1
       do is =1,idegd
       iv(is,ir) = 0
       end do
       iv(k,ir) =1
       do is =1,idegd
       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is,ir) =iqi(k,jv)
42     end do
50     end do
       do i=1,idegd
       print *,'iqis',iqi(i,1),iqi(i,2),iqi(i,3),iqi(i,4),iqi(i,5),&
       iqi(i,6),iqi(i,7),iqi(i,8)
       end do
       do jr =1,ir
       print *,'ivs',iv(1,jr),iv(2,jr),iv(3,jr),iv(4,jr),iv(5,jr),&
       iv(6,jr),iv(7,jr),iv(8,jr)
       print *,'ir',ir
       
       end do
       

58     if (ngcd.eq.300)goto 399      
       ibarray(1) =1
       do i =1,idegd
       ibarray(i+1) =ip(i)
       end do
       do i =1,idegd
       isarray(i) =0
       end do
       
       idegb = idegd
       do k = 1,ir
       iran(k) = mod((797 *nmul),131072)
       nmul =iran(k)
       end do
       do k=1,ir
       iran(k) =mod(iran(k),ipd)
       end do


       
       do k =1,ir
       do i =1,idegd
       
       isarray(i)=isarray(i)+iv(idegd +1-i,k) * iran(k)
       end do
       
       end do
       do i =1,idegd
       isarray(i) =mod(isarray(i),ipd)
       if (isarray(i).ge.0)goto 60
       isarray(i) =isarray(i) +ipd
60     end do
       print *,'izar',isarray(1),isarray(2),isarray(3),isarray(4),&
       isarray(5),isarray(6)
       
       idegs =idegb-1
       do i =1,idegd
       if (isarray(i).ne.0)goto 229
       idegs = idegs-1
       end do
       
       
       
       
       
       
       
229    do i = 1,idegs+1
       mult1(i)=isarray(i)
       mult5(i)=mult1(i)
       end do
       
       idegm=idegs
       
       
230    if(ie.eq.1)goto 250       
       ie =ie/2
       do i=1,idegm +1
       mult2(i) =mult1(i)
       end do
       idegn = idegm 
       call multy(idegm,idegn,ipd)
       idegm =idegm +idegn
       do i = 1,idegm +1
       mult1(i) =mult3(i)
       end do


       
       if(nn.ge.ie)goto 240
       goto 230
240    nn = nn-ie   
       do i = 1,idegs
       mult2(i) = mult5(i)
       end do
       
       idegn =idegs
       call multy(idegm,idegn,ipd)
       idegm =idegm +idegn
       do i = 1,idegm +1
       mult1(i)=mult3(i)
       end do
       goto 230
       
250    do i =1,idegb+1   
       isarray(i) =ibarray(i)
       end do
       idegs = idegb
       print *,'idegs',idegs
       print *,'isarray',isarray(1),isarray(2),isarray(3),isarray(4),&
       isarray(5),isarray(6),isarray(7) 
       
       do i=1,idegm+1
       ibarray(i) =mult1(i)
       end do
       idegb =idegm
       ibarray(idegb+1) =ibarray(idegb +1) - 1
       if(ibarray(idegb+1).ge.0)goto 255
       ibarray(idegb +1) =ibarray(idegb+1) +ipd
255    call subgcd(idegs,idegb,idegg,ipd)
       print *,'idegg',idegg,idegs,idegb,isarray(1),ibarray(1)
       print *,igarray(1),igarray(2),igarray(3),igarray(4),&
       igarray(5),igarray(6),igarray(7),igarray(8),igarray(9)
       ngcd=ngcd+1
       
       if (idegg.ne.idegd/2)goto 58
       
       goto 300
399    stop       
       
       
       
       
       
       
       
       
       
       
       call subgcd(6,8,idegg,ipd)
       print *,igarray(1),igarray(2),igarray(3),igarray(4),igarray(5),&
       igarray(6),igarray(7),igarray(8),igarray(9)
       print *,'idegg',idegg
300    idegr =idegg       
       return
       end
       subroutine subgcd(idegs,idegb,idegg,ipd)
       
       common ibarray(2000),isarray(1000),igarray(1000),inv(25)
       common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
       common ipowar(400)
       dimension iws(1000),iwb(2000),itempb(1000)
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
       iwb(i) =0
       do j =2,iwsd +1
       iwb(i+j-1)=iwb(i+j-1)-iws(j) * iqt
       iwb(i+j-1) =mod(iwb(i+j-1),ipd)
       print *,'imodd',iwb(i+j-1)
       
       if(iwb(i+j-1).ge.0)goto 12
       iwb(i+j-1)=iwb(i+j-1) +ipd
12     end do
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8),iwb(9)
       
       end do
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8),iwb(9)
       
       

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
       common ibarray(2000),isarray(1000),igarray(1000),inv(25)
       common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common kara(50),kbarb(50),kard(50),karp(50),karv(50)
       common iarq(2)
       common kmpol(6,20)
       common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
       common ipowar(400)
       
       
       do i=1,2000
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
       
       print *,ib
       
       if(ib.ge.0)goto 10
       ib =ia + ib
10     iu =1
       print *,'iaib',ia,ib
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
870    print *,'gcd=',id
       if (id.gt.1)goto 890
       print *,'inverse=',iv
       goto 890
888    iv =ib
890    return
       end
       subroutine subbw4
       dimension ibarray(20),isarray(10),irarray(10),inv(25)
      dimension iqt(10),iwb(10),iws(10)
      ipd =23
      inv(1) =1
      



      ibarray(1) =1
      ibarray(2) =3
      ibarray(3)=3
      ibarray(4)=7
      idegb =3
      isarray(1) =1
      isarray(2) =2
      isarray(3) =4
      idegs=2


      do i =1,idegs+1
      iws(i) =isarray(i)
      end do
      do i=1,idegb+1
      iwb(i) =ibarray(i)
      end do
      iwsd = idegs
      iwbd = idegb
      loopl = iwbd -iwsd +1
      mul = inv(iws(1))
      do i =1,loopl
      iqt(i) =iwb(i) *mul
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
      irarray(jj) =iwb(i+jj-1)
      end do
      idegr = itempbd -1
      goto 110
100   do i=1,10
      irarray(1) = 0
      end do
      idegr =0
110   print *,'quos',iqt(1),iqt(2),iqt(3)
      print *,'rem',irarray(1),irarray(2)
      return
      end
      subroutine subbw3(idd,ipd)
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
      do i=1,idd
      print *,'firmatrix',(mat2(i,jf),jf=1,idd),'irhs2',irhs2(i)
      end do
      
      idd2 =idd 
      goto 1
      mat2(1,1)=19
      mat2(1,2)=29
      mat2(1,3)=129
      mat2(1,4)=89
      mat2(1,5)=37
      mat2(2,1)=60
      mat2(2,2)=102
      mat2(2,3)=18
      mat2(2,4)=96
      mat2(2,5)=87
      mat2(3,1)=70
      mat2(3,2)=12
      mat2(3,3)=94
      mat2(3,4)=143
      mat2(3,5)=81
      mat2(4,1)=119
      mat2(4,2)=13
      mat2(4,3)=77
      mat2(4,4)=140
      mat2(4,5)=60
      mat2(5,1)=117
      mat2(5,2)=94
      mat2(5,3)=133
      mat2(5,4)=139
      mat2(5,5)=67
      irhs2(1)=0
      irhs2(2)=0
      irhs2(3)=0
      irhs2(4)=0
      irhs2(5)=1
1     do i = 1,idd2
      irhs1(i) =irhs2(i)
      end do
      do i =1,idd2
      do j = 1,idd2
      mat1(i,j) =mat2(i,j)
      end do
      end do
      do j = 1,idd2
      jmark(j) =0
      markr(j) =0
      jx(j) =0
      end do
      do j =1,idd2
      do i =1,idd2
      
      if(mat1(i,j) .eq.0)goto 137 
      if(markr(i).eq.1)goto 137
      jmark(j) =i
      markr(i) =1
      
      ib=mat1(i,j)
      call subbw6(ipd,ib,iv)
      mul = iv
      do j1 = j,idd2
      itt = mat1(i,j1) *mul
      mat1(i,j1) = mod(itt,ipd)
      if (mat1(i,j1).ge.0)goto 10
      mat1(i,j1)=mat1(i,j1) +ipd
10    print *,'mully',mul,i,j,j1
      end do
      irhs1(i) =irhs1(i) *mul
      itt =irhs1(i)
      irhs1(i) = mod(itt,ipd)
      if(irhs1(i).ge.0)goto 12
      irhs1(i) =irhs1(i) +ipd
      if(i.eq.idd2)goto 150
12    do ik = i+1,idd2
      if(markr(ik).eq.1)goto 134
      if(mat1(ik,j).eq.0)goto 134
      mul = mat1(ik,j)
      do jj =j,idd2
      
      
      itt = mat1(ik,jj) -mat1(i,jj) *mul
      mat1(ik,jj) = mod(itt,ipd)
      if(mat1(ik,jj).ge.0)goto 132
      mat1(ik,jj) =mat1(ik,jj) +ipd
      print *,'mulijikjj',mul,i,j,ik,jj

132   end do  
      irhs1(ik) =irhs1(ik) - mul * irhs1(i)
      itt = irhs1(ik)
      irhs1(ik) = mod(itt,ipd)
      if(irhs1(ik).ge.0)goto 134
      irhs1(ik) = irhs1(ik) + ipd
134   end do
      if (j.ne.4)goto 150
      do iz=1,5
      print *,'intmat',(mat1(iz,jz),jz=1,5),irhs1(iz)
      end do
      
      goto 150
137   end do

      
      jx(j) =1 
      
      

150   end do
      do i =1,idd2
      print *,'matrix',mat1(i,1),mat1(i,2),mat1(i,3),mat1(i,4),&
      mat1(i,5),irhs1(i)
      print *,'marks',jmark(i),markr(i)
      end do
      
      do i = 1,idd2
      print *,jx(i)
      end do
      
      do i =1,idd2
      if (jx(i).eq.0)goto 152 
      print *,'matrix singular'
      stop
      goto 200
152   end do
      

      isol(idd2) =irhs1(idd2)
      print *,'jmarks',(jmark(iz),iz=1,5)
      print *,'irhs1',(irhs1(iz),iz=1,5)
      print *,'idd=',idd,'idd2=',idd2
      do j =1,idd-1
      ind =idd2 -j+1
      ind2 =ind-1
      
      im=jmark(ind2)
      isol(ind2)=irhs1(im)
      print *,'ind=',ind,'ind2=',ind2,'im',im,'mt',(mat1(im,jz),jz=1,5)
      do jj =ind,idd2
      isol(ind2)=isol(ind2)-isol(jj) *mat1(im,jj)
      isol(ind2)=mod(isol(ind2),ipd)
      if (ind.ne.2)goto 1521
      print *,'ind2',ind2,'sol',isol(ind2),'jj',jj,'im',im
1521  end do
      end do
      
      do i =1,idd2
      itt = isol(i)
      isol(i) =mod(itt,ipd)
      if(isol(i).ge.0)goto 160
      isol(i) =isol(i) +ipd
160   end do
      print *,'solutions',isol(1),isol(2),isol(3)
      
200   return
      end

       
       
      
      subroutine bmid2
      
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
      dimension idrv5(5),idrv4(5),idrv3(5),idrv2(5),idrv1(5)
      dimension match(1000),iv(1000)
      dimension ksm(20,100),ktemp(5,20),ktempl(5)
      m1(1)=0
      m1(2)=1
      m1(3)=13
      m2(1)=0
      m2(2)=1
      m2(3)=18
      ia5(1)=1
      ia5(2)=1
      ia5(3)=43
      
      idrv5(1)=1
      idrv5(2)=3
      idrv5(3)=1
      idrv5(4)=7215
      idrv5(5)=7435
      idrv4(1)=0
      idrv4(2)=1
      idrv4(3)=32
      idrv4(4)=8440
      idrv3(1)=0
      idrv3(2)=1
      idrv3(3) =27
      idrv3(4)=2529
      idrv2(1)=0
      idrv2(2)=1
      idrv2(3)=6
      idrv2(4)=784
      idrv1(1)=0
      idrv1(2)=1
      idrv1(3)=10
      idrv1(4)=282
      mm=200
      
1     format(200i5)
2     format(400i5)
      do i=1,mm
      if (iv(i).ne.0)goto 3
      end do
      print *,'file error'
      stop

3     ktempl(1)=ia5(2)
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
      
      kmpol(2,1)=mod(ia4(1)+1,2)
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
      kmpol(3,1)=mod(ia3(1)+ia5(1)+1,2)
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
      kmpol(4,1)=mod(ia2(1)+ia5(1)*2+1,2)
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
      kmpol(5,1)=mod(ia1(1)+ia5(1)*3+1,2)
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
      kmpol(6,1)=mod(ia0(1)+ia5(1)*4+1,2)
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
      
      end do
      call doubmul(idegm,idegn,5,idprod)
      do ii=1,idprod+1
      do jj=1,kdoub3(ii,2)+2
      kdoub2(ii,jj)=kdoub3(ii,jj)
      end do
      print *,'sqdrv',(kdoub2(ii,jf),jf=1,kdoub2(ii,2)+2)
      end do
      
      do ii=1,idegsm+1
      do jj=1,ksm(ii,2)+2
      kdoub1(ii,jj)=ksm(ii,jj)
      end do
      end do
      idegm=idegsm
      idegn=idprod
      call doubmul(idegm,idegn,5,idprod)





      
      


500   return
      end


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
      print *,'lenzy',ilen,ilen2
      do i=1,800
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
      
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
      dimension kdum(400),isub(400)
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
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
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


      
      
      

      subroutine subgcd2(ibig,little,igcd2)
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2)
      
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
      dimension mmul(800)
      goto 1
      kmpol(1,1)=0
      kmpol(1,2)=1
      kmpol(1,3)=8
      kmpol(2,1)=0
      kmpol(2,2)=2
      kmpol(2,3)=1
      kmpol(2,4)=1494
      kmpol(3,1)=0
      kmpol(3,2)=1
      kmpol(3,3)=5547
      kmpol(4,1)=0
      kmpol(4,2)=1
      kmpol(4,3)=957
      kmpol(5,1)=0
      kmpol(5,2)=1
      kmpol(5,3)=6373





      
1     do i=1,12
      do j=1,800
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
      common ibarray(2000),isarray(1000),igarray(1000),inv(25)
      common mult1(1000),mult2(1000),mult3(2000),mult4(2000),mult5(1000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      common karr(400),kbarr(400),kcarr(800),ipqt(400),irrr(400)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2)
      common kmpol(6,20)
      common kdoub1(6,400),kdoub2(6,400),kdoub3(12,800)
      common ipowar(400)
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
