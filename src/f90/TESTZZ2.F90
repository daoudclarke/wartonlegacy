
       program testzz2
! program for precomputations for speed multiplication of      
! ordinary numbers using Chinese Remainder Theorem. Superseded.
       
      
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)

       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)

       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
       
       
       
       dimension n(100),igg1(2600),igg2(2600)
       dimension nn(50),ie(50),igcd(100,50),ibpd(50),ibpd2(50)
       dimension kbigp(2600),kbigw(300,2600),kbigpl(2600),numb1(2000)
       dimension numb2(2000),numb3(4000),igcdar(1000,30),iprodd(50)
       dimension jprodd(1000,30),moj1(30),moj2(30),moj3(30),nne(50)
       dimension iee(50),idivcon(1000),idivcon2(1000),iear(200,20)
       dimension jeet(20),jeetm(20)
       nnv=15
       ipard=4
       marr(1)=0
       marr(2)=16
       marr(3)=137
       marr(4)=1879
       marr(5)=1742
       marr(6)=274
       marr(7)=3346
       marr(8)=8669
       marr(9)=4951
       marr(10)=7064
       marr(11)=6282
       marr(12)=1701
       marr(13)=4677
       marr(14)=3306
       marr(15)=3072
       marr(16)=4911
       marr(17)=1467
       marr(18)=6516
       marr(1)=0
       marr(2)=8
       marr(3)=1234
       marr(4)=5678
       marr(5)=1234
       marr(6)=5678
       marr(7)=1234
       marr(8)=5678
       marr(9)=1234
       marr(10)=5678
       mbarr(1)=0
       mbarr(2)=7
       mbarr(3)=1
       mbarr(4)=0
       mbarr(5)=1
       mbarr(6)=0
       mbarr(7)=1
       mbarr(8)=1
       mbarr(9)=1
!       call mpdiv2
       ipard=5
       call mpdiv3(ipard)
       print *,'mcarr new',(mcarr(jf),jf=1,mcarr(2)+2)
       call mendiv
       print *,'mcarr con',(mcarr(jf),jf=1,mcarr(2)+2)
!       stop 
       
       
       
       
       iear(1,1)=0
       iear(1,2)=1
       iear(1,3)=1
       
       do i=2,200
       do jf=1,iear(i-1,2)+2
       karr(jf)=iear(i-1,jf)
       kbarr(jf)=iear(i-1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iear(i,jf)=kcarr(jf)
       end do
       end do
       do i=186,200
       print *,'i',i,'iear',(iear(i,jf),jf=1,iear(i,2)+2)
       end do
       jeetm(1)=1
       do i=1,15
       nlen=i+1
       jeet(1)=0
       jeet(2)=nlen
       jeet(3)=1
       do jf=4,3+i
       jeet(jf)=0
       end do
       do ibig=1,200
       do jf=2,jeet(2)+2
       if (jeet(jf).gt.iear(ibig,jf))goto 183
       if (jeet(jf).lt.iear(ibig,jf))goto 182
       end do
183    end do
       print *,'ibig',ibig,'i',i
       
182    jeetm(i)=ibig-1
       end do
       print *,'jeetm',(jeetm(jf),jf=1,15)


       
       
!       stop

       numblen=0
       igcdar(1,1)=0
       igcdar(1,2)=4
       igcdar(1,3)=1
       igcdar(1,4)=0
       igcdar(1,5)=0
       igcdar(1,6)=1
       igcdar(1,7)=0
       igcdar(1,8)=1
       igcdcon=1
       numb1(1)=0
       numb1(2)=5
       numb1(3)=1
       do jf=4,numb1(2)+1
       numb1(jf)=0
       end do
       injj=numb1(2)+2
       numb1(injj)=1
!       do jbig=1,7
       do jbig=1,7
       loopl=2**(numb1(2)-2)
       iconl=0
       jk=3
167    do ibig=1,igcdcon
       do jf=1,numb1(2)+2
       karb(jf)=numb1(jf)
       end do
       do jf=1,igcdar(ibig,2)+2
       kara(jf)=igcdar(ibig,jf)
       end do
       call subgcd2
       if ((kard(2).eq.1).and.(kard(3).eq.1))goto 164
       goto 165
164    end do
       igcdcon=igcdcon+1
       numblen=numblen+numb1(2)
       print *,'gcd count',igcdcon,'numblen',numb1(2),'jk',jk,'tot',&
       numblen,'co-prime',(numb1(jf),jf=1,numb1(2)+2)
       do jf=1,numb1(2)+2
       igcdar(igcdcon,jf)=numb1(jf)
       end do
165    iconl=iconl+1
      if (iconl.eq.loopl)goto 170
       
171    jk=numb1(2)+1
168    numb1(jk)=numb1(jk)+1
!       print *,'numm1',(numb1(jf),jf=1,numb1(2)+2)
       if (numb1(jk).eq.1)goto 167
       numb1(jk)=0
       jk=jk-1
       if (jk.eq.3)goto 171
       goto 168
170    numb1(2)=numb1(2)+1
       
       numb1(3)=1
       do jf=4,numb1(2)+1
       numb1(jf)=0
       end do
       injj=numb1(2)+2
       numb1(injj)=1
       
       end do
!       igcdcon=19
       do i=1,igcdcon
       idivcon(i)=0
       do jf=4,igcdar(i,2)+1
       if (igcdar(i,jf).eq.0)goto 173
       idivcon(i)=idivcon(i)+1
173    end do       
       end do
       do i=1,igcdcon
       do jf=4,igcdar(i,2)+1
       if (igcdar(i,jf).eq.0)goto 174
       idivcon2(i)=igcdar(i,2)+1-jf
       goto 175
174    end do
       idivcon2(i)=-1
175    end do

       print *,'idivcon',(idivcon(jf),jf=1,10)
      print *,'idivcon2',(idivcon2(jf),jf=1,10)
!       stop
       goto 142
       stop
       do i=1,1
       do j=1,1
       mmat(i,j)=i*j
       end do
       end do
       
       print *,'end of precomputations'
!       goto 163
       do ibig=1,100
       do jf=3,3002
       marr(jf)=1111
       mbarr(jf)=6666
       end do
       
       marr(1)=0
       mbarr(1)=0
       marr(2)=3000
       mbarr(2)=3000
       call menmul
       end do
       print *,'end conv. mult. mcarr(2)',mcarr(2)
       stop
163    do ibig=1,100
       do jf=1,4000
       marr(jf+2)=111
       mbarr(jf+2)=666
       end do
       marr(1)=0
       mbarr(1)=0
       marr(2)=4000
       mbarr(2)=4000
       call mpmul2
       end do
       print *,'end new mult. mcarr2',mcarr(2)
       stop

!       marr(1)=0
!       marr(2)=9
!       do jf=3,11
!       marr(jf)=1111*(jf-2)
!       end do
!       mbarr(1)=0
!       mbarr(2)=4
!       mbarr(3)=1
!       mbarr(4)=0
!       mbarr(5)=0
!       mbarr(6)=1
!       call mpdiv2
!       
!       print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
!       stop
!       goto 162
       do ii=1,10000
       marr(1)=0
       marr(2)=5650
       
       do jf=3,5651,2
       marr(jf)=5678
       end do
       do jf=4,5652,2
       marr(jf)=1234
       end do

       
       mbarr(1)=0
       mbarr(2)=4
       mbarr(3)=1234
       mbarr(4)=5678
       mbarr(5)=9876
       mbarr(6)=5432
       call mendiv
       
       end do
       print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
       stop
162    do ii=1,100000
       marr(1)=0
       marr(2)=5650
       do jf=3,5651,2
       marr(jf)=5678
       end do
       do jf=4,5652,2
       marr(jf)=1234
       end do

       
       mbarr(1)=0
       mbarr(2)=4
       mbarr(3)=1
       mbarr(4)=0
       mbarr(5)=0
       mbarr(6)=1
       call mpdiv2
       end do
       print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
       stop

       
       
       open(unit=2,file='e:\testz',access='direct',form=&
       'formatted',recl=4000,status='old')
4992   format (1000(i4))       
!       read (2,4992,rec=2000)(kbigw(jf),jf=1,10000)
!       print *,'kbigw1-100',(kbigw(jf),jf=1,200)
       
       
       open (unit=3,file='recl.dat',access='direct',form=&
       'formatted',recl=390000,status='old')
       read (3,5,rec=1)(ipr(jf),jf=1,65000)
       print *,'pr 122-302',(ipr(jf),jf=122,302)
5      format (65000i6)       
       print *,'ipr1000',ipr(1000),'ipr1100',ipr(1100),'ipr1200',ipr(1200)
       
       goto 142
150    do i=1,2600
       numb3(i)=0
       end do
       numb1(1)=0
       numb1(2)=1200
!       numb1(2)=3
       do jf=3,1201,2
       numb1(jf)=1234
       end do
       do jf=4,1202,2
       numb1(jf)=5678
       end do
       
       numb2(1)=0
       numb2(2)=1200
!       numb2(2)=2
       do jf=3,1201,2
       numb2(jf)=1111
       end do
       do jf=4,1202,2
       numb2(jf)=2222
       end do

       do jf=1,numb1(2)+2
       marr(jf)=numb1(jf)
       end do
       do jf=1,numb2(2)+2
       mbarr(jf)=numb2(jf)
       end do
       call menmul
       print *,'end conventional multiplication'
       do ibig=1,igcdcon
       do jf=1,numb1(2)+2
       marr(jf)=numb1(jf)
       end do
       
       do jf=1,igcdar(ibig,2)+2
       mbarr(jf)=igcdar(ibig,jf)
       end do
       if ((idivcon2(ibig).eq.0).and.(idivcon(ibig).eq.0))goto 190
       ipard=idivcon2(ibig)
       call mpdiv3(ipard)
!       call mendiv
       goto 191
190    call mpdiv2
!190    call mendiv
191    if (mcarr(2).eq.0)goto 140
!       print *,'first div',(mcarr(jf),jf=1,mcarr(2)+2)
       do jf=1,mcarr(2)+2
       moj1(jf)=mcarr(jf)
       end do
       do jf=1,numb2(2)+2
       marr(jf)=numb2(jf)
       end do
       if ((idivcon2(ibig).eq.0).and.(idivcon(ibig).eq.0))goto 180
       ipard=idivcon2(ibig)
       call mpdiv3(ipard)
!       call mendiv
       goto 181
180    call mpdiv2
!180    call mendiv
181    if (mcarr(2).eq.0)goto 140
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,moj1(2)+2
       mbarr(jf)=moj1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,igcdar(ibig,2)+2
       mbarr(jf)=igcdar(ibig,jf)
       end do
       if ((idivcon2(ibig).eq.0).and.(idivcon(ibig).eq.0))goto 210
       ipard=idivcon2(ibig)
       call mpdiv3(ipard)
!       call mendiv
       goto 211
210    call mpdiv2
!210    call mendiv
211    if (mcarr(2).eq.0)goto 140
       do jf=1,mcarr(2)+2
       nne(jf)=mcarr(jf)
       end do
       injj=nne(2)
       if (injj.gt.1)goto 194
       ieen=1
       goto 195
194     ieen=jeetm(injj-1)
195    do jf=1,iear(ieen,2)+2
       iee(jf)=iear(ieen,jf)
       end do
!       print *,'ibig',ibig

!       print *,'nne',(nne(jf),jf=1,nne(2)+2)
!       print *,'ieen',ieen,'iee',(iee(jf),jf=1,iee(2)+2)
!       print *,'igcdar',(igcdar(ibig,jf),jf=1,igcdar(ibig,2)+2)
!       stop
152    do jf=2,iee(2)+2
       if (iee(jf).gt.nne(jf))goto 153
       if (iee(jf).lt.nne(jf))goto 1531
       end do

1531   ieen=ieen+1
       do jf=1,iear(ieen,2)+2
       iee(jf)=iear(ieen,jf)
       end do
       goto 152

       
153    ieen=ieen-1
       do jf=1,iear(ieen,2)+2
       iee(jf)=iear(ieen,jf)
       end do


       do jf=1,nne(2)+2
       karr(jf)=nne(jf)
       end do
       do jf=1,iee(2)+2
       kbarr(jf)=iee(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nne(jf)=kcarr(jf)
       end do
       ibrec=ibig


!       do ibig=7,1100
!       do jf=1,numb1(2)+2
!       marr(jf)=numb1(jf)
!       end do
!       mbarr(1)=0
!       mbarr(2)=1
!       mbarr(3)=ipr(ibig)
!       call mendiv
!       if (mcarr(2).eq.0)goto 140
!       moj1=mcarr(3)
!       do jf=1,numb2(2)+2
!       marr(jf)=numb2(jf)
!       end do
!       call mendiv
!       if (mcarr(2).eq.0)goto 140
!       moj2=mcarr(3)
!       moj3=mod(moj1*moj2,ipr(ibig))
!       nne=moj3
!       iee=1
!152    if (iee.gt.nne)goto 153
!       iee=iee+iee
!       goto 152
!153    iee=iee/2
!       ibrec=ibig-6
!       print *,'nne',nne,'iee',iee,'kbigw2',kbigw(ibrec,2)
!       nne=nne-iee
       
       do jf=1,kbigw(ibrec,2)+2
       igg1(jf)=kbigw(ibrec,jf)
       igg2(jf)=kbigw(ibrec,jf)
       end do
154    if ((iee(2).eq.1).and.(iee(3).eq.1))goto 160
!       iee=iee/2
       ieen=ieen-1
       do jf=1,iear(ieen,2)+2
       iee(jf)=iear(ieen,jf)
       end do
       

       do jf=1,igg2(2)+2
       karr(jf)=igg2(jf)
       kbarr(jf)=igg2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       igg2(jf)=kcarr(jf)
       end do
       do jf=2,nne(2)+2
       if (nne(jf).lt.iee(jf))goto 154
       if (nne(jf).gt.iee(jf))goto 1541
       end do
1541   do jf=1,nne(2)+2
       karr(jf)=nne(jf)
       end do
       do jf=1,iee(2)+2
       kbarr(jf)=iee(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nne(jf)=kcarr(jf)
       end do
!       if (nne.lt.iee)goto 154
!       nne=nne-iee
       do jf=1,igg1(2)+2
       karr(jf)=igg1(jf)
       end do
       do jf=1,igg2(2)+2
       kbarr(jf)=igg2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       igg2(jf)=kcarr(jf)
       end do
       goto 154
160    do jf=1,igg2(2)+2
       karr(jf)=igg2(jf)
       end do
!       print *,'kbigw',(kbigw(ibig,jf),jf=1,kbigw(ibig,2)+2)
!       print *,'igg2',(igg2(jf),jf=1,igg2(2)+2)
       do jf=1,numb3(2)+2
       kbarr(jf)=numb3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       numb3(jf)=kcarr(jf)
       end do





140    end do
!       stop
!       print *,'product',(numb3(jf),jf=1,numb3(2)+2)
!       print *,'product reached'
       do jf=1,numb3(2)+2
       marr(jf)=numb3(jf)
       end do
       do jf=1,kbigp(2)+2
       mbarr(jf)=kbigp(jf)
       end do
       call mendiv
       print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
       stop
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=6
       mbarr(3)=1
       mbarr(4)=1
       mbarr(5)=0
       mbarr(6)=1
       mbarr(7)=1
       mbarr(8)=1
       mbarr(9)=1
       call mendiv
!       print *,'mcarr con',(mcarr(jf),jf=1,mcarr(2)+2)
       ipard=3
       call mpdiv3(ipard)
!       print *,'mcarr new',(mcarr(jf),jf=1,mcarr(2)+2)
       mbarr(1)=0
       mbarr(2)=6
       mbarr(3)=1
       mbarr(4)=1
       mbarr(5)=0
       mbarr(6)=1
       mbarr(7)=0
       mbarr(8)=1
       call mendiv
!       print *,'mcarr con',(mcarr(jf),jf=1,mcarr(2)+2)
       ipard=3
       call mpdiv3(ipard)
!       print *,'mcarr new',(mcarr(jf),jf=1,mcarr(2)+2)
       marr(1)=0
       marr(2)=8
       marr(3)=1234
       marr(4)=5678
       marr(5)=1234
       marr(6)=5678
       marr(7)=1234
       marr(8)=5678
       marr(9)=1234
       marr(10)=5678
       do jf=1,igcdar(nnv,2)+2
       mbarr(jf)=igcdar(nnv,jf)
       end do
       print *,'divisor',(mbarr(jf),jf=1,mbarr(2)+2)
       call mendiv
       print *,'mcarr con',(mcarr(jf),jf=1,mcarr(2)+2)
       if (idivcon(nnv).eq.0)goto 197
       ipard=idivcon2(nnv)
       call mpdiv3(ipard)
       print *,'mcarr new',(mcarr(jf),jf=1,mcarr(2)+2),'divcon',&
       idivcon(nnv)
       if (idivcon(nnv).ne.0)goto 198
197    call mpdiv2
       print *,'mcarr di2',(mcarr(jf),jf=1,mcarr(2)+2)
198    marr(1)=0
       marr(2)=12
       marr(3)=137
       marr(4)=1879
       marr(5)=1742
       marr(6)=274
       marr(7)=3346
       marr(8)=8669
       marr(9)=4677
       marr(10)=3306
       marr(11)=3072
       marr(12)=4911
       marr(13)=1467
       marr(14)=6516
       call mendiv
!       print *,'mcarr ori',(mcarr(jf),jf=1,mcarr(2)+2)
       


       stop
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=6
       kbarr(3)=137
       kbarr(4)=1879
       kbarr(5)=1741
       kbarr(6)=9369
       kbarr(7)=0343
       kbarr(8)=0974
       call mpadd(1)
       print *,'sum',(kcarr(jf),jf=1,kcarr(2)+2)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=4
       kbarr(3)=9863
       kbarr(4)=7628
       kbarr(5)=205
       kbarr(6)=9233
       call mpadd(1)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=5
       do jf=3,7
       mbarr(jf)=1
       end do
       mbarr(6)=0
       call mendiv
       print *,'sec mcarr',(mcarr(jf),jf=1,mcarr(2)+2)
       
       print *,'kbigp2',kbigp(2),'mcarr2',mcarr(2),'mcarr3',mcarr(3)
       print *,'kbigp',(kbigp(jf),jf=1,kbigp(2)+2)
       stop
142    a=a
       do i=1,igcdcon
       iprodd(1)=0
       iprodd(2)=1
       iprodd(3)=1
       

       do j=1,igcdcon
       if (i.eq.j)goto 130
       do jf=1,iprodd(2)+2
       marr(jf)=iprodd(jf)
       end do
       do jf=1,igcdar(j,2)+2
       mbarr(jf)=igcdar(j,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,igcdar(i,2)+2
       mbarr(jf)=igcdar(i,jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       iprodd(jf)=mcarr(jf)
       end do
130    end do
       do jf=1,iprodd(2)+2
       karb(jf)=iprodd(jf)
       end do
       do jf=1,igcdar(i,2)+2
       kara(jf)=igcdar(i,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       jprodd(i,jf)=karv(jf)
       end do
       end do
       marr(1)=0
       marr(2)=1
       marr(3)=1
       do ibig=1,igcdcon
       do jf=1,igcdar(ibig,2)+2
       mbarr(jf)=igcdar(ibig,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       end do
       do jf=1,marr(2)+2
       kbigp(jf)=marr(jf)
       end do
       
       print *,'kbigp2',kbigp(2)
       do ibig=1,igcdcon
       do i=1,2500
       kbigw(ibig,i)=0
       end do
       
       do jf=1,kbigp(2)+2
       marr(jf)=kbigp(jf)
       end do
       do jf=1,igcdar(ibig,2)+2
       mbarr(jf)=igcdar(ibig,jf)
       end do
       
       call mendiv
       if (mcarr(2).ne.0)goto 199
       do jf=1,mdarr(2)+2
       marr(jf)=mdarr(jf)
       end do
       do jf=1,jprodd(ibig,2)+2
       mbarr(jf)=jprodd(ibig,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       kbigw(ibig,jf)=mcarr(jf)
       end do
       end do


!       mbarr(1)=0
!       mbarr(2)=1
!       mbarr(3)=jprodd(ibig)
!       call menmul
!       do jf=1,mcarr(2)+2
!       kbigpl(jf)=mcarr(jf)
!       end do
!       ibrec=ibig-6
       
!        do i=1,kbigpl(2)+2
!       kbigw(ibig-6,i)=kbigpl(i)
!       end do
       
!       end do
       print *,'end memory creation phase'
       print *,'igcdcon',igcdcon
       do i=2,igcdcon
       do jf=1,igcdar(i,2)+2
       kara(jf)=igcdar(i,jf)
       end do
       do j=1,i-1
       do jf=1,igcdar(j,2)+2
       karb(jf)=igcdar(j,jf)
       end do
       call subgcd2
       if ((kard(2).eq.1).and.(kard(3).eq.1))goto 200 
       goto 199
200    end do 
       print *,'i ok',i
       end do
!       stop
       goto 150
199    print *,'big big graunch'
       stop



       
       stop

       
       end
       
       subroutine menmul
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)

       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)

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
9      mcarr(1)=0
       mcarr(2)=0
10     return
       end
       subroutine mendiv
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)

       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
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
6      if (icont.eq.0)goto 12
       do jf=1,icont
       mdarr(jf+2)=ipqt(jf)
       end do
       mdarr(2)=icont
       mdarr(1)=mod(marr(1)+mbarr(1),2)
       goto 15
11     print *,'halted : attempted division by zero'
       stop
10     mcarr(1)=0
       mcarr(2)=0
       goto 6
9      mcarr(1)=0 
       mcarr(2)=0
12     mdarr(1)=0
       mdarr(2)=0
15     return
       end



       
       subroutine subbw6(ia,ib,iv)
       
       
       ib =mod(ib,ia)
       
       
       
       if(ib.ge.0)goto 10
       ib =ia + ib
10     iu =1
       
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
       if (id.gt.1)goto 890
       
       goto 890
888    iv =ib
890    return
       end
       
   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
      
      
      
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

      subroutine mpmul2
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
       ilen3=marr(2)+mbarr(2)
      
      
      do i=1,ilen3
      mcarr(i+2) =0
      
      end do
      do i =1,marr(2)
      do j=1,mbarr(2)
!      itemp =karr(i) *kbarr(j)
      ii=marr(i+2)
      jj=mbarr(j+2)


      itemp=mmat(ii,jj)
      itemp2=int(itemp/1000)
      irem1 =itemp-itemp2 *1000
      
      
      mcarr(i+j+2)=mcarr(i+j+2) +irem1
      mcarr(i+j+1)=mcarr(i+j+1) +itemp2
      
      if (mcarr(i+j+2).lt.1000)goto 20
      mcarr(i+j+2)=mcarr(i+j+2)-1000
      mcarr(i+j+1) =mcarr(i+j+1)+1
20    do k =1,i+j-2
      if (mcarr(i+j-k+2).lt.1000)goto 22
      mcarr(i+j-k+2)=mcarr(i+j-k+2) -1000
      mcarr(i+j-k+1)=mcarr(i+j-k+1)+1
      end do
22    end do
      end do
      
      if(mcarr(3).ne.0)goto 100
      ilen3 =ilen3-1
      do i=1,ilen3
      mcarr(i+2)=mcarr(i+3)
      end do
100   mcarr(1)=mod(marr(1)+mbarr(1),2)   
      mcarr(2)=ilen3
      
      
      return
      end

      subroutine mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
      
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
      print *,'ilen',ilen,'karr',(karr(jf),jf=1,ilen)
      print *,'ilen2',ilen2,'kbarr',(kbarr(jf),jf=1,ilen2)
      stop
      
910   return      
      end
      subroutine mpdiv2
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
      
      
      
      if (marr(2).lt.mbarr(2))goto 35
      do jf=3,marr(2)+2
      kdum(jf-2)=marr(jf)
      end do
      do jf=3,mbarr(2)+2
      isub(jf-2)=mbarr(jf)
      end do
      
      icnt1=1
      icnt2=mbarr(2)
      icnte=marr(2)-mbarr(2)+1
      icntb=mbarr(2)
      isub(1)=kdum(1)
      isub(icntb)=kdum(1)
      do jf=icnt1+1,icnt2-1
      if (kdum(jf).ne.0)goto 1
      end do
      if (kdum(icnt2).eq.isub(icntb))goto 2
      if (kdum(icnt2).gt.isub(icntb))goto 1
      if (kdum(icnt1).eq.1)goto 3
      isub(icnt1)=isub(icnt1)-1
      isub(icnt2)=isub(icnt2)-1
      do jf=icnt2,icnt1+1,-1
      kdum(jf)=kdum(jf)-isub(icntb-icnt2+jf)
      if (kdum(jf).ge.0)goto 2
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
      end do
      print *,'error stop 1'
      stop
1     do jf=icnt2,icnt1+1,-1
!      print *,'jf',jf,'icnt1',icnt1,'icnt2',icnt2,'icntb',icntb
      kdum(jf)=kdum(jf)-isub(icntb-icnt2+jf)
!      print *,'isub',isub(icntb-icnt2+jf)
      if (kdum(jf).ge.0)goto 2
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
      end do
!      print *,'error stop 3'
      stop
3     icnt2=icnt2+1     
      if (icnt2.gt.marr(2))goto 30
      kdum(icnt2)=kdum(icnt2)-1
      if (kdum(icnt2).ge.0)goto 31
      kdum(icnt2)=kdum(icnt2)+10000
31    kdum(icnt1)=0      
      kdum(icnt1+1)=9998
      if (icnt1+2.eq.icnt2)goto 21
      do jf=icnt1+2,icnt2-1
      kdum(jf)=9999
      end do
21    icnt1=icnt1+1
      goto 22


2     icnt1=icnt1+1
      if (icnt1.gt.icnte)goto 30
      icnt2=icnt2+1
!      print *,'at 2 icnt1',icnt1,'icnt2',icnt2
!      print *,'kdum at 2',(kdum(jf),jf=1,marr(2)),'icnt1',icnt1,&
!      'icnt2',icnt2
22    isub(icntb)=kdum(icnt1)
      do jf=icnt1,marr(2)
      if (kdum(jf).eq.0)goto 17
      icnt1=jf
      if (icnt1.gt.icnte)goto 30
      goto 18
17    end do      
18    icnt2=icntb+jf-1
!      print *,'at 18 icnt2',icnt2,'icnt1',icnt1
      isub(icntb)=kdum(icnt1)
!      print *,'kdum at 18',(kdum(jf),jf=icnt1,icnt2-1)
      do jf=icnt1+1,icnt2-1
      if (kdum(jf).ne.0)goto 1
      end do
      if (kdum(icnt2).eq.isub(icntb))goto 1
      if (kdum(icnt2).gt.isub(icntb))goto 1
      if (kdum(icnt1).eq.1)goto 3
!      print *,'it came in here'
      isub(icntb)=isub(icntb)-1
      do jf=icnt2,icnt1+1,-1
      kdum(jf)=kdum(jf)-isub(icntb-icnt2+jf)
      if (kdum(jf).ge.0)goto 2
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
      end do
!      print *,'error stop 2'
!      print *,'kdum at error 2',(kdum(jf),jf=1,8)
      goto 2
35    do jf=1,marr(2)+2      
      mcarr(jf)=marr(jf)
      end do
      goto 36
23    kdum(icnt2)=0      
      goto 2
      
      stop
30    mcarr(1)=0  
      do jf=icnt1,marr(2)
      if (kdum(jf).eq.0)goto 34
      mcarr(2)=marr(2)-jf+1
      goto 37
34    end do
      mcarr(1)=0
      mcarr(2)=0
      goto 36
37    do jk=3,mcarr(2)+2
      mcarr(jk)=kdum(jf)
      jf=jf+1
      end do
      
36    return
      end 
      subroutine mpdiv3(ipard)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
     
      

      if (marr(2).lt.mbarr(2))goto 35
      do jf=3,marr(2)+2
      kdum(jf-2)=marr(jf)
      end do
      do jf=3,mbarr(2)+2
      isub(jf-2)=mbarr(jf)
      end do
      
      icnt1=1
      icnt2=mbarr(2)
      icnte=marr(2)-mbarr(2)+1
      icntb=mbarr(2)
      do jf=3,mbarr(2)+2
      if (isub(jf-2).eq.0)goto 100
      isub(jf-2)=kdum(1)
100   end do      
      do jf=icnt1+1,icnt2
      if (kdum(jf).lt.isub(jf))goto 101
      if (kdum(jf).gt.isub(jf))goto 1
      end do
      goto 2
101   if (kdum(icnt1).eq.1)goto 3      
107   do jf=1,icntb
      if (isub(jf).eq.0)goto 102
      isub(jf)=isub(jf)-1
102   end do
      jk=icntb
      do jf=icnt2,icnt1,-1
      kdum(jf)=kdum(jf)-isub(jk)
      jk=jk-1
      if (kdum(jf).ge.0)goto 110
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
110   end do      
!      print *,'kdum at 110',kdum(icnt1)
      if (kdum(icnt1).eq.1)goto 3
      
      goto 2
1     do jf=icnt2,icnt2-1-ipard,-1
      kdum(jf)=kdum(jf)-isub(icntb-icnt2+jf)
!      jazz=isub(icntb-icnt2+jf)
!      print *,'sub',jazz
      if (kdum(jf).ge.0)goto 103
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
103   end do
!      print *,'kdum',(kdum(jf),jf=1,marr(2))
      goto 2
      
3     icnt2=icnt2+1     
      if (icnt2.gt.marr(2))goto 30
      do jf=3,mbarr(2)+2
      isub(jf-2)=mbarr(jf)
      end do
108   do jf=icnt2,icnt2-icntb+1,-1
      kdum(jf)=kdum(jf)-isub(icntb-icnt2+jf)

      if (kdum(jf).ge.0)goto 104
      kdum(jf)=kdum(jf)+10000
      kdum(jf-1)=kdum(jf-1)-1
104   end do
      
21    icnt1=icnt1+1
      goto 22


2     icnt1=icnt1+1
      if (icnt1.gt.icnte)goto 30
      icnt2=icnt2+1
!      print *,'at 2 icnt1',icnt1,'icnt2',icnt2,'icntb',icntb
!      print *,'kdum at 2',(kdum(jf),jf=1,marr(2))
22    a=a
      do jf=icnt1,marr(2)
      if (kdum(jf).eq.0)goto 17
      icnt1=jf
      if (icnt1.gt.icnte)goto 30
      goto 18
17    end do      
18    icnt2=icntb+jf-1
!      print *,'at 18 icnt2',icnt2,'icnt1',icnt1,'icntb',icntb
!      print *,'isub before',(isub(jf),jf=1,mbarr(2))
      do jf=3,mbarr(2)+2
      if (isub(jf-2).eq.0)goto 105
      isub(jf-2)=kdum(icnt1)
105   end do      
!      print *,'icnt2',icnt2,'isub',(isub(jf),jf=1,mbarr(2)),&
!      'icntb',icntb
      jk=2
      do jf=icnt1+1,icnt2
      if (kdum(jf).lt.isub(jk))goto 106
      if (kdum(jf).gt.isub(jk))goto 1
      jk=jk+1
      end do
      goto 1
106   if (kdum(icnt1).eq.1)goto 3      
      goto 107
      
35    do jf=1,marr(2)+2      
      mcarr(jf)=marr(jf)
      end do
      goto 36
      
      
      
30    mcarr(1)=0  
      do jf=icnt1,marr(2)
      if (kdum(jf).eq.0)goto 34
      mcarr(2)=marr(2)-jf+1
      goto 37
34    end do
      mcarr(1)=0
      mcarr(2)=0
      goto 36
37    do jk=3,mcarr(2)+2
      mcarr(jk)=kdum(jf)
      jf=jf+1
      end do
      
36    return
      end 
      subroutine mpadd(isora)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
      
      
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
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
       do jf=3,karb(2)+2
       karr(jf-2)=karb(jf)
       end do
       ilen=karb(2)
       do jf=3,kara(2)+2
       kbarr(jf-2)=kara(jf)
       end do
       ilen2=kara(2)
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
       do jf=1,ilen2
       kard(jf+2)=kbarr(jf)
       end do
       return
       end







      subroutine mpgcd
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
      
      
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

      subroutine mpkron(k)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(12000),kbarr(12000),kcarr(12000),ipqt(12000),irrr(12000)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(200,50),jpol(100,50),jpol2(100,50)
       common iqt(200,50),ipd(50)
       common marr(12000),mbarr(12000),mcarr(12000),mdarr(12000)
       common nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100),mmat(1,1)
       common kdum(12000),isub(12000)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      
      
      
      dimension karae(50),karde(50),karpe(50),karh(50),ktemp(50)
      do jf=2,kard(2)+2
      if (kard(jf).ne.karp(jf))goto 2
      end do
      goto 3
2     a=a      
      if ((kard(2).eq.1).and.(kard(3).eq.0))goto 3
      
      if (kard(2).ne.0)goto 1
      
3     print *,'kard',(kard(jf),jf=1,kard(2)+2)
      print *,'number is pseudo prime'
      stop
1     a=a
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
513   a=a
!513   print *,'k',k
      return
      end
       
      
          
       

       


      
      
