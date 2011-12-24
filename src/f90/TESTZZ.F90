
       program testzz
! program for precomputations for speed multiplication of      
! ordinary integers using Chinese Remainder Theorem.Superseded
       
      
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
       
       
       
       dimension n(100),ipol(100,50),igg1(200,50),igg2(200,50),nfsol(100,50)
       dimension nn(50),ie(50),igcd(100,50),jprodd(1200),ibpd(50),ibpd2(50)
       dimension kbigp(12000),kbigw(12000),kbigpl(12000),numb1(1000)
       dimension numb2(1000),numb3(1000)
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
       
!       goto 142
       do i=1,1000
       numb3(i)=0
       end do
       numb1(1)=0
       numb1(2)=450
       do jf=3,452
       numb1(jf)=1234
       end do
       numb2(1)=0
       numb2(2)=450
       do jf=3,452
       numb2(jf)=5678
       end do
       do jf=1,numb1(2)+2
       marr(jf)=numb1(jf)
       end do
       do jf=1,numb2(2)+2
       mbarr(jf)=numb2(jf)
       end do
       call menmul
       

       
       do ibig=7,1100
       do jf=1,numb1(2)+2
       marr(jf)=numb1(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(ibig)
       call mendiv
       if (mcarr(2).eq.0)goto 140
       moj1=mcarr(3)
       do jf=1,numb2(2)+2
       marr(jf)=numb2(jf)
       end do
       call mendiv
       if (mcarr(2).eq.0)goto 140
       moj2=mcarr(3)
       moj3=mod(moj1*moj2,ipr(ibig))
       if (ibig.ne.200)goto 146
       print *,'moj3',moj3
!       stop
146    moj4=mod(moj3,10)
       if (moj4.eq.0)goto 143
       ibrec=(ibig-7)*9+moj4
       read (2,4992,rec=ibrec)(kbigw(jf),jf=1,1000)
       do jf=1,numb3(2)+2
       karr(jf)=numb3(jf)
       end do
       do jf=1,kbigw(2)+2
       kbarr(jf)=kbigw(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       numb3(jf)=kcarr(jf)
       end do
143    if (moj3.lt.10)goto 140
       moj5=mod(moj3,100)-moj4
       moj5=moj5/10
       if (moj5.eq.0)goto 144
       ibrec=(ibig-7)*9+moj5
       read (2,4992,rec=ibrec)(kbigw(jf),jf=1,1000)
       do jf=1,kbigw(2)+2
       marr(jf)=kbigw(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=10
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,numb3(2)+2
       kbarr(jf)=numb3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       numb3(jf)=kcarr(jf)
       end do
144    if (moj3.lt.100)goto 140
       moj6=mod(moj3,1000)-10*moj5-moj4
       moj6=moj6/100
       if (moj6.eq.0)goto 145
       ibrec=(ibig-7)*9+moj6
       read (2,4992,rec=ibrec)(kbigw(jf),jf=1,1000)
       do jf=1,kbigw(2)+2
       marr(jf)=kbigw(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=100
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,numb3(2)+2
       kbarr(jf)=numb3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       numb3(jf)=kcarr(jf)
       end do
145    if (moj3.lt.1000)goto 140
       moj7=moj3-100*moj6-10*moj5-moj4
       moj7=moj7/1000
       if (ibig.ne.200)goto 147
       print *,'moj7',moj7,'moj6',moj6,'moj5',moj5,'moj4',moj4
!       stop
147    ibrec=(ibig-7)*9+moj7
       read (2,4992,rec=ibrec)(kbigw(jf),jf=1,1000)
       do jf=1,kbigw(2)+2
       marr(jf)=kbigw(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=1000
       call menmul
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,numb3(2)+2
       kbarr(jf)=numb3(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       numb3(jf)=kcarr(jf)
       end do
140    end do
!       print *,'product',(numb3(jf),jf=1,numb3(2)+2)
       print *,'product reached'
       marr(1)=0
       marr(2)=1
       marr(3)=1
       do i=7,1100
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(i)
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       end do
       do jf=1,marr(2)+2
       kbigp(jf)=marr(jf)
       end do
       do jf=1,numb3(2)+2
       marr(jf)=numb3(jf)
       end do
       do jf=1,kbigp(2)+2
       mbarr(jf)=kbigp(jf)
       end do
       call mendiv
       print *,'mcarr',(mcarr(jf),jf=1,mcarr(2)+2)


       
       print *,'kbigp2',kbigp(2),'mcarr2',mcarr(2),'mcarr3',mcarr(3)

       stop
142    a=a















       do i=7,1100
       iprodd=1
       do j=7,1100
       if (i.eq.j)goto 130
       ia=ipr(i)
       iprodd=iprodd*ipr(j)
       iprodd=mod(iprodd,ia)
130    end do
       ib=iprodd
       call subbw6(ia,ib,iv)
       jprodd(i)=iv
       end do
       print *,'jprodd',(jprodd(jf),jf=1000,1100)
       marr(1)=0
       marr(2)=1
       marr(3)=1
       do i=7,1100
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(i)
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       end do
       do jf=1,marr(2)+2
       kbigp(jf)=marr(jf)
       end do
       print *,'kbigp2',kbigp(2)
       do ibig=7,1100
       
       do jf=1,kbigp(2)+2
       marr(jf)=kbigp(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(ibig)
       call mendiv
       do jf=1,mdarr(2)+2
       marr(jf)=mdarr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=jprodd(ibig)
       call menmul
       do jf=1,mcarr(2)+2
       kbigpl(jf)=mcarr(jf)
       end do
       do ii=1,9
       do i=1,1000
       kbigw(i)=0
       end do
       do jf=1,kbigpl(2)+2
       marr(jf)=kbigpl(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ii
       call menmul
       do jf=1,mcarr(2)+2
       kbigw(jf)=mcarr(jf)
       end do
       ibrec=(ibig-7)*9+ii
       write(2,4992,rec=ibrec)(kbigw(jf),jf=1,1000) 
       
       print *,'rec=',ibrec,'ibig',ibig,'kbigw2',kbigw(2)
       end do
       end do







       
       stop

       
       marr(1)=0
       marr(2)=1
       marr(3)=ipr(305)
       do i=7,1200
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(i)
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       end do
       print *,'prod',(marr(jf),jf=1,marr(2)+2)
       do jf=1,marr(2)+2
       ibpd(jf)=marr(jf)
       ipd(jf)=marr(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(330)
       call mendiv
       do jf=1,mdarr(2)+2
       ibpd2(jf)=mdarr(jf)
       end do
       mult1(1,1)=0
       mult1(1,2)=1
       mult1(1,3)=1
       mult1(2,1)=0
       mult1(2,2)=1
       mult1(2,3)=1
       idegm=1
       idegn=1
       do i=1,10
       mult2(1,1)=0
       mult2(1,2)=1
       mult2(1,3)=1
       mult2(2,1)=0
       mult2(2,2)=1
       mult2(2,3)=i+1
       call multy(idegm,idegn)
       idegm=idegm+idegn 
       do ii=1,idegm+1
       do jf=1,mult3(ii,2)+2
       mult1(ii,jf)=mult3(ii,jf)
       end do
       end do
       end do
       do i=1,idegm+1
       print *,'pol1',(mult1(i,jf),jf=1,mult1(i,2)+2)
       end do
       do i=1,idegm+1
       do jk=1,mult1(i,2)+2
       ibarray(i,jk)=mult1(i,jk)
       end do
       end do
       idegb=idegm
       isarray(1,1)=0
       isarray(1,2)=1
       isarray(1,3)=1
       isarray(2,1)=0
       isarray(2,2)=1
       isarray(2,3)=5
       idegs=1
       call subbw4(idegs,idegb,idegg)
       print *,'idegg',idegg,'idegb',idegb,'idegs',idegs
       do i=1,idegb-idegs+1
       print *,'i',i,'iqt',(iqt(i,jf),jf=1,iqt(i,2)+2)
       end do
!       ia=1009
       do i=1,1100
       
       iprodd=1
       do j=1,1100
       if (i.eq.j)goto 120
       iprodd=iprodd*(i-j)
       if (iprodd.gt.0)goto 119
       iprodd=ia+iprodd
119    iprodd=mod(iprodd,ia)
120    end do
       ib=iprodd
       if (i.gt.24)goto 118
       ia=ipr(305+i)
       goto 117
118    ia =ipr(330)      
117    call subbw6(ia,ib,iv)
       jprodd(i)=iv
       end do
       
       print *,'jprodd',(jprodd(jf),jf=1,100)
       mult1(1,1)=0
       mult1(1,2)=1
       mult1(1,3)=1
       mult1(2,1)=0
       mult1(2,2)=1
       mult1(2,3)=1
       idegm=1
       idegn=1
       do i=1,1001
       mult2(1,1)=0
       mult2(1,2)=1
       mult2(1,3)=1
       mult2(2,1)=0
       mult2(2,2)=1
       mult2(2,3)=i+1
       call multy(idegm,idegn)
       idegm=idegm+idegn 
       do ii=1,idegm+1
       do jf=1,mult3(ii,2)+2
       mult1(ii,jf)=mult3(ii,jf)
       end do
       end do
       end do
       do i=1,idegm+1
!       print *,'pol1',(mult1(i,jf),jf=1,mult1(i,2)+2)
       end do
       do ii=1,idegm+1
       do jf=1,mult1(ii,2)+2
       mult5(ii,jf)=mult1(ii,jf)
       end do
       end do
       
       
       
       
       
       
       
       do ibig=1,1001
       if (ibig.gt.24)goto 121
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=ipr(305+ibig)
       do jf=1,ibpd(2)+2
       marr(jf)=ibpd(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       ipd(jf)=mdarr(jf)
       end do
       goto 122
121    mbarr(1)=0       
       mbarr(2)=1
       mbarr(3)=ipr(330)
       do jf=1,ibpd(2)+2
       marr(jf)=ibpd(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       ipd(jf)=mdarr(jf)
       end do
122    a=a
       end do


       
       stop
       issz=0
       ncony=0
       nfcon=0
       ibbsw=0
       
       print *,'type 1 for root-finding,2 for other equations'
       read *,ityp
       if (ityp.eq.2)goto 22
       print *,'modulus length in 4 digit words'
       ipd(1)=0
       read *,ipd(2)
       print *,'modulus'
       read *,(ipd(jf),jf=3,ipd(2)+2)
       print *,'degree of root'
       read *,ipold
       print *,'length of positive number to be rooted'
       read *,marr(2)
       print *,'number'
       read *,(marr(jf),jf=3,marr(2)+2)
       marr(1)=0
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
        
       
       
       call mendiv
       if (ipold.eq.1)goto 110
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       ipol(ipold+1,jf)=kcarr(jf)
       end do
       
       ipol(1,1)=0
       ipol(1,2)=1
       ipol(1,3)=1
       do jf=2,ipold
       ipol(jf,1)=0
       ipol(jf,2)=0
       end do
       goto 41


22     ipd(1)=0
       ipd(2)=8
       ipd(3)=1234
       ipd(4)=5678
       ipd(5)=9012
       ipd(6)=3456
       ipd(7)=7890
       ipd(8)=1234
       ipd(9)=5678
       ipd(10)=9097




!       ipd=40031
!       ipd=11677
!       call subbw6(11677,10000,iv)
!       print *,'iv',iv
!       stop


       goto 41
41     if (ityp.eq.1)goto 301       
       ipd(1)=0
       ipd(2)=1
       ipd(3)=3
301    do jf=1,ipd(2)+2
       nn(jf)=ipd(jf)
       end do
       ie(1)=0
       ie(2)=1
       ie(3)=1
2      do jf=2,ie(2)+2
       if (ie(jf).lt.nn(jf))goto 73
       if (ie(jf).gt.nn(jf))goto 74
       end do
73     do jf=1,ie(2)+2       
       karr(jf)=ie(jf)
       kbarr(jf)=ie(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ie(jf)=kcarr(jf)
       end do
       goto 2
74     do jf=1,ie(2)+2
       marr(jf)=ie(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ie(jf)=mdarr(jf)
       end do
       do jf=1,nn(2)+2
       karr(jf)=nn(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nn(jf)=kcarr(jf)
       end do
!       print *,'ie',(ie(jf),jf=1,ie(2)+2)
!       print *,'nn',(nn(jf),jf=1,nn(2)+2)
       
       if (ityp.eq.1)goto 21
       ipold=10
       ipol(1,1)=0
       ipol(1,2)=1
       ipol(1,3)=1
       ipol(2,1)=0
       ipol(2,2)=8
       ipol(2,3)=1234
       ipol(2,4)=5678
       ipol(2,5)=9012
       ipol(2,6)=3456
       ipol(2,7)=7890
       ipol(2,8)=1234
       ipol(2,9)=5678
       ipol(2,10)=9082-4
       
       ipol(3,1)=0
       ipol(3,2)=1
       ipol(3,3)=85
       ipol(4,1)=0
       ipol(4,2)=8
       ipol(4,3)=1234
       ipol(4,4)=5678
       ipol(4,5)=9012
       ipol(4,6)=3456
       ipol(4,7)=7890
       ipol(4,8)=1234
       ipol(4,9)=5678
       ipol(4,10)=8872
       
       ipol(5,1)=0
       ipol(5,2)=1
       ipol(5,3)=274
       ipol(6,1)=0
       ipol(6,2)=8
       ipol(6,3)=1234
       ipol(6,4)=5678
       ipol(6,5)=9012
       ipol(6,6)=3456
       ipol(6,7)=7890
       ipol(6,8)=1234
       ipol(6,9)=5678
       
       ipol(6,10)=8977
       ipol(7,1)=0
       ipol(7,2)=1
       ipol(7,3)=1
       ipol(8,1)=0
       ipol(8,2)=1
       ipol(8,3)=2
       ipol(9,1)=0
       ipol(9,2)=1
       ipol(9,3)=3
       ipol(10,1)=0
       ipol(10,2)=1
       ipol(10,3)=4
       ipol(11,1)=0
       ipol(11,2)=1
       ipol(11,3)=5

       ipold=31
       ipol(1,1)=0
       ipol(2,1)=1
       ipol(3,1)=1
       do i=2,31
       ipol(i,1)=0
       ipol(i,2)=0
       end do
       ipol(32,1)=0
       ipol(32,2)=8
       ipol(32,3)=884
       ipol(32,4)=3759
       ipol(32,5)=8411
       ipol(32,6)=9020
       ipol(32,7)=8411
       ipol(32,8)=6493
       ipol(32,9)=1086
       ipol(32,10)=6889
       ipold=2
       
       ipol(2,1)=0
       ipol(2,2)=1
       ipol(2,3)=2
       ipol(3,1)=0
       ipol(3,2)=1
       ipol(3,3)=1
21     idegp=1
       
       igg1(1,1)=0
       igg1(1,2)=1
       igg1(1,3)=1
       igg1(2,1)=0
       igg1(2,2)=0
       
       igg2(1,1)=0
       igg2(1,2)=1
       igg2(1,3)=1
       igg2(2,1)=0
       igg2(2,2)=0
       
       idegn=1
4      if ((ie(2).eq.1).and.(ie(3).eq.1))goto 30
!       print *,'ie',ie(3),'nn',nn(3)
       do jf=1,ie(2)+2
       marr(jf)=ie(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ie(jf)=mdarr(jf)
       end do

       
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult1(jf,jk)=igg2(jf,jk)
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       idegm=idegn
       call multy(idegm,idegn)
!       print *,'mult3',((mult3(jf,jk),jk=1,mult3(jf,2)+2),jf=1,&
!       idegm+idegn+1)
       
       idegn=idegm+idegn
       if (idegn.lt.ipold)goto 10
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,ipold+1
       do jk=1,ipol(jf,2)+2
       isarray(jf,jk)=ipol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=ipold
       ncony=ncony+1
       ncong=1
       call subbw4(idegs,idegb,idegg)
!       if (ncony.eq.3)goto 6666
!      print *,'big','idegg',idegg
!       if (idegg.eq.0)goto 6666

       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
6      do jf=2,nn(2)+2
       if (nn(jf).lt.ie(jf))goto 4
       if (nn(jf).gt.ie(jf))goto 12
       end do
       goto 12

10     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 6
12     do jf=1,nn(2)+2
       karr(jf)=nn(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nn(jf)=kcarr(jf)
       end do





!       print *,'coming thru','idegn',idegn
       do jf=1,idegp+1
       do jk=1,igg1(jf,2)+2
       mult1(jf,jk)=igg1(jf,jk)
       end do
       end do
       idegm=idegp
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       call multy(idegm,idegn)
!       print *,'mult3',((mult3(jf,jk),jk=1,mult3(jf,2)+2),jf=1,&
!       idegm+idegn+1)
       
       idegn=idegm+idegn
       
       

       if (idegn.lt.ipold)goto 20
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,ipold+1
       do jk=1,ipol(jf,2)+2
       isarray(jf,jk)=ipol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=ipold
       ncony=ncony+1
       ncong=2
       call subbw4(idegs,idegb,idegg)
!       print *,'small','idegg',idegg
!       if (ncony.eq.3)goto 6666
!       if (idegg.eq.0)goto 6667
       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
       goto 4
6666   a=a
!6666   print *,'bigstop','idegm',idegm,'idegn',idegn,'ncong',ncong       
!       print *,'mult1',((mult1(jf,jk),jk=1,mult1(jf,2)+2),jf=1,&
!       idegm+1)
!       print *,'mult2',((mult2(jf,jk),jk=1,mult2(jf,2)+2),jf=1,&
!       idegn+1)
       
       
       stop
6667   print *,'smallstop'
       stop
20     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 4
30     a=a
!30     print *,'igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
       if (idegn.eq.0)goto 1300
       do jf=1,igg2(idegn,2)+2
       karr(jf)=igg2(idegn,jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       if (kcarr(1).eq.0)goto 32
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       igg2(idegn,jf)=kcarr(jf)
       end do
       goto 78
1300   idegn=1
       do jf=1,igg2(1,2)+2
       igg2(2,jf)=igg2(1,jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       igg2(1,jf)=kcarr(jf)
       end do
32     do jf=1,kcarr(2)+2 
       igg2(idegn,jf)=kcarr(jf)
       end do
78     a=a
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
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
       nnh(jf)=mdarr(jf)
       end do

       ihalf(1)=0
       ihalf(2)=1
       
       ihalf(3) =1
37     do jf=2,ihalf(2)+2
       if (ihalf(jf).lt.nnh(jf))goto 79
       if (ihalf(jf).gt.nnh(jf))goto 38
       end do
79     do jf=1,ihalf(2)+2 
       karr(jf)=ihalf(jf)
       kbarr(jf)=ihalf(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ihalf(jf)=kcarr(jf)
       end do
       goto 37
38     do jf=1,ihalf(2)+2
       marr(jf)=ihalf(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ihalf(jf)=mdarr(jf)
       end do
       do jf=1,nnh(2)+2
       karr(jf)=nnh(jf)
       end do
       
       
       do jf=1,ihalf(2)+2
       
       kbarr(jf)=ihalf(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nnh(jf)=kcarr(jf)
       end do

       





       
       
       do i=1,idegn+1
       if (igg2(i,2).eq.0)goto 33
       goto 34
33     end do
!       print *,'gcd',((ipol(jf,jk),jk=1,ipol(jf,2)+2),jf=1,ipold+1)
!       print *,'stopping where expected'
       
       jpold=ipold
       do jf=1,ipold+1
       do jk=1,ipol(jf,2)+2
       jpol(jf,jk)=ipol(jf,jk)
       igcd(jf,jk)=ipol(jf,jk)
       igarray(jf,jk)=ipol(jf,jk)
       end do
       end do
       igcdd=ipold
       

!       if (ipold.eq.1)goto 100
       if (ipold.eq.2)goto 200
       if (ipold.ge.3)goto 300
!       if (ipold.eq.4)goto 400
!       if (ipold.eq.5)goto 500
       stop
34     idegs=idegn+1-i       
       print *,'idegs',idegs
       do jf=1,idegs+1
       do jk=1,igg2(jf+idegn-idegs,2)+2
       isarray(jf,jk)=igg2(jf+idegn-idegs,jk)
       end do
       end do
       do jf=1,ipold+1
       do jk=1,ipol(jf,2)+2
       ibarray(jf,jk)=ipol(jf,jk)
       end do
       end do
       idegb=ipold
       call subgcd(idegs,idegb,idegg)
!       print *,'idegg from subgcd',idegg
       
!       print *,'gcd',((igarray(jf,jk),jk=1,igarray(jf,2)+2),jf=1,idegg+1)
       
       
       
       jpold=idegg
       do jf=1,idegg+1
       do jk=1,igarray(jf,2)+2
       jpol(jf,jk)=igarray(jf,jk)
       igcd(jf,jk)=igarray(jf,jk)
       end do
       end do
       igcdd=idegg
!       call split(jpold,jpold2)
       if (idegg.eq.1)goto 100
       if (idegg.eq.2)goto 200
       if (idegg.ge.3)goto 300
!       if (idegg.eq.4)goto 400
!       if (idegg.eq.5)goto 500
       print *,'polynomial no solution in integers'


       stop
110    print *,'1 solution',(mcarr(jf),jf=1,mcarr(2)+2)
       stop
100    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       
       do jf=1,igarray(1,2)+2
       karb(jf)=igarray(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,igarray(2,2)+2
       mbarr(jf)=igarray(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
!       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
       goto 1000








200    call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
!       print *,'nfcon',nfcon,'sols',isol(1),isol(2)
       goto 1000
300    a=a
!       print *,'here we are at 300 jpold',jpold
       
       call split(jpold,jpold2)
!       print *,'jpold',jpold,'jpold2',jpold2
!       stop
       if (jpold2.ne.1)goto 606
!       if (jpold-jpold2.gt.2)goto 606
!       print *,'jpold2',jpold2,'jpold',jpold
       do i=1,jpold+1
!       print *,'i',i,'jpol',(jpol(i,jf),jf=1,jpol(i,2)+2)
       end do
       do i=1,jpold2+1
!       print *,'i',i,'jpol2',(jpol2(i,jf),jf=1,jpol2(i,2)+2)
       end do
       do i=1,3
!       print *,'i',i,'iqt',(iqt(i,jf),jf=1,iqt(i,2)+2)
       end do
       
606    a=a       
!       print *,'jpold',jpold,'jpold2',jpold2
       
       if ((jpold.eq.4).and.(jpold2.eq.2))goto 420
       if ((jpold.eq.3).and.(jpold2.eq.2))goto 350
       if ((jpold.eq.3).and.(jpold2.eq.1))goto 610
       if ((jpold.gt.4).and.(jpold2.eq.2))goto 620
       if ((jpold.ge.4).and.(jpold2.eq.1))goto 630
       if (jpold-jpold2.eq.2)goto 640
       if (jpold-jpold2.eq.1)goto 650
       if (jpold2.gt.jpold-jpold2)goto 660
       jpold=jpold2
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       goto 300
660    do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       jpold=jpold-jpold2
       goto 300
640    a=a
       do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       idegg=jpold-jpold2
       call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       nfcon=nfcon+2
       idegb=igcdd
       do jf=1,idegb+1
       do jk=1,igcd(jf,2)+2
       ibarray(jf,jk)=igcd(jf,jk)
       end do
       end do
       idegs=jpold-jpold2
       do jf=1,idegs+1
       do jk=1,iqt(jf,2)+2
       isarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       call subbw4(idegs,idegb,idegg)
       igcdd=idegb-idegs
       do jf=1,igcdd+1
       do jk=1,iqt(jf,2)+2
       igcd(jf,jk)=iqt(jf,jk)
       end do
       end do
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       jpold=jpold2
       goto 300
620    do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       jpold=jpold-jpold2
       idegb=igcdd
       do jf=1,igcdd+1
       do jk=1,igcd(jf,2)+2
       ibarray(jf,jk)=igcd(jf,jk)
       end do
       end do
       idegs=jpold2
       do jf=1,idegs+1
       do jk=1,jpol2(jf,2)+2
       isarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw4(idegs,idegb,idegg)
       igcdd=idegb-idegs
       do jf=1,igcdd+1
       do jk=1,iqt(jf,2)+2
       igcd(jf,jk)=iqt(jf,jk)
       end do
       end do
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       igarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       nfcon=nfcon+2
       goto 300
630    a=a
       jpold=jpold-jpold2
       do jf=1,jpold+1
       do jk=1,iqt(jf,2)+2
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       karb(jf)=igarray(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,igarray(2,2)+2
       mbarr(jf)=igarray(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
!       print *,'igarray1',(igarray(1,jf),jf=1,igarray(1,2)+2)
!       print *,'igarray2',(igarray(2,jf),jf=1,igarray(2,2)+2)
!       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
!       print *,'stop630'
       
       idegb=igcdd
       do jf=1,igcdd+1
       do jk=1,igcd(jf,2)+2
       ibarray(jf,jk)=igcd(jf,jk)
       end do
       end do
       idegs=jpold2
       do jf=1,idegs+1
       do jk=1,jpol2(jf,2)+2
       isarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw4(idegs,idegb,idegg)
       igcdd=idegb-idegs
       do jf=1,igcdd+1
       do jk=1,iqt(jf,2)+2
       igcd(jf,jk)=iqt(jf,jk)
       end do
       end do
       goto 300


650    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,iqt(1,2)+2
       karb(jf)=iqt(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,iqt(2,2)+2
       mbarr(jf)=iqt(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
!       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
!       print *,'stop650'
       
       idegb=igcdd
       do jf=1,igcdd+1
       do jk=1,igcd(jf,2)+2
       ibarray(jf,jk)=igcd(jf,jk)
       end do
       end do
       idegs=jpold-jpold2
       do jf=1,idegs+1
       do jk=1,iqt(jf,2)+2
       isarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       call subbw4(idegs,idegb,idegg)
       igcdd=idegb-idegs
       do jf=1,igcdd+1
       do jk=1,iqt(jf,2)+2
       igcd(jf,jk)=iqt(jf,jk)
       end do
       end do
       jpold=jpold2
       do jf=1,jpold+1
       do jk=1,jpol2(jf,2)+2
       jpol(jf,jk)=jpol2(jf,jk)
       end do
       end do
       do jf=1,jpold
!       print *,'jpold',jpold,'jpol',(jpol(jf,jk),jk=1,jpol(jf,2)+2)
       end do
!       stop

       goto 300
610    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       karb(jf)=igarray(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,igarray(2,2)+2
       mbarr(jf)=igarray(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
!       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
!       print *,'nfcon1','sol',(nfsol(1,jf),jf=1,nfsol(1,2)+2)
!       stop
       do jf=1,jpold-jpold2+1
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       idegg=jpold-jpold2
       call subbw2(idegd,idegr)
324    do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
       if (igcdd.eq.jpold)goto 1000
       idegb=igcdd
       do jf=1,idegb+1
       do jk=1,igcd(jf,2)+2
       ibarray(jf,jk)=igcd(jf,jk)
       end do
       end do
       idegs=jpold
       do jf=1,idegs+1
       do jk=1,jpol(jf,2)+2
       isarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       call subbw4(idegs,idegb,idegg)
       jpold=idegb-idegs
       igcdd=jpold
       do jf=1,idegb+1
       do jk=1,iqt(jf,2)+2
       igcd(jf,jk)=iqt(jf,jk)
       jpol(jf,jk)=iqt(jf,jk)
       end do
       end do
       goto 300

       
350    a=a
       do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,iqt(1,2)+2
       karb(jf)=iqt(1,jf)
       end do
       call mpgcd
       do jf=1,karv(2)+2
       marr(jf)=karv(jf)
       end do
       do jf=1,iqt(2,2)+2
       mbarr(jf)=iqt(2,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       nfcon=nfcon+1
       do jf=1,kcarr(2)+2
       nfsol(nfcon,jf)=kcarr(jf)
       end do
!       print *,'nfcon',nfcon,'sol',(nfsol(nfcon,jf),jf=1,nfsol(nfcon,2)+2)
!       if (issz.eq.1)goto 309
       
       
       
!       print *,'nfcon',nfcon,'nfsol',nfsol(nfcon),'ibbsw',ibbsw
       
       
       do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       igarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       call subbw2(idegd,idegr)
       goto 324


420    do jf=1,jpold2+1
       do jk=1,jpol2(jf,2)+2
       igarray(jf,jk)=jpol2(jf,jk)
       end do
       end do
       do i=1,3
!       print *,'iqt',(iqt(i,jf),jf=1,iqt(i,2)+2)
       end do
!       stop
       call subbw2(idegd,idegr)
       do jf=1,2
       do jk=1,isol(jf,2)+2
       nfsol(nfcon+jf,jk)=isol(jf,jk)
       end do
       end do
       
       nfcon=nfcon+2
!       print *,'iqt1',(iqt(1,jf),jf=1,iqt(1,2)+2)
       do jf=1,3
       do jk=1,iqt(jf,2)+2
       igarray(jf,jk)=iqt(jf,jk)
       end do
       end do
       idegg=2
       call subbw2(idegd,idegr)
!       print *,'jpold',jpold,'jpold2',jpold2
       
       
       goto 324

1000   a=a
      print *,'nfcon',nfcon,'sols',((nfsol(jf,jk),jk=1,nfsol(jf,2)+2),&       
       jf=1,nfcon)
       
              
       end
       subroutine split(jpold,jpold2)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)

       
       
       
       dimension n(100),ipol(100,50),igg1(200,50),igg2(200,50),iapar(50)
       dimension ie(50),nn(50)
       
       
       idegp=1
       iapar(1)=0
       iapar(2)=1
       iapar(3)=1
61     igg1(1,1)=0
       igg1(1,2)=1
       igg1(1,3)=1
       do jf=1,iapar(2)+2
       igg1(2,jf)=iapar(jf)
       igg2(2,jf)=iapar(jf)
       end do
       
       igg2(1,1)=0
       igg2(1,2)=1
       igg2(1,3)=1
       
       idegn=1
       do jf=1,ihalf(2)+2
       ie(jf)=ihalf(jf)
       end do
       do jf=1,nnh(2)+2
       nn(jf)=nnh(jf)
       end do
       
       
!       print *,'ihalf',(ihalf(jf),jf=1,ihalf(2)+2)
!       print *,'nnh',(nnh(jf),jf=1,nnh(2)+2)
!       print *,'ie',(ie(jf),jf=1,ie(2)+2)
!       if (iapar(3).ne.2)goto 4
       
       
4      if ((ie(2).eq.1).and.(ie(3).eq.1))goto 30
!       print *,'ie',ie(3),'nn',nn(3)
       do jf=1,ie(2)+2
       marr(jf)=ie(jf)
       end do
       mbarr(1)=0
       mbarr(2)=1
       mbarr(3)=2
       call mendiv
       do jf=1,mdarr(2)+2
       ie(jf)=mdarr(jf)
       end do

       
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult1(jf,jk)=igg2(jf,jk)
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       idegm=idegn
       call multy(idegm,idegn)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 10
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       isarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg)
!       print *,'big'
       
       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
6      do jf=2,nn(2)+2
       if (nn(jf).lt.ie(jf))goto 4
       if (nn(jf).gt.ie(jf))goto 12
       end do
       goto 12

10     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 6
12     do jf=1,nn(2)+2
       karr(jf)=nn(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       nn(jf)=kcarr(jf)
       end do





!       print *,'coming thru'
       do jf=1,idegp+1
       do jk=1,igg1(jf,2)+2
       mult1(jf,jk)=igg1(jf,jk)
       end do
       end do
       idegm=idegp
       do jf=1,idegn+1
       do jk=1,igg2(jf,2)+2
       mult2(jf,jk)=igg2(jf,jk)
       end do
       end do
       call multy(idegm,idegn)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 20
       do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       ibarray(jf,jk)=mult3(jf,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       isarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg)
!       print *,'small'
       do jf=1,idegg+1
       do jk=1,irarray(jf,2)+2
       igg2(jf,jk)=irarray(jf,jk)
       end do
       end do
       idegn=idegg
       goto 4
20     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 4
30     a=a
!30     print *,'igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
       do jf=1,igg2(idegn+1,2)+2
       karr(jf)=igg2(idegn+1,jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       if (kcarr(1).eq.0)goto 32
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
32     do jf=1,kcarr(2)+2        
       igg2(idegn+1,jf)=kcarr(jf)
       end do







       
       do i=1,idegn+1
       
       if (igg2(i,2).eq.0)goto 33
       goto 34
33     end do
!       print *,'split gcd',(jpol(jf),jf=1,jpold+1)
!       iapar=iapar+1
!       if (iapar.gt.ipd)goto 71
       do jf=1,iapar(2)+2
       karr(jf)=iapar(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iapar(jf)=kcarr(jf)
       end do
       
       goto 61
71     print *,'problem in split'
       stop
34     idegs=idegn+1-i       
!       print *,'idegs',idegs
       do jf=1,idegs+1
       do jk=1,igg2(jf+idegn-idegs,2)+2
       isarray(jf,jk)=igg2(jf+idegn-idegs,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       ibarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=jpold
       call subgcd(idegs,idegb,idegg)
!       print *,'split gcd',((igarray(jf,jk),jk=1,igarray(jf,2)+2),&
!       jf=1,idegg+1),'idegg',idegg,'idegb',idegb
       
!       iapar=iapar+1
!       if (iapar.gt.ipd)goto 41
       do jf=1,iapar(2)+2
       karr(jf)=iapar(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(0)
       do jf=1,kcarr(2)+2
       iapar(jf)=kcarr(jf)
       end do
       if (iapar(3).lt.500)goto 42
       print *,'isind=1'
       stop 
42     if ((idegg.eq.0).or.(idegg.eq.jpold))goto 61       
       do jf=1,idegg+1
       do jk=1,igarray(jf,2)+2
       isarray(jf,jk)=igarray(jf,jk)
       jpol2(jf,jk)=igarray(jf,jk)
       end do
       end do
       jpold2=idegg
       
       idegs=idegg

       call subbw4(idegs,idegb,idegg)
       
!       print *,'iquo',(iqt(jf),jf=1,idegb-idegs+1)
       
       
       return
       end
       subroutine menmul
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)

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
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
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



       subroutine subbw2(idegd,idegr)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
       
       
       
       dimension nnumb1(50),nnumb2(50),nnumb3(50)
       dimension nnumb4(50)
       dimension ninv(50),iprar(50)
       dimension ipol1(6,25)
       


19     ia0(1)=0
       ia0(2)=1
       ia0(3)=201
       ia0(4)=7448
       ia0(5)=5699
       ia1(1)=0
       ia1(2)=1
       ia1(3)=301
       ia1(4)=3728
       ia1(5)=8717
       ia2(1)=0
       ia2(2)=1
       ia2(3)=401
       ia2(4)=9533
       ia2(5)=7210
       ia3(1)=0
       ia3(2)=1
       ia3(3)=501
       ia3(4)=8587
       ia3(5)=2630
       ia4(1)=0
       ia4(2)=1
       ia4(3)=600
       ia4(4)=5901
       ia4(5)=9090
       ia5(1)=0
       ia5(2)=1
       ia5(3)=1
       ia5(4)=8984
       ia5(5)=4842
       do jf=1,ipd(2)+2
       iprar(jf)=ipd(jf)
       end do
       do jf=1,igarray(1,2)+2
       nnumb1(jf)=igarray(1,jf)
       end do
       do jf=1,igarray(2,2)+2
       nnumb2(jf)=igarray(2,jf)
       end do
       do jf=1,igarray(3,2)+2
       nnumb3(jf)=igarray(3,jf)
       end do
!       print *,'nnumb1',(nnumb1(jk),jk=1,nnumb1(2)+2)
!       print *,'nnumb2',(nnumb2(jk),jk=1,nnumb2(2)+2)
!       print *,'nnumb3',(nnumb3(jk),jk=1,nnumb3(2)+2)






312    do jf=1,iprar(2)+2
       kara(jf)=iprar(jf)
       end do
       do jf=1,igarray(1,2)+2
       karr(jf)=igarray(1,jf)
       kbarr(jf)=igarray(1,jf)
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
       karb(jf)=mcarr(jf)
       end do






       
       
!       itemp=2*igarray(1)
316    call mpgcd
       do jf=1,karv(2)+2
       ninv(jf)=karv(jf)
       end do
       do jf=1,nnumb2(2)+2
       marr(jf)=nnumb2(jf)
       mbarr(jf)=nnumb2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       nnumb4(jf)=mcarr(jf)
       end do
       
       
!       print *,'nnumb4',(nnumb4(jk),jk=1,nnumb4(2)+2)
!       print *,'nnumb1',(nnumb1(jk),jk=1,nnumb1(2)+2)
       do jf=3,nnumb1(2)+2
       karr(jf-2)=nnumb1(jf)
       end do
       ilen=nnumb1(2)
       do jf=3,nnumb3(2)+2
       kbarr(jf-2)=nnumb3(jf)
       end do
       ilen2=nnumb3(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf)=kcarr(jf)
       end do
       ilen2=ilen3
       karr(1)=4
       ilen=1
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(2)=ilen3
       kbarr(1)=0

       
       do jf=1,nnumb4(2)+2
       karr(jf)=nnumb4(jf)
       end do
       call mpadd(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       isgn=kcarr(1)
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)


       
       
       do jf=1,irlen
       ncom(1,jf+2)=irrr(jf)
       end do
       ncom(1,1)=0
       ncom(1,2)=irlen
       if (isgn.eq.0)goto 318
       do jf=2,ncom(1,2)+2
       karr(jf)=ncom(1,jf)
       end do
       karr(1)=1


       do jf=1,iprar(2)+2
       kbarr(jf)=iprar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ncom(1,jf)=kcarr(jf)
       end do
318    ncom(2,1)=0
       ncom(2,2)=0
       a=a
!318    print *,'ncom',(ncom(1,jk),jk=1,ncom(1,2)+2)
       call cornsq
       
       do jf=1,igarray(2,2)+2
       kbarr(jf)=igarray(2,jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,isqurar(1,1,2)+2
       kbarr(jf)=isqurar(1,1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       if (mcarr(1).eq.0)goto 100
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       goto 102
100    do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
102    do jf=1,ninv(2)+2
       mbarr(jf)=ninv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       isol(1,jf)=mcarr(jf)
       end do
!       print *,'isol1',(isol(1,jf),jf=1,isol(1,2)+2)
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,igarray(2,2)+2
       kbarr(jf)=igarray(2,jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,isqurar(1,1,2)+2
       kbarr(jf)=isqurar(1,1,jf)
       end do
       call mpadd(1)
       
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       
       if (mcarr(1).eq.0)goto 150
       do jf=1,mcarr(2)+2
       karr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       kbarr(jf)=ipd(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       goto 152
150    do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
152    do jf=1,ninv(2)+2
       mbarr(jf)=ninv(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       isol(2,jf)=mcarr(jf)
       end do
!       print *,'isol2',(isol(2,jf),jf=1,isol(2,2)+2)
       
       
       

       
       
400    a=a              
       return
       end
       subroutine subgcd(idegs,idegb,idegg)
       
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       dimension iws(200,50),iwb(200,50),itempb(200,50),mul(50),ib(50)
!       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
!       jf=1,idegb+1)
!       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
!       jf=1,idegs+1)

       
       do i= 1,idegs +1
       do jf=1,isarray(i,2)+2
       
       iws(i,jf) = isarray(i,jf)
       end do
       end do
       do i=1,idegb+1
       do jf=1,ibarray(i,2)+2
       iwb(i,jf) = ibarray(i,jf)
       end do
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
!       print *,'iwbd',iwbd,'iwsd',iwsd
       do jf=1,iws(1,2)+2
       ib(jf)=iws(1,jf)

       end do
       if (ib(2).ne.0)goto 101
       print *,'worries'
       stop
101    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,ib(2)+2
       karb(jf)=ib(jf)
       end do
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       
       

501    do i =1,loopl
       
       if (iwb(i,2).eq.0)goto 60
       if (mul(2).ne.0)goto 5019
       print *,'worries'
       stop
5019   do jf=3,iwb(i,2)+2
       karr(jf-2)=iwb(i,jf)
       end do
       ilen=iwb(i,2)
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ipd(2)+2
       kbarr(jf-2)=ipd(jf)
       end do
       ilen2=ipd(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       iqt(i,jf+2)=irrr(jf)
       end do
       iqt(i,1)=0
       iqt(i,2)=irlen
       
       
509    iwb(i,1)=0
       iwb(i,2)=0
       
       
       
       

       do j =2,iwsd +1
       do jf=1,iws(j,2)+2
       marr(jf)=iws(j,jf)
       end do
       do jf=1,iqt(i,2)+2
       mbarr(jf)=iqt(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,iwb(i+j-1,2)+2
       karr(jf)=iwb(i+j-1,jf)
       end do
       call mpadd(1)
       if (kcarr(1).eq.0)goto 20
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
       
20     do jf=1,kcarr(2)+2
       iwb(i+j-1,jf)=kcarr(jf)
       end do
       end do
60     end do
       do i=1,iwbd+1
       if (iwb(i,2).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       do jj=1,itempbd
       do jf=1,iwb(i+jj-1,2)+2
       itempb(jj,jf)=iwb(i+jj-1,jf)
       end do
       end do
       do i=1,iwsd+1
       do jf=1,iws(i,2)+2
       iwb(i,jf)=iws(i,jf)
       end do
       end do
       iwbd=iwsd
       do jj=1,itempbd
       do jf=1,itempb(jj,2)+2
       iws(jj,jf)=itempb(jj,jf)
       end do
       end do
       iwsd=itempbd-1
       goto 10
100    idegg=iwsd
       do i=1,iwsd+1
       do jf=1,iws(i,2)+2
       igarray(i,jf)=iws(i,jf)
       end do
       end do
       return
       end

       subroutine multy(idegm,idegn)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)

       
       do i=1,idegm+idegn+1
       do jf=1,50
       mult3(i,jf) =0
       end do
       end do
        
       do i =1,idegm +1
       do j =1,idegn+1
       
       do jf=1,mult1(i,2)+2
       marr(jf)=mult1(i,jf)
       end do
       do jf=1,mult2(j,2)+2
       mbarr(jf)=mult2(j,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,mult3(i+j-1,2)+2
       karr(jf)=mult3(i+j-1,jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       marr(jf)=kcarr(jf)
       end do
       
       call mendiv
       do jf=1,mcarr(2)+2
       mult3(i+j-1,jf)=mcarr(jf)
       end do
       end do
       end do


!       mult3(i+j-1)=mult3(i+j-1) +mult1(i) *mult2(j)
!       mult3(i+j-1) =mod(mult3(i+j-1),ipd)
!       print *,'ok','idegm',idegm,'idegn',idegn
       return
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
       subroutine subbw4(idegs,idegb,idegg)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
       dimension iws(200,50),iwb(200,50),mul(50),ib(50)
!       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
!       jf=1,idegb+1)
!       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
!       jf=1,idegs+1)
       
       
       
       do i= 1,idegs +1
       do jf=1,isarray(i,2)+2
       
       iws(i,jf) = isarray(i,jf)
       end do
       end do
       do i=1,idegb+1
       do jf=1,ibarray(i,2)+2
       iwb(i,jf) = ibarray(i,jf)
       end do
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       do jf=1,iws(1,2)+2
       ib(jf)=iws(1,jf)

       end do
       if (ib(2).ne.0)goto 101
       print *,'worries'
       stop
101    do jf=1,ipd(2)+2
       kara(jf)=ipd(jf)
       end do
       do jf=1,ib(2)+2
       karb(jf)=ib(jf)
       end do
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       
       

501    do i =1,loopl
       
       if (iwb(i,2).eq.0)goto 60
       if (mul(2).ne.0)goto 5019
       print *,'worries'
       stop
5019   do jf=3,iwb(i,2)+2
       karr(jf-2)=iwb(i,jf)
       end do
       ilen=iwb(i,2)
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,ipd(2)+2
       kbarr(jf-2)=ipd(jf)
       end do
       ilen2=ipd(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,irlen
       iqt(i,jf+2)=irrr(jf)
       end do
       iqt(i,1)=0
       iqt(i,2)=irlen
!       print *,'iqti',i,(iqt(i,jf),jf=1,iqt(i,2)+2)
       
509    iwb(i,1)=0
       iwb(i,2)=0
       
       
       
       

       do j =2,iwsd +1
       do jf=1,iws(j,2)+2
       marr(jf)=iws(j,jf)
       end do
       do jf=1,iqt(i,2)+2
       mbarr(jf)=iqt(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       marr(jf)=mcarr(jf)
       end do
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
       end do
       call mendiv
       do jf=1,mcarr(2)+2
       kbarr(jf)=mcarr(jf)
       end do
       do jf=1,iwb(i+j-1,2)+2
       karr(jf)=iwb(i+j-1,jf)
       end do
       call mpadd(1)
       if (kcarr(1).eq.0)goto 20
       do jf=1,kcarr(2)+2
       kbarr(jf)=kcarr(jf)
       end do
       do jf=1,ipd(2)+2
       karr(jf)=ipd(jf)
       end do
       call mpadd(0)
       
20     do jf=1,kcarr(2)+2
       iwb(i+j-1,jf)=kcarr(jf)
       end do
       end do
       goto 601
60     iqt(i,1)=0
       iqt(i,2)=0
601    end do
       do i=1,iwbd+1
       if (iwb(i,2).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       
      do jj =1,itempbd
      do jf=1,iwb(jj+i-1,2)+2
      irarray(jj,jf) =iwb(jj+i-1,jf)
      end do
      end do
      idegr = itempbd -1
      goto 110

100   idegr =0
110   a=a
!     print *,'quos',(iqt(jf),jf=1,loopl)
!      print *,'rem',((irarray(jf,jk),jk=1,irarray(jf,2)+2),jf=1,idegr+1)
      idegg=idegr
      return
      end
      
       

       
       
      
      


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
      
      
      
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
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
      
      
      dimension kdum(12000),isub(12000)
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





      subroutine mpadd(isora)
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
      
      
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
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
      
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
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
      
      
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
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
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
       
      
          
       

       subroutine cornsq
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       
     
       
       dimension iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iprecod(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       do jf=1,ipd(2)+2
       ip(jf)=ipd(jf)
       end do
       
!       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
!       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
!       print *,'ip',(ip(jf),jf=1,ip(2)+2)
!       stop
       if ((ncom(1,2).eq.1).and.(ncom(1,3).eq.0))goto 232
       if (ncom(1,2).eq.0)goto 232
       if (ncom(1,3).eq.0)goto 232
       
       iconz=1
       nn=2
       do jbig=1,2
       if (ncom(jbig,2).eq.0)goto 1104
       do jf=3,ncom(jbig,2)+2
       karr(jf-2)=ncom(jbig,jf)
       end do
       ilen=ncom(jbig,2)
       
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       do jf=1,irlen
       iaa(jbig,jf+2)=irrr(jf)
       end do
       iaa(jbig,1)=0
       iaa(jbig,2)=irlen
       goto 1103
1104   iaa(jbig,1)=0       
       iaa(jbig,2)=0
1103   end do
       

       
1108   do i=1,2
       do j=1,3
       ix(i,j)=0
       end do
       end do
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       
       
       
!       do jf=3,ip(2)+2
!       karr(jf-2)=ip(jf)
!       kbarr(jf-2)=ip(jf)
!       end do
!       ilen=ip(2)
!       ilen2=ilen
!       call mpmul(ilen,ilen2,ilen3)
!       do jf=1,ilen3
!       karr(jf+2)=kcarr(jf)
!       end do
!       karr(1)=0
!       karr(2)=ilen3
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=1
       call mpadd(1)
       do jf=1,kcarr(2)+2
       itot(jf)=kcarr(jf)
       iprecod(jf)=kcarr(jf)
       end do
!       print *,'firprecod',(iprecod(jf),jf=1,iprecod(2)+2)
       
       do ibig=1,200
       do jf=3,iprecod(2)+2
       karr(jf-2)=iprecod(jf)
       end do
       ilen=iprecod(2)
       
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (icont.eq.0)goto 84
       if (irlen.ne.0)goto 84
       do jf=1,icont
       iprecod(jf+2)=ipqt(jf)
       end do
       iprecod(1)=0
       iprecod(2)=icont
       end do
       print *,'loop too short' 
       stop
84     ivan=iprecod(2)+2
       irem=mod(iprecod(ivan),2)
!       print *,'precod',(iprecod(jf),jf=1,iprecod(2)+2),'irem',irem
       
       
       do jf=1,iprecod(2)+2
       karr(jf)=iprecod(jf)
       end do
       kbarr(1)=0
       kbarr(2)=1
       kbarr(3)=2-irem
       call mpadd(0)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ipn(jf+2)=ipqt(jf)
       end do
       ipn(1)=0
       ipn(2)=icont
!       print *,'ok1',' ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 1600
98     do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ix(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
!       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2),'ix2',(ix(2,jf),jf=1,&
!       ix(2,2)+2)
!       print *,'ipn',(ipn(jf),jf=1,ipn(2)+2)
       
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       do ibig=1,2
       do jf=1,iaa(ibig,2)+2
       iaas(ibig,jf)=iaa(ibig,jf)
       end do
       end do
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       ib(ibig,jf)=ibprod(ibig,jf)
       ibperm(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ixperm(ibig,jf)=ix(ibig,jf)
       end do
       end do
132    nn=nn+1        
       nn=mod(nn,5003)
       
       if (nn.eq.0)goto 232
!       if (ip(2).gt.1)goto 403
403    iaas(1,1)=0
       iaas(1,2)=1
       iaas(1,3)=nn
       kard(1)=0
       kard(2)=1
       kard(3)=nn
       do jf=1,ip(2)+2
       karp(jf)=ip(jf)
       end do
       call mpkron(k)
       if (k.ne.-1)goto 132
       
!       nn=nn*607
!       nn=mod(nn,1000)
!       iaas(2,1)=0
!       iaas(2,2)=1
!       iaas(2,3)=nn
       iaas(2,1)=0
       iaas(2,2)=0
       do jf=1,iprecod(2)+2
       ipn(jf)=iprecod(jf)
       end do
       if (ibsw.eq.1)goto 1620
112    call sub516
!       print *,'ok4'
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iz(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 172
       goto 174
172    if (ib(2,2).eq.0)goto 228
174    icon=1
173    do ibig=1,2
       do jf=1,ib(ibig,2)+2
       iacn(ibig,jf)=ib(ibig,jf)
       end do
       do jf=1,iz(ibig,2)+2
       ibcn(ibig,jf)=iz(ibig,jf)
       end do
       end do
!       print *,'ok5',' icon',icon,'ib1',(ib(1,jf),jf=1,ib(1,2)+2)
!       print *,'ok5','ib2',(ib(2,jf),jf=1,ib(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ib(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if ((ib(1,2).eq.1).and.(ib(1,3).eq.1))goto 192
       goto 194
192    if (ib(2,2).eq.0)goto 200
194    icon=icon+1
       if (icon.eq.4500)goto 280
       goto 173
200    if (mod(icon,2).eq.1)goto 280
       do ibig=1,2
       do jf=1,iz(ibig,2)+2
       iaas(ibig,jf)=iz(ibig,jf)
       end do
       end do
       ipn(1)=0
       ipn(2)=1
       ipn(3)=icon/2
       call sub516
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       do jf=1,ix(ibig,2)+2
       ibcn(ibig,jf)=ix(ibig,jf)
       end do
       end do
!       print *,'iacn',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',iacn(2,3)
!       print *,'ibcn',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',(ibcn(2,jf&
!       ),jf=1,ibcn(2,2)+2)
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ix(ibig,jf)=icprod(ibig,jf)
       end do
       end do
228    a=a
!228    print *,'square root=',(ix(1,jf),jf=1,ix(1,2)+2),'ipn',ipn(3)
!       print *,'square root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       isqurar(1,ibig,jf)=ix(ibig,jf)
       end do
       end do
       nsq=1
       
       goto 300
280    a=a
!280    print *,'no square root exists',' iconz',iconz
       
       do ibig=1,2
       do jf=1,ixperm(ibig,2)+2
       ix(ibig,jf)=ixperm(ibig,jf)
       end do
!       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2)
       
       
       do jf=1,ibperm(ibig,2)+2
       ib(ibig,jf)=ibperm(ibig,jf)
       end do
       end do
       iconz=iconz+1
       if (iconz.eq.6000)goto 232
       goto 132
1600   ipn(1)=0
       ipn(2)=1
       ipn(3)=2
       ibsw=1
       goto 98
1620   ipn(1)=0
       ipn(2)=1
       ipn(3)=3
       goto 112
232    nsq=0
       do i=1,3
!       print *,'i',i,'igarray',(igarray(i,jf),jf=1,igarray(i,2)+2)
       end do
       print *,'ncom1',(ncom(1,jf),jf=1,ncom(1,2)+2),'nn',nn
       print *,'probably not prime',(ip(jf),jf=1,ip(2)+2)
       stop
300    return
       end








       
       subroutine sub400(id,ipz,k)
       ip=ipz
       ide =int(id/2)
       ipe= int(ip/2)
       ide2 =id -ide *2
       ipe2 =ip -ipe*2
       if(ide2.eq.0)goto 414
       goto 416
414    if(ipe2.eq.0)goto 512
416    iv = 0
       ipe = ip
       ii = 0
419    a=a
!419    print *,'pefirst',ipe
       ipe2 =int(ipe/2)
       ipe3 =ipe -ipe2 *2
       if(ipe3.eq.1)goto 432
       ipe =ipe2
       iv = iv+1
       ii =ii +1
       if(ii.lt.50)goto 419
432    ive = int(iv/2)
       ive2 =iv -ive*2
       if(ive2.eq.0)goto 450
       iae =(id **2 -1)/8
       iae1 =int(iae/2)
       iae2 =iae -iae1 *2
       iae3 =iae2 +2
       k=(-1)**iae3
       goto 451
450    k =1
451    ide = id
452    if(ide.eq.0)goto 510 
       ive = 0
       ii = 0
456    ide1 =int(ide/2)
       ide2 =ide -ide1 *2
       if(ide2.eq.1)goto 470
       iv = iv+1
       ide =ide1
       ii =ii+1
       if(ii.lt.50)goto 456
470    ive =int(iv/2)
       ive2=iv -ive *2
       if(ive2.eq.0)goto 486
       iae=(ipe **2 -1)/8
       iae2 =int(iae/2)
       iae3 =iae -iae2 *2
       iae3 = iae3 +2
       
       
       
       
       k = (-1)**iae3*k
486    iae2 =((ide -1)*(ipe-1))/4       
       iae3 =int(iae2/2)
       iae4 =iae2-iae3 *2
       iae4 =iae4 +2
       k =(-1) **iae4 *k
       ir=abs(ide)
       
       itemp=int(ipe/ir)
       ide =ipe -itemp*ir
       ipe = ir
       
       goto 452
510    if(ipe.eq.1)goto 513
512    k =0
513    a=a
!513    print *,k,ipe
       return
       end
       
       
       
       
       
       
       
                   
       
       
       
          
       
       
       
       subroutine sub516
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)

       
       

       
       dimension ipn2(100),ie(100),igg(2,100)


       if ((ipn(2).eq.1).and.(ipn(3).eq.1))goto 300
       goto 301
300    do ibig=1,2
       do jf=1,iaas(ibig,2)+2
       icprod(ibig,jf)=iaas(ibig,jf)
       ibprod(ibig,jf)=iaas(ibig,jf)
       end do
       end do
       goto 620

301    ie(1)=0
       ilen=1
       karr(1)=2
       ilen2=1
       kbarr(1)=2
4      call mpmul(ilen,ilen2,ilen3)
       if (ilen3.lt.ipn(2))goto 22
       if (ilen3.gt.ipn(2))goto 20
       do i=1,ilen3
       if (ipn(i+2).lt.kcarr(i))goto 20
       if (ipn(i+2).gt.kcarr(i))goto 22
       end do
       goto 22
20     ie(2)=ilen3       
       do i=1,ilen3
       ie(i+2)=kcarr(i)
       end do
       goto 30
22     do i=1,ilen3       
       kbarr(i)=kcarr(i)
       end do
       ilen2=ilen3
       goto 4
30     do i=1,ipn(2)+2
       ipn2(i)=ipn(i)
       end do
       do i=1,ie(2)
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
!       print *,'ie',(ie(jf),jf=1,ie(2)+2)
       
       do i=1,ipn2(2)+2
       karr(i)=ipn2(i)
       end do
       do i=1,ie(2)+2
       kbarr(i)=ie(i)
       end do
       call mpadd(1)
       do i=1,kcarr(2)+2
       ipn2(i)=kcarr(i)
       end do
       do ibig=1,2
       do jf=1,iaas(ibig,2)+2
       igg(ibig,jf)=iaas(ibig,jf)
       ibprod(ibig,jf)=iaas(ibig,jf)
       end do
       end do
31     if ((ie(2).eq.1).and.(ie(3).eq.1))goto 620
       do jf=3,ie(2)+2
       karr(jf-2)=ie(jf)
       end do
       ilen=ie(2)
       ilen2=1
       kbarr(1)=2
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,icont
       ie(jf+2)=ipqt(jf)
       end do
       ie(2)=icont
!       print *,'ok2'
       
!       if (ipn2(2).lt.ie(2))goto 31
!       if (ipn2(2).gt.ie(2))goto 42
!       do jf=3,ipn2(2)+2
!       if (ipn2(jf).lt.ie(jf))goto 31
!       if (ipn2(jf).gt.ie(jf))goto 42
!       end do
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       ibcn(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ibprod(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       if (ipn2(2).lt.ie(2))goto 31
       if (ipn2(2).gt.ie(2))goto 42
       do jf=3,ipn2(2)+2
       if (ipn2(jf).lt.ie(jf))goto 31
       if (ipn2(jf).gt.ie(jf))goto 42
       end do
42     a=a       
       
       do jf=1,ipn2(2)+2
       karr(jf)=ipn2(jf)
       end do
       do jf=1,ie(2)+2
       kbarr(jf)=ie(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       ipn2(jf)=kcarr(jf)
       end do
       do ibig=1,2
       do jf=1,ibprod(ibig,2)+2
       iacn(ibig,jf)=ibprod(ibig,jf)
       end do
       end do
       do ibig=1,2
       do jf=1,igg(ibig,2)+2
       ibcn(ibig,jf)=igg(ibig,jf)
       end do
       end do
       call sub1100
       do ibig=1,2
       do jf=1,icprod(ibig,2)+2
       ibprod(ibig,jf)=icprod(ibig,jf)
       end do
       end do
       goto 31
620    a=a
! 620    print *,'ok3'
       return
       end


       subroutine sub1100
       common ibarray(200,50),isarray(200,50),igarray(200,50),inv(25)
       common mult1(1200,50),mult2(1200,50),mult3(1200,50),mult4(1200,50)
       common mult5(1200,50)
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
       common nsq,isqurar(2,2,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       do jf=1,ipd(2)+2
       ip(jf)=ipd(jf)
       end do
       
       if ((iacn(1,2).eq.0).or.(ibcn(1,2).eq.0))goto 4
       do jf=3,iacn(1,2)+2
       karr(jf-2)=iacn(1,jf)
       end do
       ilen=iacn(1,2)
       do jf=3,ibcn(1,2)+2
       kbarr(jf-2)=ibcn(1,jf)
       end do
       ilen2=ibcn(1,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       icprod(1,jf+2)=kcarr(jf)
       end do
       icprod(1,1)=0
       icprod(1,2)=ilen3
       goto 41
4      icprod(1,1)=0 
       icprod(1,2)=0
41     if ((iacn(2,2).eq.0).or.(ibcn(2,2).eq.0))goto 42

       do jf=3,iacn(2,2)+2
       karr(jf-2)=iacn(2,jf)
       end do
       ilen=iacn(2,2)
       do jf=3,ibcn(2,2)+2
       kbarr(jf-2)=ibcn(2,jf)
       end do
       ilen2=ibcn(2,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=ilen3
       do jf=1,icprod(1,2)+2
       karr(jf)=icprod(1,jf)
       end do
       goto 51
42     do jf=1,icprod(1,2)+2       
       kcarr(jf)=icprod(1,jf)
       end do
       goto 52

51     call mpadd(1)
52     if (kcarr(2).eq.0)goto 1
       isgn=kcarr(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       if (ip(2).eq.0)goto 999
       if (ip(3).eq.0)goto 999
       
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1
       if (isgn.eq.0)goto 2
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=isgn
       kbarr(2)=irlen
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       icprod(1,jf)=kcarr(jf)
       end do
       goto 3
1      icprod(1,1)=0
       icprod(1,2)=0
       goto 3
2      do jf=1,irlen
       icprod(1,jf+2)=irrr(jf)
       end do
       icprod(1,1)=isgn
       icprod(1,2)=irlen
       
3      if ((iacn(1,2).eq.0).or.(ibcn(2,2).eq.0))goto 104
       do jf=3,iacn(1,2)+2
       karr(jf-2)=iacn(1,jf)
       end do
       ilen=iacn(1,2)
       do jf=3,ibcn(2,2)+2
       kbarr(jf-2)=ibcn(2,jf)
       end do
       ilen2=ibcn(2,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       icprod(2,jf+2)=kcarr(jf)
       end do
       icprod(2,1)=0
       icprod(2,2)=ilen3
       goto 141
104    icprod(2,1)=0
       icprod(2,2)=0
141    if ((iacn(2,2).eq.0).or.(ibcn(1,2).eq.0))goto 142       
       do jf=3,iacn(2,2)+2
       karr(jf-2)=iacn(2,jf)
       end do
       ilen=iacn(2,2)
       do jf=3,ibcn(1,2)+2
       kbarr(jf-2)=ibcn(1,jf)
       end do
       ilen2=ibcn(1,2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       kbarr(jf+2)=kcarr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=ilen3
       do jf=1,icprod(2,2)+2
       karr(jf)=icprod(2,jf)
       end do
       goto 151
142    do jf=1,icprod(2,2)+2
       kcarr(jf)=icprod(2,jf)
       end do
       goto 152
151    call mpadd(0)
152    if (kcarr(2).eq.0)goto 101
       isgn=kcarr(1)
       do jf=3,kcarr(2)+2
       karr(jf-2)=kcarr(jf)
       end do
       ilen=kcarr(2)
       do jf=3,ip(2)+2
       kbarr(jf-2)=ip(jf)
       end do
       ilen2=ip(2)
       if (ip(3).eq.0)goto 999
       if (ip(2).eq.0)goto 999
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 101
       if (isgn.eq.0)goto 102
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=isgn
       kbarr(2)=irlen
       do jf=1,ip(2)+2
       karr(jf)=ip(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       icprod(2,jf)=kcarr(jf)
       end do
       goto 200
999    print *,'ip',(ip(jf),jf=1,ip(2)+2)
       stop
101    icprod(2,1)=0
       icprod(2,2)=0
       goto 200
102    do jf=1,irlen
       icprod(2,jf+2)=irrr(jf)
       end do
       icprod(2,1)=isgn
       icprod(2,2)=irlen
200    a=a
! 200    print *,'icprod1',(icprod(1,jf),jf=1,icprod(1,2)+2)
!       print *,'icprod2',(icprod(2,jf),jf=1,icprod(2,2)+2)
!       print *,'iacn1',(iacn(1,jf),jf=1,iacn(1,2)+2),'iacn2',&
!       (iacn(2,jf),jf=1,iacn(2,2)+2)
!       print *,'ibcn1',(ibcn(1,jf),jf=1,ibcn(1,2)+2),'ibcn2',&
!       (ibcn(2,jf),jf=1,ibcn(2,2)+2)
       
       return
       end
       
       
       
       



      
      
