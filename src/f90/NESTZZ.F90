
       program nestzz
! recursive speed multiplication technique      

      
       
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
       
       
       
       dimension n(100),numb1(5000),numb2(5000),numb3(10000)
       dimension nn(50),ie(50),nnn(50),nnl(50)
       dimension iarray1(510000),iarray2(510000)
       dimension iprodd(50),iprod1(5000),iprod2(5000),iprod3(5000)
       dimension iee(50),iear(200,20)
       dimension jeet(20),jeetm(20)
      goto 160
      ilen=16
      ilen2=16
      do jf=1,ilen
      karr(jf)=jf+2
      kbarr(jf)=jf+3
      end do
      do i=1,10000
      call mpmul(ilen,ilen2,ilen3)
      end do
      print *,'ilen3',ilen3
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
      print *,'ilen3',ilen3
      stop


      
      
      goto 160
!       do jbig=1,1000       
       do ibig=1,1
       numb1(1)=0
       numb1(2)=2048
       numb1(2)=16
       do jf=3,2050
       numb1(jf)=1234
       end do
       numb2(1)=0
       numb2(2)=2048
       numb2(2)=16
       do jf=3,2050
       numb2(jf)=5678
       end do
      do jf=1,numb1(2)+2
      marr(jf)=numb1(jf)
      end do
      do jf=1,numb2(2)+2
      mbarr(jf)=numb2(jf)
      end do
      call menmul
      end do
      print *,'old mul',(mcarr(jf),jf=1,100)
      print *,'there now old'
      stop
      iaa=31
      ilen=298
      ilen2=298
      
      

       
       
       
       
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
!       print *,'i',i,'iear',(iear(i,jf),jf=1,iear(i,2)+2)
       end do
160    print *,'new mult. start'
       do jbig=1,1000       
!       do jbig=1,1
       numb1(1)=0
       numb1(2)=2048
!       numb1(2)=16
       do jf=3,2050
       numb1(jf)=jf
       end do
       numb2(1)=0
       numb2(2)=2048
!       numb2(2)=16
       do jf=3,2050
       numb2(jf)=jf+1
       end do
!       limb=2*(numb1(2)+2)+500 
       limb=2*(numb1(2)+2)
!       print *,'limb',limb
       do jf=1,limb
!        do jf=1,3000
        iarray1(jf)=0
        iarray2(jf)=0
        end do
        do jf=1,numb1(2)+2
        iarray1(jf)=numb1(jf)
        end do
        do jf=1,numb2(2)+2
        iarray1(numb1(2)+2+jf)=numb2(jf)
        end do
        isgn1=iarray1(1)
        ilen=iarray1(2)
        do jf=1,ilen
        karr(jf)=iarray1(2+jf)
        end do
        isgn2=iarray1(numb1(2)+3)
        ilen2=iarray1(numb1(2)+4)
        do jf=1,ilen2
        kbarr(jf)=iarray1(numb1(2)+4+jf)
        end do
       minj=(numb1(2))/2+2
!        print *,'ilen',ilen,'ilen2',ilen2
        call mpmul0(ilen,ilen2,ilen3,minj)
!       print *,'lengths',jnn1(2),jnn2(2),jnn3(2),jnn4(2),jnn5(2),jnn6(2)
!       print *,'jnn1',(jnn1(jf),jf=1,jnn1(2)+2)
!       print *,'jnn2',(jnn2(jf),jf=1,jnn2(2)+2)
!       print *,'jnn3',(jnn3(jf),jf=1,jnn3(2)+2)
!       print *,'jnn4',(jnn4(jf),jf=1,jnn4(2)+2)
!       print *,'jnn5',(jnn5(jf),jf=1,jnn5(2)+2)
!       print *,'jnn6',(jnn6(jf),jf=1,jnn6(2)+2)
       
       ind2=0
!       minj=(numb1(2))/2+2
!       minj=514
!      minj=4
       do jf=1,jnn1(2)+2
       iarray2(ind2+jf)=jnn1(jf)
       end do
       iarray2(ind2+1)=isgn1
       ind2=ind2+minj
       do jf=1,jnn2(2)+2
       iarray2(ind2+jf)=jnn2(jf)
       end do
       iarray2(ind2+1)=isgn2
       ind2=ind2+minj
       do jf=1,jnn3(2)+2
       iarray2(ind2+jf)=jnn3(jf)
       end do
       iarray2(ind2+1)=mod(isgn1+jnn3(1),2)
       ind2=ind2+minj
       do jf=1,jnn4(2)+2
       iarray2(ind2+jf)=jnn4(jf)
       end do
       iarray2(ind2+1)=mod(isgn2+jnn4(1),2)
       ind2=ind2+minj
       do jf=1,jnn5(2)+2
       iarray2(ind2+jf)=jnn5(jf)
       end do
       iarray2(ind2+1)=isgn1
       ind2=ind2+minj
       do jf=1,jnn6(2)+2
       iarray2(ind2+jf)=jnn6(jf)
       end do
       iarray2(ind2+1)=isgn2
!       print *,'first iarray2',(iarray2(jf),jf=1,200)
       
!       limb=(3*limb)/2
       limb=(8*limb)/5
!       print *,'limb',limb
       do jf=1,limb
!       do jf=1,3000
       iarray1(jf)=iarray2(jf)
       end do
       do jf=1,limb
!       do jf=1,3000
       iarray2(jf)=0
       end do


!  lots of 514
       kinj=3
       linj=minj
!       linj=514
!       linj=4
       minj=(linj-2)/2+2
90     do ibig=1,kinj
       ind1=(ibig-1)*linj*2
       if (iarray1(ind1+2).eq.0)goto 100
       if (iarray1(ind1+linj+2).eq.0)goto 100
       isgn1=iarray1(ind1+1)
       ilen=iarray1(ind1+2)
       do jf=1,ilen
       karr(jf)=iarray1(ind1+2+jf)
       end do
       isgn2=iarray1(ind1+linj+1)
       ilen2=iarray1(ind1+linj+2)
       do jf=1,ilen2
       kbarr(jf)=iarray1(ind1+linj+2+jf)
       end do
!       print *,'ilen',ilen,'ilen2',ilen2,'karr',(karr(jf),jf=1,ilen),&
!       'kbarr',(kbarr(jf),jf=1,ilen2)
       
       call mpmul0(ilen,ilen2,ilen3,minj)
!       print *,'jnn1',(jnn1(jf),jf=1,jnn1(2)+2),'linj',linj
!       print *,'jnn2',(jnn2(jf),jf=1,jnn2(2)+2),'ind1',ind1
!       print *,'jnn3',(jnn3(jf),jf=1,jnn3(2)+2),'ibig',ibig
!       print *,'jnn4',(jnn4(jf),jf=1,jnn4(2)+2)
!       print *,'jnn5',(jnn5(jf),jf=1,jnn5(2)+2)
!       print *,'jnn6',(jnn6(jf),jf=1,jnn6(2)+2)
!       print *,'lengths',jnn1(2),jnn2(2),jnn3(2),jnn4(2),jnn5(2),jnn6(2),&
!       'ibig',ibig
!       stop
       
       ind2=(ibig-1)*minj*6
       do jf=1,jnn1(2)+2
       iarray2(ind2+jf)=jnn1(jf)
       end do
       iarray2(ind2+1)=isgn1
!       print *,'iarray2 l',iarray2(ind2+2)
       ind2=ind2+minj
       do jf=1,jnn2(2)+2
       iarray2(ind2+jf)=jnn2(jf)
       end do
       iarray2(ind2+1)=isgn2
!       print *,'iarray2 l',iarray2(ind2+2)
       ind2=ind2+minj
       do jf=1,jnn3(2)+2
       iarray2(ind2+jf)=jnn3(jf)
       end do
       iarray2(ind2+1)=mod(isgn1+jnn3(1),2)
!       print *,'iarray2 l',iarray2(ind2+2)
       ind2=ind2+minj
       
       do jf=1,jnn4(2)+2
       iarray2(ind2+jf)=jnn4(jf)
       end do
       iarray2(ind2+1)=mod(isgn2+jnn4(1),2)
!       print *,'iarray2 l',iarray2(ind2+2)
       ind2=ind2+minj
       do jf=1,jnn5(2)+2
       iarray2(ind2+jf)=jnn5(jf)
       end do
       iarray2(ind2+1)=isgn1
!       print *,'iarray2 l',iarray2(ind2+2)
       ind2=ind2+minj
       do jf=1,jnn6(2)+2
       iarray2(ind2+jf)=jnn6(jf)
       end do
       iarray2(ind2+1)=isgn2
!       print *,'iarray2 l',iarray2(ind2+2)
!       stop
       goto 101
100    ind2=(ibig-1)*minj*6
       do i=1,6*minj
       iarray2(ind2+i)=0
       
       end do
!       print *,'zero on ibig',ibig
101    end do
!       print *,'jnn1',(jnn1(jf),jf=1,jnn1(2)+2)
!       print *,'jnn2',(jnn2(jf),jf=1,jnn2(2)+2)
!       print *,'jnn3',(jnn3(jf),jf=1,jnn3(2)+2)
!       print *,'jnn4',(jnn4(jf),jf=1,jnn4(2)+2)
!       print *,'jnn5',(jnn5(jf),jf=1,jnn5(2)+2)
!       print *,'jnn6',(jnn6(jf),jf=1,jnn6(2)+2)
!       print *,'second iarray2',(iarray2(jf),jf=1,300),'kinj',kinj,&
!       'linj',linj
       
       if (kinj.gt.200)goto 110
!       if (kinj.gt.3)goto 110
       kinj=kinj*3       
       linj=(linj-2)/2+2
       minj=(linj-2)/2+2
!       print *,'pre linj',linj,'minj',minj
!       if (kinj.gt.1000)goto 110
!       do jf=1,3000
!       iarray1(jf)=iarray2(jf)
!       end do
       limb=(8*limb)/5
!       print *,'limb',limb,'kinj',kinj
!       limb=(3*limb)/2
       do jf=1,limb
!       do jf=1,3000
       iarray1(jf)=0
       end do
       do jf=1,limb
       iarray1(jf)=iarray2(jf)
       end do
       do jf=1,limb
       iarray2(jf)=0
       end do
!       print *,'has been a funny'
       goto 90
110    a=a
       
       limb=(8*limb)/5
!       limb=(3*limb)/2
!       print *,'limb',limb,'kinj',kinj
       
       do jf=1,limb
!       do jf=1,3000
       iarray1(jf)=0
       end do
       
       ind1=0
       ind2=0

       linj=minj+minj-2
       kinj=kinj*3
       leng1=minj-2
       leng2=linj-2
       linj=(minj-2)*4+2
!       print *,'minj',minj,'linj',linj,'kinj',kinj,'leng1',leng1,&
!       'leng2',leng2
       
       do ibig=1,kinj
       do jf=1,iarray2(ind1+2)+2
       marr(jf)=iarray2(ind1+jf)
       end do
       do jf=1,iarray2(ind1+minj+2)+2
       mbarr(jf)=iarray2(ind1+minj+jf)
       end do
       
       call menmul
       if (mcarr(2).eq.0)goto 139
!       print *,'marr2',marr(2),'mbarr2',mbarr(2),'leng1',leng1,'leng2',&
!       leng2
       
139    do jf=1,mcarr(2)+2
       iprod1(jf)=mcarr(jf)
       end do
       ind1=ind1+minj+minj
       do jf=1,iarray2(ind1+2)+2
       marr(jf)=iarray2(ind1+jf)
       end do
       do jf=1,iarray2(ind1+minj+2)+2
       mbarr(jf)=iarray2(ind1+minj+jf)
       end do
       
       call menmul
       if (mcarr(2).eq.0)goto 138
!       print *,'marr2',marr(2),'mbarr2',mbarr(2),'leng1',leng1,'leng2',&
!       leng2
138    do jf=1,mcarr(2)+2
       iprod2(jf)=mcarr(jf)
       end do
       ind1=ind1+minj+minj
       do jf=1,iarray2(ind1+2)+2
       marr(jf)=iarray2(ind1+jf)
       end do
       do jf=1,iarray2(ind1+minj+2)+2
       mbarr(jf)=iarray2(ind1+minj+jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       iprod3(jf)=mcarr(jf)
       end do
       do jf=1,iprod1(2)+2
       karr(jf)=iprod1(jf)
       end do
       do jf=1,iprod2(2)+2
       kbarr(jf)=iprod2(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,iprod3(2)+2
       kbarr(jf)=iprod3(jf)
       end do
       call mpadd(0)
!       print *,'very first sum',(kcarr(jf),jf=1,kcarr(2)+2)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       
       if (karr(2).eq.0)goto 130
       lenj=karr(2)
       karr(2)=karr(2)+leng1
       
       do jf=1,leng1
       karr(lenj+2+jf)=0
       end do
       
130    do jf=1,iprod1(2)+2
       kbarr(jf)=iprod1(jf)
       end do
       if (kbarr(2).eq.0)goto 131
       lenj2=kbarr(2)
       kbarr(2)=kbarr(2)+leng2
       
       do jf=1,leng2
       kbarr(lenj2+2+jf)=0
       end do
131    call mpadd(0)
!       print *,'first karr',(karr(jf),jf=1,karr(2)+2)
!       print *,'first kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
       
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,iprod3(2)+2
       kbarr(jf)=iprod3(jf)
       end do
       call mpadd(0)
!       print *,'first sum',(kcarr(jf),jf=1,kcarr(2)+2)
       do jf=1,kcarr(2)+2
       iarray1(ind2+jf)=kcarr(jf)
       end do
       ind1=ind1+minj+minj
       ind2=ind2+linj
!       print *,'at end inds','ind1',ind1,'ind2',ind2
       end do
!       print *,'iarray1 after mult.',(iarray1(jf),jf=1,300),'ind2',ind2
!       print *,'prod1',(iprod1(jf),jf=1,iprod1(2)+2)
!       print *,'prod2',(iprod2(jf),jf=1,iprod2(2)+2)
!       print *,'prod3',(iprod3(jf),jf=1,iprod3(2)+2)
       
       limb=(5*limb)/7
!       print *,'limb',limb
!       stop
!       limb=(2*limb)/3
       do jf=1,limb
       iarray2(jf)=iarray1(jf)
       end do
       do jf=1,limb
       iarray1(jf)=0
       end do
!       print *,'prod1',(iprod1(jf),jf=1,iprod1(2)+2)
!       print *,'prod2',(iprod2(jf),jf=1,iprod2(2)+2)
!       print *,'prod3',(iprod3(jf),jf=1,iprod3(2)+2)
       
       
       
       
!       print *,'iarray1',(iarray1(jf),jf=1,100)
       
       
120    minj=linj
       linj=(minj-2)*2+2
       kinj=kinj/3
       leng1=2*leng1
       leng2=2*leng2
!       print *,'next ints minj',minj,'linj',linj,'kinj',kinj,'leng1',leng1,&
!       'leng2',leng2
!       stop
       ind1=0
       ind2=0
       do ibig=1,kinj/3
       do jf=1,iarray2(ind1+2)+2
       karr(jf)=iarray2(ind1+jf)
       end do
       do jf=1,iarray2(ind1+minj+2)+2
       kbarr(jf)=iarray2(ind1+minj+jf)
       end do
       call mpadd(0)
!       print *,' fir karr',(karr(jf),jf=1,karr(2)+2)
!       print *,'fir kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       do jf=1,iarray2(ind1+minj+minj+2)+2
       kbarr(jf)=iarray2(ind1+minj+minj+jf)
       end do
       call mpadd(0)
!       print *,' sec karr',(karr(jf),jf=1,karr(2)+2)
!       print *,'sec kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       if (karr(2).eq.0)goto 140
       lenj=karr(2)


       karr(2)=karr(2)+leng1
       
       do jf=1,leng1
       karr(lenj+2+jf)=0
       end do
140    ind1=ind1
       do jf=1,iarray2(ind1+2)+2
       kbarr(jf)=iarray2(ind1+jf)
       end do
       if (kbarr(2).eq.0)goto 141
       lenj2=kbarr(2)
       kbarr(2)=kbarr(2)+leng2
       
       do jf=1,leng2
       kbarr(lenj2+2+jf)=0
       end do
141    call mpadd(0)
!       print *,'th karr',(karr(jf),jf=1,karr(2)+2)
!       print *,'th kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
       
       
       do jf=1,kcarr(2)+2
       karr(jf)=kcarr(jf)
       end do
       
       ind1=ind1+minj+minj
       do jf=1,iarray2(ind1+2)+2
       kbarr(jf)=iarray2(ind1+jf)
       end do
       call mpadd(0)
!       print *,'karr',(karr(jf),jf=1,karr(2)+2)
!       print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
!       print *,'ssum',(kcarr(jf),jf=1,kcarr(2)+2)
       do jf=1,kcarr(2)+2
       iarray1(ind2+jf)=kcarr(jf)
       end do
       
       ind1=ind1+minj
       ind2=ind2+linj
       end do
!       print *,'next iarray1',(iarray1(jf),jf=1,200)
       
       limb=(limb*5)/7
!       print *,'limb',limb
       do jf=1,limb
       iarray2(jf)=iarray1(jf)
       end do
       do jf=1,limb
       iarray1(jf)=0
       end do
       if (kinj.gt.3)goto 120
       
       end do
!       print *,'answer',(iarray2(jf),jf=1,iarray2(2)+2)
       print *,'there now new'
       print *,'answer',(iarray2(jf),jf=1,100),'limb',limb
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)

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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      
      
      
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      
      
      
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
     
      

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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      
      
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      
      
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
      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
       
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
       
      
          
       

       


      
      





                    
      
      
      
      subroutine mpmul0(ilen,ilen2,ilen3,minj)
! input karr,ilen,kbarr,ilen2 ; output jnn1,jnn2,jnn3,jnn4,jnn5,jnn6      
! all arrays except ilen,ilen2      
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

      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
      dimension karrp(8000),kbarrp(8000),karra(2,8000),kbarra(2,8000)
      dimension ilena(2),ilena2(2),ibdf(8000)
      dimension iarr(8000),iadf(8000)
      ilenp=ilen
      ilenp2=ilen2
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      ilenas=minj-2
      goto 2
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
      
!      print *,'ilenp2',ilenp2,'ilenas',ilenas,'ilena21',ilena2(1)
      
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
!      print *,'jf',jf
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      jazzy=jf-1+jk
!      print *,'jk',jk,'jazzy',jazzy,'kbarr',kbarr(jazzy),'ilen2',ilen2
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
!      print *,'kbarr bef',(kbarr(jz),jz=1,ilen2)
!      print *,'found kbarra2',(kbarra(2,jz),jz=1,ilena2(2))
      
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
!      print *,'kbarr imp',(kbarr(jf),jf=1,ilenp2)


      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i)),&
!      'karra',(karra(i,jf),jf=1,ilena(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      jnn1(jf+2)=karra(1,jf)
      end do
      jnn1(1)=0
      jnn1(2)=ilena(1)
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      jnn2(jf+2)=kbarra(1,jf)
      end do
      jnn2(1)=0
      jnn2(2)=ilena2(1)
      ilen2=ilena2(1)
      isub=2
!      call mpmul2(ilen,ilen2,ilen3)
!      print *,'ilen3',ilen3
!      goto 540
      goto 545
535   iarr(1)=0
      iarr(2)=0
      jnn1(1)=0
      jnn1(2)=0
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
      if (kcarr(2).eq.0)goto 5521
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
5522  if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5521  jnn3(1)=0
      jnn3(2)=0
      
      jnn4(1)=0 
      jnn4(2)=0
      goto 552
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
      if (kcarr(2).eq.0)goto 5521
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      do jf=3,kcarr(2)+2
      jnn4(jf)=kcarr(jf)
      kbarr(jf-2)=kcarr(jf)
      end do 
      jnn4(1)=kcarr(1)
      jnn4(2)=kcarr(2)
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      jnn3(jf)=iadf(jf)
      end do
      jnn3(1)=iadf(1)
      jnn3(2)=iadf(2)
      ilen=iadf(2)
      isub=3
!      call mpmul2(ilen,ilen2,ilen3)
!      print *,'2ilen3',ilen3
      goto 552
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
      jnn5(jf+2)=karra(2,jf)
      end do
      jnn5(1)=0
      jnn5(2)=ilena(2)
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      jnn6(jf+2)=kbarra(2,jf)
      end do
      jnn6(1)=0
      jnn6(2)=ilena2(2)
      ilen2=ilena2(2)
      isub=4 
!      call mpmul2(ilen,ilen2,ilen3)
!      print *,'3ilen3',ilen3
!      print *,'kbarra',(kbarra(2,jf),jf=1,ilen2)
!      print *,'kbarrp',(kbarrp(jf),jf=1,ilen2)
      
      goto 110
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
 
      
      subroutine mpmul2(ilen,ilen2,ilen3)
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

      common jnn1(5000),jnn2(5000),jnn3(5000),jnn4(5000),jnn5(5000)
      common jnn6(5000)
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
      print *,'after 2'
      
      
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
      
      print *,'ilenp2',ilenp2,'ilenas',ilenas
      
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
      isubb=2
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
      isubb=3
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
      isubb=4 
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
      print *,'three ilen3',ilen3
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

100   if (isubb.eq.2)goto 540
      if (isubb.eq.3)goto 550
      if (isubb.eq.4)goto 560
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
      
