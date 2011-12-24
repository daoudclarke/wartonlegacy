       
       program bobhrxn4
!    single precision version of bobhrx4
!      program for rigorous testing of polynomials for primality
! used for very high degree polynomials modulo small p
!      involves powering to compute  GCD's of     
!      polynomials       
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
       
       
       
       dimension n(100),ipol(3000),igg1(3000),igg2(3000)
       dimension nfsol(10,50)
       dimension nn(200),ie(200),inpol(2000),ifact(100),ipow1(200)
       dimension indegs(2000),ipg(200)
       print *,'enter degree of polynomial and modulus and number of non'
       print *,'-zero coefficients'
       read *,ipold,ipd,nnz
       do jf=1,ipold+1
       inpol(jf)=0
       indegs(jf)=0
       end do
       igcon=1

       print *,'list degrees of non-zero terms and corresponding' 
       print *,' coefficients in pairs'
       read *,(indegs(jf),inpol(jf),jf=1,nnz)
       print *,'degs and inpols',(indegs(jf),inpol(jf),jf=1,nnz)
       open (unit=3,file='recl.dat',access='direct',form=&
       'formatted',recl=390000,status='old')
       read (3,5,rec=1)(ipr(jf),jf=1,65000)
5      format(65000i6)       
       print *,'primes',(ipr(jf),jf=2,100)
       nfact=0
       do jf=2,ipold
       if (ipr(jf).gt.ipold)goto 801
       idiv=mod(ipold,ipr(jf))
       if (idiv.ne.0)goto 80
       nfact=nfact+1
       ifact(nfact)=ipr(jf)
80     end do
801    print *,'nfact',nfact,'factors of degree',(ifact(jf),jf=1,nfact)
       
       ipg(1)=0
       ipg(2)=1
       ipg(3)=ipd
       do jf=1,ipold+1
       ipol(jf)=0
       
       end do
       do jf=1,nnz
       ind=indegs(jf)
       ipol(ind)=inpol(jf)
       
       end do
       print *,'ipol nz',(ipol(jf),jf=1,ipold+1)
       do ii=1,igcon
       ipol(ii)=1
       
       end do

93     a=a       
       do ibigl=1,nfact+1
       print *,'ibigl',ibigl
       iee=1
       if (ibigl.ne.1)goto 85
       nnpow=ipold
       
       goto 81
85     a=a     
       iee=1
       ibbbsw=2
       nnpow=ipold/ifact(ibigl-1)

81     if (iee.gt.nnpow)goto 82
       iee=iee*2
       goto 81
82     iee=iee/2
       npnz=nnpow-iee
       do jf=1,ipg(2)+2
       ipow1(jf)=ipg(jf)
       end do
83     if (iee.eq.1)goto 84         
       iee=iee/2
       do jf=1,ipow1(2)+2
       marr(jf)=ipow1(jf)
       mbarr(jf)=ipow1(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       ipow1(jf)=mcarr(jf)
       end do
       if (npnz.lt.iee)goto 83
       npnz=npnz-iee
       do jf=1,ipow1(2)+2
       marr(jf)=ipow1(jf)
       end do
       do jf=1,ipg(2)+2
       mbarr(jf)=ipg(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       ipow1(jf)=mcarr(jf)
       end do
       goto 83
84     print *,'ipow1',(ipow1(jf),jf=1,ipow1(2)+2)
       
       
       

       issz=0
       ncony=0
       nfcon=0
       ibbsw=0
       







       
41     a=a       
       do jf=1,ipow1(2)+2
       nn(jf)=ipow1(jf)
!       do jf=1,ipd(2)+2
!       nn(jf)=ipd(jf)
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
       print *,'ie',(ie(jf),jf=1,ie(2)+2)
       print *,'nn',(nn(jf),jf=1,nn(2)+2)
       
       
       
       idegp=1
       
       igg1(1)=1
       
       igg1(2)=0
       
       
       igg2(1)=1
       
       igg2(2)=0
       
       idegn=1
4      if ((ie(2).eq.1).and.(ie(3).eq.1))goto 30
       print *,'ie',(ie(jf),jf=1,ie(2)+2),'nn',(nn(jf),jf=1,nn(2)+2),&
       'igcon',igcon
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
       
       mult1(jf)=igg2(jf)
       mult2(jf)=igg2(jf)
       end do
       
       idegm=idegn
       call multy(idegm,idegn)
!       print *,'mult3',((mult3(jf,jk),jk=1,mult3(jf,2)+2),jf=1,&
!       idegm+idegn+1)
       
       idegn=idegm+idegn
       if (idegn.lt.ipold)goto 10
       do jf=1,idegn+1
       
       ibarray(jf)=mult3(jf)
       end do
       
       do jf=1,ipold+1
       
       isarray(jf)=ipol(jf)
       
       end do
       idegb=idegn
       idegs=ipold
       ncony=ncony+1
       ncong=1
!       print *,'ok1 ncony',ncony,'idegs',idegs,'idegb',idegb
       do jl=1,idegs+1
!       print *,'isarray',(isarray(jl,jf),jf=1,isarray(jl,2)+2)
       end do
       do jl=1,idegb+1
!       print *,'ibarray',(ibarray(jl,jf),jf=1,ibarray(jl,2)+2)
       end do
       call subbw4(idegs,idegb,idegg)
!       print *,'ok2 ncony',ncony
!       if (ncony.eq.3)goto 6666
!      print *,'big','idegg',idegg
!       if (idegg.eq.0)goto 6666

       do jf=1,idegg+1
       
       igg2(jf)=irarray(jf)
       end do
       
       idegn=idegg
6      do jf=2,nn(2)+2
       if (nn(jf).lt.ie(jf))goto 4
       if (nn(jf).gt.ie(jf))goto 12
       end do
       goto 12

10     do jf=1,idegn+1
       
       igg2(jf)=mult3(jf)
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
       
       mult1(jf)=igg1(jf)
       end do
       
       idegm=idegp
       do jf=1,idegn+1
       
       mult2(jf)=igg2(jf)
       end do
       
       call multy(idegm,idegn)
!       print *,'mult3',((mult3(jf,jk),jk=1,mult3(jf,2)+2),jf=1,&
!       idegm+idegn+1)
       
       idegn=idegm+idegn
       
       

       if (idegn.lt.ipold)goto 20
       do jf=1,idegn+1
       
       ibarray(jf)=mult3(jf)
       
       end do
       do jf=1,ipold+1
       
       isarray(jf)=ipol(jf)
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
       
       igg2(jf)=irarray(jf)
       
       end do
       idegn=idegg
       goto 4
6666   print *,'bigstop','idegm',idegm,'idegn',idegn,'ncong',ncong       
       print *,'mult1',(mult1(jf),jf=1,&
       idegm+1)
       print *,'mult2',(mult2(jf),jf=1,&
       idegn+1)
       
       
       stop
6667   print *,'smallstop'
       stop
20     do jf=1,idegn+1
       
       igg2(jf)=mult3(jf)
       end do
       
       goto 4
30     print *,'igg2',(igg2(jf),jf=1,idegn+1)
       igg2(idegn)=igg2(idegn)-1
       if (igg2(idegn).ge.0)goto 32
       igg2(idegn)=igg2(idegn)+ipd
       goto 78
32     a=a
78     a=a
       





       print *,'adj igg2',(igg2(jf),jf=1,idegn+1)
       
       
       do i=1,idegn+1
       if (igg2(i).eq.0)goto 33
       goto 91
33     end do
       print *,'gcd',(ipol(jf),jf=1,ipold+1)
       print *,'stopping where expected','ibigl',ibigl,'nfact',nfact
       if (ibigl.ne.1)goto 90
       
       goto 92
       

       


90     print *,'polynomial is composite',(ipol(jf),jf=1,ipold+1),&
       'mod',ipd,'degee',ipold
       
       igcon=igcon+1
       ipol(igcon)=1
       if(igcon.lt.ipold)goto 93
       
       
       
       
       
       stop
91     idegs=idegn+1-i       
       if (ibigl.eq.1)goto 90
       print *,'idegs',idegs
       do jf=1,idegs+1
       
       isarray(jf)=igg2(jf+idegn-idegs)
       
       end do
       print *,'isarray before',(isarray(jf),jf=1,idegs+1)
       
       do jf=1,ipold+1
       
       ibarray(jf)=ipol(jf)
       end do
       
       idegb=ipold
       call subgcd(idegs,idegb,idegg)
       print *,'idegg from subgcd',idegg
       if (idegg.ne.0)goto 90
!       print *,'gcd',(igarray(jf),jf=1,idegg+1)
92      end do 
      print *,'polynomial is prime',(ipol(jf),jf=1,ipold+1),& 
      'mod',ipd,'count',igcon,'degree',ipold 
       
1000   a=a
      
       
              
       end
       subroutine menmul
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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



       subroutine subgcd(idegs,idegb,idegg)
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
       
       dimension iws(3000),iwb(3000),itempb(3000)
       print *,'ibar',(ibarray(jf),&
       jf=1,idegb+1)
       print *,'isar',(isarray(jf),&
       jf=1,idegs+1)

       
       do i= 1,idegs +1
       
       
       iws(i) = isarray(i)
       end do
       
       do i=1,idegb+1
       
       iwb(i) = ibarray(i)
       
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       print *,'iwbd',iwbd,'iwsd',iwsd
       print *,'iws',(iws(k1),k1=1,iwsd+1)
       
       ib=iws(1)

       
       if (ib.ne.0)goto 101
       print *,'worries'
       stop
101    ia=ipd
       call subbw6(ia,ib,iv)
       
       

       mul=iv
       
       
       

501    do i =1,loopl
       
       if (iwb(i).eq.0)goto 60
       if (mul.ne.0)goto 5019
       print *,'worries'
       stop
5019   iqt(i)=iwb(i)*mul
       iqt(i)=mod(iqt(i),ipd)
       iwb(i)=0

       
       
       
       

       do j =2,iwsd +1
       iwb(i+j-1)=iwb(i+j-1)-iws(j)*iqt(i)
       iwb(i+j-1)=mod(iwb(i+j-1),ipd)
       if (iwb(i+j-1).ge.0)goto 5020
       iwb(i+j-1)=iwb(i+j-1)+ipd



       
5020   end do       
       

       
       
60     end do
       do i=1,iwbd+1
       if (iwb(i).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       do jj=1,itempbd
       
       itempb(jj)=iwb(i+jj-1)
       
       end do
       do i=1,iwsd+1
       
       iwb(i)=iws(i)
       
       end do
       iwbd=iwsd
       do jj=1,itempbd
       
       iws(jj)=itempb(jj)
       
       end do
       iwsd=itempbd-1
       goto 10
100    idegg=iwsd
       do i=1,iwsd+1
       
       igarray(i)=iws(i)
       
       end do
       return
       end

       subroutine multy(idegm,idegn)
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       

       
       do i=1,idegm+idegn+1
       
       mult3(i) =0
       end do
       
        
       do i =1,idegm +1
       do j =1,idegn+1
       mult3(i+j-1)=mult3(i+j-1)+mult1(i)*mult2(j)
       mult3(i+j-1)=mod(mult3(i+j-1),ipd)
       if (mult3(i+j-1).ge.0)goto 10
       mult3(i+j-1)=mult3(i+j-1)+ipd
10     end do
       end do
       



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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
       
       dimension iws(3000),iwb(3000),itempb(3000)
!       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
!       jf=1,idegb+1)
!       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
!       jf=1,idegs+1)
!       print *,'ibar',(ibarray(jf),&
!       jf=1,idegb+1)
!       print *,'isar',(isarray(jf),&
!       jf=1,idegs+1)

       
       do i= 1,idegs +1
       
       
       iws(i) = isarray(i)
       end do
       
       do i=1,idegb+1
       
       iwb(i) = ibarray(i)
       
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       print *,'iwbd',iwbd,'iwsd',iwsd
!       print *,'iws',(iws(k1),k1=1,iwsd+1)
       
       ib=iws(1)

       
       if (ib.ne.0)goto 101
       print *,'worries'
       stop
101    ia=ipd
       call subbw6(ia,ib,iv)
       
       

       mul=iv
       
       
       

501    do i =1,loopl
       
       if (iwb(i).eq.0)goto 60
       if (mul.ne.0)goto 5019
       print *,'worries'
       stop
5019   iqt(i)=iwb(i)*mul
       iqt(i)=mod(iqt(i),ipd)
       iwb(i)=0

       
       
       
       

       do j =2,iwsd +1
       iwb(i+j-1)=iwb(i+j-1)-iws(j)*iqt(i)
       iwb(i+j-1)=mod(iwb(i+j-1),ipd)
       if (iwb(i+j-1).ge.0)goto 5020
       iwb(i+j-1)=iwb(i+j-1)+ipd



       
5020   end do       
       

       
       
60     end do
       do i=1,iwbd+1
       if (iwb(i).ne.0)goto 30
       end do
       goto 100
30     itempbd=iwbd+2-i
       do jj=1,itempbd
       
       itempb(jj)=iwb(i+jj-1)
       
       end do
       iwbd=iwsd
       do jj=1,itempbd
       
       irarray(jj)=itempb(jj)
       
       end do
       idegr=itempbd-1
       goto 110
100    idegr=0
110    a=a
      
      idegg=idegr
      return
      end
      
       

       
       
      
      


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
      
      
      dimension kdum(600),isub(600)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
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
       common ibarray(3000),isarray(3000),igarray(3000),inv(25)
       common mult1(3000),mult2(3000),mult3(3000),mult4(3000)
       common mult5(3000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(3000),jpol(10,50),jpol2(10,50)
       common iqt(3000),ipd
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       

      
      
      
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
       
      
          
       

        
       
      
      
