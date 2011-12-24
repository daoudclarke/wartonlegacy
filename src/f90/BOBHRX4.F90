       
       program bobhrx4
!      program for rigorous testing of polynomials for primality
! used for very high degree polynomials modulo small p
!      involves powering to compute  GCD's of     
!      polynomials       
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
       
       
       
       dimension n(100),ipol(1100,3),igg1(1100,3),igg2(1100,3)
       dimension nfsol(10,50)
       dimension nn(200),ie(200),inpol(2000),ifact(100),ipow1(200)
       dimension indegs(2000)
       print *,'enter degree of polynomial and modulus and number of non'
       print *,'-zero coefficients'
       read *,ipold,ipd(3),nnz
       do jf=1,ipold+1
       inpol(jf)=0
       indegs(jf)=0
       end do
       igcon=nnz

       print *,'list degrees + 1 of non-zero terms and corresponding' 
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
       idiv=mod(ipold,ipr(jf))
       if (idiv.ne.0)goto 80
       nfact=nfact+1
       ifact(nfact)=ipr(jf)
80     end do
       print *,'nfact',nfact,'factors of degree',(ifact(jf),jf=1,nfact)
       
       ipd(1)=0
       ipd(2)=1
       do jf=1,ipold+1
       ipol(jf,1)=0
       ipol(jf,2)=0
       ipol(jf,3)=0
       end do
       do jf=1,nnz
       ind=indegs(jf)
       ipol(ind,3)=inpol(jf)
       ipol(ind,2)=1
       end do
       print *,'ipol nz',(ipol(jf,3),jf=1,ipold+1)
!       do ii=1,igcon
!       ipol(ii,2)=1
!       ipol(ii,3)=1
!       end do

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
       do jf=1,ipd(2)+2
       ipow1(jf)=ipd(jf)
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
       do jf=1,ipd(2)+2
       mbarr(jf)=ipd(jf)
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
       print *,'ok1 ncony',ncony,'idegs',idegs,'idegb',idegb
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
6666   print *,'bigstop','idegm',idegm,'idegn',idegn,'ncong',ncong       
       print *,'mult1',((mult1(jf,jk),jk=1,mult1(jf,2)+2),jf=1,&
       idegm+1)
       print *,'mult2',((mult2(jf,jk),jk=1,mult2(jf,2)+2),jf=1,&
       idegn+1)
       
       
       stop
6667   print *,'smallstop'
       stop
20     do jf=1,idegn+1
       do jk=1,mult3(jf,2)+2
       igg2(jf,jk)=mult3(jf,jk)
       end do
       end do
       goto 4
30     print *,'igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
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
32     do jf=1,kcarr(2)+2 
       igg2(idegn,jf)=kcarr(jf)
       end do
78     a=a
       





       print *,'adj igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
       
       do i=1,idegn+1
       if (igg2(i,2).eq.0)goto 33
       goto 91
33     end do
       print *,'gcd',((ipol(jf,jk),jk=1,ipol(jf,2)+2),jf=1,ipold+1)
       print *,'stopping where expected','ibigl',ibigl,'nfact',nfact
       if (ibigl.ne.1)goto 90
       
       goto 92
       

       


90     print *,'polynomial is composite',(ipol(jf,3),jf=1,ipold+1),&
       'mod',ipd(3),'degee',ipold
       igcon=igcon+1
       ipol(igcon,2)=1
       ipol(igcon,3)=1
       if (igcon.lt.ipold)goto 93
       
       
       
       
       
       stop
91     idegs=idegn+1-i       
       if (ibigl.eq.1)goto 90
       print *,'idegs',idegs
       do jf=1,idegs+1
       do jk=1,igg2(jf+idegn-idegs,2)+2
       isarray(jf,jk)=igg2(jf+idegn-idegs,jk)
       
       end do
       print *,jf,'isarray before',(isarray(jf,jr),jr=1,isarray(jf,2)+2)
       end do
       do jf=1,ipold+1
       do jk=1,ipol(jf,2)+2
       ibarray(jf,jk)=ipol(jf,jk)
       end do
       end do
       idegb=ipold
       call subgcd(idegs,idegb,idegg)
       print *,'idegg from subgcd',idegg
       if (idegg.ne.0)goto 90
!       print *,'gcd',(igarray(jf),jf=1,idegg+1)
92      end do 
      print *,'polynomial is prime',(ipol(jf,3),jf=1,ipold+1),& 
      'mod',ipd(3),'count',igcon,'degree',ipold 
       
1000   a=a
      
       
              
       end
       subroutine split(jpold,jpold2)
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)

       
       
       
       dimension n(100),ipol(10,50),igg1(1100,3),igg2(1100,3),iapar(50)
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
       
       
       print *,'ihalf',(ihalf(jf),jf=1,ihalf(2)+2)
       print *,'nnh',(nnh(jf),jf=1,nnh(2)+2)
       print *,'ie',(ie(jf),jf=1,ie(2)+2)
!       if (iapar(3).ne.2)goto 4
       
       
4      if ((ie(2).eq.1).and.(ie(3).eq.1))goto 30
       print *,'ie',ie(3),'nn',nn(3)
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
       print *,'big'
       
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





       print *,'coming thru'
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
       print *,'small'
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
30     print *,'igg2',((igg2(jf,jk),jk=1,igg2(jf,2)+2),jf=1,idegn+1)
       
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
       print *,'idegs',idegs
       do jf=1,idegs+1
       do jk=1,igg2(jf,2)+2
       isarray(jf,jk)=igg2(jf,jk)
       end do
       end do
       do jf=1,jpold+1
       do jk=1,jpol(jf,2)+2
       ibarray(jf,jk)=jpol(jf,jk)
       end do
       end do
       idegb=jpold
       call subgcd(idegs,idegb,idegg)
       print *,'split gcd',((igarray(jf,jk),jk=1,igarray(jf,2)+2),&
       jf=1,idegg+1),'idegg',idegg,'idegb',idegb
       
       iapar=iapar+1
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       dimension iws(1100,3),iwb(1100,3),itempb(1100,3),mul(100),ib(100)
       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
       jf=1,idegb+1)
       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
       jf=1,idegs+1)

       
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
       print *,'iwbd',iwbd,'iwsd',iwsd
       print *,'iws',((iws(k1,k2),k2=1,iws(k1,2)+2),k1=1,iwsd+1)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)

       
       do i=1,idegm+idegn+1
       do jf=1,3
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
       dimension iws(1100,3),iwb(1100,3),itempb(1100,3),mul(100),ib(100)
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
60     end do
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
      
      
      
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
       
      
          
       

       subroutine cornsq
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
       
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
       
!       common ix1(1,100),iy1(1,100),ix2(1,100),iy2(1,100),icc(1,100)
!       common n(100),jaa(100),jbb(100),kbcorn(200),kycorn(120)
       
     
       
       dimension itempz(2),iaa(2,100)
       dimension ix(2,100),itot(200),iz(2,100),ixperm(2,100)
       dimension iq(100),iprecod(100),iprem(100),ib(2,100),ibperm(2,100)                               
       ibsw=0
       do jf=1,ipd(2)+2
       ip(jf)=ipd(jf)
       end do
       
       print *,'ncom',(ncom(1,jf),jf=1,ncom(1,2)+2)
       print *,'ncom2',(ncom(2,jf),jf=1,ncom(2,2)+2)
       print *,'ip',(ip(jf),jf=1,ip(2)+2)
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
       iaas(1,1)=0
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
       if (icon.eq.50)goto 280
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
228    print *,'square root=',(ix(1,jf),jf=1,ix(1,2)+2),'ipn',ipn(3)
       print *,'square root comp=',(ix(2,jf),jf=1,ix(2,2)+2)
       do ibig=1,2
       do jf=1,ix(ibig,2)+2
       isqurar(1,ibig,jf)=ix(ibig,jf)
       end do
       end do
       nsq=1
       
       goto 300
280    print *,'no square root exists',' iconz',iconz
       do ibig=1,2
       do jf=1,ixperm(ibig,2)+2
       ix(ibig,jf)=ixperm(ibig,jf)
       end do
       print *,'ix',(ix(1,jf),jf=1,ix(1,2)+2)
       
       
       do jf=1,ibperm(ibig,2)+2
       ib(ibig,jf)=ibperm(ibig,jf)
       end do
       end do
       iconz=iconz+1
       if (iconz.eq.300)goto 232
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
419    print *,'pefirst',ipe
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
513    print *,k,ipe
       return
       end
       
       
       
       
       
       
       
                   
       
       
       
          
       
       
       
       subroutine sub516
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
       
       
       
       
       
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
       common ibarray(1100,3),isarray(1100,3),igarray(1100,3),inv(25)
       common mult1(1100,3),mult2(1100,3),mult3(1100,3),mult4(1100,3)
       common mult5(1100,3)
       
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(1100,3),jpol(10,50),jpol2(10,50)
       common iqt(1100,3),ipd(50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
       
       
       
       
       
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
       
       
       
       



      
      
