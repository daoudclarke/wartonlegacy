       program bobfacp
!      program for trying out different methods of poly. factorization
!      involves powering to compute  GCD's of     
!      polynomials. finds all roots of equations of degree 5 and less with       
!      modulus less than 44000
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       
       
       dimension n(100),ipol(100),igg1(100),igg2(100),nfsol(100)
       dimension inpos(100),innos(100),jnpos(100),jnnos(100)
       dimension inarray(5000),imarray(5000)
       print *,'enter type of run parameter,1 for inverse'
       print *,'    2 for integer solutions of equation'
       read *,intype
       if (intype.ne.1)goto 41
       print *,'input degree of polynomial to be inverted '
       print *,' then modulus,then no. of non-zero coefficients'
       read *,idegin,ipd,nnz
       print *,'input in pairs coeff. no. and cefficient'
       read *,(inpos(jf),innos(jf),jf=1,nnz)
       do i=1,idegin+1
       inarray(i)=0
       end do
       
       do i=1,nnz
       inz=inpos(i)
       inarray(inz)=innos(i)
       end do
       print *,'input degree of polynomial modulus'
       print *,'then no. of non-zero coefficients'
       read *,idegmd,nnz2
       print *,'input in pairs coeff. no. and coefficients'
       read *,(jnpos(jf),jnnos(jf),jf=1,nnz2)
       do i=1,idegmd+1
       imarray(i)=0
       end do
       do i=1,nnz2
       jnz=jnpos(i)
       
       imarray(jnz)=jnnos(i)
       end do
       if (idegin.le.idegmd)goto 800
       do i=1,idegin+1
       ibarray(i)=inarray(i)
       end do
       idegb=idegin
       do i=1,idegmd+1
       isarray(i)=imarray(i)
       end do
       idegs=idegmd
       call subbw4(idegs,idegb,idegg,ipd)
       if (idegg.eq.-1)goto 801
       do jf=1,idegg+1
       isarray(jf)=igarray(jf)
       end do
       idegs=idegg
       if (idegs.eq.0)goto 803
802    do jf=1,idegmd+1
       ibarray(jf)=imarray(jf)
       end do
       idegb=idegmd
       call subinv(idegs,idegb,idegg,ipd)
       stop
800    do jf=1,idegin+1
       isarray(jf)=inarray(jf)
       end do
       idegs=idegin
       if (idegs.gt.0)goto 802
803    call subbw6(ipd,isarray(1),iv)
       print *,'inverse=',iv
       stop
       


801    print *,'polynomial divides mod polynomial exactly'
       stop


       
       ipd=5
       idegb=3
       ibarray(1)=8
       ibarray(2)=8
       ibarray(3)=2
       ibarray(4)=1
       idegs=2
       isarray(1)=1
       isarray(2)=2
       isarray(3)=3
       call subinv(idegs,idegb,idegg,ipd)
       stop


       
       ipd=2
       idegb=6
       ibarray(1)=1
       ibarray(2)=1
       ibarray(3)=1
       ibarray(4)=0
       ibarray(5)=1
       ibarray(6)=0
       ibarray(7)=1
       idegs=4
       isarray(1)=1
       isarray(2)=0
       isarray(3)=0
       isarray(4)=1
       isarray(5)=1
       call subinv(idegs,idegb,idegg,ipd)
       stop

       
       
       ncony=0
       nfcon=0
       ibbsw=0
       ipd=701
       

!       ipd=40031
!       ipd=11677
!       call subbw6(11677,10000,iv)
!       print *,'iv',iv
!       stop


       goto 41
       goto 40
       idegs=3
       isarray(1)=657
       isarray(2)=247
       isarray(3)=436
       isarray(4)=346
       isarray(5)=353
       idegb=4
       ibarray(1)=297
       ibarray(2)=246
       ibarray(3)=269
       ibarray(4)=421
       ibarray(5)=353

       call subbw4(idegs,idegb,idegg,ipd)
       stop
       call subgcd(idegs,idegb,idegg,ipd)
       print *,'gcd',(igarray(jf),jf=1,idegg+1)
       stop


!       goto 40
41     a=a       
       print *,'enter degree of equation (less than 6)and greater than 1)'
       print *,' and modulus less than 10000' 
       read *,ipold,ipd
       print *,'enter all coefficients'
       
       
       
       read *,(ipol(jf),jf=1,ipold+1)
       do ij=1,ipold+1
       ipol(ij)=mod(ipol(ij),ipd)
       if (ipol(ij).ge.0)goto 805
       ipol(ij)=ipol(ij)+ipd
805    end do
       ncony=0
       nfcon=0
       ibbsw=0
       
       
       nn=ipd
       ie =1
2      if (ie.gt.nn)goto 3
       ie = ie *2
       goto 2
3      ie = ie/2
       
       nn = nn-ie
!       ipold=5
!       ipol(1)=1
!       ipol(2)=ipd-15
!       ipol(3)=85
!       ipol(4)=ipd-225
!       ipol(5)=274
!       ipol(6)=ipd-120
       idegp=1
       
       igg1(1)=1
       igg1(2)=0
       igg2(1)=1
       igg2(2)=0
       idegn=1
4      if (ie.eq.1)goto 30
!       print *,'ie',ie,'nn',nn
       ie=ie/2
       do jf=1,idegn+1
       mult1(jf)=igg2(jf)
       mult2(jf)=igg2(jf)
       end do
       idegm=idegn
       call multy(idegm,idegn,ipd)
       print *,'mult3',(mult3(jf),jf=1,idegm+idegn+1)
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
       call subbw4(idegs,idegb,idegg,ipd)
!       if (ncony.eq.2)goto 6666
!       print *,'big'
       do jf=1,idegg+1
       igg2(jf)=irarray(jf)
       end do
       idegn=idegg
6      if (nn.ge.ie)goto 12
       goto 4
10     do jf=1,idegn+1
       igg2(jf)=mult3(jf)
       end do
       goto 6
12     nn=nn-ie
!       print *,'coming thru'
       do jf=1,idegp+1
       mult1(jf)=igg1(jf)
       end do
       idegm=idegp
       do jf=1,idegn+1
       mult2(jf)=igg2(jf)
       end do
       call multy(idegm,idegn,ipd)
       print *,'mult3',(mult3(jf),jf=1,idegm+idegn+1)
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
       call subbw4(idegs,idegb,idegg,ipd)
!       if (ncony.eq.2)goto 6666
!       print *,'small'
       do jf=1,idegg+1
       igg2(jf)=irarray(jf)
       end do
       idegn=idegg
       goto 4
6666  print *,'idegm',idegm,'idegn',idegn,'ncong',ncong 
       print *,'mult1',(mult1(jf),jf=1,idegm+1)
       print *,'mult2',(mult2(jf),jf=1,idegn+1)
       stop
20     do jf=1,idegn+1
       igg2(jf)=mult3(jf)
       end do
       goto 4
30     print *,'igg2',(igg2(jf),jf=1,idegn+1)
       if (idegn.eq.0)goto 1300
       igg2(idegn)=igg2(idegn)-1
       igg2(idegn)=mod(igg2(idegn),ipd)
       if (igg2(idegn).ge.0)goto 32
       igg2(idegn)=igg2(idegn)+ipd
32     print *,'igg2el',igg2(idegn)
       
       do i=1,idegn+1
       if (igg2(i).eq.0)goto 33
       goto 34
33     end do
       goto 1301
1300   idegn=1
       igg2(2)=igg2(1)
       igg2(1)=ipd-1
       goto 32
1301   a=a
       print *,'gcd',(ipol(jf),jf=1,ipold+1)
       
       jpold=ipold
       do jf=1,ipold+1
       jpol(jf)=ipol(jf)
       end do
       nnh=(ipd-1)/2
       ihalf =1
37     if (ihalf.gt.nnh)goto 38
       ihalf = ihalf *2
       goto 37
38     ihalf = ihalf/2
       if (ipold.eq.1)goto 100
       if (ipold.eq.3)goto 300
       if (ipold.eq.4)goto 400
       if (ipold.eq.5)goto 500
       do jf=1,ipold+1
       igarray(jf)=ipol(jf)
       end do
       if (ipold.eq.2)goto 200
       stop
34     idegs=idegn+1-i       
       print *,'idegs',idegs
       do jf=1,idegs+1
       isarray(jf)=igg2(jf+idegn-idegs)
       end do
       do jf=1,ipold+1
       ibarray(jf)=ipol(jf)
       end do
       idegb=ipold
       call subgcd(idegs,idegb,idegg,ipd)
       print *,'gcd',(igarray(jf),jf=1,idegg+1)
       
       nnh=(ipd-1)/2
       ihalf =1
52      if (ihalf.gt.nnh)goto 53
       ihalf = ihalf *2
       goto 52
53      ihalf = ihalf/2
       jpold=idegg
       do jf=1,idegg+1
       jpol(jf)=igarray(jf)
       end do
       
       
!       call split(jpold,jpold2,ihalf,ipd)
       if (idegg.eq.1)goto 100
       if (idegg.eq.2)goto 200
       if (idegg.eq.3)goto 300
       if (idegg.eq.4)goto 400
       if (idegg.eq.5)goto 500
       print *,'polynomial irreducible'


       stop
100    if (ipd.lt.10000)goto 102 
       kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
       goto 104
102    kara(1)=0
       kara(2)=1
       kara(3)=ipd
104    if (igarray(1).lt.10000)goto 110
       karb(1)=0
       karb(2)=1
       karb(3)=igarray(1)/10000
       karb(4)=igarray(1)-karb(3)*10000
       goto 112
110    karb(1)=0
       karb(2)=1
       karb(3)=igarray(1)
112    call mpgcd
       if (igarray(2).lt.10000)goto 120
       ilen=2
       karr(1)=igarray(2)/10000
       karr(2)=igarray(2)-karr(1)*10000
       goto 122
120    ilen=1
       karr(1)=igarray(2)
122    do jf=3,karv(2)+2
       kbarr(jf-2)=karv(jf)
       end do
       ilen2=karv(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,kara(2)+2
       kbarr(jf-2)=kara(jf)
       end do
       ilen2=kara(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,kara(2)+2
       karr(jf)=kara(jf)
       end do
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=irlen
       call mpadd(1)
       nfcon=nfcon+1
       isum=0
       do jf=kcarr(2)+2,3,-1
       isum=isum+kcarr(jf)*10000**(jf-kcarr(2)-2)
       end do
       print *,'isum',isum,'kcarr',(kcarr(jf),jf=1,kcarr(2)+2)
       nfsol(nfcon)=isum
       
       print *,'nfcon',nfcon,'nfsol',nfsol(1)
       goto 1000
200    call subbw2(idegd,ipd,idegr)
       nfsol(nfcon+1)=isol(1)
       nfsol(nfcon+2)=isol(2)
       nfcon=nfcon+2
       print *,'nfcon',nfcon,'sols',isol(1),isol(2)
       goto 1000
300    a=a
       
       call split(jpold,jpold2,ihalf,ipd)
       print *,'jpold2',jpold2
       if (jpold2.eq.2)goto 350
301    if (jpold2.eq.3)goto 350
303    if (jpold2.eq.4)goto 350       
       if (ipd.lt.10000)goto 302 
       kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
       goto 304
302    kara(1)=0
       kara(2)=1
       kara(3)=ipd
304    if (igarray(1).lt.10000)goto 310
       karb(1)=0
       karb(2)=2
       karb(3)=igarray(1)/10000
       karb(4)=igarray(1)-karb(3)*10000
       goto 312
310    karb(1)=0
       karb(2)=1
       karb(3)=igarray(1)
312    call mpgcd
       if (igarray(2).lt.10000)goto 320
       ilen=2
       karr(1)=igarray(2)/10000
       karr(2)=igarray(2)-karr(1)*10000
       goto 322
320    ilen=1
       karr(1)=igarray(2)
322    do jf=3,karv(2)+2
       kbarr(jf-2)=karv(jf)
       end do
       print *,'karv',(karv(jf),jf=1,karv(2)+2)
       ilen2=karv(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,kara(2)+2
       kbarr(jf-2)=kara(jf)
       end do
       ilen2=kara(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       do jf=1,kara(2)+2
       karr(jf)=kara(jf)
       end do
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=irlen
       call mpadd(1)
       nfcon=nfcon+1
       isum=0
       do jf=kcarr(2)+2,3,-1
       isum=isum+kcarr(jf)*10000**(kcarr(2)+2-jf)
       
       end do
       print *,'isum',isum,'kcarr',(kcarr(jf),jf=1,kcarr(2)+2)
       nfsol(nfcon)=isum
       print *,'igarray1',igarray(1),'igarray2',igarray(2)
       print *,'nfcon',nfcon,'nfsol',nfsol(nfcon)
       
       if (ibbsw.eq.1)goto 410
       if (ibbsw.eq.2)goto 510
       do jf=1,jpold-jpold2+1
       igarray(jf)=iqt(jf)
       end do
       idegg=jpold-jpold2
       call subbw2(idegd,ipd,idegr)
324    nfsol(nfcon+1)=isol(1)
       nfsol(nfcon+2)=isol(2)
       nfcon=nfcon+2
       goto 1000
       
350    a=a
       if (ipd.lt.10000)goto 352 
       kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
       goto 354
352    kara(1)=0
       kara(2)=1
       kara(3)=ipd
354    if (iqt(1).lt.10000)goto 360
       karb(1)=0
       karb(2)=2
       karb(3)=iqt(1)/10000
       karb(4)=iqt(1)-karb(3)*10000
       goto 362
360    karb(1)=0
       karb(2)=1
       karb(3)=iqt(1)
362    call mpgcd
       if (iqt(2).lt.10000)goto 370
       ilen=2
       karr(1)=iqt(2)/10000
       karr(2)=iqt(2)-karr(1)*10000
       goto 372
370    ilen=1
       karr(1)=iqt(2)
372    do jf=3,karv(2)+2
       kbarr(jf-2)=karv(jf)
       end do
       print *,'karr',(karr(jf),jf=1,ilen)
       print *,'karv',(karv(jf),jf=1,karv(2)+2)
       ilen2=karv(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,kara(2)+2
       kbarr(jf-2)=kara(jf)
       end do
       ilen2=kara(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       print *,'kara',(kara(jf),jf=1,kara(2)+2)
       do jf=1,kara(2)+2
       karr(jf)=kara(jf)
       end do
       do jf=1,irlen
       kbarr(jf+2)=irrr(jf)
       end do
       kbarr(1)=0
       kbarr(2)=irlen
       call mpadd(1)
       print *,'kcarr',(kcarr(jf),jf=1,kcarr(2)+2)
       nfcon=nfcon+1
       isum=0
       do jf=kcarr(2)+2,3,-1
       isum=isum+kcarr(jf)*10000**(jf-kcarr(2)-2)
       end do
       nfsol(nfcon)=isum
       
       print *,'nfcon',nfcon,'nfsol',nfsol(nfcon),'ibbsw',ibbsw
       
       if (ibbsw.eq.1)goto 430
       if (ibbsw.eq.2)goto 530
       do jf=1,jpold2+1
       igarray(jf)=jpol2(jf)
       end do
       call subbw2(idegd,ipd,idegr)
       goto 324
400    call split(jpold,jpold2,ihalf,ipd)
       print *,'jpold2',jpold2
       
       if (jpold2.eq.2)goto 420
       ibbsw=1
       goto 301
410    ibbsw=0
       do jf=1,4
       jpol(jf)=iqt(jf)
       end do
       jpold=3
       goto 300

420    do jf=1,jpold2+1
       igarray(jf)=jpol2(jf)
       end do
       call subbw2(idegd,ipd,idegr)
       nfsol(nfcon+1)=isol(1)
       nfsol(nfcon+2)=isol(2)
       nfcon=nfcon+2
       do jf=1,3
       igarray(jf)=iqt(jf)
       end do
       idegg=2
       call subbw2(idegd,ipd,idegr)
       goto 324
430    do jf=1,jpold2+1
       jpol(jf)=jpol2(jf)
       end do
       jpold=jpold2
       ibbsw=0
       goto 300
500    call split(jpold,jpold2,ihalf,ipd)
       print *,'jpold2',jpold2
       
       if (jpold2.eq.2)goto 520
       if (jpold2.eq.3)goto 522
       ibbsw=2
       goto 301
510    ibbsw=0
       do jf=1,jpold2+1
       jpol(jf)=jpol2(jf)
       end do
       jpold=jpold2
       goto 400
520    call subbw2(idegd,ipd,idegr)
       nfsol(nfcon+1)=isol(1)
       nfsol(nfcon+2)=isol(2)
       nfcon=nfcon+2
       do jf=1,jpold-jpold2+1
       jpol(jf)=iqt(jf)
       end do
       jpold=3
       ibbsw=0
       goto 300
522    do jf=1,3
       igarray(jf)=iqt(jf)
       end do
       call subbw2(idegd,ipd,idegr)
       nfsol(nfcon+1)=isol(1) 
       nfsol(nfcon+2)=isol(2)
       nfcon=nfcon+2
       do jf=1,4
       jpol(jf)=jpol2(jf)
       end do
       jpold=jpold2
       ibbsw=0
       goto 300
530    ibbsw=0
       do jf=1,jpold2+1
       jpol(jf)=jpol2(jf)
       end do
       jpold=jpold2
       goto 400
40     a=a
       idegs=5
       isarray(1)=601
       isarray(2)=501
       isarray(3)=401
       isarray(4)=301
       isarray(5)=201
       isarray(6)=101
       ibarray(1)=1
       do jf=2,ipd+1
       ibarray(jf)=0
       end do
       ibarray(ipd)=ipd-1
       idegb=ipd
       call subgcd(idegs,idegb,idegg,ipd)
       print *,'gcd',(igarray(jf),jf=1,idegg+1)
!       call subbw4(idegs,ipd,idegg,ipd)
       stop
       ipd=5
       idegb=3
       ibarray(1)=1
       ibarray(2)=3
       ibarray(3)=3
       ibarray(4)=1
       idegs=2
       isarray(1)=1
       isarray(2)=2
       isarray(3)=1
       call subgcd(idegs,idegb,idegg,ipd)
       print *,'gcd',(igarray(jf),jf=1,idegg+1)
       
       stop
       call subbw2(5,ipd,idegr)
1000   print *,'nfcon',nfcon,'sols',(nfsol(jf),jf=1,nfcon)       
       
       
              
       end
       subroutine split(jpold,jpold2,ihalf,ipd)
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       
       
       dimension n(100),ipol(100),igg1(100),igg2(100)
       idegp=1
       iapar=1
61     igg1(1)=1
       igg1(2)=iapar
       igg2(1)=1
       igg2(2)=iapar
       idegn=1
       ie=ihalf
       print *,'ihalf',ihalf
       nn=(ipd-1)/2-ie
4      if (ie.eq.1)goto 30
       print *,'ie',ie,'nn',nn
       
       ie=ie/2
       do jf=1,idegn+1
       mult1(jf)=igg2(jf)
       mult2(jf)=igg2(jf)
       end do
       idegm=idegn
       call multy(idegm,idegn,ipd)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 10
       do jf=1,idegn+1
       ibarray(jf)=mult3(jf)
       end do
       do jf=1,jpold+1
       isarray(jf)=jpol(jf)
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg,ipd)
       print *,'big'
       
       do jf=1,idegg+1
       igg2(jf)=irarray(jf)
       end do
       idegn=idegg
6      if (nn.ge.ie)goto 12
       goto 4
10     do jf=1,idegn+1
       igg2(jf)=mult3(jf)
       end do
       goto 6
12     nn=nn-ie
       print *,'coming thru'
       do jf=1,idegp+1
       mult1(jf)=igg1(jf)
       end do
       idegm=idegp
       do jf=1,idegn+1
       mult2(jf)=igg2(jf)
       end do
       call multy(idegm,idegn,ipd)
       idegn=idegm+idegn
       if (idegn.lt.jpold)goto 20
       do jf=1,idegn+1
       ibarray(jf)=mult3(jf)
       end do
       do jf=1,jpold+1
       isarray(jf)=jpol(jf)
       end do
       idegb=idegn
       idegs=jpold
       call subbw4(idegs,idegb,idegg,ipd)
       print *,'small'
       do jf=1,idegg+1
       igg2(jf)=irarray(jf)
       end do
       idegn=idegg
       goto 4
20     do jf=1,idegn+1
       igg2(jf)=mult3(jf)
       end do
       goto 4
30     print *,'split igg2',(igg2(jf),jf=1,idegn+1)
       
       igg2(idegn+1)=igg2(idegn+1)-1
       igg2(idegn+1)=mod(igg2(idegn+1),ipd)
       if (igg2(idegn+1).ge.0)goto 32
       igg2(idegn+1)=igg2(idegn+1)+ipd
32     print *,'igg2el',igg2(idegn+1)
       do i=1,idegn+1
       if (igg2(i).eq.0)goto 33
       goto 34
33     end do
       print *,'split gcd',(jpol(jf),jf=1,jpold+1)
       iapar=iapar+1
       if (iapar.gt.ipd)goto 71
       goto 61
71     print *,'problem in split'
       stop
34     idegs=idegn+1-i       
       print *,'idegs',idegs
       do jf=1,idegs+1
       isarray(jf)=igg2(jf+idegn-idegs)
       end do
       do jf=1,jpold+1
       ibarray(jf)=jpol(jf)
       end do
       idegb=jpold
       call subgcd(idegs,idegb,idegg,ipd)
       print *,'split gcd',(igarray(jf),jf=1,idegg+1),'iapar',iapar
       
       iapar=iapar+1
       if (iapar.gt.ipd)goto 41
       goto 42
41     iapar=1 
42     if ((idegg.eq.0).or.(idegg.eq.jpold))goto 61       
       do jf=1,idegg+1
       isarray(jf)=igarray(jf)
       jpol2(jf)=igarray(jf)
       end do
       jpold2=idegg
       idegs=idegg
       call subbw4(idegs,idegb,idegg,ipd)
       
       print *,'iquo',(iqt(jf),jf=1,idegb-idegs+1)
       
       
       return
       end
       
       
       subroutine subbw2(idegd,ipd,idegr)
        
       
       common  ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common  mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       dimension indic(20),iq(20,20),iqi(20,20),iv(20,20)
       dimension irowp(20),irown(20),ic(20)
       dimension iran(20),nnumb1(4),nnumb2(4),nnumb3(4),nnumb4(4),nnumb5(4)
       dimension ninv(4),iprar(4),nsol(2),ixold(6,25),ixnew(6,25)
       dimension ipol1(6,25),ipol2(6,25),ipowr(5,6,3),jxold(6),jxnew(6)
       dimension kr(10000,10)
       


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




        if (ipd.lt.10000)goto 200
        iprar(1)=0
        iprar(2)=2
        iprar(3)=ipd/10000
        iprar(4)=ipd-iprar(3)*10000
        goto 300
200     iprar(1)=0       
        iprar(2)=1
        iprar(3)=ipd

300    if (igarray(1).lt.10000)goto 302
       nnumb1(1)=0
       nnumb1(2)=2
       nnumb1(3)=igarray(1)/10000
       nnumb1(4)=igarray(1)-nnumb1(3)*10000
       goto 304
302    nnumb1(1)=0
       nnumb1(2)=1
       nnumb1(3)=igarray(1)
304    if (igarray(2).lt.10000)goto 306
       nnumb2(1)=0
       nnumb2(2)=2
       nnumb2(3)=igarray(2)/10000
       nnumb2(4)=igarray(2)-nnumb2(3)*10000
       goto 308
306    nnumb2(1)=0
       nnumb2(2)=1
       nnumb2(3)=igarray(2)
308    if (igarray(3).lt.10000)goto 310
       nnumb3(1)=0
       nnumb3(2)=2
       nnumb3(3)=igarray(3)/10000
       nnumb3(4)=igarray(3)-nnumb3(3)*10000
       goto 312
310    nnumb3(1)=0
       nnumb3(2)=1
       nnumb3(3)=igarray(3)
312    do jf=1,iprar(2)+2
       kara(jf)=iprar(jf)
       end do
       
       itemp=2*igarray(1)
       itemp=mod(itemp,ipd)
       if (itemp.lt.10000)goto 314
       karb(3)=itemp/10000
       karb(4)=itemp-karb(3)*10000
       karb(2)=2
       karb(1)=0
       goto 316
314    karb(1)=0
       karb(2)=1
       karb(3)=itemp
316    call mpgcd
       do jf=1,karv(2)+2
       ninv(jf)=karv(jf)
       end do
       do jf=3,nnumb2(2)+2
       karr(jf-2)=nnumb2(jf)
       kbarr(jf-2)=nnumb2(jf)
       end do
       ilen=nnumb2(2)
       ilen2=nnumb2(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       nnumb4(jf+2)=kcarr(jf)
       end do
       nnumb4(2)=ilen3
       nnumb4(1)=0
       print *,'nnumb4',(nnumb4(jk),jk=1,nnumb4(2)+2)
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
       ncom(jf+2)=irrr(jf)
       end do
       ncom(1)=0
       ncom(2)=irlen
       if (isgn.eq.0)goto 318
       do jf=2,ncom(2)+2
       karr(jf)=ncom(jf)
       end do
       karr(1)=1


       do jf=1,iprar(2)+2
       kbarr(jf)=iprar(jf)
       end do
       call mpadd(0)
       do jf=1,kcarr(2)+2
       ncom(jf)=kcarr(jf)
       end do
318    print *,'ncom',(ncom(jk),jk=1,ncom(2)+2)
       
       call bwq5(ipd,ians)
       item=ipd-igarray(2)+ians
       item=mod(item,ipd)
       do ind1=1,2
       if (item.lt.10000)goto 320
       karr(1)=item/10000
       karr(2)=item-karr(1)*10000
       ilen=2
       goto 322
320    karr(1)=item 
       ilen=1
322    do jf=3,ninv(2)+2
       kbarr(jf-2)=ninv(jf)
       end do
       ilen2=ninv(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       do jf=3,iprar(2)+2
       kbarr(jf-2)=iprar(jf)
       end do
       ilen2=iprar(2)
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 324
       nsol(ind1)=0
       do jf=irlen,1,-1
       nsol(ind1)=nsol(ind1)+irrr(jf)*10000**(irlen-jf)
       end do
       goto 326
324    nsol(ind1)=0
326    item=2*ipd-igarray(2)-ians
       end do
       print *,'no. of sols. =2',nsol(1),nsol(2)
       isol(1)=nsol(1)
       isol(2)=nsol(2)
       
       goto 400









      

















       
       stop

       stop
       idegs=5
       ipd=43037
       
       nmul =1
       ngcd=0
       
       
       
       
58     nn=(ipd -1)/2
       
       ie =1
2      if (ie.gt.nn)goto 3
       ie = ie *2
       goto 2
3      ie = ie/2
       
       nn = nn-ie
       
       
       
       iran(1) = mod((797 *nmul),131072)
       nmul =mod(iran(1),ipd)
       
       

       
       
400    a=a              
       return
       end
       subroutine subgcd(idegs,idegb,idegg,ipd)
       
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       dimension iws(5000),iwb(5000),itempb(5000),mul(4)
       icon=0
       ipdc=2000000000/ipd
       
       do i= 1,idegs +1
       iws(i) = isarray(i)
       end do
       do i=1,idegb+1
       iwb(i) = ibarray(i)
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       ib=iws(1)
       if (ib.ne.0)goto 101
!      call subbw6(ipd,ib,iv)
       print *,'big big stop',(iws(jk),jk=1,iwsd)
       stop
       iv=1
       
       mulg=iv
       
       icon = icon+1
       ib=iws(1)
101    if (ipd.ge.45500)goto 15
       call subbw6(ipd,ib,iv)
       mulg=iv
       goto 501
       kara(1)=0
       kara(2)=1
       kara(3)=ipd
       goto 16
15     kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
16     if (ib.lt.10000)goto 50
       karb(1)=0
       karb(2)=2
       karb(3)=ib/10000
       karb(4)=ib-karb(3)*10000
       goto 52
50     karb(1)=0       
       karb(2)=1
       karb(3)=ib
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       if (mul(2).eq.1)goto 521
       mulg=mul(3)*10000+mul(4)
       
       goto 501
521    mulg=mul(3)
       

501    do i =1,loopl
       
       if (iwb(i).eq.0)goto 60
       if (mulg.ne.0)goto 5019
       print *,'worries'
       stop

5019   iqf=iwb(i)*mulg
       iqf=mod(iqf,ipd)
       if (iwb(i).gt.2100000000/mulg)goto 70
       
       
       if (ipd.lt.45500)goto 509
!      if (iwb(i).lt.ipdc)goto 502
!      iwb(i)=0
!      goto 502
       
70     if (iwb(i).lt.10000)goto 701
       karr(1)=iwb(i)/10000
       karr(2)=iwb(i)-karr(1)*10000
       ilen=2
       
       goto 502
701    karr(1)=iwb(i)
       ilen=1

502    iwb(i) =0
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       iqt(jf+2)=kcarr(jf)
       end do
       iqt(2)=ilen3
       iqt(1)=0
       
       if ((iqt(2).gt.2).and.(iqt(3).gt.20))goto 509
       if (iqt(2).eq.1)goto 72
       if (iqt(2).eq.3)goto 721
       iqf=iqt(3)*10000+iqt(4)
       goto 73
721    iqf=iqt(3)*100000000+iqt(4)*10000+iqt(5)
       iqf=mod(iqf,ipd)
       goto 73

72     iqf=iqt(3)
73     iqf=mod(iqf,ipd)
       if (iqf.ne.0)goto 509
       print *,'big trouble',(iqt(jk),jk=1,iqt(2)+2),'ipd=',ipd
       stop

509    iwb(i)=0
       do j =2,iwsd +1
!      goto 503
       if (ipd.lt.45500)goto 503
!      if (iwb(i).lt.ipdc)goto 503
       if (iws(j).eq.0)goto 503
       if ((iqt(2).gt.2).and.(iqt(3).gt.20))goto 78
       if (iqf.eq.0)goto 112
       if (iws(j).lt.2100000000/iqf)goto 503
       
       
78     if (iws(j).lt.10000)goto 711
       karr(1)=iws(j)/10000
       karr(2)=iws(j)-karr(1)*10000
       ilen=2 
       goto 712
112    print *,'big big big stop'       
       stop
711    karr(1)=iws(j)
       ilen=1
712    do jf=3,iqt(2)+2
       kbarr(jf-2)=iqt(jf)
       end do
       ilen2=iqt(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       kbarr(1)=ipd/10000
       kbarr(2)=ipd-kbarr(1)*10000
       ilen2=2
       
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       isum=0
       if(irlen.eq.0)goto 503
       if (irlen.eq.1)goto 601
       isum=irrr(1)*10000 +irrr(2)
       goto 603
601    isum=irrr(1)       
       goto 603
       


503    isum=iws(j)*iqf
       isum=mod(isum,ipd) 
       
       
603    iwb(i+j-1)=iwb(i+j-1)-isum
       iwb(i+j-1)=mod(iwb(i+j-1),ipd)
       
       
       if(iwb(i+j-1).ge.0)goto 12
       iwb(i+j-1)=iwb(i+j-1) +ipd
12     end do
       
       
60     end do
       
       
       

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
       print *,'iws',(iws(jf),jf=1,iwsd+1)
       
       
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
       subroutine subinv(idegs,idegb,idegg,ipd)
       
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       dimension iws(5000),iwb(5000),itempb(5000),mul(4)
       dimension iarrv1(5000),iarrt1(5000),iarrq(5000),iarrv(5000)
       dimension iarru(5000),isaveb(5000),isaves(5000),isaveg(5000)
       icon=0
       ipdc=2000000000/ipd
       idegu=0
       iarru(1)=1
       idegv1=-1
       iswdeg=0



       
       do i= 1,idegs +1
       iws(i) = isarray(i)
       end do
       do i=1,idegb+1
       iwb(i) = ibarray(i)
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       ib=iws(1)
       if (ib.ne.0)goto 101
!      call subbw6(ipd,ib,iv)
       print *,'big big stop',(iws(jk),jk=1,iwsd)
       stop
       iv=1
       
       mulg=iv
       
       icon = icon+1
       ib=iws(1)
101    if (ipd.ge.45500)goto 15
       call subbw6(ipd,ib,iv)
       mulg=iv
       goto 501
       
       
       
       
       
       kara(1)=0
       kara(2)=1
       kara(3)=ipd
       goto 16
15     kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
16     if (ib.lt.10000)goto 50
       karb(1)=0
       karb(2)=2
       karb(3)=ib/10000
       karb(4)=ib-karb(3)*10000
       goto 52
50     karb(1)=0       
       karb(2)=1
       karb(3)=ib
52     call mpgcd
       do jf=1,karv(2)+2
       mul(jf)=karv(jf)
       end do
       
       if (mul(2).eq.1)goto 521
       mulg=mul(3)*10000+mul(4)
       
       goto 501
521    mulg=mul(3)
       

501    do i =1,loopl
       iarrq(i)=0
       if (iwb(i).eq.0)goto 60
       if (mulg.ne.0)goto 5019
       print *,'worries'
       stop

5019   iqf=iwb(i)*mulg
       iqf=mod(iqf,ipd)
       iarrq(i)=iqf
!       print *,'mulg',mulg,'iwb',iwb(i),'iqf',iqf
       if (iwb(i).gt.2100000000/mulg)goto 70
       
       
       if (ipd.lt.45500)goto 509
!      if (iwb(i).lt.ipdc)goto 502
!      iwb(i)=0
!      goto 502
       
70     if (iwb(i).lt.10000)goto 701
       karr(1)=iwb(i)/10000
       karr(2)=iwb(i)-karr(1)*10000
       ilen=2
       
       goto 502
701    karr(1)=iwb(i)
       ilen=1

502    iwb(i) =0
       do jf=3,mul(2)+2
       kbarr(jf-2)=mul(jf)
       end do
       ilen2=mul(2)
       
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       iqt(jf+2)=kcarr(jf)
       end do
       iqt(2)=ilen3
       iqt(1)=0
       
       if ((iqt(2).gt.2).and.(iqt(3).gt.20))goto 509
       if (iqt(2).eq.1)goto 72
       if (iqt(2).eq.3)goto 721
       iqf=iqt(3)*10000+iqt(4)
       goto 73
721    iqf=iqt(3)*100000000+iqt(4)*10000+iqt(5)
       iqf=mod(iqf,ipd)
       goto 73

72     iqf=iqt(3)
73     iqf=mod(iqf,ipd)
       if (iqf.ne.0)goto 509
       print *,'big trouble',(iqt(jk),jk=1,iqt(2)+2),'ipd=',ipd
       stop

509    iwb(i)=0
       do j =2,iwsd +1
!      goto 503
       if (ipd.lt.45500)goto 503
!      if (iwb(i).lt.ipdc)goto 503
       if (iws(j).eq.0)goto 503
       if ((iqt(2).gt.2).and.(iqt(3).gt.20))goto 78
       if (iqf.eq.0)goto 112
       if (iws(j).lt.2100000000/iqf)goto 503
       
       
78     if (iws(j).lt.10000)goto 711
       karr(1)=iws(j)/10000
       karr(2)=iws(j)-karr(1)*10000
       ilen=2 
       goto 712
112    print *,'big big big stop'       
       stop
711    karr(1)=iws(j)
       ilen=1
712    do jf=3,iqt(2)+2
       kbarr(jf-2)=iqt(jf)
       end do
       ilen2=iqt(2)
       call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       kbarr(1)=ipd/10000
       kbarr(2)=ipd-kbarr(1)*10000
       ilen2=2
       
       call mpdiv(ilen,ilen2,irlen,icont,iswq)
       
       isum=0
       if(irlen.eq.0)goto 504
       if (irlen.eq.1)goto 601
       isum=irrr(1)*10000 +irrr(2)
       
       goto 603
601    isum=irrr(1)       
       
       goto 603
504    print *,'b up stop'   
       stop

503    isum=iws(j)*iqf
       isum=mod(isum,ipd) 
       
       
603    iwb(i+j-1)=iwb(i+j-1)-isum
       iwb(i+j-1)=mod(iwb(i+j-1),ipd)
       
       
       if(iwb(i+j-1).ge.0)goto 12
       iwb(i+j-1)=iwb(i+j-1) +ipd
12     end do
       
       
60     end do
       
       
       

       do i= 1,iwbd +1

       if(iwb(i).ne.0)goto 20
       end do
       iswdeg=1
       goto 290
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
!       print *,'iws',(iws(jf),jf=1,iwsd+1)
290    if (idegv1.eq.-1)goto 200
       do jf=1,idegv1+1
       mult1(jf)=iarrv1(jf)
       end do
       idegm=idegv1
       do jf=1,loopl
       mult2(jf)=iarrq(jf)
       end do
       idegn=loopl-1
       call multy(idegm,idegn,ipd)
       idegm=idegm+idegn
       do jf=1,idegm+1
       iarrt1(jf)=0
       end do
       
       do jf=1,idegm+1
       if (mult3(jf).eq.0)goto 201
       iarrt1(jf)=ipd-mult3(jf)
201    end do
       if (idegm.ge.idegu)goto 202
       do jf=idegm+1,1,-1
       iarrt1(jf+idegu-idegm)=iarrt1(jf)
       end do
       do jf=1,idegu-idegm
       iarrt1(jf)=0
       end do
       idegm=idegu
202    do jf=idegm+1,1,-1
       iarrt1(jf)=iarrt1(jf)+iarru(jf+idegu-idegm)
       iarrt1(jf)=mod(iarrt1(jf),ipd)
       if (iarrt1(jf).ge.0)goto 203
       iarrt1(jf)=iarrt1(jf)+ipd
203    end do
       do jf=1,idegm+1
       if (iarrt1(jf).ne.0)goto 204
       end do
       idegt1=-1
       goto 206
200    idegt1=idegu
       if (idegt1.eq.-1)goto 206
       
       do jf=1,idegt1+1
       iarrt1(jf)=iarru(jf)
       end do
       goto 206






204    jin=jf-1
       do jf=1+jin,idegm+1
       iarrt1(jf-jin)=iarrt1(jf)
       end do
       idegt1=idegm-jin
206    idegu=idegv1
       if (idegv1.eq.-1)goto 210
       do jf=1,idegv1+1
       iarru(jf)=iarrv1(jf)
       end do
210    idegv1=idegt1
       if (idegt1.eq.-1)goto 10
       do jf=1,idegt1+1
       iarrv1(jf)=iarrt1(jf)
       end do
!       print *,'q',(iarrq(jf),jf=1,loopl)
!       print *,'t1',(iarrt1(jf),jf=1,idegt1+1)
       if (idegu.eq.-1)goto 211
!       print *,'iarru',(iarru(jf),jf=1,idegu+1)
211    a=a       
       
       if (iswdeg.eq.1)goto 100
!       print *,'iwb',(iwb(jf),jf=1,iwbd+1)
       
       
       goto 10
100    idegg=iwsd
       do i=1,10
       igarray(i)=0
       end do
       do i = 1,iwsd +1
       igarray(i) = iws(i)
       end do
       if (idegg.gt.0)goto 400
!       print *,'igarray',(igarray(jf),jf=1,idegg+1)
       
       
       if (idegu.eq.-1)goto 220
       do jf=1,idegb+1
       mult1(jf)=ibarray(jf)
       end do
       idegm=idegb
       do jf=1,idegu+1
       mult2(jf)=iarru(jf)
       end do
       idegn=idegu
       call multy(idegm,idegn,ipd)
       idegm=idegm+idegn
       do jf=1,idegm
       iarrv(jf)=0
       end do
       iarrv(idegm+1)=igarray(idegg+1)
       do jf=1,idegm+1
       if (mult3(jf).eq.0)goto 300
       iarrv(jf)=iarrv(jf)+ipd-mult3(jf)
300    end do
       idegv=idegm
       iarrv(idegm+1)=mod(iarrv(idegm+1),ipd)
       goto 303
400    print *,'polynomials not co-prime gcd=',(igarray(jf),jf=1&
       ,idegg+1)
       stop
220    print *,'error to be investigated'
       stop

303    isavedb=idegb
       isaveds=idegs
       isavedg=idegg
       do jf=1,idegb+1
       isaveb(jf)=ibarray(jf)
       end do
       do jf=1,idegs+1
       isaves(jf)=isarray(jf)
       end do
       do jf=1,idegg+1
       isaveg(jf)=igarray(jf)
       end do
       do jf=1,idegv+1
       ibarray(jf)=iarrv(jf)
       end do
       idegb=idegv
       do jf=1,isaveds+1
       isarray(jf)=isaves(jf)
       end do
       idegs=isaveds
       call subbw4(idegs,idegb,idegg,ipd)
       do jf=1,idegb-idegs+1
       iarrv(jf)=iqt(jf)
       end do
       idegv=idegb-idegs
       print *,'iarrv',(iarrv(jf),jf=1,idegv+1)
!      reconstitute arrays        
       do jf=1,isaveds+1
       isarray(jf)=isaves(jf)
       end do
       idegs=isaveds
       do jf=1,isavedb+1
       ibarray(jf)=isaveb(jf)
       end do
       idegb=isavedb
       do jf=1,isavedg+1
       igarray(jf)=isaveg(jf)
       end do
       idegg=isavedg
       indiv=iarru(idegu+1)*ibarray(idegb+1)+iarrv(idegv+1)*&
       isarray(idegs+1)
       
       
       indiv=mod(indiv,ipd)
       print *,'indiv',indiv,'big',ibarray(idegb+1),'small',&
       isarray(idegs+1),'ucon',iarru(idegu+1),'vcon',iarrv(idegv+1)
       
       call subbw6(ipd,indiv,iv)
       indiv=iv
       print *,'indiv',indiv,'idegv',idegv,'idegu',idegu
       
       do jf=1,idegv+1
       iarrv(jf)=iarrv(jf)*indiv
       iarrv(jf)=mod(iarrv(jf),ipd)
       end do
       print *,'inverse',(iarrv(jf),jf=1,idegv+1)
       
       do jf=1,idegu+1
       iarru(jf)=iarru(jf)*indiv
       iarru(jf)=mod(iarru(jf),ipd)
       end do
       print *,'multiple of mod',(iarru(jf),jf=1,idegu+1)
       stop





       
       
       
       
       
       
       return
       end
       subroutine multy(idegm,idegn,ipd)
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),kbarb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       
       do i=1,idegm+idegn+1
       mult3(i) =0
       
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
       
       subroutine subbw4(idegs,idegb,idegg,ipd)
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ip(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
       dimension iws(5000),iwb(5000),itempb(5000),mmul(10)
!       print *,'ibar',(ibarray(jf),jf=1,idegb+1)
!       print *,'idegs',(isarray(jf),jf=1,idegs+1)
!      ipd =5
!      inv(1) =1
       


!      ibarray(1) =1
!      ibarray(2) =601
!      ibarray(3)=501
!      ibarray(4)=401
!      ibarray(5)=301
!      ibarray(6)=201
!      idegb =5
!      isarray(1) =4
!      isarray(2) =386
!      isarray(3) =494
!      isarray(4)=802
!      isarray(5)=301
!      idegs=4
      
      do i =1,idegs+1
      iws(i) =isarray(i)
      end do
      do i=1,idegb+1
      iwb(i) =ibarray(i)
      end do
      ib=iws(1)
      iwsd = idegs
      iwbd = idegb
      loopl = iwbd -iwsd +1
!      call subbw6(ipd,iws(1),iv)
!      mul = iv
101    if (ipd.ge.45500)goto 15
       call subbw6(ipd,ib,iv)
       mul=iv
       goto 501
       kara(1)=0
       kara(2)=1
       kara(3)=ipd
       goto 16
15     kara(1)=0
       kara(2)=2
       kara(3)=ipd/10000
       kara(4)=ipd-kara(3)*10000
16     if (ib.lt.10000)goto 50
       karb(1)=0
       karb(2)=2
       karb(3)=ib/10000
       karb(4)=ib-karb(3)*10000
       goto 52
50     karb(1)=0       
       karb(2)=1
       karb(3)=ib
52     call mpgcd
       do jf=1,karv(2)+2
       mmul(jf)=karv(jf)
       end do
       
       if (mmul(2).eq.1)goto 521
       mul=mmul(3)*10000+mmul(4)
       
       goto 501
521    mul=mmul(3)
501    a=a      
      do i =1,loopl
      iqt(i) =iwb(i) *mul
!      print *,'i',i,'iqt',iqt(i)
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
      irarray(jj) =iwb(i+jj-1)
      end do
      idegr = itempbd -1
      goto 110
100   do i=1,10
      irarray(1) = 0
      end do
      idegr =0
110   a=a
!     print *,'quos',(iqt(jf),jf=1,loopl)
!      print *,'rem',(irarray(jf),jf=1,idegr+1)
      idegg=idegr
      return
      end
      subroutine subbw3(idd,ipd)
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
      do i=1,idd
      print *,'firmatrix',(mat2(i,jf),jf=1,idd),'irhs2',irhs2(i)
      end do
      



      
      idd2 =idd 
      




      do i = 1,idd2
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
      goto 150
137   end do

      
      jx(j) =1 
      
      
      
      
      
      



150   end do
      do i =1,idd2
      print *,'matrix',mat1(i,1),mat1(i,2),mat1(i,3),mat1(i,4),&
      irhs1(i)
      print *,'marks',jmark(i),markr(i)
      end do
      
      do i = 1,idd2
      print *,jx(i)
      end do
      
      do i =1,idd2
      if (jx(i).eq.0)goto 152 
      print *,'matrix singular'
      ipd=1
      goto 200
152   end do
      
      isol(idd2) =irhs1(idd2)
      do j =1,idd-1
      ind =idd2 -j+1
      ind2 =ind-1
      im=jmark(ind2)
      isol(ind2)=irhs1(im)
      do jj =ind,idd2
      isol(ind2)=isol(ind2)-isol(jj) *mat1(im,jj)
      isol(ind2)=mod(isol(ind2),ipd)
      end do
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

       
       
      
      


   
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
      
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
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
      dimension kdum(1600),isub(1600)
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
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
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
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
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
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
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
      common ibarray(5000),isarray(5000),igarray(5000),inv(25)
      common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
      common mult5(5000)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      common ip(20)
      
      
      common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
      common mnum(50)
      common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
      common ipr(65000),norma(20),gpr(15000)
      common kara(50),karb(50),kard(50),karp(50),karv(50)
      
      common iarq(2),ncom(20),irarray(5000),jpol(100),jpol2(100),iqt(5000)
      
      
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





      
    

      
      
      
      
       subroutine bwq5(ip,ians)
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(10),jpol2(10),iqt(5000)
       dimension itempz(2)
       
       
       
       do jf=3,ncom(2)+2
       karr(jf-2)=ncom(jf)
       end do
       ilen=ncom(2)
       if (ip.lt.10000)goto 1100
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 1102
1100   kbarr(1)=ip
       ilen2=1
1102   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 1104
       if (irlen.eq.1)goto 1106
       iaa=irrr(1)*10000+irrr(2)
       goto 1108
1104   print *,'factor found early',ip
       stop
1106   iaa=irrr(1)


       
       
1108   ix = 0
       
       
       
       
       
10     iprecod =ip -1
       i=0
26     itemp = int(iprecod/2)
       irem1 =iprecod -itemp*2
       if(itemp.eq.0)goto 46
       if(irem1.gt.0)goto 40
       iprecod = itemp
       i =i+1
       if(i.lt.200)goto 26
40     iq =iprecod
       ie = i
       goto 48
46     iq =1
       ie = i
48     i =1
       
       n = 1
52     n = n*607
54     itemp=int(n/1000)
56     irem1 =n-1000 *itemp
       n = irem1
       id = n
       if (ip.lt.10000)goto 5901
       karp(1)=0
       karp(2)=2
       karp(3)=int(ip/10000)
       karp(4)=ip-karp(3)*10000
       goto 5902
5901   karp(1)=0       
       karp(2)=1
       karp(3)=ip
5902   kard(1)=0
       kard(2)=1
       kard(3)=id

       call mpkron(k)
       
       
       if(k.eq.-1)goto 68
       
       
       i =i +1
       if(i.lt.1000)goto 52
68     ipn =iq
       
       
       
5601   iaas = n
       call sub516(ibprod,iaas,ipn,ip)
       iz = ibprod
       
       
       iy = iz
       ir = ie
       ipn = (iq-1)/2
       iaas = iaa
       call sub516(ibprod,iaas,ipn,ip)
       ix = ibprod
       print *,ix
       if (ix.lt.10000)goto 140
       karr(1)=int(ix/10000)
       karr(2)=ix-karr(1)*10000
       ilen=2
       goto 142
140    karr(1)=ix       
       ilen=1
142    if (iaa.lt.10000)goto 150       
       kbarr(1)=int(iaa/10000)
       kbarr(2)=iaa-kbarr(1)*10000
       ilen2=2
       goto 152
150    kbarr(1)=iaa
       ilen2=1
152    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 160
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 162
160    kbarr(1)=ip
       ilen2=1
162    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 169
       do jf=1,irlen
       karr(jf)=irrr(jf)
       itempz(jf)=irrr(jf)
       end do
       itempzl=irlen
       ilen=irlen
       goto 170
169    print *,'error type 2'
       stop
170    if (ix.lt.10000)goto 172
       kbarr(1)=int(ix/10000)
       kbarr(2)=ix-kbarr(1)*10000
       ilen2=2
       goto 174
172    kbarr(1)=ix
       ilen2=1
174    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 181
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 182
181    kbarr(1)=ip
       ilen2=1
182    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 169
       if (irlen.eq.1)goto 190
       ib=irrr(1)*10000+irrr(2)
       goto 192
190    ib=irrr(1)
192    if (itempzl.eq.1)goto 196
       ix=itempz(1)*10000+itempz(2)
       goto 100
196    ix=itempz(1)














       
       
       
       
100    irem1 =mod(ib,ip)
       
       if(irem1.eq.1)goto 200
       i=1
       m=1
110    ipow =2**m
       goto 700
112    irem2 =ibprod
       
       if(irem2.eq.1)goto 130
       m =m+1
       goto 110
130    if(m.eq.ir)goto 180
       ipow =2**(ir-m-1)
       goto 800
134    it = ibprod
       
       if (it.lt.10000)goto 2000
       karr(1)=int(it/10000)
       karr(2)=it-karr(1)*10000
       kbarr(1)=karr(1)
       kbarr(2)=karr(2)
       ilen=2
       ilen2=2
       goto 2002
2000   karr(1)=it       
       kbarr(1)=it
       ilen=1
       ilen2=1
2002   call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 2004
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000 
       ilen2=2
       goto 2006
2004   kbarr(1)=ip
       ilen2=1
2006   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2010
       iy=irrr(1)*10000+irrr(2)
       goto 2012
2008   print *,'error type 3'
       stop
2010   iy=irrr(1)
2012   ir=mod(m,ip)
       
       
       if (ix.lt.10000)goto 2014
       karr(1)=int(ix/10000)
       karr(2)=ix-karr(1)*10000
       ilen=2
       goto 2016
2014   karr(1)=ix
       ilen=1
2016   if (it.lt.10000)goto 2018
       kbarr(1)=int(it/10000)
       kbarr(2)=it-kbarr(1)*10000
       ilen2=2
       goto 2020
2018   kbarr(1)=it
       ilen2=1

2020   call mpmul(ilen,ilen2,ilen3)
       
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       
       if (ip.lt.10000)goto 2022
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 2024
2022   kbarr(1)=ip
       ilen2=1
2024   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2026
       ix=irrr(1)*10000+irrr(2)
       goto 2028
2026   ix=irrr(1)
2028   if (ib.lt.10000)goto 2030
       
       karr(1)=int(ib/10000)
       karr(2)=ib-karr(1)*10000
       ilen=2
       goto 2032
2030   karr(1)=ib
       ilen=1
2032   if (iy.lt.10000)goto 2034
       kbarr(1)=int(iy/10000)
       kbarr(2)=iy-kbarr(1)*10000
       ilen2=2
       goto 2036
2034   kbarr(1)=iy
       ilen2=1
2036   call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 2038
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 2040
2038   kbarr(1)=ip
       ilen2=1
   

2040   call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 2008
       if (irlen.eq.1)goto 2042
       ib=irrr(1)*10000+irrr(2)
       
       goto 100
2042   ib=irrr(1)
       
       goto 100



       
       
       
180    print *,'no square root exists'
       print *,'m',m,'ir',ir,'ib',ib,'ip',ip
       stop
       goto 900
200    print *,'square root=',ix
       ians=ix
       
       
       goto 900
700    ipn =ipow
       if (ipn.gt.10000000)goto 950
       
       iaas =ib
       call sub516(ibprod,iaas,ipn,ip)
       goto 112
800    ipn =ipow
       iaas=iy
       call sub516(ibprod,iaas,ipn,ip)
       goto 134
900    goto 1000       
       print *,'loop too short'
       stop
950    ians=999999
1000   return       
       end
       subroutine sub400(id,ip,k)
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
       
       subroutine sub516(ibprod,iaas,ipn,ip)
       common ibarray(5000),isarray(5000),igarray(5000),inv(25)
       common mult1(5000),mult2(5000),mult3(5000),mult4(5000)
       common mult5(5000)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       common ipp(20)
       
       
       
       common karr(1600),kbarr(1600),kcarr(3200),ipqt(1600),irrr(1600)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(20),irarray(5000),jpol(10),jpol2(10),iqt(5000)
       
       iy=iaas
       n=ipn
       ipow=2
14     if (ipow.gt.n)goto 20
       ipow=ipow*2
       goto 14
20     ie=int(ipow/2)
       n1=n
       n1=n1-ie
26     if (ie.eq.1)goto 100
       ie=int(ie/2)
       if (iy.lt.10000)goto 200
       karr(1)=int(iy/10000)
       karr(2)=iy-karr(1)*10000
       ilen=2
       kbarr(1)=karr(1)
       kbarr(2)=karr(2)
       ilen2=2
       goto 202
200    karr(1)=iy       
       ilen=1
       kbarr(1)=iy
       ilen2=1

202    call mpmul(ilen,ilen2,ilen3)
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 220
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 222
220    kbarr(1)=ip
       ilen2=1
222    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 230
       if (irlen.eq.1)goto 240
       iy=irrr(1)*10000+irrr(2)
       goto 242
230    print *,'error number 2'
       iy=0
240    iy=irrr(1)
242    if (n1.lt.ie)goto 26
       n1=n1-ie
       if (iy.lt.10000)goto 250
       karr(1)=int(iy/10000)
       karr(2)=iy-karr(1)*10000
       ilen=2
       goto 252
250    karr(1)=iy
       ilen=1
252    if (iaas.lt.10000)goto 260       
       kbarr(1)=int(iaas/10000)
       kbarr(2)=iaas-kbarr(1)*10000
       ilen2=2
       goto 262
260    kbarr(1)=iaas       
       ilen2=1
262    call mpmul(ilen,ilen2,ilen3)        
       do jf=1,ilen3
       karr(jf)=kcarr(jf)
       end do
       ilen=ilen3
       if (ip.lt.10000)goto 270
       kbarr(1)=int(ip/10000)
       kbarr(2)=ip-kbarr(1)*10000
       ilen2=2
       goto 272
270    kbarr(1)=ip
       ilen2=1
272    call mpdiv(ilen,ilen2,irlen,icont,iswq)
       if (irlen.eq.0)goto 230
       if (irlen.eq.1)goto 280
       iy=irrr(1)*10000+irrr(2)
       goto 26
280    iy=irrr(1)
       goto 26

       
       
100    ibprod=iy
       
620    return
       end
