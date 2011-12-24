      

      program bobalgm
      character*200 istring,ostring
      character(200)gstring
      character*1 loc(200)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      dimension idrv5(5),idrv4(5),idrv3(5),idrv2(5),idrv1(5)
      dimension match(1000),iv(1000)
      dimension ksm(20,1000),ktemp(5,20),ktempl(5),kres(100),kdenres(100)
      dimension iyew(6,500),nar(100),nar2(200),nartot(200),ivee(6,500)
      dimension kktem(2,100)
      iyes=0
      ino=1
      ivdeg=0
      ivee(1,1)=0
      ivee(1,2)=1
      ivee(1,3)=1
      open (unit=5,file='compdata',access='direct',form=&
      'formatted',recl=212,status='old')
      read (5,398,rec=1)kia,kib,isgn,(nar(i),i=1,70),&
      (nar2(j),j=1,119)
398   format (i6,i6,i1,70i1,119i1,10i1)      
      goto 66


      iyew(1,1)=0
      iyew(1,2)=1
      iyew(1,3)=7
      iyew(2,1)=0
      iyew(2,2)=1
      iyew(2,3)=8
      kpoldeg=3
      kmpol(1,1)=0
      kmpol(1,2)=1
      kmpol(1,3)=1
      kmpol(2,1)=0
      kmpol(2,2)=0
      kmpol(3,1)=0
      kmpol(3,2)=0
      kmpol(4,1)=1
      kmpol(4,2)=1
      kmpol(4,3)=2

! for dividing algebraic integers      
      idegs=1
      isarray(1,1)=0
      isarray(1,2)=1
      isarray(1,3)=2
      isarray(2,1)=0
      isarray(2,2)=1
      isarray(2,3)=3
      isden(1,1)=0
      isden(1,2)=1
      isden(1,3)=1
      isden(2,1)=0
      isden(2,2)=1
      isden(2,3)=1
      idegb=2
      ibarray(1,1)=0
      ibarray(1,2)=1
      ibarray(1,3)=6
      ibarray(2,1)=0
      ibarray(2,2)=1
      ibarray(2,3)=7
      ibarray(3,1)=0
      ibarray(3,2)=1
      ibarray(3,3)=8
      ibden(1,1)=0
      ibden(1,2)=1
      ibden(1,3)=9
      ibden(2,1)=0
      ibden(2,2)=1
      ibden(2,3)=10
      ibden(3,1)=0
      ibden(3,2)=1
      ibden(3,3)=11
      do i=1,kpoldeg+1
      do jf=1,kmpol(i,2)+2
      ibarray(i,jf)=kmpol(i,jf)
      end do
      end do
      do i=1,kpoldeg+1
     ibden(i,1)=0
     ibden(i,2)=1
     ibden(i,3)=1
      end do
     idegb=kpoldeg 
      
      call subalgd(idegs,idegb,idegg)
      print *,'iqt1',(iqt(1,jf),jf=1,iqt(1,2)+2)
      print *,'iqtden1',(iqtden(1,jf),jf=1,iqtden(1,2)+2)
      print *,'iqt2',(iqt(2,jf),jf=1,iqt(2,2)+2)
      print *,'iqtden2',(iqtden(2,jf),jf=1,iqtden(2,2)+2)
      print *,'idegg',idegg,'irarray',((irarray(i,jf),jf=1,irarray(i,2)+2)&
      ,i=1,idegg+1),'irden',((irden(i,jf),jf=1,irden(i,2)+2),i=1,idegg+1)
      do jf=1,irarray(1,2)+2
      kres(jf)=irarray(1,jf)
      end do
      do jf=1,irden(1,2)+2
      kdenres(jf)=irden(1,jf)
      end do
      
      kdegs=idegb-idegs
      
      do i=1,kdegs+1
      do jf=1,iqt(i,2)+2
      kdoub1(i,jf)=iqt(i,jf)
      end do
      do jf=1,iqtden(i,2)+2
      kden1(i,jf)=iqtden(i,jf)
      end do
      end do
      kdegb=1
      do i=1,2
      do jf=1,iyew(i,2)+2
      kdoub2(i,jf)=iyew(i,jf)
      end do
      end do
      kdoub2(1,1)=mod(iyew(1,1)+1,2)
      kdoub2(2,1)=mod(iyew(2,1)+1,2)
      kden2(1,1)=0
      kden2(1,2)=1
      kden2(1,3)=1
      kden2(2,1)=0
      kden2(2,2)=1
      kden2(2,3)=1


      
      
      
      
      call subalgm(kdegs,kdegb,kpoldeg,kdegg)
      idegg=kdegg
      print *,'iqt1',(iqt(1,jf),jf=1,iqt(1,2)+2)
      print *,'iqtden1',(iqtden(1,jf),jf=1,iqtden(1,2)+2)
      print *,'iqt2',(iqt(2,jf),jf=1,iqt(2,2)+2)
      print *,'iqtden2',(iqtden(2,jf),jf=1,iqtden(2,2)+2)
      print *,'idegg',idegg,'irarray',((irarray(i,jf),jf=1,irarray(i,2)+2)&
      ,i=1,idegg+1),'irden',((irden(i,jf),jf=1,irden(i,2)+2),i=1,idegg+1)
      print *,'kres',(kres(jf),jf=1,kres(2)+2),'kresden',&
      (kdenres(jf),jf=1,kdenres(2)+2)
      do i=1,idegg+1
      do jf=1,irarray(i,2)+2
      marr(jf)=irarray(i,jf)
      end do
      do jf=1,kdenres(2)+2
      mbarr(jf)=kdenres(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      irarray(i,jf)=mcarr(jf)
      end do
      do jf=1,irden(i,2)+2
      marr(jf)=irden(i,jf)
      end do
      do jf=1,kres(2)+2
      mbarr(jf)=kres(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      irden(i,jf)=mcarr(jf)
      end do
      irarray(i,1)=mod(irarray(i,1)+irden(i,1),2)
      irden(i,1)=0
      do jf=1,irarray(i,2)+2
      kara(jf)=irarray(i,jf)
      end do
      do jf=1,irden(i,2)+2
      karb(jf)=irden(i,jf)
      end do
      kara(1)=0
      karb(1)=0
      call subgcd2
      do jf=1,irarray(i,2)+2
      marr(jf)=irarray(i,jf)
      end do
      do jf=1,kard(2)+2
      mbarr(jf)=kard(jf)
      end do
      call mendiv
      do jf=2,mdarr(2)+2
      irarray(i,jf)=mdarr(jf)
      end do
      do jf=1,irden(i,2)+2
      marr(jf)=irden(i,jf)
      end do
      call mendiv
      do jf=2,mdarr(2)+2
      irden(i,jf)=mdarr(jf)
      end do
      end do
      print *,'idegg',idegg,'irarray',((irarray(i,jf),jf=1,irarray(i,2)+2)&
      ,i=1,idegg+1),'irden',((irden(i,jf),jf=1,irden(i,2)+2),i=1,idegg+1)

      
      
      
      stop
      
      
      
66    a=a       
      
      irecnn=0
      ipoin=0
      open(unit=1,file='kernel',access='direct',form=&
      'formatted',recl=1000,status='old')
      mm=200
      read(1,1,rec=3)(iv(ii),ii=1,mm)
      close(unit=1)
      
      open(unit=2,file='matches',access='direct',form=&
      'formatted',recl=2000,status='old')
      open(unit=3,file='square',access='direct',form=&
      'formatted',recl=2000,status='old')
      open(unit=4,file='minpol',access='direct',form=&
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
      ia0(3)=10
      ia0(4)=2813
      ia1(1)=0
      ia1(2)=1
      ia1(3)=10
      ia1(4)=282
      ia2(1)=0
      ia2(2)=1
      ia2(3)=3
      ia2(4)=392
      ia3(1)=0
      ia3(2)=1
      ia3(3)=9
      ia3(4)=843
      ia4(1)=0
      ia4(2)=1
      ia4(3)=8
      ia4(4)=7110
      ia5(1)=1
      ia5(2)=1
      ia5(3)=43
      ia5(4)=1487
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
      
1     format(200i5)
2     format(400i5)
      do i=1,mm
      if (iv(i).ne.0)goto 3
      end do
      print *,'file error'
      stop
3     read(2,2,rec=1)(match(ii),ii=1,mm*2)
      read (5,398,rec=i)kia,kib,isgn,(nar(ij),ij=1,70),&
      (nartot(j),j=1,119)
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
      do ix=1,2
      do jf=1,kdoub1(ix,2)+2
      iyew(ix,jf)=kdoub1(ix,jf)
      end do
      end do
      print *,'idegm',idegm,kdoub1(1,1),kdoub1(1,2),kdoub1(1,3)
      icc=1
      do i=ione+1,mm
      if (iv(i).eq.0)goto 43
      icc=icc+1
      read (5,398,rec=i)kia,kib,isgn,(nar(ij),ij=1,70),&
      (nar2(j),j=1,119)
          
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
      kpoldeg=5
!      goto 60
      do ix=20,119
      if (nartot(ix)-nar2(ix).lt.0)goto 60
      end do
      print *,'i',i,'icc',icc
      
      kpoldeg=5
      do ix=1,kpoldeg+1
      do jf=1,kmpol(ix,2)+2
      ibarray(ix,jf)=kmpol(ix,jf)
      end do
      end do
      do ix=1,kpoldeg+1
     ibden(ix,1)=0
     ibden(ix,2)=1
     ibden(ix,3)=1
      end do
     idegb=kpoldeg 
     idegs=1
     do ix=1,2
     do jf=1,kdoub2(ix,2)+2
     isarray(ix,jf)=kdoub2(ix,jf)
     end do 
     end do 
     do ix=1,2
     isden(ix,1)=0
     isden(ix,2)=1
     isden(ix,3)=1
     end do
     call subalgd(idegs,idegb,idegg) 
      do jf=1,irarray(1,2)+2
      kres(jf)=irarray(1,jf)
      end do
      do jf=1,irden(1,2)+2
      kdenres(jf)=irden(1,jf)
      end do
      
      kdegs=idegb-idegs
      
      do ix=1,kdegs+1
      do jf=1,iqt(ix,2)+2
      kdoub1(ix,jf)=iqt(ix,jf)
      end do
      do jf=1,iqtden(ix,2)+2
      kden1(ix,jf)=iqtden(ix,jf)
      end do
      
      end do
      
      kdegb=idegm
      do ix=1,2
      do jf=1,kdoub2(ix,2)+2
      kktem(ix,jf)=kdoub2(ix,jf)
      end do
      end do
      do ix=1,idegm+1
      do jf=1,iyew(ix,2)+2
      kdoub2(ix,jf)=iyew(ix,jf)
      end do
      end do
      do ix=1,idegm+1
      
      kdoub2(ix,1)=mod(iyew(ix,1)+1,2)
      
      kden2(ix,1)=0
      kden2(ix,2)=1
      kden2(ix,3)=1
      end do


      
      
      
      
      call subalgm(kdegs,kdegb,kpoldeg,kdegg)
      idegg=kdegg
      print *,'iqt1',(iqt(1,jf),jf=1,iqt(1,2)+2)
      print *,'iqtden1',(iqtden(1,jf),jf=1,iqtden(1,2)+2)
      print *,'iqt2',(iqt(2,jf),jf=1,iqt(2,2)+2)
      print *,'iqtden2',(iqtden(2,jf),jf=1,iqtden(2,2)+2)
      print *,'idegg',idegg,'irarray',((irarray(ix,jf),jf=1,irarray(ix,2)+2)&
      ,ix=1,idegg+1),'irden',((irden(ix,jf),jf=1,irden(ix,2)+2),ix=1,idegg+1)
      print *,'kres',(kres(jf),jf=1,kres(2)+2),'kresden',&
      (kdenres(jf),jf=1,kdenres(2)+2)
      do ix=1,idegg+1
      do jf=1,irarray(ix,2)+2
      marr(jf)=irarray(ix,jf)
      end do
      do jf=1,kdenres(2)+2
      mbarr(jf)=kdenres(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      irarray(ix,jf)=mcarr(jf)
      end do
      do jf=1,irden(ix,2)+2
      marr(jf)=irden(ix,jf)
      end do
      do jf=1,kres(2)+2
      mbarr(jf)=kres(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      irden(ix,jf)=mcarr(jf)
      end do
      irarray(ix,1)=mod(irarray(ix,1)+irden(ix,1),2)
      irden(ix,1)=0
      do jf=1,irarray(ix,2)+2
      kara(jf)=irarray(ix,jf)
      end do
      do jf=1,irden(ix,2)+2
      karb(jf)=irden(ix,jf)
      end do
      kara(1)=0
      karb(1)=0
      call subgcd2
      do jf=1,irarray(ix,2)+2
      marr(jf)=irarray(ix,jf)
      end do
      do jf=1,kard(2)+2
      mbarr(jf)=kard(jf)
      end do
      call mendiv
      do jf=2,mdarr(2)+2
      irarray(ix,jf)=mdarr(jf)
      end do
      do jf=1,irden(ix,2)+2
      marr(jf)=irden(ix,jf)
      end do
      call mendiv
      do jf=2,mdarr(2)+2
      irden(ix,jf)=mdarr(jf)
      end do
      end do
      print *,'idegg',idegg,'irarray',((irarray(ix,jf),jf=1,irarray(ix,2)+2)&
      ,ix=1,idegg+1),'irden',((irden(ix,jf),jf=1,irden(ix,2)+2),ix=1,idegg+1)
      
      do ix=1,idegg+1
      if ((irden(ix,2).eq.1).and.(irden(ix,3).eq.1))goto 61 
      goto 62
61    end do       
      idegm=idegg
      itemdeg=idegm
      do ix=1,idegg+1
      do jf=1,irarray(ix,2)+2
      iyew(ix,jf)=irarray(ix,jf)
      end do
      end do
      do ix=1,ivdeg+1
      do jf=1,ivee(ix,2)+2
      kdoub1(ix,jf)=ivee(ix,jf)
      end do
      end do
      do ix=1,2
      do jf=1,kktem(ix,2)+2
      kdoub2(ix,jf)=kktem(ix,jf)
      end do
      end do
      idegm=ivdeg
      idegn=1
      idegp=5
      call doubmul(idegm,idegn,idegp,idprod)
      ivdeg=idprod
      idegm=itemdeg
      do ii=1,ivdeg+1
      do jj=1,kdoub3(ii,2)+2
      ivee(ii,jj)=kdoub3(ii,jj)
      end do
      end do
      do kf=1,ivdeg+1
      print *,'ivee',(ivee(kf,jf),jf=1,ivee(kf,2)+2)
      end do
      do jf=1,119
      nartot(jf)=nartot(jf)-nar2(jf)
      end do
      print *,'division','idegm',idegm
      iyes=iyes+1
      goto 43
62    do ix=1,2
      do jf=1,kktem(ix,2)+2
      kdoub2(ix,jf)=kktem(ix,jf)
      end do
      end do
60    a=a      
      do ix=1,idegm+1
      do jf=1,iyew(ix,2)+2
      kdoub1(ix,jf)=iyew(ix,jf)
      end do
      end do
      idegn=1
      idegp=5
      
      call doubmul(idegm,idegn,idegp,idprod)
      idegm=idprod
      do ii=1,idegm+1
      do jj=1,kdoub3(ii,2)+2
      kdoub1(ii,jj)=kdoub3(ii,jj)
      iyew(ii,jj)=kdoub3(ii,jj)
      end do
      end do
      do kf=1,idegm+1
      print *,'k',(kdoub1(kf,jf),jf=1,kdoub1(kf,2)+2)
      end do
      do jf=1,119
      nartot(jf)=nartot(jf)+nar2(jf)
      end do
      print *,'no division','idegm',idegm
      ino=ino+1
      
43    end do
      
      print *,'icc',icc,'idprod',idprod
      do ix=1,idegm+1
      print *,'ix',ix,'iyew',(iyew(ix,jf),jf=1,iyew(ix,2)+2)
      end do
      
      do ix=1,ivdeg+1
      print *,'ix',ix,'ivee',(ivee(ix,jf),jf=1,ivee(ix,2)+2)
      end do
      print *,'yesses',iyes,'noes',ino
      print *,'nartot',(nartot(jf),jf=1,119)
      
      itemdeg=idegm
      idegm=ivdeg
      idegn=ivdeg
      idegp=kpoldeg
      do ix=1,kpoldeg
      do jf=1,ivee(ix,2)+2
      kdoub1(ix,jf)=ivee(ix,jf)
      kdoub2(ix,jf)=ivee(ix,jf)
      end do
      end do
      call doubmul(idegm,idegn,idegp,idprod)
      idegm=idprod
      do ii=1,idegm+1
      do jj=1,kdoub3(ii,2)+2
      kdoub1(ii,jj)=kdoub3(ii,jj)
      end do
      end do
      idegn=itemdeg
      do ix=1,idegn+1
      do jf=1,iyew(ix,2)+2
      kdoub2(ix,jf)=iyew(ix,jf)
      end do
      end do
      idegp=kpoldeg
      call doubmul(idegm,idegn,idegp,idprod)
      idegm=idprod
      print *,'kdoub3 1',(kdoub3(1,jf),jf=1,kdoub3(1,2)+2)
      idegm=itemdeg





      
      
      idegsm=idprod
      idegm=itemdeg
      do ii=1,idegm+1
      do jj=1,iyew(ii,2)+2
      ksm(ii,jj)=iyew(ii,jj)
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


      idegsm=idegm
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
      print *,'icc',icc,'yesses',iyes,'noes',ino
      print *,'all written'

500   close(unit=2)
      close(unit=1)
      open (unit=1,file='veelit',access='direct',form=&
      'formatted',recl=2000,status='old')
      do i=1,idprod+1
      do jf=1,500
      kdoub3(i,jf)=0
      end do
      end do
      do i=1,idprod+1
      do jf=1,ivee(i,2)+2
      kdoub3(i,jf)=ivee(i,jf)
      end do
      end do
      do i=1,idprod+1
      write(1,1003,rec=i)(kdoub3(i,jf),jf=1,500)
      end do
      close (unit=1)



      
      
      
      end
      


      subroutine mpprime(icorp,inlen)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
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
      goto 22
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

      
      
      subroutine doubmul(idegm,idegn,idegp,idprod)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      
      
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
     
      
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
      if(kdoub3(i,2).ne.0)goto 70
      end do
      print *,'error algebraic nos'
      stop
70    idprod=idegg+1-i
      igap=i-1
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
       
       

       
       subroutine menmul
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)

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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
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
       
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
       dimension iws(20,50),iwb(20,50),itempb(20,50),mul(50),ib(50)
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)

       
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
       
       dimension iws(20,50),iwb(20,50),itempb(20,50),mul(50),ib(50)
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
       subroutine subalgd(idegs,idegb,idegg)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
       
       dimension iws(20,50),iwb(20,50),itempb(20,50),mul(50),ib(50)
       dimension iwsden(20,50),iwbden(20,50),mulden(50)
       dimension itemp1(100),itemp2(100),itemp3(100),itemp4(100),itemp5(100)

!       print *,'ibar',((ibarray(jf,jk),jk=1,ibarray(jf,2)+2),&
!       jf=1,idegb+1)
!       print *,'isar',((isarray(jf,jk),jk=1,isarray(jf,2)+2),&
!       jf=1,idegs+1)
       
       
       
       do i= 1,idegs +1
       do jf=1,isarray(i,2)+2
       
       iws(i,jf) = isarray(i,jf)
       end do
       do jf=1,isden(i,2)+2
       iwsden(i,jf)=isden(i,jf)
       end do
       
       
       
       end do
       do i=1,idegb+1
       do jf=1,ibarray(i,2)+2
       iwb(i,jf) = ibarray(i,jf)
       end do
       do jf=1,ibden(i,2)+2
       iwbden(i,jf)=ibden(i,jf)
       end do
       
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       if (loopl.lt.1)goto 220
       do jf=1,iws(1,2)+2
       mulden(jf)=iws(1,jf)

       end do
       do jf=1,iwsden(1,2)+2
       mul(jf)=iwsden(1,jf)
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
       kara(jf+2)=kcarr(jf)
       end do
       isgnn=mod(mul(1)+iwb(i,1),2)
       kara(1)=0
       kara(2)=ilen3
       do jf=1,iwbden(i,2)+2
       marr(jf)=iwbden(i,jf)
       end do
       do jf=1,mulden(2)+2
       mbarr(jf)=mulden(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karb(jf)=mcarr(jf)
       end do
       isgnd=mod(mulden(1)+iwbden(i,1),2)
       karb(1)=0
       call subgcd2
       print *,'kard',(kard(jf),jf=1,kard(2)+2)
       do jf=1,kara(2)+2
       marr(jf)=kara(jf)
       end do
       do jf=1,kard(2)+2
       mbarr(jf)=kard(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       iqt(i,jf)=mdarr(jf)
       end do
       iqt(i,1)=isgnn
       do jf=1,karb(2)+2
       marr(jf)=karb(jf)
       end do
       do jf=1,kard(2)+2
       mbarr(jf)=kard(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       iqtden(i,jf)=mdarr(jf)
       end do
       iqtden(i,1)=isgnd
       iqt(i,1)=mod(iqt(i,1)+iqtden(i,1),2)
       iqtden(i,1)=0





!       print *,'iqti',i,(iqt(i,jf),jf=1,iqt(i,2)+2)
       
509    iwb(i,1)=0
       iwb(i,2)=0
       iwbden(i,1)=0
       iwbden(i,2)=0
       
       

       do j =2,iwsd +1
       do jf=1,iws(j,2)+2
       marr(jf)=iws(j,jf)
       end do
       do jf=1,iqt(i,2)+2
       mbarr(jf)=iqt(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       itemp1(jf)=mcarr(jf)
       end do
       do jf=1,iwsden(j,2)+2
       marr(jf)=iwsden(j,jf)
       end do
       do jf=1,iqtden(i,2)+2
       mbarr(jf)=iqtden(i,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       itemp2(jf)=mcarr(jf)
       end do
       do jf=1,iwb(i+j-1,2)+2
       marr(jf)=iwb(i+j-1,jf)
       end do
       do jf=1,itemp2(2)+2
       mbarr(jf)=itemp2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       itemp3(jf)=mcarr(jf)
       end do
       do jf=1,itemp1(2)+2
       marr(jf)=itemp1(jf)
       end do
       do jf=1,iwbden(i+j-1,2)+2
       mbarr(jf)=iwbden(i+j-1,jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       itemp4(jf)=mcarr(jf)
       end do
       do jf=1,itemp3(2)+2
       karr(jf)=itemp3(jf)
       end do
       do jf=1,itemp4(2)+2
       kbarr(jf)=itemp4(jf)
       end do
       call mpadd(1)
       do jf=1,kcarr(2)+2
       itemp5(jf)=kcarr(jf)
       kara(jf)=kcarr(jf)
       end do
       if (itemp5(2).eq.0)goto 200
!  numerator unreduced in itemp5
       do jf=1,iwbden(i+j-1,2)+2
       marr(jf)=iwbden(i+j-1,jf)
       end do
       do jf=1,itemp2(2)+2
       mbarr(jf)=itemp2(jf)
       end do
       call menmul
       do jf=1,mcarr(2)+2
       karb(jf)=mcarr(jf)
       end do
       kara(1)=0
       
       call subgcd2
       print *,'kara2',(kara(jf),jf=1,kara(2)+2)
       print *,'karb2',(karb(jf),jf=1,karb(2)+2)
       print *,'kard2',(kard(jf),jf=1,kard(2)+2)
       do jf=1,kara(2)+2
       marr(jf)=kara(jf)
       end do
       do jf=1,kard(2)+2
       mbarr(jf)=kard(jf)
       end do
       call mendiv
! putting reduced result back into iwb and iwbden
       print *,'mdarr for iwb',(mdarr(jf),jf=1,mdarr(2)+2)
       do jf=1,mdarr(2)+2
       iwb(i+j-1,jf)=mdarr(jf)
       end do
       iwb(i+j-1,1)=itemp5(1)
       do jf=1,karb(2)+2
       marr(jf)=karb(jf)
       end do
       do jf=1,kard(2)+2
       mbarr(jf)=kard(jf)
       end do
       call mendiv
       do jf=1,mdarr(2)+2
       iwbden(i+j-1,jf)=mdarr(jf)
       end do
       print *,'i',i,'j',j,'iwb',(iwb(i+j-1,jk),jk=1,iwb(i+j-1,2)+2)
       print *,'iwbden',(iwbden(i+j-1,jk),jk=1,iwbden(i+j-1,2)+2)
       
       
      


      goto 201
200   iwb(i+j-1,1)=0
      iwb(i+j-1,2)=0
      iwbden(i+j-1,1)=0
      iwbden(i+j-1,2)=0
       
       
201    end do
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
      do jf=1,iwbden(jj+i-1,2)+2
      irden(jj,jf)=iwbden(jj+i-1,jf)
      end do
      
      
      end do
      idegr = itempbd -1
      goto 110
220   do i=1,idegb+1
      do jf=1,ibarray(i,2)+2
      irarray(i,jf)=ibarray(i,jf)
      end do
      do jf=1,ibden(i,2)+2
      irden(i,jf)=ibden(i,jf)
      end do
      end do
      idegr=idegb
      iqt(1,1)=0
      iqt(1,2)=0
      iqtden(1,1)=0
      iqtden(1,2)=1
      iqtden(1,3)=1
      goto 110



100   idegr =0
110   a=a
!     print *,'quos',(iqt(jf),jf=1,loopl)
      print *,'rem',((irarray(jf,jk),jk=1,irarray(jf,2)+2),jf=1,idegr+1)
      idegg=idegr
      
      return
      end
      subroutine subalgm(kdegs,kdegb,kpoldeg,kdegg)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      dimension itemp1(100),itemp2(100),itemp3(100),itemp4(100),itemp5(100)
      
      idegs=kdegs
      idegb=kdegb
      print *,'idegs',idegs,'idegb',idegb
      do i=1,kdegs+1
      do jf=1,kdoub1(i,2)+2
      isarray(i,jf)=kdoub1(i,jf)
      end do
      do jf=1,kden1(i,2)+2
      isden(i,jf)=kden1(i,jf)
      end do
      
      
      end do
      do i=1,kdegb+1
      do jf=1,kdoub2(i,2)+2
      ibarray(i,jf)=kdoub2(i,jf)
      end do
      do jf=1,kden2(i,2)+2
      ibden(i,jf)=kden2(i,jf)
      end do
      end do



      
      
      do i=1,12 
      kdoub3(i,1)=0
      kdoub3(i,2)=0 
      kden3(i,1)=0 
      kden3(i,2)=1
      kden3(i,3)=1
      end do
      do i=1,idegs+1
      do j=1,idegb+1
      do jf=1,isarray(i,2)+2
      marr(jf)=isarray(i,jf)
      end do
      do jf=1,ibarray(j,2)+2
      mbarr(jf)=ibarray(j,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      itemp1(jf)=mcarr(jf)
      end do
      
      do jf=1,isden(i,2)+2
      marr(jf)=isden(i,jf)
      end do
      do jf=1,ibden(j,2)+2
      mbarr(jf)=ibden(j,jf)
      end do
      
      call menmul
      do jf=1,mcarr(2)+2
      itemp2(jf)=mcarr(jf)
      end do
      print *,'itemp2',(itemp2(jf),jf=1,itemp2(2)+2)
      do jf=1,itemp2(2)+2
      marr(jf)=itemp2(jf)
      end do
      do jf=1,kdoub3(i+j-1,2)+2
      mbarr(jf)=kdoub3(i+j-1,jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      itemp3(jf)=mcarr(jf)
      end do
      do jf=1,itemp1(2)+2
      marr(jf)=itemp1(jf)
      end do
      do jf=1,kden3(i+j-1,2)+2
      mbarr(jf)=kden3(i+j-1,jf)
      end do

      call menmul
      do jf=1,mcarr(2)+2
      itemp4(jf)=mcarr(jf)
      end do

      do jf=1,itemp3(2)+2
      karr(jf)=itemp3(jf)
      end do
      do jf=1,itemp4(2)+2
      kbarr(jf)=itemp4(jf)
      end do
      call mpadd(0)
      if (kcarr(2).eq.0)goto 1
      do jf=1,kcarr(2)+2
      itemp5(jf)=kcarr(jf)
      kara(jf)=kcarr(jf)
      end do

      kara(1)=0
      do jf=1,kden3(i+j-1,2)+2
      marr(jf)=kden3(i+j-1,jf)
      end do
      do jf=1,itemp2(2)+2
      mbarr(jf)=itemp2(jf)
      end do
      call menmul
      do jf=1,mcarr(2)+2
      karb(jf)=mcarr(jf)
      end do
      print *,'here we go','karb',(karb(jf),jf=1,karb(2)+2)
      call subgcd2
      print *,'ok1',' kard',(kard(jf),jf=1,kard(2)+2)
      do jf=1,kara(2)+2
      marr(jf)=kara(jf)
      end do
      do jf=1,kard(2)+2
      mbarr(jf)=kard(jf)
      end do
      print *,'mbarr2',mbarr(2)
      call mendiv
      
      print *,'ok2'
      do jf=1,mdarr(2)+2
      kdoub3(i+j-1,jf)=mdarr(jf)
      end do
      
      kdoub3(i+j-1,1)=itemp5(1)
      do jf=1,karb(2)+2
      marr(jf)=karb(jf)
      end do
      print *,'mbarr2',mbarr(2)
      call mendiv
      print *,'ok3'
      do jf=1,mdarr(2)+2
      kden3(i+j-1,jf)=mdarr(jf)
      end do
      goto 2
1     kdoub3(i+j-1,1)=0 
      kdoub3(i+j-1,2)=0
      kden3(i+j-1,1)=0
      kden3(i+j-1,2)=1
      kden3(i+j-1,3)=1
2     end do
      end do
      print *,'ok4'
      do i=1,kdegs+kdegb+1
      do jf=1,kdoub3(i,2)+2
      ibarray(i,jf)=kdoub3(i,jf)
      end do
      do jf=1,kden3(i,2)+2
      ibden(i,jf)=kden3(i,jf)
      end do
      end do
      do i=1,kpoldeg+1
      do jf=1,kmpol(i,2)+2
      isarray(i,jf)=kmpol(i,jf)
      end do
      end do
      do i=1,kpoldeg+1
      isden(i,1)=0
      isden(i,2)=1
      isden(i,3)=1
      end do
      idegs=kpoldeg
      idegb=kdegs+kdegb
      print *,'mul idegs',idegs,'idegb',idegb
      do i=1,idegb+1
!      print *,'i',i,'ibarray',(ibarray(i,jf),jf=1,ibarray(i,2)+2)
!      print *,'i',i,'ibden',(ibden(i,jf),jf=1,ibden(i,2)+2)
!      print *,'i',i,'isarray',(isarray(i,jf),jf=1,isarray(i,2)+2)
!      print *,'i',i,'isden',(isden(i,jf),jf=1,isden(i,2)+2)
      end do
      
      
      
      
      
      
      
      call subalgd(idegs,idegb,idegg)
      kdegg=idegg
      
      
      
      return
      end














      
      
      subroutine mpmul(ilen,ilen2,ilen3)
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      
      
      
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      
      
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      
      
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
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
       do jf=1,ilen2
       kard(jf+2)= kbarr(jf)
       end do
       return
       end
       

       subroutine mpgcd
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
      
      
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
       common ibarray(20,50),isarray(20,50),igarray(20,50),inv(25)
       common ibden(20,50),isden(20,50),igden(20,50),irden(20,50)
       common mult1(20,50),mult2(20,50),mult3(20,50),mult4(20,50)
       common mult5(20,50)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(10,50)
       common ipp(20)
       common karr(200),kbarr(200),kcarr(400),ipqt(200),irrr(200)
       common mnum(50)
       common ia5(5),ia4(5),ia3(5),ia2(5),ia1(5),ia0(5),m1(5),m2(5)
       common ipr(65000),norma(20),gpr(15000)
       common kara(50),karb(50),kard(50),karp(50),karv(50)
       common iarq(2),ncom(2,100),irarray(20,50),jpol(10,50),jpol2(10,50)
       common iqt(10,50),ipd(50),iqtden(20,50)
       common marr(200),mbarr(200),mcarr(400),mdarr(200),nnh(50),ihalf(50)
!       common ncom(2,100),iprar(50),iansar(50),narc(50)
       common ip(100)
       common iaas(2,100),ipn(100),iacn(2,100),ibcn(2,100)
       common ibprod(2,100),icprod(2,100)
       common nsq,isqurar(2,2,100)
      common kmpol(6,20)
      common kdoub1(6,500),kdoub2(6,500),kdoub3(12,1000)
      common kden1(20,50),kden2(20,50),kden3(20,50)
       
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
       
       
       
       
       
       
       
       
       
       
                   
       
       
       
          
       
       
       
       



      
      
