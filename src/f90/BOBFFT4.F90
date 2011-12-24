       program fftmp4
!      Version of FFT multiplication which prevents errors from   
!      propogating

       double precision uni(66000,2),fta(66000,2),ftb(66000,2)
       double precision kcarr(66000,2),ftc(66000,2) 
       double precision triar(2,66000,2),wrt(60),wrt2(60)
       double precision argr,pi,x1,tempr,tempc
!       dimension uni(40000,2),fta(40000,2),ftb(40000,2)
!       dimension kcarr(40000,2),ftc(40000,2)
!       dimension triar(2,40000,2),wrt(60),wrt2(60)
       
       
       dimension karr(66000),kbarr(66000),numb(66000)
       dimension inn1(66000),inn2(66000),ibarr(60),ibarr2(60),ibarr3(60)
       dimension ibarr4(60),indar(60)
       
       
!       pi = 3.14159265389324d0
       pi=3.14159
!       print *,pi
!      5123*1473       
       do i=1,66000
       karr(i) =0
       kbarr(i)=0
       kcarr(i,1)=0
       kcarr(i,2)=0
       end do
       inn1(1)=0
       inn1(2)=8
       inn1(3)=1
       inn1(4)=2
       inn1(5)=3
       inn1(6)=4
       inn1(7)=4
       inn1(8)=5
       inn1(9)=6
       inn1(10)=7
       inn1(11)=8
       inn1(12)=9
       inn1(13)=1
       inn1(14)=2
       inn1(15)=3
       inn1(16)=4
       inn1(17)=6
       inn1(18)=3
       inn1(19)=1


       inn2(1)=0
       inn2(2)=8
       inn2(3)=2
       inn2(4)=3
       inn2(5)=4
       inn2(6)=5
       inn2(7)=2
       inn2(8)=2
       inn2(9)=6
       inn2(10)=5
       inn2(11)=3
       inn2(12)=2
       inn2(13)=1
       inn2(14)=1
       inn2(15)=4
       inn2(16)=5
       inn2(17)=3
       inn2(18)=1
       inn2(19)=7
       inn2(20)=3
       inn2(21)=3
       inn2(22)=1
       inn2(23)=3
       inn2(24)=5
       inn2(25)=3
       inn2(26)=2
       inn2(27)=8
       inn2(28)=2
       inn2(29)=3
       inn2(30)=8
       inn2(31)=3
       do ibbb=1,10
       
!       goto 2

       
       inn1(1)=0
       inn1(2)=32768
       inn2(1)=0
       inn2(2)=32768
       do jf=3,inn1(2)+2
!       inn1(jf)=1
!       inn2(jf)=2
       inn1(jf)=111
       inn2(jf)=222
       end do
2      ilens=inn1(2)
       do jf=3,ilens+2
       karr(ilens+3-jf)=inn1(jf)
       end do
       ilens=inn2(2)
       do jf=3,ilens+2
       kbarr(ilens+3-jf)=inn2(jf)
       end do
       goto 1

       
       
       
       
       
       
       
1      ilen =32768*2
       klen=16
       wrt(1)=-1.0d0
       wrt2(1)=0.0d0
       wrt(2)=0.0d0
       wrt2(2)=1.0d0
       do i=3,klen
       wrt(i)=((1.0+wrt(i-1))/2.0)**0.5
       wrt2(i)=((1-wrt(i-1))/2.0)**0.5
       end do
       do i=1,60
       indar(i)=0
       end do
       do ibig=1,ilen-1
       do j=1,60
       if (indar(j).eq.0)goto 400
       indar(j)=0
       end do
       goto 402
400    indar(j)=1
402    jfar=j
       tempr=1.0d0
       tempc=0.0d0
!       print *,'indar',(indar(jf),jf=1,klen)
       

       do jj=1,klen
       if (indar(jj).eq.0)goto 404
       uni(ibig,1)=tempr*wrt(klen+1-jj)-tempc*wrt2(klen+1-jj)
       uni(ibig,2)=tempr*wrt2(klen+1-jj)+tempc*wrt(klen+1-jj)
       tempr=uni(ibig,1)
       tempc=uni(ibig,2)
404    end do
       end do
       uni(ilen,1)=1.0d0
       uni(ilen,2)=0.0d0
       do ii=1,ilen
       
       
!       print *,'uni ii ',uni(ii,1),uni(ii,2)
       end do
       
       goto 3




       
       












       
       
       goto 3
7      uni(1,1)=0.707106781562044d0
       
       


3      do i =1,ilen
       fta(i,1)= 0.0d0
       ftb(i,1)= 0.0d0
       ftc(i,1) =0.0d0
       fta(i,2)=0.0d0
       ftb(i,2)=0.0d0
       ftc(i,2)=0.0d0
       
       
       end do
       
       isw1=1
       do jf=1,ilen
       triar(1,jf,1)=karr(jf)
       triar(1,jf,2)=0.0d0
       end do
       goto 99
       do i=1,ilen
       do j=1,ilen
       ipow=mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 10
       ipow=ipow+ilen
10     fta(i,1)=fta(i,1)+uni(ipow,1)*karr(j)
       fta(i,2)=fta(i,2)+uni(ipow,2)*karr(j)
       ftb(i,1)=ftb(i,1)+uni(ipow,1)*kbarr(j)
       ftb(i,2)=ftb(i,2)+uni(ipow,2)*kbarr(j)
       end do
!       print *,'fta',fta(i,1),fta(i,2)
       end do
       isw1=2
       goto 250


99     do k=1,klen
       do jjk=1,60
       ibarr(jjk)=0
       end do
       do jf=1,ilen
!       jz=jf-1
       
       
       
       do jm=1,klen
!      ibarr(jm)=0 
       ibarr2(jm)=0
       ibarr3(jm)=0
       ibarr4(jm)=0
       end do
       ind1=1
!       jq=jz
!  136    jp=jq/2
!       itemp=jq-jp*2
!       ibarr(ind1)=itemp
!       jq=jp
!       if (jp.eq.0)goto 138 
!       ind1=ind1+1
!       goto 136
        if (jf.eq.1)goto 140
        do jjr=1,60
        if (ibarr(jjr).eq.0)goto 136
        ibarr(jjr)=0
        end do
        goto 138
136     ibarr(jjr)=1
        
138     do jxx=klen,1,-1
        if (ibarr(jxx).eq.0)goto 139
        ind1=jxx
        goto 140
139     end do
140     do jv=1,ind1       
       ibarr2(jv)=ibarr(jv)
       ibarr3(jv)=ibarr(jv)
       ibarr4(jv)=0
       end do
       inj=klen-k+1
       ibarr2(inj)=0
       ibarr3(inj)=1
       do jv=1,k
       ibarr4(klen-jv+1)=ibarr(klen-k+jv)
       end do
       ind31=0
       ind32=0
       indn=0
       ipow=1
       do jv=1,klen
       
       ind31=ind31+ibarr2(jv)*ipow
       ind32=ind32+ibarr3(jv)*ipow
       indn=indn+ibarr4(jv)*ipow
       ipow=ipow*2
       end do
!       print *,'ibarr4',(ibarr4(jy),jy=1,klen)
!       print *,'ibarr',(ibarr(jy),jy=1,klen),'jjr',jjr,'jf',jf,'ind1',ind1
       
       if (mod(indn,ilen).eq.0)goto 2000
       indp=mod(indn,ilen)
       goto 2002
2000   indp=ilen

2002   triar(2,jf,1)=triar(1,ind31+1,1)+uni(indp,1)*triar(1,ind32+1,1)&
       -uni(indp,2)*triar(1,ind32+1,2)
       triar(2,jf,2)=triar(1,ind31+1,2)+uni(indp,1)*triar(1,ind32+1,2)&
       +uni(indp,2)*triar(1,ind32+1,1)
       i31=ind31+1
       i32=ind32+1
!       print *,'jf',jf,'i31',i31,'indp',indp,'i32',i32,'k',k
!       print *,'ff4 fir triar',triar(2,jf,1),triar(2,jf,2)
       
       if (isw1.ne.3)goto 2005
       if (jf.ne.8)goto 2005
       if (k.ne.klen)goto 2005
       
       
       
2005   if (k.ne.klen)goto 991
       kcarr(indn+1,1)=triar(2,jf,1)
       kcarr(indn+1,2)=triar(2,jf,2)
!       print *,'triar',triar(2,jf,1),triar(2,jf,2)
       



991    end do
       
       
       do jf=1,ilen
       do kf=1,2
       triar(1,jf,kf)=triar(2,jf,kf)
       end do
       end do
       
       
       end do
!       print *,'fft4 kcarrs',(kcarr(jf,1),jf=1,ilen)
        
       
       
       
       if (isw1.eq.2)goto 200
       if (isw1.eq.3)goto 300
       
       isw1=2
       
       
       do jf=1,ilen
       do kf=1,2
       fta(jf,kf)=kcarr(jf,kf)
       end do
       end do
       do jf=1,ilen
       
       triar(1,jf,1)=kbarr(jf)
       triar(1,jf,2)=0.0d0
       end do
       goto 99

200    isw1=3       
!       print *,'kcars',(kcarr(jf,1),jf=1,ilen)
!       print *,'kcar2',(kcarr(jf,2),jf=1,ilen)
       
       do jf=1,ilen
       do kf=1,2
       ftb(jf,kf)=kcarr(jf,kf)
       end do
       end do
250    isw1=3
       do i=1,ilen
       ftc(i,1) =fta(i,1) *ftb(i,1)- fta(i,2)*ftb(i,2)
       ftc(i,2) =fta(i,1)*ftb(i,2)+fta(i,2)*ftb(i,1)
       end do
!       print *,'ftc1',(ftc(jf,1),jf=1,ilen)
!       print *,'ftc2',(ftc(jf,2),jf=1,ilen)
       do jf=1,ilen
       triar(1,jf,1)=ftc(jf,1)
       triar(1,jf,2)=ftc(jf,2)
       end do
       
       goto 99
       do i=1,ilen
       kcarr(i,1)=0.0
       end do
       
       do i=1,ilen
       do j=1,ilen
       ipow=mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 20
       ipow=ilen
20     kcarr(i,1)=kcarr(i,1)+uni(ipow,1) *ftc(j,1)-uni(ipow,2)&
       *ftc(j,2)
       end do
       end do

300    do i=1,ilen-1
       if (i.ne.1)goto 32
       compt =kcarr(1,1)
32     kcarr(i,1)=kcarr(i+1,1)/ilen
       end do
       kcarr(ilen,1) =compt/ilen
!       print *,'kcarrs',(kcarr(jf,1),jf=1,ilen)
       
       
       rum=kcarr(16,1) +0.5
       num=int(rum)
!       print *,'num',num
       do jf=1,ilen
       numb(jf)=0
       
       end do
       do jf=ilen,2,-1
       rum=kcarr(jf,1)+0.5
       numb(jf)=int(rum)
       end do
       
       
       do jf=ilen,2,-1
       itemp=numb(jf)
!       numb(jf)=mod(itemp,10)
!       numb(jf-1)=numb(jf-1)+int(itemp/10)
       numb(jf)=mod(itemp,1000)
       numb(jf-1)=numb(jf-1)+int(itemp/1000)
       end do
       
       do jf=ilen,1,-1
       numb(jf+2)=numb(jf)
       end do
       numb(1)=0
       do jf=3,ilen+2
       if (numb(jf).ne.0)goto 40
       end do
       print *,'number zero'
       stop
40     ilen2=jf-2        
       do jk=ilen2+2,ilen+2
       numb(jk-ilen2+1)=numb(jk)
       end do
       numb(2)=ilen+1-ilen2
       end do
       print *,'numb',(numb(jk),jk=1,numb(2)+2)
!       print *,'fft4 karr',(karr(jf),jf=1,ilen)
!       print *,'fft4 kbarr',(kbarr(jf),jf=1,ilen)
       end
